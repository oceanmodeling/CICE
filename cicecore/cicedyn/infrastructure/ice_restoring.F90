!=======================================================================
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! authors: Elizabeth C. Hunke, LANL

      module ice_restoring

      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_constants, only: c0, c1, c2, p2
      use ice_domain_size, only: ncat, max_blocks, n_aero
      use ice_forcing, only: trestore, trest, &
          aicen_bry, vicen_bry, vsnon_bry,    &
          Tsfc_bry, Tinz_bry, Sinz_bry, alvln_bry,vlvln_bry,&
          apondn_bry, hpondn_bry, ipondn_bry,iage_bry, Tsnz_bry, &
          uvel_bry, vvel_bry !pedrocice 
      use ice_state, only: aicen, vicen, vsnon, trcrn, bound_state, &
          aice_init, aice0, aice, vice, vsno, trcr, trcr_depend, &
          uvel, vvel,& 
          divu,shear,strength! !pedrocice  
          
          
      use icepack_tracers, only: ntrcr,tr_pond_lvl,nbtrcr  
      
      use ice_timers, only: ice_timer_start, ice_timer_stop, timer_bound
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_diag
      use icepack_intfc, only: icepack_warnings_flush, icepack_warnings_aborted
      use icepack_intfc, only: icepack_init_trcr
      use icepack_intfc, only: icepack_query_parameters, &
          icepack_query_tracer_sizes, icepack_query_tracer_flags, &
          icepack_query_tracer_indices
      
      use ice_domain, only:sea_ice_time_bry
      
      implicit none
      private
      public :: ice_HaloRestore_init, ice_HaloRestore

      logical (kind=log_kind), public :: &
         restore_ice                 ! restore ice state if true
      
      real (kind=dbl_kind), dimension (:,:,:), allocatable :: & !pedrocice
         uvel_rest , & ! ice velocity
         vvel_rest 
      !-----------------------------------------------------------------
      ! state of the ice for each category
      !-----------------------------------------------------------------
    
      real (kind=dbl_kind), dimension (:,:,:,:), allocatable, public :: &
         aicen_rest , & ! concentration of ice
         vicen_rest , & ! volume per unit area of ice          (m)
         vsnon_rest     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (:,:,:,:,:), allocatable, public :: &
         trcrn_rest     ! tracers

!=======================================================================

      contains

!=======================================================================

!  Allocates and initializes arrays needed for restoring the ice state
!  in cells surrounding the grid.


 subroutine ice_HaloRestore_init

      use ice_blocks, only: block, get_block, nblocks_x, nblocks_y
      use ice_communicate, only: my_task, master_task
      use ice_domain, only: ew_boundary_type, ns_boundary_type, &
          nblocks, blocks_ice
      use ice_grid, only: tmask, hm
      use ice_flux, only: Tf, Tair, salinz, Tmltz
      use ice_restart_shared, only: restart_ext

   integer (int_kind) :: &
     i,j,iblk,nt,n,      &! dummy loop indices
     ilo,ihi,jlo,jhi,    &! beginning and end of physical domain
     iglob(nx_block),    &! global indices
     jglob(ny_block),    &! global indices
     iblock, jblock,     &! block indices
     ibc,                &! ghost cell column or row
     ntrcr,              &!
     npad                 ! padding column/row counter

   character (len=7), parameter :: &
!     restore_ic = 'defined' ! otherwise restore to initial ice state
     restore_ic = 'initial' ! restore to initial ice state

   type (block) :: &
     this_block  ! block info for current block

   character(len=*), parameter :: subname = '(ice_HaloRestore_init)'

   if (.not. restore_ice) return

   call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
   call icepack_warnings_flush(nu_diag)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
      file=__FILE__, line=__LINE__)

   if ((ew_boundary_type == 'open' .or. &
        ns_boundary_type == 'open') .and. .not.(restart_ext)) then
      if (my_task == master_task) write (nu_diag,*) 'ERROR: restart_ext=F and open boundaries'
      call abort_ice(error_message=subname//'open boundary and restart_ext=F', &
         file=__FILE__, line=__LINE__)
   endif

   allocate (aicen_rest(nx_block,ny_block,ncat,max_blocks), &
             vicen_rest(nx_block,ny_block,ncat,max_blocks), &
             vsnon_rest(nx_block,ny_block,ncat,max_blocks), &
             trcrn_rest(nx_block,ny_block,ntrcr,ncat,max_blocks),&
             uvel_rest(nx_block,ny_block,max_blocks),&  !pedrocice
             vvel_rest(nx_block,ny_block,max_blocks))

!    aicen_rest(:,:,:,:) = c0
!    vicen_rest(:,:,:,:) = c0
!    vsnon_rest(:,:,:,:) = c0
!    trcrn_rest(:,:,:,:,:) = c0
   aicen_rest(:,:,:,:) = c0!aicen(:,:,:,:)
   vicen_rest(:,:,:,:) = c0!vicen(:,:,:,:)
   vsnon_rest(:,:,:,:) = c0!vsnon(:,:,:,:)
   trcrn_rest(:,:,:,:,:) = c0!trcrn(:,:,:,:,:)
   uvel_rest(:,:,:)  = c0
   vvel_rest(:,:,:)  = c0

!-----------------------------------------------------------------------
! initialize
! halo cells have to be filled manually at this stage
! these arrays could be set to values read from a file...
!-----------------------------------------------------------------------

   if (trim(restore_ic) == 'defined') then

      ! restore to defined ice state
      !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
      !$OMP                     iglob,jglob,iblock,jblock)
      do iblk = 1, nblocks
         this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi
         iglob = this_block%i_glob
         jglob = this_block%j_glob
         iblock = this_block%iblock
         jblock = this_block%jblock

         call set_restore_var (nx_block,            ny_block,            &
                               ilo, ihi,            jlo, jhi,            &
                               iglob,               jglob,               &
                               iblock,              jblock,              &
                               Tair (:,:,    iblk), &
                               Tf   (:,:,    iblk),                      &
                               salinz(:,:,:, iblk), Tmltz(:,:,:,  iblk), &
                               tmask(:,:,    iblk),                      &
                               aicen_rest(:,:,  :,iblk), &
                               trcrn_rest(:,:,:,:,iblk), ntrcr,         &
                               vicen_rest(:,:,  :,iblk), &
                               vsnon_rest(:,:,  :,iblk))
      enddo ! iblk
      !$OMP END PARALLEL DO

   else  ! restore_ic

   ! restore to initial ice state

! the easy way
!   aicen_rest(:,:,:,:) = c0!aicen(:,:,:,:)
!   vicen_rest(:,:,:,:) = c0!vicen(:,:,:,:)
!   vsnon_rest(:,:,:,:) = c0!vsnon(:,:,:,:)
!   trcrn_rest(:,:,:,:,:) = c0!trcrn(:,:,:,:,:)
!   uvel_rest(:,:,:,:)  = c0
!   vvel_rest(:,:,:,:)  = c0
! the more precise way
   !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
   !$OMP                     i,j,n,nt,ibc,npad)
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, ny_block
            do i = 1, ilo
               aicen_rest(i,j,n,iblk) = aicen(ilo,j,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(ilo,j,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(ilo,j,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ilo,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = 1, ny_block
            do i = ihi, ibc
               aicen_rest(i,j,n,iblk) = aicen(ihi,j,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(ihi,j,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(ihi,j,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ihi,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, jlo
            do i = 1, nx_block
               aicen_rest(i,j,n,iblk) = aicen(i,jlo,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(i,jlo,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(i,jlo,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ilo,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = jhi, ibc
            do i = 1, nx_block
               aicen_rest(i,j,n,iblk) = aicen(i,jhi,n,iblk)
               vicen_rest(i,j,n,iblk) = vicen(i,jhi,n,iblk)
               vsnon_rest(i,j,n,iblk) = vsnon(i,jhi,n,iblk)
               do nt = 1, ntrcr
                  trcrn_rest(i,j,nt,n,iblk) = trcrn(ihi,j,nt,n,iblk)
               enddo
            enddo
            enddo
            enddo
         endif
      endif

   enddo ! iblk
   !$OMP END PARALLEL DO

   endif ! restore_ic

      !-----------------------------------------------------------------
      ! Impose land mask
      !-----------------------------------------------------------------

   do iblk = 1, nblocks
      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen_rest(i,j,n,iblk) = aicen_rest(i,j,n,iblk) * hm(i,j,iblk)
            vicen_rest(i,j,n,iblk) = vicen_rest(i,j,n,iblk) * hm(i,j,iblk)
            vsnon_rest(i,j,n,iblk) = vsnon_rest(i,j,n,iblk) * hm(i,j,iblk)
            do nt = 1, ntrcr
               trcrn_rest(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk) &
                                                            * hm(i,j,iblk)
            enddo
         enddo
         enddo
      enddo
   enddo

   if (my_task == master_task) &
      write (nu_diag,*) 'ice restoring timescale = ',trestore,' days'

 end subroutine ice_HaloRestore_init

!=======================================================================

! initialize restoring variables, based on set_state_var
! this routine assumes boundaries are not cyclic

    subroutine set_restore_var (nx_block, ny_block, &
                                ilo, ihi, jlo, jhi, &
                                iglob,    jglob,    &
                                iblock,   jblock,   &
                                Tair, &
                                Tf,                 &
                                salinz,   Tmltz,    &
                                tmask,    aicen,    &
                                trcrn,    ntrcr,    &
                                vicen,    vsnon)

! authors: E. C. Hunke, LANL

      use ice_arrays_column, only: hin_max
      use ice_blocks, only: nblocks_x, nblocks_y
      use icepack_intfc, only: icepack_init_trcr
      use icepack_parameters, only: c0, c1, c2, p2, p5, rhoi, rhos, Lfresh, &
           cp_ice, cp_ocn, Tsmelt, Tffresh
      use ice_domain_size, only: nilyr, nslyr, ncat

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo, ihi          , & ! physical domain indices
         jlo, jhi          , & !
         iglob(nx_block)   , & ! global indices
         jglob(ny_block)   , & !
         iblock            , & ! block indices
         jblock            , & !
         ntrcr                 ! number of tracers in use

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         Tair    , & ! air temperature  (K)
         Tf          ! freezing temperature (C)

      real (kind=dbl_kind), dimension (nx_block,ny_block,nilyr), intent(in) :: &
         salinz  , & ! initial salinity profile
         Tmltz       ! initial melting temperature profile

      logical (kind=log_kind), dimension (nx_block,ny_block), intent(in) :: &
         tmask      ! true for ice/ocean cells

      real (kind=dbl_kind), dimension (nx_block,ny_block,ncat), intent(out) :: &
         aicen , & ! concentration of ice
         vicen , & ! volume per unit area of ice          (m)
         vsnon     ! volume per unit area of snow         (m)

      real (kind=dbl_kind), dimension (nx_block,ny_block,ntrcr,ncat), intent(out) :: &
         trcrn     ! ice tracers
                   ! 1: surface temperature of ice/snow (C)

      ! local variables

      integer (kind=int_kind) :: &
         i, j        , & ! horizontal indices
         ij          , & ! horizontal index, combines i and j loops
         ibc         , & ! ghost cell column or row
         npad        , & ! padding column/row counter
         k           , & ! ice layer index
         n           , & ! thickness category index
         it          , & ! tracer index
         nt_Tsfc     , & !
         nt_fbri     , & !
         nt_qice     , & !
         nt_sice     , & !
         nt_qsno     , & !
         icells          ! number of cells initialized with ice

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with restoring
   
         
      logical (kind=log_kind) :: &
         tr_brine

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! compressed indices for cells with restoring

      real (kind=dbl_kind) :: &
         Tsfc, hbar, &
         hsno_init       ! initial snow thickness

      real (kind=dbl_kind), dimension(ncat) :: &
         ainit, hinit    ! initial area, thickness

      real (kind=dbl_kind), dimension(nilyr) :: &
         qin             ! ice enthalpy (J/m3)

      real (kind=dbl_kind), dimension(nslyr) :: &
         qsn             ! snow enthalpy (J/m3)

      character(len=*), parameter :: subname = '(set_restore_var)'
      
      
      
      
     
      
      
      call icepack_query_tracer_flags(tr_brine_out=tr_brine)
      call icepack_query_tracer_indices(nt_Tsfc_out=nt_Tsfc, nt_fbri_out=nt_fbri, &
           nt_qice_out=nt_qice, nt_sice_out=nt_sice, nt_qsno_out=nt_qsno)
      call icepack_warnings_flush(nu_diag)
      if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
         file=__FILE__, line=__LINE__)

      indxi(:) = 0
      indxj(:) = 0

      !-----------------------------------------------------------------
      ! Initialize restoring variables everywhere on grid
      !-----------------------------------------------------------------

      do n = 1, ncat
         do j = 1, ny_block
         do i = 1, nx_block
            aicen(i,j,n) = c0
            vicen(i,j,n) = c0
            vsnon(i,j,n) = c0
            if (tmask(i,j)) then
               trcrn(i,j,nt_Tsfc,n) = Tf(i,j)  ! surface temperature
            else
               trcrn(i,j,nt_Tsfc,n) = c0  ! on land gridcells
            endif
            if (ntrcr >= 2) then
               do it = 2, ntrcr
                  trcrn(i,j,it,n) = c0
               enddo
            endif
            if (tr_brine) trcrn(i,j,nt_fbri,n) = c1
         enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      ! initial area and thickness in ice-occupied restoring cells
      !-----------------------------------------------------------------

      hbar = c2  ! initial ice thickness
      hsno_init = 0.20_dbl_kind ! initial snow thickness (m)
      do n = 1, ncat
         hinit(n) = c0
         ainit(n) = c0
         if (hbar > hin_max(n-1) .and. hbar < hin_max(n)) then
            hinit(n) = hbar
            ainit(n) = 0.95_dbl_kind ! initial ice concentration
         endif
      enddo

      !-----------------------------------------------------------------
      ! Define cells where ice is placed (or other values are used)
      ! Edges using initial values (zero, above) are commented out
      !-----------------------------------------------------------------

      icells = 0
      if (iblock == 1) then              ! west edge
            do j = 1, ny_block
            do i = 1, ilo
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (iblock == nblocks_x) then      ! east edge
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (iglob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + jglob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do j = 1, ny_block
            do i = ihi, ibc
               if (tmask(i,j)) then
               icells = icells + 1
               indxi(icells) = i
               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (jblock == 1) then              ! south edge
            do j = 1, jlo
            do i = 1, nx_block
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      if (jblock == nblocks_y) then      ! north edge
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (jglob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + iglob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do j = jhi, ibc
            do i = 1, nx_block
               if (tmask(i,j)) then
!               icells = icells + 1
!               indxi(icells) = i
!               indxj(icells) = j
               endif
            enddo
            enddo
      endif

      !-----------------------------------------------------------------
      ! Set restoring variables
      !-----------------------------------------------------------------

         do n = 1, ncat

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)

               ! ice volume, snow volume
               aicen(i,j,n) = ainit(n)
               vicen(i,j,n) = hinit(n) * ainit(n) ! m
               vsnon(i,j,n) = min(aicen(i,j,n)*hsno_init,p2*vicen(i,j,n))

               call icepack_init_trcr(Tair=Tair(i,j),    Tf=Tf(i,j),  &
                                      Sprofile=salinz(i,j,:),         &
                                      Tprofile=Tmltz(i,j,:),          &
                                      Tsfc=Tsfc,                      &
                                      qin=qin(:),        qsn=qsn(:))

               ! surface temperature
               trcrn(i,j,nt_Tsfc,n) = Tsfc ! deg C
               ! ice enthalpy, salinity
               do k = 1, nilyr
                  trcrn(i,j,nt_qice+k-1,n) = qin(k)
                  trcrn(i,j,nt_sice+k-1,n) = salinz(i,j,k)
               enddo
               ! snow enthalpy
               do k = 1, nslyr
                  trcrn(i,j,nt_qsno+k-1,n) = qsn(k)
               enddo               ! nslyr

            enddo               ! ij
         enddo                  ! ncat

         call icepack_warnings_flush(nu_diag)
         if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
            file=__FILE__, line=__LINE__)

   end subroutine set_restore_var

!=======================================================================

!  This subroutine is intended for restoring the ice state to desired
!  values in cells surrounding the grid.
!  Note: This routine will need to be modified for nghost > 1.
!        We assume padding occurs only on east and north edges.

 subroutine ice_HaloRestore

      use ice_blocks, only: block, get_block, nblocks_x, nblocks_y
      use ice_calendar, only: dt
      use ice_domain, only: ew_boundary_type, ns_boundary_type, &
          nblocks, blocks_ice
          
      use ice_domain_size, only: nblyr,nilyr, nslyr, ncat
      use ice_communicate, only: my_task, master_task
      use ice_fileunits, only: nu_diag
      
       use icepack_parameters, only: c0, c1, c2, p2, p5, rhoi, rhos, Lfresh, &
           cp_ice, cp_ocn, Tsmelt, Tffresh,hs_min, p01
      
      use icepack_tracers, only: nt_Tsfc, nt_qice, nt_qsno, nt_sice, &  
          nt_fbri, tr_brine, nt_vlvl, nt_alvl, nt_iage, &
          nt_apnd, nt_hpnd, nt_ipnd, tr_aero, tr_pond_topo, nbtrcr
          
      
          
       

      use icepack_parameters, only: ktherm
      use ice_flux, only: Tmltz,Tf
     
      use icepack_mushy_physics, only: icepack_enthalpy_snow, icepack_enthalpy_mush
      use ice_dyn_shared, only: a_min
      
      
      use ice_flux, only: fpond, fresh, fhocn, fsalt,Tf
      use ice_flux_bgc, only: flux_bio, faero_ocn,fiso_ocn
      
      use icepack_itd, only: cleanup_itd
      use ice_state, only: trcr_base, nt_strata, n_trcr_strata 
          
          
!       use ice_therm_shared, only: heat_capacity
      use ice_arrays_column, only: hin_max, first_ice
      
      use ice_exit, only: abort_ice
     
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     i,j,k,iblk,nt,n,      &! dummy loop indices
     ilo,ihi,jlo,jhi,    &! beginning and end of physical domain
     ibc,                &! ghost cell column or row
     ntrcr,              &!
     npad                 ! padding column/row counter

   type (block) :: &
     this_block  ! block info for current block

   real (dbl_kind) :: &
     secday,             &!
     ctime,Ti,slope,cslope,trest_ice               ! dt/trest
   
   logical(kind=log_kind) :: &
         lsnow, &          ! snow presence: T: has snow, F: no snow
         lice              ! ice presence: T: has ice, F: no ice
        
   
   character(len=*), parameter :: subname = '(ice_HaloRestore)'
   
   
   logical (kind=log_kind) :: &
         l_stop          ! if true, abort model

   integer (kind=int_kind) :: &
      istop, jstop    ! indices of grid cell where model aborts

   l_stop = .false.

   call ice_timer_start(timer_bound)
   call icepack_query_parameters(secday_out=secday)
   call icepack_query_tracer_sizes(ntrcr_out=ntrcr)
   call icepack_warnings_flush(nu_diag)
   if (icepack_warnings_aborted()) call abort_ice(error_message=subname, &
      file=__FILE__, line=__LINE__)

!-----------------------------------------------------------------------
!
!  Initialize
!
!-----------------------------------------------------------------------

      ! for now, use same restoring constant as for SST
!       if (trestore == 0) then
!          trest = dt          ! use data instantaneously
!       else
          trest_ice = real(.5,kind=dbl_kind) * 3600!secday ! seconds
!       endif
!       trest = real(5,kind=dbl_kind) * secday
	
        ctime = dt/trest_ice
        
        
!         print *, 'trest_ice,ctime: ',trest_ice,ctime
!       if (my_task == master_task ) then
!            write (nu_diag,*) 'Restoring ice'
!       end if
!-----------------------------------------------------------------------
!
!  Restore values in cells surrounding the grid
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblk,ilo,ihi,jlo,jhi,this_block, &
   !$OMP                     i,j,n,nt,ibc,npad)
   do iblk = 1, nblocks
      this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

      if (this_block%iblock == 1) then              ! west edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, ny_block
            do i = 1, ilo
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%iblock == nblocks_x) then  ! east edge
         if (trim(ew_boundary_type) /= 'cyclic') then
            ! locate ghost cell column (avoid padding)
            ibc = nx_block
            do i = nx_block, 1, -1
               npad = 0
               if (this_block%i_glob(i) == 0) then
                  do j = 1, ny_block
                     npad = npad + this_block%j_glob(j)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do n = 1, ncat
            do j = 1, ny_block
            do i = ihi, ibc
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == 1) then              ! south edge
         if (trim(ns_boundary_type) /= 'cyclic') then
            do n = 1, ncat
            do j = 1, jlo
            do i = 1, nx_block
               aicen(i,j,n,iblk) = aicen(i,j,n,iblk) &
                  + (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime
               vicen(i,j,n,iblk) = vicen(i,j,n,iblk) &
                  + (vicen_rest(i,j,n,iblk)-vicen(i,j,n,iblk))*ctime
               vsnon(i,j,n,iblk) = vsnon(i,j,n,iblk) &
                  + (vsnon_rest(i,j,n,iblk)-vsnon(i,j,n,iblk))*ctime
               do nt = 1, ntrcr
                  trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk) &
                     + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,j,nt,n,iblk))*ctime
               enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (this_block%jblock == nblocks_y) then  ! north edge
         if (trim(ns_boundary_type) /= 'cyclic' .and. &
             trim(ns_boundary_type) /= 'tripole' .and. &
             trim(ns_boundary_type) /= 'tripoleT') then
            ! locate ghost cell row (avoid padding)
            ibc = ny_block
            do j = ny_block, 1, -1
               npad = 0
               if (this_block%j_glob(j) == 0) then
                  do i = 1, nx_block
                     npad = npad + this_block%i_glob(i)
                  enddo
               endif
               if (npad /= 0) ibc = ibc - 1
            enddo

            do j = jhi,ibc  
             cslope  = real(1) !- real(ibc-j)/real(3)
            
	    !write(nu_diag,*) "#####################################################"
	    !write(nu_diag,*) ,'jhi,j,jhi-j,cslope = ',jhi,j,jhi+1-j,cslope	    
	    
            do i = 1, nx_block
            do n = 1, ncat
               if (sea_ice_time_bry) then 
		
! 		  write (nu_diag,*) 'aice_bry  : ', aicen_bry(i,ibc,n,iblk)
! 		  write (nu_diag,*) 'vice_bry  : ', vicen_bry(i,ibc,n,iblk)
! 		  write (nu_diag,*) 'hice_bry  : ', (vicen_bry(i,j,n,iblk)/aicen_bry(i,j,n,iblk))
		  
	    
		    aicen_rest(i,j,n,iblk) 	    = aicen_bry(i,ibc,n,iblk)
		    vicen_rest(i,j,n,iblk) 	    = vicen_bry(i,ibc,n,iblk)
		    vsnon_rest(i,j,n,iblk)          = vsnon_bry(i,ibc,n,iblk)
		    trcrn_rest(i,j,nt_Tsfc,n,iblk)  = Tsfc_bry(i,ibc,n,iblk)   
		    trcrn_rest(i,j,nt_alvl,n,iblk)  = alvln_bry(i,ibc,n,iblk)
		    trcrn_rest(i,j,nt_vlvl,n,iblk)  = vlvln_bry(i,ibc,n,iblk) 
		    trcrn_rest(i,j,nt_iage,n,iblk)  = iage_bry(i,ibc,n,iblk) 
		    if (tr_pond_lvl) then
		      trcrn_rest(i,j,nt_apnd,n,iblk) = apondn_bry(i,ibc,n,iblk) 
		      trcrn_rest(i,j,nt_hpnd,n,iblk) = hpondn_bry(i,ibc,n,iblk)
		      trcrn_rest(i,j,nt_ipnd,n,iblk) = ipondn_bry(i,ibc,n,iblk)  
		    endif
                  endif
                  
!                   if (aicen_bry(i,ibc,n,iblk) > c0 ) then
! 			write(nu_diag,*) "#####################################################"
! 			write(nu_diag,*) "aicen_rest(i,j,n,iblk) = ",aicen_rest(i,j,n,iblk)
! 			write(nu_diag,*) "aicen_bry(i,j,n,iblk) = ",aicen_bry(i,j,n,iblk)
! 			write(nu_diag,*) "vicen_rest(i,j,n,iblk) = ",vicen_rest(i,j,n,iblk)
! 			write(nu_diag,*) "aicen_bry(i,j,n,iblk) = ",vicen_bry(i,j,n,iblk)
! 			write(nu_diag,*) "vsnon_rest(i,j,n,iblk) = ",vsnon_rest(i,j,n,iblk)
! 		        write(nu_diag,*) "vicen_rest(i,j,n,iblk)/aicen_rest(i,j,n,iblk) = ",vicen_rest(i,j,n,iblk)/aicen_rest(i,j,n,iblk)
 
!                   endif
                  
                  do k = 1,nilyr
		     trcrn_rest(i,j,nt_sice+k-1,n,iblk) = 19.539*((real(k)/real(nilyr))**2) - 19.93*(real(k)/real(nilyr)) + 8.913
!                    trcrn_rest(i,j,nt_sice+k-1,n,iblk) = Sinz_bry(i,ibc,k,n,iblk)                 
                     if (ktherm == 2) then
                        trcrn_rest(i,j,nt_qice+k-1,n,iblk) = icepack_enthalpy_mush(Tinz_bry(i,ibc,k,n,iblk),Sinz_bry(i,ibc,k,n,iblk))
                     else             
			trcrn_rest(i,j,nt_qice+k-1,n,iblk) = icepack_enthalpy_mush(c0-p01,c0)
		     endif
                  enddo  
		  do k = 1,nslyr  
                     trcrn_rest(i,j,nt_qsno+k-1,n,iblk) = -rhos*(Lfresh - cp_ice * min(c0,trcrn_rest(i,j,nt_Tsfc,n,iblk))) 
                  enddo
               
                  
               !if ((vicen_bry(i,ibc,n,iblk)/aicen_bry(i,ibc,n,iblk)) > real(.1)) then
               
!                if (SUM(aicen(i,ibc,n,iblk),3) > 1) then
!                   write(nu_diag,*) ' aicen exceeds 1: aicen= ', aicen(i,ibc,n,iblk),3)
!                endif 

! 	       aicen(i,j,n,iblk) = aicen_rest(i,j,n,iblk)  
!                vicen(i,j,n,iblk) = vicen_rest(i,j,n,iblk) 
!                vsnon(i,j,n,iblk) = vsnon_rest(i,j,n,iblk)  
! 	      if ((c0 < aicen(i,j,n,iblk)) .or. (c0 < aicen_rest(i,j,n,iblk))) then
! 	      write(nu_diag,*) "############################################################################"
! 	      write(nu_diag,*) "i,j,n,iblk = ", i, j, n, iblk 
! 	      write(nu_diag,*) "Reseting ice values, aicen(i,jhi,n,iblk) , aicen_rest(i,j,n,iblk) = " , aicen(i,jhi,n,iblk), aicen_rest(i,j,n,iblk)
!               write(nu_diag,*) "cslope = ", cslope
!               endif
                !+ (aicen_rest(i,j,n,iblk)-aicen(i,j,n,iblk))*ctime*cslope
! 	      if ((c0 < aicen(i,j,n,iblk))) then
! 		write(nu_diag,*) "Value post reset, aicen(i,j,n,iblk) = " , aicen(i,j,n,iblk)
!               endif
               
               
! 		
               if    (vvel(i,jhi,iblk) >=0) then 
               
		aicen(i,j,n,iblk)  = aicen(i,jhi,n,iblk)  + (c0-aicen(i,jhi,n,iblk))*ctime
                vicen(i,j,n,iblk)  = vicen(i,jhi,n,iblk)  + (c0-vicen(i,jhi,n,iblk))*ctime
                vsnon(i,j,n,iblk)  = vsnon(i,jhi,n,iblk)  + (c0-vsnon(i,jhi,n,iblk))*ctime
                
		uvel(i,j,iblk)     =  uvel(i,j-1,iblk)                
	        vvel(i,j,iblk)     =  vvel(i,j-1,iblk) 
	        
		divu(i,j,iblk)     = c0!real(0.1)
		shear(i,j,iblk)    = c0!real(0.1)
		strength(i,j,iblk) = c0!real(0.1)
              else
              
		aicen(i,j,n,iblk)  = aicen(i,jhi,n,iblk)  + (aicen_rest(i,j,n,iblk)-aicen(i,jhi,n,iblk))*ctime
                vicen(i,j,n,iblk)  = vicen(i,jhi,n,iblk)  + (vicen_rest(i,j,n,iblk)-vicen(i,jhi,n,iblk))*ctime
                vsnon(i,j,n,iblk)  = vsnon(i,jhi,n,iblk)  + (vsnon_rest(i,j,n,iblk)-vsnon(i,jhi,n,iblk))*ctime
                
		uvel(i,j,iblk)     = uvel_rest(i,jhi,iblk)                
	        vvel(i,j,iblk)     = vvel_rest(i,jhi,iblk) 
	        
! 	        divu(i,j,iblk)     = c0
! 		shear(i,j,iblk)    = c0
! 		strength(i,j,iblk) = c0
		
		divu(i,j,iblk)     = divu(i,j-1,iblk)
		shear(i,j,iblk)    = shear(i,j-1,iblk)
		strength(i,j,iblk) = strength(i,j-1,iblk)
              endif
              
             ! if ((vicen(i,j,n,iblk)/aicen(i,j,n,iblk) < real(1e-2)) .or. (aicen(i,j,n,iblk) < real(1e-2))) then
		
            !   aicen(i,j,n,iblk)  = c0
	!	vicen(i,j,n,iblk)  = c0  
        !       vsnon(i,j,n,iblk)  = c0
        !      endif
              
              
              do nt = 1, ntrcr-nbtrcr
                  if  ((sea_ice_time_bry).and.((nt == nt_qice).or. &
                      (nt == nt_sice))) then
                     do k = 1,nilyr
                    trcrn(i,j,nt+k-1,n,iblk) = trcrn_rest(i,j,nt+k-1,n,iblk)!+(trcrn_rest(i,j,nt+k-1,n,iblk)- trcrn(i,j,nt+k-1,n,iblk))*ctime*cslope
!                     trcrn(i,j,nt+k-1,n,iblk) = trcrn(i,jhi,nt+k-1,n,iblk) +(trcrn_rest(i,j,nt+k-1,n,iblk)- trcrn(i,jhi,nt+k-1,n,iblk))*ctime
                     enddo
                  else if ((sea_ice_time_bry).and.(nt == nt_qsno)) then
                     do k = 1,nslyr
		      trcrn(i,j,nt+k-1,n,iblk) = trcrn_rest(i,j,nt+k-1,n,iblk) !+(trcrn_rest(i,j,nt+k-1,n,iblk)- trcrn(i,j,nt+k-1,n,iblk))*ctime*cslope
!                       trcrn(i,j,nt+k-1,n,iblk) = trcrn(i,j,nt+k-1,n,iblk) +(trcrn_rest(i,j,nt+k-1,n,iblk)- trcrn(i,j,nt+k-1,n,iblk))*cslope!                            &
!                      trcrn(i,j,nt+k-1,n,iblk) = trcrn(i,jhi,nt+k-1,n,iblk)+ (trcrn_rest(i,j,nt+k-1,n,iblk)-trcrn(i,jhi,nt+k-1,n,iblk))*ctime
                     enddo 
                  else 
                     trcrn(i,j,nt,n,iblk) = trcrn_rest(i,j,nt,n,iblk)!+(trcrn_rest(i,j,nt,n,iblk)- trcrn(i,j,nt,n,iblk))*ctime*cslope
!                       trcrn(i,j,nt,n,iblk) = trcrn(i,j,nt,n,iblk)&
! 					       +(trcrn_rest(i,j,nt,n,iblk)- trcrn(i,j,nt,n,iblk))*cslope
!                       trcrn(i,j,nt,n,iblk) = trcrn(i,jhi,nt,n,iblk) + (trcrn_rest(i,j,nt,n,iblk)-trcrn(i,jhi,nt,n,iblk))*ctime
                  endif 
               enddo
           
            if (aicen(i,j,n,iblk) >real(0)) then
                if ((vicen(i,j,n,iblk)/aicen(i,j,n,iblk) < real(1e-2)) .or. (aicen(i,j,n,iblk) < real(1e-2))) then
                
                        aicen(i,j,n,iblk)  = c0
                        vicen(i,j,n,iblk)  = c0
                        vsnon(i,j,n,iblk)  = c0
                        trcrn(i,j,:,n,iblk) = c0
                endif
            endif
            enddo !n
            call cleanup_itd (dt,         ntrcr,            &
                        nilyr,                nslyr,            &
                        ncat,                 hin_max,          &
                        aicen(i,j,:,iblk),    trcrn(i,j,1:ntrcr,:,iblk),                        &
                        vicen(i,j,:,iblk),    vsnon(i,j,:,iblk),                        &
                        aice0(i,j,iblk),      aice(i,j,iblk),                       &
                        n_aero,                                 &
                        nbtrcr,               nblyr,            &
                        tr_aero,                                &
                        tr_pond_topo,                           &
                        first_ice(i,j,:,iblk),                              &
                        trcr_depend(1:ntrcr), trcr_base,                 &
                        n_trcr_strata,        nt_strata,        &
                        fpond(i,j,iblk),      fresh(i,j,iblk),                      &
                        fsalt(i,j,iblk),      fhocn(i,j,iblk),                      &
                        faero_ocn(i,j,:,iblk),fiso_ocn(i,j,:,iblk),                     &
                        flux_bio(i,j,1:nbtrcr,iblk),Tf(i,j,iblk))
            enddo !i
            enddo !j
            
!             call cleanup_itd (dt,                  ntrcr,            &
!                         nilyr,                nslyr,            &
!                         ncat,                 hin_max,          &
!                         aicen(i,j,:,iblk),   	  &
!                         trcrn(i,j,1:ntrcr,:,iblk),&
!                         vicen   (i,j,:,iblk),	  &
!                         vsnon (i,j,  :,iblk),	  &
!                         aice0   (i,j,  iblk),	  & 
!                         aice      (i,j,iblk),     &
!                         n_aero,           	  &
!                         nbtrcr,           	  &   
!                         nblyr,            	  &
!                         tr_aero,                  &
!                         tr_pond_topo,             &
!                         first_ice(i,j,:,iblk),    &
!                         trcr_depend(1:ntrcr),     &
!                         trcr_base,        	  &
!                         n_trcr_strata,            &
!                         nt_strata,                &
!                         fpond(i,j,iblk),          &
!                         fresh(i,j,  iblk),	  &
!                         fsalt(i,j,iblk),   	  &
!                         fhocn(i,j,  iblk),	  &
!                         faero_ocn(i,j,:,iblk),	  &
!                         fiso_ocn(i,j,:,iblk), 	  &
!                         fzsal(i,j,  iblk),    	  &
!                         flux_bio(i,j,1:nbtrcr,iblk))
!             



            

                               
                               
!             if (l_stop) then
!                 write (nu_diag,*) ' my_task, iblk =', &
!                                    my_task, iblk
!                 write (nu_diag,*) 'Global block:', this_block%block_id
!                 if (istop > 0 .and. jstop > 0) &
!                      write(nu_diag,*) 'Global i and j:', &
!                                       this_block%i_glob(istop), &
!                                       this_block%j_glob(jstop) 
!                 call abort_ice ('ice: ITD cleanup error in ice_HaloRestore')
!             endif
            
         endif
!       endif
      
      

   enddo ! iblk
   !$OMP END PARALLEL DO
!     write (nu_diag,*) 'Restoring ice north edge...RESTORED!'
   
   call bound_state (aicen, vicen, vsnon, ntrcr, trcrn)
   
   
   call ice_timer_stop(timer_bound)

 end subroutine ice_HaloRestore

!=======================================================================

      end module ice_restoring

!=======================================================================
