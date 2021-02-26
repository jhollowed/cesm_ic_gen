    Program main

!=======================================================================
!     CAM3 grid and resolutions
!=======================================================================
      use cam_grid, only: model_version, tc_name, tc_nr, c_alpha, &
                          choice_testcase,                        &
                          hyam, hybm
!=======================================================================
!     CAM3 variables
!=======================================================================
      use cam_variables
!=======================================================================
!     test cases
!=======================================================================
      use dcmip_initial_conditions_test_123
      use dcmip_initial_conditions_test_4
      use baroclinic_wave

      implicit none

!=======================================================================
!     local variables
!=======================================================================
      integer i, j, k, ai, bi, imh                          ! loop indeces
      integer dry                                           ! 0: dry, 1: moist
      integer deep                                          ! 0: shallow, 1: deep
      integer zcoords                                       ! 0: p, 1: z-levels
      integer pertt                                         ! 0: exponential, 1: streamfunction
      logical lperturb                                      ! .FALSE. : steady-state, .TRUE. : baroclinic wave
      real(8)    d1, d2, d3, d4, d5, d6, d7, d8, d9         ! dummy vars
      real(8)    d10, d11, d12, d13                         ! dummy vars
      real(8)    zdummy, p100, bigX                         ! dummy vars
      integer    t_profile                                  ! 0: constant lapse rate, 1: isothermal
      
!-----------------------------------------------------------------------
!     choose a test case
!-----------------------------------------------------------------------
 10   write (*,*)
      write (*,'(A)') 'Choose the test case:'
      write (*,'(A)') '(11) ad: 3D Deformation Advection'
!      write (*,'(A)') '(12) ah: 3D Hadley Cell Advection'
!      write (*,'(A)') '(13) ah: 3D Advection in Presence of Orography'
!      write (*,'(A)') '(21) cm: Small Earth Schar Mountain Constant Flow'
!      write (*,'(A)') '(22) sm: Small Earth Schar Mountain Shear Flow'
      write (*,'(A)') '(200) sr: Dry state at rest with mountain, constant lapse rate'
      write (*,'(A)') '(31)  gw: Dry gravity waves'
      write (*,'(A)') '(32)  hs: Balanced initial state driven by heat source'
      write (*,'(A)') '(400) ss: Regular-Earth dry steady-state initial conditions'
      write (*,'(A)') '(410) bw: Regular-Earth dry baroclinic wave'
!      write (*,'(A)') '(411) bw: Small-Earth baroclinic wave X=10'
!      write (*,'(A)') '(412) bw: Small-Earth baroclinic wave X=100'
!      write (*,'(A)') '(413) bw: Small-Earth baroclinic wave X=1000'
      write (*,'(A)') '(42)  mw: Moist baroclinic wave X=1, large-scale condensation, same initial data as (43)'
      write (*,'(A)') '(43)  mw: Moist baroclinic wave X=1, large-scale condensation, sensible and latent heat, turb.'
      write (*,'(A)') '(60)  md: Dry mountain-generated Rossby waves'
      write (*,'(A)') '(61)  mm: Moist mountain-generated Rossby waves'
      write (*,'(A)') '(70)  ud: Ullrich et al. dry baroclinic wave X=1 with exponential u perturb'
      write (*,'(A)') '(71)  um: Ullrich et al. moist baroclinic wave X=1 with exponential u perturb'

      read (*,*) choice_testcase
      
!-----------------------------------------------------------------------
!     choose model version
!-----------------------------------------------------------------------
 20   write (*,*) 'Choose a CESM dynamical core:'
      write (*,*) '(1) Eulerian spectral transform model'
      write (*,*) '(2) Finite volume model'
      write (*,*)
      write (*,*) 'Choose model version:'
      read  (*,*) model_version

      select case (model_version)
      case (1)
         model = 'eul'
      case (2)
         model = 'fv'
      case default
         goto 20
      end select
      
!-----------------------------------------------------------------------
!     choose the horizontal resolution and set resolution dependent parameters
!-----------------------------------------------------------------------
      call select_resolution
!-----------------------------------------------------------------------
!     allocate data arrays
!-----------------------------------------------------------------------
      call allocate_data_arrays
!=======================================================================
!     allocate and initialize the grid
!=======================================================================
      call initialize_grid
!=======================================================================
!     initialize the model variables
!=======================================================================
      zdummy = 0.0
      dry = 0           ! default is set to dry conditions for test 4XX
      lperturb = .TRUE. ! default is set to perturbed initial state for test 4XX

      select case (choice_testcase)
      case (11)
          tc_name = 'ad'
          tc_nr   = '11'
          print*,'Creating Initial Data for Test 11'
          if (model_version == 1) then
            do i=1, nlon
              do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_deformation(lon(i),lat(j),p100,zdummy,0,u(i,k,j),v(i,k,j),d4,t(i,k,j), & 
                        phis(i,j),ps(i,j),d5,q(i,k,j),tt_lw(i,k,j),& 
                        tt_md(i,k,j),tt_hi(i,k,j),ttrmd(i,k,j))
                enddo
               enddo
            enddo
          else
            do i=1, nlon
              do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_deformation(lon(i),lat(j),p100,zdummy,0,d2,d3,d4,t(i,k,j), & 
                        phis(i,j),ps(i,j),d5,q(i,k,j),tt_lw(i,k,j),& 
                        tt_md(i,k,j),tt_hi(i,k,j),ttrmd(i,k,j))
                enddo
               enddo
            enddo
            do i=1, nlon
              do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_deformation(slon(i),lat(j),p100,zdummy,0,d3,vs(i,k,j), & 
                        d4,d5,d6,d7,d8,d9,d10,d11,d12,d13)
                enddo
              enddo
            enddo
            do i=1, nlon
              do j=1, nslat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_deformation(lon(i),slat(j),p100,zdummy,0,us(i,k,j), & 
                        d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13)
                enddo
              enddo
            enddo
         endif
      case (12)
          tc_name = 'ah'
          tc_nr   = '012'
          print*,'Creating Initial Data for Test 12'
          if (model_version == 1) then
            do i=1, nlon
              do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_hadley(lon(i),lat(j),p100,zdummy,0,u(i,k,j),v(i,k,j),d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j))
                enddo
               enddo
            enddo
          else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_hadley(lon(i),lat(j),p100,zdummy,0,d1,d2,d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev 
                   p100 = 100.0*lev(k)
                   call test1_advection_hadley(slon(i),lat(j),p100,zdummy,0,d1,vs(i,k,j), & 
                        d3,d4,d5,d6,d7,d8,d9)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_hadley(lon(i),slat(j),p100,zdummy,0,us(i,k,j),d2, & 
                        d3,d4,d5,d6,d7,d8,d9)
                enddo
             enddo
          enddo
          endif

      case (13)
          tc_name = 'ca'
          tc_nr   = '013'
          print*,'Creating Initial Data for Test 13'
          if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_orography(lon(i),lat(j),p100,zdummy,0,u(i,k,j),v(i,k,j),d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j),& 
                        tt_md(i,k,j),tt_hi(i,k,j),ttrmd(i,k,j))
                enddo
             enddo
          enddo
          else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_orography(lon(i),lat(j),p100,zdummy,0,d1,d2,d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j),& 
                        tt_md(i,k,j),tt_hi(i,k,j),ttrmd(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev 
                   p100 = 100.0*lev(k)
                   call test1_advection_orography(slon(i),lat(j),p100,zdummy,0,d1,vs(i,k,j), & 
                        d3,d4,d5,d6,d7,d8,d9,d10,d11,d12)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test1_advection_orography(lon(i),slat(j),p100,zdummy,0,us(i,k,j),d2, & 
                        d3,d4,d5,d6,d7,d8,d9,d10,d11,d12)
                enddo
             enddo
          enddo
          endif

      case (200)
          tc_name = 'sr'
          tc_nr   = '200'
          t_profile = 0
          print*,'Creating Initial Data for Test 200'
          if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   call test2_rest_mountain(hyam(k),hybm(k),lon(i),lat(j),t_profile,u(i,k,j),v(i,k,j),d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   call test2_rest_mountain(hyam(k),hybm(k),lon(i),lat(j),t_profile,d1,d2,d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   call test2_rest_mountain(hyam(k),hybm(k),slon(i),lat(j),t_profile,d1, & 
                        vs(i,k,j),d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   call test2_rest_mountain(hyam(k),hybm(k),lon(i),slat(j),t_profile,us(i,k,j), & 
                        d1,d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          endif

      case (21)
          tc_name = 'cm'
          tc_nr   = '021'
          print*,'Creating Initial Data for Test 21'
          if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(lon(i),lat(j),p100,zdummy,0,1,hyam(k),hybm(k),0,u(i,k,j),v(i,k,j),d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(lon(i),lat(j),p100,zdummy,0,1,hyam(k),hybm(k),0,d1,d2,d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(slon(i),lat(j),p100,zdummy,0,1,hyam(k),hybm(k),0,d1, & 
                        vs(i,k,j),d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(lon(i),slat(j),p100,zdummy,0,1,hyam(k),hybm(k),0,us(i,k,j), & 
                        d1,d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          endif

      case (22)
          tc_name = 'sm'
          tc_nr   = '022'
          print*,'Creating Initial Data for Test 22'
          if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(lon(i),lat(j),p100,zdummy,0,1,hyam(k),hybm(k),1,u(i,k,j),v(i,k,j),d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(lon(i),lat(j),p100,zdummy,0,1,hyam(k),hybm(k),1,d1,d2,d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(slon(i),lat(j),p100,zdummy,0,1,hyam(k),hybm(k),1,d1, & 
                        vs(i,k,j),d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test2_schaer_mountain(lon(i),slat(j),p100,zdummy,0,1,hyam(k),hybm(k),1,us(i,k,j), & 
                        d1,d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          endif

      case (31,32)
          if (choice_testcase .eq.31) then
            tc_name = 'gw'
            tc_nr   = '31'
            lperturb = .true.
            print*,'Creating Initial Data for Test 31'
          else
            tc_name = 'hs'
            tc_nr   = '32'
            lperturb = .false.
            print*,'Creating Initial Data for Test 32 (run with heat source)'
          endif
          if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   call test3_gravity_wave_new(lperturb,lon(i),lat(j),hyam(k),hybm(k),u(i,k,j),v(i,k,j),d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   !p100 = 100.0*lev(k)
                   !call test3_gravity_wave(lon(i),lat(j),p100,zdummy,0,d1,d2,d3,t(i,k,j),& 
                   !     phis(i,j),ps(i,j),d4,q(i,k,j))

                   call test3_gravity_wave_new(lperturb,lon(i),lat(j),hyam(k),hybm(k),d1,d2,d3,t(i,k,j), & 
                        phis(i,j),ps(i,j),d4,q(i,k,j))

                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   !p100 = 100.0*lev(k)
                   !call test3_gravity_wave(slon(i),lat(j),p100,zdummy,0,d1,vs(i,k,j),d2,d3,d4,d5,d6,d7)

                   call test3_gravity_wave_new(lperturb,slon(i),lat(j),hyam(k),hybm(k),d1,vs(i,k,j),d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   !p100 = 100.0*lev(k)
                   !call test3_gravity_wave(lon(i),slat(j),p100,zdummy,0,us(i,k,j),d1,d2,d3,d4,d5,d6,d7)

                   call test3_gravity_wave_new(lperturb,lon(i),slat(j),hyam(k),hybm(k),us(i,k,j),d1,d2,d3,d4,d5,d6,d7)

                enddo
             enddo
          enddo
          endif

      case (400)
          tc_name = 'ss'
          tc_nr   = '400'
          bigx = 1.0
          lperturb = .FALSE.  ! no perturbation is added to the steady-state
          print*,'Creating Initial Data for Test 400'

      case (410)
          tc_name = 'bw'
          tc_nr   = '410'
          bigx = 1.0
          print*,'Creating Initial Data for Test 410'

      case (411)
          tc_name = 'bw'
          tc_nr   = '411'
          bigx = 10.0
          print*,'Creating Initial Data for Test 411'

      case (412)
          tc_name = 'bw'
          tc_nr   = '412'
          bigx = 100.0
          print*,'Creating Initial Data for Test 412'

      case (413)
          tc_name = 'bw'
          tc_nr   = '413'
          bigx = 1000.0
          print*,'Creating Initial Data for Test 413'

      case (42)
          tc_name = 'mw'
          tc_nr   = '042'
          bigx = 1.0
          dry = 1
          print*,'Creating Initial Data for Test 42'

      case (43)
          tc_name = 'mw'
          tc_nr   = '043'
          bigx = 1.0
          dry = 1
          print*,'Creating Initial Data for Test 43'

      case (70)
          tc_name = 'ud'
          tc_nr   = '070'
          dry     = 0    ! 0: dry, 1:moist
          zcoords = 0    ! pressure is specified as the vertical coordinate
          bigx    = 1.0  ! factor for a reduced size earth
          pertt   = 0    ! 0: exponential, 1: streamfunction 
          deep    = 0    ! 0: shallow, 1: deep 
          print*,'Creating Initial Data for Test 70'

      case (71)
          tc_name = 'um'
          tc_nr   = '071'
          dry     = 1    ! 0: dry, 1:moist
          zcoords = 0    ! pressure is specified as the vertical coordinate
          bigx    = 1.0  ! factor for a reduced size earth
          pertt   = 0    ! 0: exponential, 1: streamfunction 
          deep    = 0    ! 0: shallow, 1: deep 
          print*,'Creating Initial Data for Test 71'

      end select

      if (choice_testcase == 400 .or. choice_testcase == 410 .or. choice_testcase == 411 .or. choice_testcase == 412 .or.  &
          choice_testcase == 413 .or. choice_testcase == 42  .or. choice_testcase == 43) then
        if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test4_baroclinic_wave (lperturb,dry,bigx,hyam(k),hybm(k),lon(i),lat(j),p100,zdummy,0,u(i,k,j), &
                        v(i,k,j),d3,t(i,k,j),phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j),tt_md(i,k,j))
                enddo
             enddo
          enddo
        else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test4_baroclinic_wave (lperturb,dry,bigx,hyam(k),hybm(k),lon(i),lat(j),p100,zdummy,0,d1,d2,d3,t(i,k,j), &
                        phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j),tt_md(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test4_baroclinic_wave (lperturb,dry,bigx,hyam(k),hybm(k),slon(i),lat(j),p100,zdummy,0,d1,vs(i,k,j),d2, &
                                               d3,d4,d5,d6,d7,d8,d9)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call test4_baroclinic_wave (lperturb,dry,bigx,hyam(k),hybm(k),lon(i),slat(j),p100,zdummy,0,us(i,k,j),d1,d2, &
                                               d3,d4,d5,d6,d7,d8,d9)

                enddo
             enddo
          enddo
        endif
      endif

      if (choice_testcase == 70 .or. choice_testcase == 71) then
        if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call baroclinic_wave_test(deep, dry, pertt, bigx, lon(i),lat(j),p100,zdummy,zcoords,&
                        u(i,k,j),v(i,k,j),t(i,k,j),d3,phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
        else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call baroclinic_wave_test(deep, dry, pertt, bigx, lon(i),lat(j),p100,zdummy,zcoords,&
                        d1,d2,t(i,k,j),d3,phis(i,j),ps(i,j),d4,q(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call baroclinic_wave_test(deep, dry, pertt, bigx, slon(i),lat(j),p100,zdummy,zcoords,&
                        d1,vs(i,k,j),d2,d3,d4,d5,d6,d7)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   p100 = 100.0*lev(k)
                   call baroclinic_wave_test(deep, dry, pertt, bigx, lon(i),slat(j),p100,zdummy,zcoords,&
                        us(i,k,j),d1,d2,d3,d4,d5,d6,d7)

                enddo
             enddo
          enddo
        endif
      endif
!---------------------------------------------------------------------------
!     Mountain-generated Rossby waves
!---------------------------------------------------------------------------
      select case (choice_testcase)
      case (60)
          tc_name = 'md'
          tc_nr   = '60'
          dry     = 0
          print*,'Creating Initial Data for Test 60'

      case (61)
          tc_name = 'mm'
          tc_nr   = '61'
          dry     = 1
          print*,'Creating Initial Data for Test 61'
      end select

      if (choice_testcase == 60 .or. choice_testcase == 61) then
        if (model_version == 1) then
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   call mountain_Rossby (dry,hyam(k),hybm(k),lon(i),lat(j),u(i,k,j), &
                        v(i,k,j),d3,t(i,k,j),phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j),tt_md(i,k,j))
                enddo
             enddo
          enddo
        else
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   call mountain_Rossby (dry,hyam(k),hybm(k),lon(i),lat(j),d1,d2,d3,t(i,k,j), &
                        phis(i,j),ps(i,j),d4,q(i,k,j),tt_lw(i,k,j), tt_md(i,k,j))
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nlat
                do k=1, nlev
                   call mountain_Rossby (dry,hyam(k),hybm(k),slon(i),lat(j),d1,vs(i,k,j),d2, &
                                               d3,d4,d5,d6,d7,d8,d9)
                enddo
             enddo
          enddo
          do i=1, nlon
             do j=1, nslat
                do k=1, nlev
                   call mountain_Rossby (dry,hyam(k),hybm(k),lon(i),slat(j),us(i,k,j),d1,d2, &
                                               d3,d4,d5,d6,d7,d8,d9)

                enddo
             enddo
          enddo
        endif
      endif


!=======================================================================
!     write the initial data to the CAM5 initialization file
!=======================================================================
      call write_to_file


!-----------------------------------------------------------------------
!     deallocate arrays
!-----------------------------------------------------------------------
      call deallocate_grid
      call deallocate_variables

  end program main

