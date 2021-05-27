      Program main
!=======================================================================
!     CAM3 grid and resolutions
!=======================================================================
      use cam_grid
!=======================================================================
!     CAM3 variables and parameters for inital vortex
!=======================================================================
      use cam_variables

!=======================================================================
!     test case
!=======================================================================
      use dcmip_initial_conditions_test_123
      use dcmip_initial_conditions_test_4
      use baroclinic_wave
!      use tc_initial_vortex
 
      implicit none
      
!=======================================================================
!     local variables
!=======================================================================
      integer i, j, k, n                                 ! loop indices                   

      integer dry                                        ! 0: dry, 1: moist
      integer deep                                       ! 0: shallow, 1: deep
      integer zcoords                                    ! 0: p, 1: z-levels
      integer pertt                                      ! 0: exponential, 1: streamfunction
      logical lperturb                                   ! .FALSE. : steady-state, .TRUE. : baroclinic wave
      real(8) d1, d2, d3, d4                             ! dummy vars
      real(8) zdummy, p100, bigX                         ! dummy vars
      integer t_profile                                  ! flag for temperature profile, 1: constant lapse rate, 0: isothermal

!-----------------------------------------------------------------------
!     choose a test case
!-----------------------------------------------------------------------
 10   write (*,*)
      write (*,'(A)') '(11) ad: 3D Deformation Advection'
!      write (*,'(A)') '(12) ah: 3D Hadley Cell Advection'
!      write (*,'(A)') '(13) ah: 3D Advection in Presence of Orography'
!      write (*,'(A)') '(21) cm: Small Earth Schar Mountain Constant Flow'
!      write (*,'(A)') '(22) sm: Small Earth Schar Mountain Shear Flow'
      write (*,'(A)') '(200) sr: Dry state at rest with mountain, constant lapse rate'
      write (*,'(A)') '(31)  gw: Dry gravity waves'
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
      write (*,'(A)') '(80)  tc: Reed and Jablonowski (2011) tropical cyclone'
      read (*,*) choice_testcase

!-----------------------------------------------------------------------
!     choose model version
!-----------------------------------------------------------------------
  5   write (*,*) 'Choose a CESM dynamical core:'
!      write (*,*) '(1) Eulerian / Semi-Lagrangian'
!      write (*,*) '(2) Finite volume model'
!      write (*,*) '(3) Spectral Element'
!      write (*,*)
!      write (*,*) 'Choose model version:'
!      read  (*,*) model_version
      model_version = 4
      select case (model_version)
      case (1)
         model = 'eul_sl'
      case (2)
         model = 'fv'
      case (3) 
         model = 'se'
      case (4) 
         write (*,*) 'Finite Volume Cubed Sphere (FV3)'
         model = 'fv3'
      case default
         goto 5
      end select

!-----------------------------------------------------------------------
!     choose horizontal resolution and set resolution dependent parameters
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
!     call print_grid

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

      case (200)
          tc_name = 'sr'
          tc_nr   = '200'
          t_profile = 0     ! 0: constant lapse rate, 1: isothermal
          print*,'Creating Initial Data for Test 200'
         
      case (31)
          tc_name = 'gw'
          tc_nr   = '31'
          print*,'Creating Initial Data for Test 31'

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

      case (60)
          tc_name = 'md'
          tc_nr   = '60'
          dry = 0
          print*,'Creating Initial Data for Test 60'

      case (61)
          tc_name = 'mm'
          tc_nr   = '61'
          dry = 1
          print*,'Creating Initial Data for Test 61'

      case (80)
          tc_name = 'tc'
          tc_nr   = '80'
          dry = 1
          print*,'Creating Initial Data for Test 80'

      end select


      select case (choice_testcase)
      case (400, 410, 411, 412, 413, 42, 43)
         !==========================================
         !         SE unstaggered grid
         !==========================================
         do k = 1, nlev
            do i = 1, ncol
                call test4_baroclinic_wave (lperturb,dry,bigx,hyam(k),hybm(k),lon(i),lat(i),zdummy,zdummy,0,u(k,i,1), &
                                            v(k,i,1),d3,t(k,i,1),phis(i,1),ps(i,1),d4,q(k,i,1),q1(k,i,1),q2(k,i,1))
            enddo
         enddo

      case (80)
         do k = 1, nlev
            do i = 1, ncol
!               call tc_initial_vortex(hyam(k),hybm(k),lon(i),lat(i),dumin,0.d0,0, &
!                    u(k,i,1),v(k,i,1),t(k,i,1), &
!                    phis(i,1),ps(i,1),rho(k,i,1),q(k,i,1))
            enddo
         enddo

      case (70,71)
         do k = 1, nlev
            p100 = 100.0*lev(k)
            do i = 1, ncol
                   call baroclinic_wave_test(deep, dry, pertt, bigx, lon(i),lat(i),p100,zdummy,zcoords,&
                        u(k,i,1),v(k,i,1),t(k,i,1),d3,phis(i,1),ps(i,1),d4,q(k,i,1))
            enddo
         enddo

       case (60, 61)
         do k = 1, nlev
            do i = 1, ncol
                call mountain_Rossby (dry,hyam(k),hybm(k),lon(i),lat(i),u(k,i,1), &
                                            v(k,i,1),d3,t(k,i,1),phis(i,1),ps(i,1),d4,q(k,i,1),q1(k,i,1),q2(k,i,1))
            enddo
         enddo

       case (11)
         do k = 1, nlev
            p100 = 100.0*lev(k)
            do i = 1, ncol
                call test1_advection_deformation (lon(i),lat(i),p100,d1,0,u(k,i,1),v(k,i,1),d3,t(k,i,1),phis(i,1), &
                                                  ps(i,1),d4,q(k,i,1),q1(k,i,1),q2(k,i,1),q3(k,i,1),q4(k,i,1))
            enddo
         enddo

       case (200)
         do k = 1, nlev
            do i = 1, ncol
                call test2_rest_mountain (hyam(k),hybm(k),lon(i),lat(i),t_profile,u(k,i,1), &
                                            v(k,i,1),d3,t(k,i,1),phis(i,1),ps(i,1),d4,q(k,i,1))
            enddo
         enddo

       case (31)
         do k = 1, nlev
            do i = 1, ncol
                   call test3_gravity_wave_new (lon(i),lat(i),hyam(k),hybm(k),u(k,i,1),&
                        v(k,i,1),d3,t(k,i,1),phis(i,1),ps(i,1),d4,q(k,i,1))
            enddo
         enddo

       end select


!=======================================================================
!     write the initial data to the CAM3 initialization file
!=======================================================================
      call write_to_file
!-----------------------------------------------------------------------
!     deallocate arrays
!-----------------------------------------------------------------------
      call deallocate_grid
      call deallocate_variables

  end

