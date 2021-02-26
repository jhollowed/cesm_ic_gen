MODULE dcmip_initial_conditions_test_123

  !=======================================================================
  !
  !  Functions for setting up initial conditions for the dynamical core tests:
  !
  !  11 - Deformational Advection Test
  !  12 - Hadley Cell Advection Test
  !  13 - Orography Advection Test
  !  200 - Steady-state at rest with constant lapse rate, mountain
  !  2X - Non-Hydrostatic Mountain Waves Over A Schar-Type Mountain
  !  31 - Non-Hydrostatic Gravity Waves
  !
  !  Given a point specified by: 
  !  	lon	longitude (radians) 
  ! 	lat	latitude (radians) 
  ! 	p/z	pressure/height
  !  the functions will return:
  !	u	zonal wind (m s^-1)
  !	v	meridional wind (m s^-1)
  !	w	vertical velocity (m s^-1)
  !	t	temperature (K)
  !	phis	surface geopotential (m^2 s^-2)
  !	ps	surface pressure (Pa)
  !	rho	density (kj m^-3)
  !	q	specific humidity (kg/kg)
  !	qi	tracers (kg/kg)
  !     p       pressure if height based (Pa)
  !
  !
  !  Authors: James Kent, Paul Ullrich, Christiane Jablonowski 
  !		(University of Michigan, dcmip@ucar.edu)
  !          version 3
  !          June/8/2012
  !
  !          corrected in version 3: 
  !          test 31: the density is now initialized with the background temperature (not the perturbed T)
  !          constants converted to double precision
  !
  !=======================================================================

!=======================================================================
! use physical constants
!=======================================================================
  use cam_physical_constants
!=======================================================================
! CAM3 grid and resolutions
!=======================================================================
  use cam_grid
!=======================================================================
! CAM3 variables
!=======================================================================
  use cam_variables

  IMPLICIT NONE

!-----------------------------------------------------------------------
!     Physical Parameters
!-----------------------------------------------------------------------

	!real(8), parameter ::	a	= 6371220.0d0,	&	! Earth's Radius (m)
	!			Rd 	= 287.0d0,	&	! Ideal gas const dry air (J kg^-1 K^1)
	!			g	= 9.80616d0,	&	! Gravity (m s^2)
	!			cp	= 1004.5d0,	&	! Specific heat capacity (J kg^-1 K^1)
	!			pi	= 4.d0*atan(1.d0)       ! pi

!-----------------------------------------------------------------------
!     Additional constants
!-----------------------------------------------------------------------

	!real(8), parameter ::	p0	= 100000.d0		! reference pressure (Pa)


CONTAINS

!==========================================================================================
! TEST CASE 11 - PURE ADVECTION - 3D DEFORMATIONAL FLOW
!==========================================================================================

! The 3D deformational flow test is based on the deformational flow test of Nair and Lauritzen (JCP 2010), 
! with a prescribed vertical wind velocity which makes the test truly 3D. An unscaled planet (with scale parameter
! X = 1) is selected. 

SUBROUTINE test1_advection_deformation (lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q,q1,q2,q3,q4)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z		! Height (m)

	real(8), intent(inout) :: p		! Pressure  (Pa)				

	integer,  intent(in) :: zcoords 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q, & 		! Specific Humidity (kg/kg)
				q1, & 		! Tracer q1 (kg/kg)
				q2, & 		! Tracer q2 (kg/kg)
				q3, & 		! Tracer q3 (kg/kg)
				q4		! Tracer q4 (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p 

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	tau     = 12.d0 * 86400.0,	&	! period of motion 12 days
			    	u0      = (2.d0*pi*a)/tau,	&	! 2 pi a / 12 days
			    	k0	= (10.d0*a)/tau,	&	! Velocity Magnitude
			    	omega0	= (23000.d0*pi)/tau,	&	! Velocity Magnitude
                            	T0      = 300.d0,		&	! temperature
                            	H       = Rd * T0 / g,		&	! scale height
                            	RR      = 1.d0/2.d0,		&	! horizontal half width divided by 'a'
                            	ZZ      = 1000.d0,		&	! vertical half width
                            	z0      = 5000.d0,		&	! center point in z
                            	lambda0 = 5.d0*pi/6.d0,		&	! center point in longitudes
                            	lambda1 = 7.d0*pi/6.d0,		&	! center point in longitudes
                            	phi0    = 0.d0,			&	! center point in latitudes
                            	phi1    = 0.d0                      
                            
      real(8) :: height							! The height of the model levels
      real(8) :: ptop							! Model top in p
      real(8) :: sin_tmp, cos_tmp, sin_tmp2, cos_tmp2			! Calculate great circle distances
      real(8) :: d1, d2, r, r2, d3, d4					! For tracer calculations 
      real(8) :: s							! Shape function
      real(8) :: lonp							! Translational longitude, depends on time
      real(8) :: time							! Initially set to zero seconds, needs
									! to be modified when used in dycore

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------
	
	! Height and pressure are aligned (p = p0 exp(-z/H))

	if (zcoords .eq. 1) then

		height = z
		p = p0 * exp(-z/H)

	else

		height = H * log(p0/p)

	endif

	! Model top in p

	ptop    = p0*exp(-12000.d0/H)	

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE TIME DEPENDENT AND THEREFORE MUST BE UPDATED
!    IN THE DYNAMICAL CORE
!-----------------------------------------------------------------------

	! These are initial conditions hence time = 0

	time = 0.d0

	! Translational longitude = longitude when time = 0

	lonp = lon - 2.d0*pi*time/tau

	! Shape function

	s = min(1.d0,2.d0*sqrt(sin(pi*(p-ptop)/(p0-ptop)) ) )

	! Zonal Velocity

	u = k0*sin(lonp)*sin(lonp)*sin(2.d0*lat)*cos(pi*time/tau) + u0*cos(lat)

	! Meridional Velocity

	v = k0*sin(2.d0*lonp)*cos(lat)*cos(pi*time/tau)

	! Vertical Velocity - can be changed to vertical pressure velocity by 
	! omega = -(g*p)/(Rd*T0)*w

	w = -((Rd*T0)/(g*p))*omega0*sin(lonp)*cos(lat)*cos(2.0*pi*time/tau)*sin(s*pi/2.d0)

!-----------------------------------------------------------------------
!    TEMPERATURE IS CONSTANT 300 K
!-----------------------------------------------------------------------

	t = T0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential) 
!-----------------------------------------------------------------------
      
	phis = 0.d0

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

!-----------------------------------------------------------------------
!     initialize tracers
!-----------------------------------------------------------------------

	! Tracer 1 - Cosine Bells

		! To calculate great circle distance

	sin_tmp = sin(lat) * sin(phi0)
	cos_tmp = cos(lat) * cos(phi0)
	sin_tmp2 = sin(lat) * sin(phi1)
	cos_tmp2 = cos(lat) * cos(phi1)

		! great circle distance without 'a'
	
	r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0)) 
	r2  = ACOS (sin_tmp2 + cos_tmp2*cos(lon-lambda1)) 
	d1 = min( 1.d0, (r/RR)**2 + ((height-z0)/ZZ)**2 )
	d2 = min( 1.d0, (r2/RR)**2 + ((height-z0)/ZZ)**2 )
	
	q1 = 0.5d0 * (1.d0 + cos(pi*d1)) + 0.5d0 * (1.d0 + cos(pi*d2))

	! Tracer 2 - Correlated Cosine Bells

	q2 = 0.9d0 - 0.8d0*q1**2

	! Tracer 3 - Slotted Ellipse

		! Make the ellipse

	!d3 = (r/RR)**2 + ((height-z0)/ZZ)**2
	!d4 = (r2/RR)**2 + ((height-z0)/ZZ)**2

	if (d1 .le. RR) then
		q3 = 1.d0
	elseif (d2 .le. RR) then
		q3 = 1.d0
	else
		q3 = 0.1d0
	endif

		! Put in the slot
	
	if (height .gt. z0 .and. abs(lat) .lt. 0.125d0) then 

		q3 = 0.1d0      

	endif

	! Tracer 4 - Sum to one

	q4 = 1.d0 - 0.3d0*(q1+q2+q3)

	! Tracer 2 3 and 4
		! Set to zero outside `buffer zone'
	
	!if (height .gt. (z0+1.25*ZZ) .or. height .lt. (z0-1.25*ZZ)) then 

	!	q2 = 0.0d0
	!	q3 = 0.0d0
	!	q4 = 0.0d0

	!endif




END SUBROUTINE test1_advection_deformation





!==========================================================================================
! TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==========================================================================================

SUBROUTINE test1_advection_hadley (lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q,q1)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z		! Height (m)

	real(8), intent(inout) :: p		! Pressure  (Pa)
				
	integer,  intent(in) :: zcoords 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q, & 		! Specific Humidity (kg/kg)
				q1 		! Tracer q1 (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	tau     = 1.d0 * 86400.d0,	&	! period of motion 1 day
			    	u0      = 40.d0,		&	! Velocity Magnitude
			    	w0	= 0.15d0,		&	! Velocity Magnitude
                            	T0      = 300.d0,		&	! temperature
                            	H       = Rd * T0 / g,		&	! scale height
                            	K       = 5.d0,			&	! number of hadley cells
                            	z1      = 2000.d0,		&	! lower tracer bound         
                            	z2      = 5000.d0,		&	! upper tracer bound         
                            	z0      = 0.5d0*(z1+z2),	&	! midpoint        
                            	ztop    = 12000.d0			! model top         
                            
      real(8) :: height							! Model level heights
      real(8) :: time							! Initially set to zero seconds, needs
									! to be modified when used in dycore

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 exp(-z/H))

	if (zcoords .eq. 1) then

		height = z
		p = p0 * exp(-z/H)

	else

		height = H * log(p0/p)

	endif

!-----------------------------------------------------------------------
!    TEMPERATURE IS CONSTANT 300 K
!-----------------------------------------------------------------------

	t = T0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
	phis = 0.d0

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE TIME DEPENDENT AND THEREFORE MUST BE UPDATED
!    IN THE DYNAMICAL CORE
!-----------------------------------------------------------------------

	time = 0.d0

	! Zonal Velocity

	u = u0*cos(lat)

	! Meridional Velocity

	v = - (a*w0*pi)/(K*ztop) *cos(lat)*sin(K*lat)*cos(pi*height/ztop)*cos(pi*time/tau)

	! Vertical Velocity - can be changed to vertical pressure velocity by 
	! omega = -(g*p)/(Rd*T0)*w

	w = (w0/K)*(-2.d0*sin(K*lat)*sin(lat) + K*cos(lat)*cos(K*lat)) &
		*sin(pi*height/ztop)*cos(pi*time/tau)


!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

!-----------------------------------------------------------------------
!     initialize tracers
!-----------------------------------------------------------------------

	! Tracer 1 - Layer

	if (height .lt. z2 .and. height .gt. z1) then
	
		q1 = 0.5d0 * (1.d0 + cos( 2.d0*pi*(height-z0)/(z2-z1) ) )

	else

		q1 = 0.d0

	endif

END SUBROUTINE test1_advection_hadley














!==========================================================================================
! TEST CASE 13 - HORIZONTAL ADVECTION OF THIN CLOUD-LIKE TRACERS IN THE PRESENCE OF OROGRAPHY
!==========================================================================================

SUBROUTINE test1_advection_orography (lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q,q1,q2,q3,q4)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z		! Height (m)

	real(8), intent(inout) :: p		! Pressure  (Pa)
				
	integer,  intent(in) :: zcoords 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q, & 		! Specific Humidity (kg/kg)
				q1, & 		! Tracer q1 (kg/kg)
				q2, &	 	! Tracer q2 (kg/kg)
				q3, &		! Tracer q3 (kg/kg)
				q4 		! Tracer q4 (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	tau     = 12.d0 * 86400.d0,	&	! period of motion 12 days
			    	u0      = 2.0*pi*a/tau,		&	! Velocity Magnitude
                            	T0      = 300.d0,		&	! temperature   
                            	H       = Rd * T0 / g,		&	! scale height  
                            	alpha   = pi/6.0,		&	! rotation angle    
                            	lambdam = 3.0*pi/2.0,		&	! rotation angle    
                            	phim    = 0.0,    		&	! rotation angle    
                            	h0      = 2000.0,		&	! rotation angle    
                            	Rm      = 3.0*pi/4.0,		&	! rotation angle    
                            	zetam   = pi/16.0,		&	! rotation angle    
                            	lambdap = pi/2.0,		&	! rotation angle    
                            	phip    = 0.0,   		&	! rotation angle    
                            	Rp      = pi/4.0,		&	! cloud deck radius    
                            	zp1     = 3050.0,		&	! midpoint            
                            	zp2     = 5050.0,		&	! midpoint            
                            	zp3     = 8200.0,		&	! midpoint          
                            	dzp1    = 1000.0,		&	! thickness            
                            	dzp2    = 1000.0,		&	! thickness            
                            	dzp3    = 400.0,		&	! thickness        
                            	ztop    = 12000.d0			! model top         
                            
      real(8) :: height							! Model level heights
      real(8) :: r							! Great circle distance
      real(8) :: rz							! height differences
      real(8) :: zs							! Surface elevation
      real(8) :: qa							! Tracer part a
      real(8) :: qb							! Tracer part b
      real(8) :: qc							! Tracer part c
      real(8) :: time							! Initially set to zero seconds, needs
									! to be modified when used in dycore

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 exp(-z/H))

	if (zcoords .eq. 1) then

		height = z
		p = p0 * exp(-z/H)

	else

		height = H * log(p0/p)

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE TIME DEPENDENT AND THEREFORE MUST BE UPDATED
!    IN THE DYNAMICAL CORE
!-----------------------------------------------------------------------

	time = 0.d0

	! Zonal Velocity

	u = u0*(cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))

	! Meridional Velocity

	v = -u0*(sin(lon)*sin(alpha))

	! Vertical Velocity - can be changed to vertical pressure velocity by 
	! omega = -(g*p)/(Rd*T0)*w

	w = 0.0

	! NOTE that if orography-following coordinates are used then the vertical 
	! velocity needs to be translated into the new coordinate system due to
	! the variation of the height along coordinate surfaces
	! See section 13 and the appendix of the test case document

!-----------------------------------------------------------------------
!    TEMPERATURE IS CONSTANT 300 K
!-----------------------------------------------------------------------

	t = T0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
	r = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

	if (r .lt. Rm) then

		zs = (h0/2.0)*(1.0+cos(pi*r/Rm))*cos(pi*r/zetam)**2.0

	else

		zs = 0.0

	endif

	phis = g*zs

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0 * exp(-zs/H)

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

!-----------------------------------------------------------------------
!     initialize tracers
!-----------------------------------------------------------------------

	! Tracer 1 - Cloud Layer
      
	r = acos( sin(phip)*sin(lat) + cos(phip)*cos(lat)*cos(lon - lambdap) )

	rz = abs(height - zp1)

	if (rz .lt. 0.5*dzp1 .and. r .lt. Rp) then

		q1 = 0.25*(1.0+cos(2.0*pi*rz/dzp1))*(1.0+cos(pi*r/Rp))

	else

		q1 = 0.0

	endif

	rz = abs(height - zp2)

	if (rz .lt. 0.5*dzp2 .and. r .lt. Rp) then

		q2 = 0.25*(1.0+cos(2.0*pi*rz/dzp2))*(1.0+cos(pi*r/Rp))

	else

		q2 = 0.0

	endif

	rz = abs(height - zp3)

	if (rz .lt. 0.5*dzp3 .and. r .lt. Rp) then

		q3 = 1.0

	else

		q3 = 0.0

	endif

	q4 = q1 + q2 + q3

END SUBROUTINE test1_advection_orography


!=========================================================================
! Test 2-0-0:  Steady-State Atmosphere at Rest in the Presence of Orography
!=========================================================================
SUBROUTINE test2_rest_mountain (hyam,hybm,lon,lat,t_profile,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
                                hyam, &		! A coefficient for hybrid-eta coordinate, at model level midpoint
				hybm		! B coefficient for hybrid-eta coordinate, at model level midpoint


	integer,  intent(in) :: t_profile 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q 		! Specific Humidity (kg/kg)

	! if t_profile = 0, then a constant lapse rate is assumed 
	! if t_profile = 1, implement different temperature profile, currently isothermal 

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
        real(8), parameter ::   T0      = 300.d0,               &       ! temperature (K)
                                gamma   = 0.0065d0,             &       ! temperature lapse rate (K/m)      
                                lambdam = 3.d0*pi/2.d0,         &       ! mountain longitude center point (radians)   
                                phim    = 0.d0,                 &       ! mountain latitude center point (radians)    
                                h0      = 2000.d0,              &       ! peak height of the mountain range (m)
!                                 zetam   = pi/16.d0,           &       ! mountain oscillation half-width (radians) 
!                                 Rm      = pi/4.d0,            &       ! mountain radius (radians)
                                Rm      = pi/24.d0                      ! mountain radius (radians)
                            
      real(8) :: height                                                 ! Model level heights (m)
      real(8) :: r                                                      ! Great circle distance (radians)
      real(8) :: zs                                                     ! Surface elevation (m)
      real(8) :: exponent                                               ! exponent: g/(Rd * gamma)
      real(8) :: exponent_rev                                           ! reversed exponent
      real(8) :: p                                                      ! Pressure  (Pa)


!-----------------------------------------------------------------------
!    compute exponents 
!-----------------------------------------------------------------------
     exponent     = g/(Rd*gamma)
     exponent_rev = 1.d0/exponent

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
        r = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

        if (r .lt. Rm) then

!               zs = (h0/2.d0)*(1.d0+cos(pi*r/Rm))*cos(pi*r/zetam)**2.d0   !  Schaer mountain 
                zs = (h0/2.d0)*(1.d0+cos(pi*r/Rm))                         !  mountain height

        else

                zs = 0.d0

        endif

        phis = g*zs

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

        ps = p0 * (1.d0 - gamma/T0*zs)**exponent

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

       p = hyam*p0 + hybm*ps       ! compute the pressure based on the surface pressure and hybrid coefficients
       height = T0/gamma * (1.d0 - (p/p0)**exponent_rev)  ! compute the height at this pressure

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE ZERO (STATE AT REST)
!-----------------------------------------------------------------------

        ! Zonal Velocity

        u = 0.d0

        ! Meridional Velocity

        v = 0.d0

        ! Vertical Velocity

        w = 0.d0

!-----------------------------------------------------------------------
!    TEMPERATURE WITH CONSTANT LAPSE RATE
!-----------------------------------------------------------------------

        if (t_profile .eq. 0) then
          t = T0 - gamma*height ! constant lapse rate
        else
          t = T0 ! isothermal
        endif

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

        rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

        q = 0.d0

END SUBROUTINE test2_rest_mountain






!==========================================================================================
! TEST CASE 2X - SCHAER MOUNTIAN - CONST AND SHEAR FLOW
!==========================================================================================

! The tests in this section examine the response of atmospheric models to flow over a 
! mountain profile, with and without wind shear. In order to ensure the simulated response
! contains both hydrostatic and non-hydrostatic features, the radius of the Earth is scaled so that the
! simulation is in the non-hydrostatic domain. We chose a non-rotating Earth with angular velocity
! = 0 s^-1 and select a reduced-size Earth with radius as = a/X. The reduction factor is set to
! X = 500. This choice leads to an Earth with a circumference at the equator of about 2pi a/X, which is 
! approximately 80 km. These underlying ideas behind the tests are based on the work of Schaer et al. (MWR 2002),
! Wedi and Smolarkiewicz (QJR 2009), and Wedi et al. (ECMWF Tech Report 2009)
! Note however that in the presence of vertical wind shear we do not recommend the isothermal conditions 
! suggested in the literature. They lead to imbalances of the initial conditions in spherical geometry

SUBROUTINE test2_schaer_mountain (lon,lat,p,z,zcoords,hybrid_eta,hyam,hybm,shear,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z, &		! Height (m)				
                                hyam, &		! A coefficient for hybrid-eta coordinate, at model level midpoint
				hybm		! B coefficient for hybrid-eta coordinate, at model level midpoint


	real(8), intent(inout) :: p		! Pressure  (Pa)
				

	integer,  intent(in) :: zcoords, &	! 0 or 1 see below
				shear, & 	! 0 or 1 see below
				hybrid_eta 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q 		! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output z
	! if zcoords = 0, then we use p

	! if shear = 1, then we use shear flow
	! if shear = 0, then we use constant u

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	X       = 500.d0,		&	! Reduced Earth reduction factor
			    	Om      = 0.d0,			&	! Rotation Rate of Earth
                            	as      = a/X,			&	! New Radius of small Earth     
			    	ueq     = 20.d0,		&	! Reference Velocity 
                            	Teq     = 300.d0,		&	! Temperature at Equator    
			    	Peq     = 100000.d0,		&	! Reference PS at Equator
                            	ztop    = 30000.d0,		&	! Model Top       
				lambdac = pi/4.d0, 		& 	! Lon of Schar Mountain Center
				phic    = 0.d0,	 		& 	! Lat of Schar Mountain Center
				h0      = 250.d0, 		& 	! Height of Mountain
				d       = 5000.d0, 		& 	! Mountain Half-Width
				xi      = 4000.d0, 		& 	! Mountain Wavelength
				cs      = 0.00025d0 		 	! Wind Shear (shear=1)
                            
      real(8) :: height							! Model level heights
      real(8) :: sin_tmp, cos_tmp					! Calculation of great circle distance
      real(8) :: r							! Great circle distance
      real(8) :: zs							! Surface height
      real(8) :: c							! Shear

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------

	sin_tmp = sin(lat) * sin(phic)
	cos_tmp = cos(lat) * cos(phic)
	
	! great circle distance with 'a/X'  

	r  = as * ACOS (sin_tmp + cos_tmp*cos(lon-lambdac))     
	zs   = h0 * exp(-(r**2)/(d**2))*(cos(pi*r/xi)**2)
	phis = g*zs

!-----------------------------------------------------------------------
!    SHEAR FLOW OR CONSTANT FLOW
!-----------------------------------------------------------------------

	if (shear .eq. 1) then

		c = cs

	else

		c = 0.d0

	endif

!-----------------------------------------------------------------------
!    TEMPERATURE 
!-----------------------------------------------------------------------

	t = Teq *(1.d0 - (c*ueq*ueq/(g))*(sin(lat)**2) )

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = peq*exp( -(ueq*ueq/(2.d0*Rd*Teq))*(sin(lat)**2) - phis/(Rd*t)    )


!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE 
!-----------------------------------------------------------------------

	if (zcoords .eq. 1) then

		height = z
		p = peq*exp( -(ueq*ueq/(2.d0*Rd*Teq))*(sin(lat)**2) - g*height/(Rd*t)    )

	else

                if (hybrid_eta .eq. 1) then
                   p = hyam*p0 + hybm*ps
                endif
		height = (Rd*t/(g))*log(peq/p) - (t*ueq*ueq/(2.d0*Teq*g))*(sin(lat)**2)

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES
!-----------------------------------------------------------------------

	! Zonal Velocity

	u = ueq * cos(lat) * sqrt( (2.d0*Teq/(t))*c*height + t/(Teq) )

	! Meridional Velocity

	v = 0.d0

	! Vertical Velocity = Vertical Pressure Velocity = 0

	w = 0.d0

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test2_schaer_mountain




SUBROUTINE test3_gravity_wave_new (lperturb,lon,lat,hyam,hybm,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	logical, intent(in)  :: lperturb        ! if true then pot. T perturbation is overlaid
        real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
                                hyam, &
                                hybm


	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q 		! Specific Humidity (kg/kg)

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
!	real(8), parameter :: 	X       = 125.d0,		&	! Reduced Earth reduction factor
        real(8)            ::  om, &                                    ! Rotation rate, depends on test case 
                               u0                                       ! Reference Velocity 

	real(8), parameter :: 	X       = 1.d0,		        &	! Regular size Earth
                            	as      = a/X,			&	! New Radius of the scaled Earth     
                            	Teq     = 300.d0,		&	! Temperature at Equator    
			    	Peq     = 100000.d0,		&	! Reference PS at Equator
                            	ztop    = 10000.d0,		&	! Model Top       
				lambdac = 2.d0*pi/3.d0, 	& 	! Lon of Pert Center
				d       = 5000.d0*125.d0, 		& 	! Width for Pert
				phic    = 0.d0,	 		& 	! Lat of Pert Center
				delta_theta = 1.d0, 		& 	! Max Amplitude of Pert
				Lz      = 20000.d0, 		& 	! Vertical Wavelength of Pert
				N       = 0.01d0, 		& 	! Brunt-Vaisala frequency
!				N       = 0.01323d0, 		& 	! Brunt-Vaisala frequency
				N2      = N*N, 		 	&	! Brunt-Vaisala frequency Squared
				bigG    = (g*g)/(N2*cp)			! Constant
                            
      real(8) :: height	, p						! Model level height
      real(8) :: sin_tmp, cos_tmp					! Calculation of great circle distance
      real(8) :: r, s							! Shape of perturbation
      real(8) :: TS 							! Surface temperature
      real(8) :: t_mean, t_pert						! Mean and pert parts of temperature
      real(8) :: theta_pert						! Pot-temp perturbation

      if (choice_testcase .eq. 31) then
        om = 0.0d0 ! no rotation
        u0 = 20.d0 ! Reference Velocity 
      else
        om = omega ! Earth's rotation
        u0 = 0.d0  ! no background flow
      endif
!-----------------------------------------------------------------------
!    THE VELOCITIES
!-----------------------------------------------------------------------

	! Zonal Velocity 

	u = u0 * cos(lat)

	! Meridional Velocity

	v = 0.d0

	! Vertical Velocity = Vertical Pressure Velocity = 0

	w = 0.d0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------

	phis = 0.d0

!-----------------------------------------------------------------------
!    SURFACE TEMPERATURE 
!-----------------------------------------------------------------------

	TS = bigG + (Teq-bigG)*exp( -(u0*N2/(4.d0*g*g))*(u0+2.d0*om*as)*(cos(2.d0*lat)-1.d0)    ) 

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = peq*exp( (u0/(4.0*bigG*Rd))*(u0+2.0*Om*as)*(cos(2.0*lat)-1.0)  ) &
		* (TS/Teq)**(cp/Rd)

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE AND MEAN TEMPERATURE
!-----------------------------------------------------------------------

        p = hyam*p0 + hybm*ps

        height = (-g/N2)*log( (TS/bigG)*( (p/ps)**(Rd/cp) - 1.d0  ) + 1.d0 )

        t_mean = bigG*(1.d0 - exp(N2*height/g))+ TS*exp(N2*height/g)

!-----------------------------------------------------------------------
!    rho (density), unperturbed using the background temperature t_mean
!-----------------------------------------------------------------------

       rho = p/(Rd*t_mean)

     if (lperturb) then ! test 31
!-----------------------------------------------------------------------
!      POTENTIAL TEMPERATURE PERTURBATION, 
!      here: converted to temperature and added to the temperature field
!      models with a prognostic potential temperature field can utilize
!      the potential temperature perturbation theta_pert directly and add it
!      to the background theta field (not included here)
!-----------------------------------------------------------------------

         sin_tmp = sin(lat) * sin(phic)
         cos_tmp = cos(lat) * cos(phic)

         ! great circle distance with 'a/X' 

          r  = as * ACOS (sin_tmp + cos_tmp*cos(lon-lambdac)) 

          s = (d**2)/(d**2 + r**2)

          theta_pert = delta_theta*s*sin(2.d0*pi*height/Lz)

          t_pert = theta_pert*(p/p0)**(Rd/cp)

          t = t_mean + t_pert
      else
!---------------------------------
!         balanced state (test 32)
!---------------------------------
          t = t_mean
      endif

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

      q = 0.d0

END SUBROUTINE test3_gravity_wave_new



!==========================================================================================
! TEST CASE 60 and 61 - Mountain generated Rossby waves, dry or moist
!==========================================================================================
  subroutine mountain_Rossby (moist,hyam,hybm,lon,lat,u,v,w,t,phis,ps,rho,q,q1,q2)

    INTEGER, INTENT(IN)  :: moist
    REAL(8), INTENT(IN)  :: &
                hyam,       & ! hybrid coefficient a at layer midpoint
                hybm,       & ! hybrid coefficient b at layer midpoint
                lon,        & ! Longitude (radians)
                lat           ! Latitude (radians)

    REAL(8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                w,          & ! Vertical Velocity (m s^-1)
                t,          & ! Temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q,          & ! Specific Humidity (kg/kg)
                q1,         & ! Tracer mixing ratio tracer 1(kg/kg)
                q2            ! Tracer mixing ratio tracer 2 (kg/kg)

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(r8),parameter :: u0      = 20.d0,                          &   ! 20 m/s
                            T0      = 288.d0,                         &   !  temperature
                            N2      = g*g/(cp*T0),                    &   !  squared Brunt Vaisala frequency N^2
                            h0      = 2000.d0,                        &   !  amplitude of the mountain, 2km
                            h1      = 4000.d0,                        &   !  amplitude of the mountain, 4km
                            d       = 1500.d3,                        &   !  half width 1500 km, single peak
                            d0      = 500.d3,                         &   !  half width 500 km, peak 1
                            d1      = 500.d3,                         &   !  half width 500 km, peak 2
                            lambda0 = 0.5d0*pi,                       &   !  center point in longitudes, peak 1
                            lambda1 = (105.d0/180.d0)*pi,             &   !  center point in longitudes, peak 2
                            phi0    = pi/6.d0,                        &   !  center point in latitudes, all peaks
                            p_sp    = 93000.d0,                       &   !  pressure at the South Pole in Pa
                            q0      = 0.8d0 * epsilon / p0 * (e0sat * exp (-L/Rv*(1/t0-1/T_0C))) ! gives around 80% RH in midlatitudes
                                                                          ! maximum specific humidity                   
!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------
      real(r8) :: sin_tmp, cos_tmp
      real(r8) :: tmp1, tmp2, tmp3
      real(r8) :: r                       ! great circle radius
      real(r8) :: pressure                ! pressure, computed with hybrid coefficients
      real(r8) :: eta                     ! hybrid level
 
!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
     u = u0 * cos(lat)
     v = 0.d0
     w = 0.d0
!-----------------------------------------------------------------------
!     initialize PHIS (surface geopotential) and the PS (surface pressure)
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2._r8 * g*g * kappa) * (u0/a + 2._r8 * omega)
      tmp2 = N2 / (g*g * kappa)
      sin_tmp = sin(lat) * sin(phi0)
      cos_tmp = cos(lat) * cos(phi0)
      tmp3 = tmp1*((sin(lat))**2 - 1._r8)
! Single peak with peak height h0 (currently set to 2000 m)
      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0)) ! great circle distance with 'a'
      phis = g*h0 * exp(-(r/d)**2)                    ! Gaussian profile of the mountain
! Double peak
!      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0)) ! great circle distance with 'a'
!      phis = g*h0 * exp(-(r/d0)**2)                    ! Gaussian profile of the mountain
!      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda1)) ! great circle distance with 'a'
!      phis = phis +  g*h1 * exp(-(r/d1)**2) ! Gaussian profile of the mountain

      if (phis .lt. 1.d-80) phis = 0.d0
      ps   = p_sp * exp( -tmp3 - tmp2*phis)

!-----------------------------------------------------------------------
!     initialize Q 
!-----------------------------------------------------------------------
      if (moist == 0) then
         q = 0.
      else
         pressure = hyam*p0 + hybm*ps
         q = q0 * pressure/ps
      endif
!-----------------------------------------------------------------------
!     consider T0 to be the virtual temperature and initialize T (temperature)
!-----------------------------------------------------------------------
      t = T0/(1+0.608d0*q)
!-----------------------------------------------------------------------
!     initialize passive tracers
!-----------------------------------------------------------------------
      eta = hyam + hybm
      if (eta.gt.0.99) then
        q1 = 1.d0                      ! lowest model level
      else
        q1 = 0.d0
      endif 

      if (eta.gt.0.83 .and. eta.lt.0.87) then
        q2 = 1.d0                      ! model level closest to the 850 hPa level away from the mountain 
      else
        q2 = 0.d0
      endif 
  end subroutine mountain_Rossby





!==========================================================================================
! TEST CASE 12 - PURE ADVECTION - 3D HADLEY-LIKE FLOW
!==========================================================================================

SUBROUTINE test1_advection_hadley_TEST (lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q,q1,q2,q3,q4)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z		! Height (m)

	real(8), intent(inout) :: p		! Pressure  (Pa)
				
	integer,  intent(in) :: zcoords 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q, & 		! Specific Humidity (kg/kg)
				q1, & 		! Tracer q1 (kg/kg)
				q2, & 		! Tracer q1 (kg/kg)
				q3, & 		! Tracer q1 (kg/kg)
				q4 		! Tracer q1 (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we use p

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	tau     = 1.d0 * 86400.d0,	&	! period of motion 1 day
			    	u0      = 40.d0,		&	! Velocity Magnitude
			    	w0	= 0.25d0,		&	! Velocity Magnitude
                            	T0      = 300.d0,		&	! temperature
                            	H       = Rd * T0 / g,		&	! scale height
                            	K       = 5.d0,			&	! number of hadley cells
                            	z1      = 3500.d0,		&	! lower tracer bound         
                            	z2      = 6500.d0,		&	! upper tracer bound         
                            	z0      = 0.5d0*(z1+z2),	&	! midpoint        
                            	z11      = 2500.d0,		&	! lower tracer bound         
                            	z21      = 5500.d0,		&	! upper tracer bound         
                            	z01      = 0.5d0*(z11+z21),	&	! midpoint        
                            	z12      = 3000.d0,		&	! lower tracer bound         
                            	z22      = 6000.d0,		&	! upper tracer bound         
                            	z02      = 0.5d0*(z12+z22),	&	! midpoint        
                            	z13      = 2000.d0,		&	! lower tracer bound         
                            	z23      = 5000.d0,		&	! upper tracer bound         
                            	z03      = 0.5d0*(z13+z23),	&	! midpoint        
                            	ztop    = 12000.d0			! model top         
                            
      real(8) :: height							! Model level heights
      real(8) :: time							! Initially set to zero seconds, needs
									! to be modified when used in dycore

!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 exp(-z/H))

	if (zcoords .eq. 1) then

		height = z
		p = p0 * exp(-z/H)

	else

		height = H * log(p0/p)

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE TIME DEPENDENT AND THEREFORE MUST BE UPDATED
!    IN THE DYNAMICAL CORE
!-----------------------------------------------------------------------

	time = 0.d0

	! Zonal Velocity

	u = u0*cos(lat)

	! Meridional Velocity

	v = - (a*w0*pi)/(K*ztop) *cos(lat)*sin(K*lat)*cos(pi*height/ztop)*cos(pi*time/tau)

	! Vertical Velocity - can be changed to vertical pressure velocity by 
	! omega = -(g*p)/(Rd*T0)*w

	w = (w0/K)*(-2.d0*sin(K*lat)*sin(lat) + K*cos(lat)*cos(K*lat)) &
		*sin(pi*height/ztop)*cos(pi*time/tau)

        ! FOR THE TEST WE DIVIDE BY W BY RHO => omega = -g w


!-----------------------------------------------------------------------
!    TEMPERATURE IS CONSTANT 300 K
!-----------------------------------------------------------------------

	t = T0

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
	phis = 0.d0

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

!-----------------------------------------------------------------------
!     initialize tracers
!-----------------------------------------------------------------------

	! Tracer 1 - Layer

	if (height .lt. z2 .and. height .gt. z1) then
	
		q1 = 0.5d0 * (1.d0 + cos( 2.d0*pi*(height-z0)/(z2-z1) ) )

	else

		q1 = 0.d0

	endif

	if (height .lt. z21 .and. height .gt. z11) then
	
		q2 = 0.5d0 * (1.d0 + cos( 2.d0*pi*(height-z01)/(z21-z11) ) )

	else

		q2 = 0.d0

	endif

	if (height .lt. z22 .and. height .gt. z12) then
	
		q3 = 0.5d0 * (1.d0 + cos( 2.d0*pi*(height-z02)/(z22-z12) ) )

	else

		q3 = 0.d0

	endif

	if (height .lt. z23 .and. height .gt. z13) then
	
		q4 = 0.5d0 * (1.d0 + cos( 2.d0*pi*(height-z03)/(z23-z13) ) )

	else

		q4 = 0.d0

	endif

END SUBROUTINE test1_advection_hadley_TEST


END MODULE dcmip_initial_conditions_test_123
