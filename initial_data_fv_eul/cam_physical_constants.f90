module cam_physical_constants
  
  implicit none

  integer, parameter :: r8 = SELECTED_REAL_KIND(12) ! 8 byte real

  real(r8), parameter ::           &
       p0         = 100000._r8,    &! reference pressure 
       Rd         = 287.0_r8,     &! gas constant J/(K kg)
       cp         = 1004.5_r8,    &! specific heat at constant pressure J/(K kg)
       kappa      = Rd/cp,         &! kappa = 2/7
       g          = 9.80616_r8,    &! gravitational acceleration (m/s^2)
       a          = 6371220._r8,   &! Earth's radius in m
       pi         = 3.14159265358979323846_r8,&  ! pi
       omega      = 2._r8*pi/86164._r8, &! Earth's angular velocity 1/s
       !omega      = 2._r8*pi/(0.05*86164._r8), &! Small Earth's angular velocity 1/s
       pih        = pi*0.5_r8,     &! pi/2
       deg2rad    = pi/180._r8,    &
       L          = 2.5e6_r8,      & ! latent heat of vaporization
       e0sat      = 610.78_r8,     & ! saturation vapor pressure at T_0C
       T_0C       = 273.16_r8,     & ! freezing point
       Rv         = 461.5_r8,      & ! gas constant for water vapor
       epsilon    = 0.622_r8,      & ! Rv/Rd -1
       Mvap       = 0.608_r8         ! Ratio of molar mass dry air/water vapor

end module cam_physical_constants
