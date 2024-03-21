module phconst_module
!
! physics constants
!
  use kind_module
  implicit none
  public
  real(kind=dp), parameter :: pi=acos(-1.0d0)
  real(kind=dp), parameter :: rad2deg=180.0d0/pi, deg2rad=pi/180.0d0
  real(kind=dp), parameter :: re=6.371d6      ! radius of Earth [m]
  real(kind=dp), parameter :: grav=9.80665d0  ! acceleration rat of the Earth [m/s^2]
  real(kind=dp), parameter :: cp=1005.7d0     ! specific heat at constant pressure [J/kg/K]
  real(kind=dp), parameter :: rd=287.04d0     ! gas constant of dry air [J/kg/K]
  real(kind=dp), parameter :: rv=461.5d0      ! gas constant of moist air [J/kg/K]
  real(kind=dp), parameter :: fvirt=rv/rd-1.0 ! factor for converting into virtual temperature: Tv=(1.0+fvirt*q)*T
  real(kind=dp), parameter :: lh=2.5104d6     ! latent heat of condensation [J/kg]
  real(kind=dp), parameter :: lapse=0.65d-2   ! lapse rate for normal atmosphere [K/m]
  real(kind=dp), parameter :: t0=273.15_dp    ! absolute temperature [K]
  ! reference temperature and pressure for energy calculation
  real(kind=dp), parameter :: tr=300.0d0, pr=800.0d2![Pa]
end module phconst_module
