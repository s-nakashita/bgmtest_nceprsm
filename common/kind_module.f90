module kind_module
!
! define kind types
!
  implicit none
  public
!
! symbolic names for kind types of 4-, 2-, and 1-byte integers:
!
  integer, parameter :: i8b = selected_int_kind(15)
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i2b = selected_int_kind(4)
!
! symbolic names for kind types of single- and double-precision reals:
!
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)
!  integer, parameter :: qp = kind(1.0q0)
!
! symbolic names for kind types of single- and double-precision complexes:
!
  integer, parameter :: spc = kind((1.0,1.0))
  integer, parameter :: dpc = kind((1.0d0,1.0d0))
!
! missing value
!
  real(kind=dp), parameter :: undef=-9.99d33
!
  integer, parameter :: filelenmax=256
!
end module kind_module
