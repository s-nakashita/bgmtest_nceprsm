module date_module
!
! date routines
! [note] use library w3_4 in sys/lib
!
! history:
! 24-03-19 create
!
  use kind_module
  implicit none
  public :: ndate, nhour
contains
!
! calculate date before or after several minutes
!
  subroutine ndate(date0,dt,date1)
    implicit none
    integer, intent(in)  :: date0(5) !year,month,day,hour,minutes
    integer, intent(in)  :: dt    !minutes
    integer, intent(out) :: date1(5) !year,month,day,hour,minutes 
    integer :: idat(8),jdat(8) !year,month,day,timezone,hour,minutes,seconds,milliseconds
    real(kind=sp) :: rinc(5) !days,hours,minutes,seconds,milliseconds

!    print *, date0
    idat = 0
    idat(1:3) = date0(1:3)
    idat(5:6) = date0(4:5)
!    print *, idat
    rinc = 0.0
    rinc(3) = real(dt,kind=sp)
!    print *, rinc

    call w3movdat(rinc,idat,jdat)
!    print *, jdat
    date1(1:3)=jdat(1:3)
    date1(4:5)=jdat(5:6)

    return
  end subroutine ndate
!
! calculate minutes between two dates
!
  subroutine nhour(date0,date1,dt)
    implicit none
    integer, intent(in)  :: date0(5) !year,month,day,hour,minutes
    integer, intent(in)  :: date1(5) !year,month,day,hour,minutes 
    integer, intent(out) :: dt    !minutes (date1 - date0)
    integer :: idat(8),jdat(8) !year,month,day,timezone,hour,minutes,seconds,milliseconds
    real(kind=sp) :: rinc(5) !days,hours,minutes,seconds,milliseconds

!    print *, date0
    idat = 0
    idat(1:3) = date0(1:3)
    idat(5:6) = date0(4:5)
!    print *, idat
!    print *, date1
    jdat = 0
    jdat(1:3) = date1(1:3)
    jdat(5:6) = date1(4:5)
!    print *, jdat

    call w3difdat(jdat,idat,3,rinc)
!    print *, rinc
    dt = nint(rinc(3))

    return
  end subroutine nhour

end module date_module
