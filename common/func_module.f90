module func_module
!
! several functions
!
! history:
! 22-12-05 SN create
!
  use kind_module
  use phconst_module
  implicit none
  private
!
  public :: calc_pfull, conv_temp, calc_rh, calc_td, calc_q, calc_q2,&
   &        calc_wd, calc_uv, prsadj, distll_1
  contains
!
! p_full
!
  subroutine calc_pfull(ix,jy,kz,sig,ps,p_full)
    implicit none
    integer, intent(in)       :: ix, jy, kz
    real(kind=dp), intent(in) :: sig(kz)
    real(kind=dp), intent(in) :: ps(ix,jy)
    real(kind=dp), intent(out):: p_full(ix,jy,kz)
    integer :: i,j,k
    do k=1,kz
      do j=1,jy
        do i=1,ix
          p_full(i,j,k) = ps(i,j) * sig(k)
        end do
      end do
    end do
    return
  end subroutine calc_pfull
!
! virtual temperature <=> temperature
!
  subroutine conv_temp(ix,jy,kz,t,q,conv)
    implicit none
    integer, intent(in) :: ix, jy, kz
    real(kind=dp), intent(inout) :: t(ix,jy,kz) 
    real(kind=dp), intent(in) :: q (ix,jy,kz) !kg/kg
    integer, intent(in) :: conv !=0:tv=>t, =1:t=>tv
    real(kind=dp) :: fact
    integer :: i,j,k

    do k=1,kz
      do j=1,jy
        do i=1,ix
          fact = 1.0_dp + fvirt*q(i,j,k)
          if(conv.eq.0) then
            t(i,j,k) = t(i,j,k)/fact
          else
            t(i,j,k) = t(i,j,k)*fact
          end if
        end do
      end do
    end do
    return
  end subroutine conv_temp
!
! relative humidity
!
  subroutine calc_rh(t,q,p,rh)
    implicit none
    real(kind=dp),parameter :: e0=6.112_dp
    real(kind=dp),parameter :: a=17.67_dp
!    real(kind=dp),parameter :: a=19.482_dp
    real(kind=dp),parameter :: b=243.5_dp
!    real(kind=dp),parameter :: c=4304.4_dp
!    real(kind=dp),parameter :: e0l=6.11_dp
!    real(kind=dp),parameter :: al=17.3_dp
!    real(kind=dp),parameter :: bl=237.3_dp
!    real(kind=dp),parameter :: e0i=6.1121_dp
!    real(kind=dp),parameter :: ai=22.587_dp
!    real(kind=dp),parameter :: bi=273.86_dp
    real(kind=dp),intent(in) :: t,q,p !t[K],q[kg/kg],p[Pa]
    real(kind=dp),intent(out):: rh
    real(kind=dp) :: e,es,tc

    if(q.lt.0.0_dp) then
      rh = 0.0_dp
      return
    end if
    e = q * p * 0.01_dp / (0.378_dp * q + 0.622_dp) ![hPa]

    tc = t - t0
!    es = exp(a-c/(tc+b)) !WMO
    es = e0*exp(a*tc/(tc+b)) !Bolton(1980)
!    if(tc.ge.0.0_dp) then
!      es = e0l * exp(al*tc/(bl+tc))
!    else if(tc.le.-15.0_dp) then
!      es = e0i * exp(ai*tc/(bi+tc))
!    else
!      es = e0l * exp(al*tc/(bl+tc)) * (15.0_dp+tc)/15.0_dp &
!       & + e0i * exp(ai*tc/(bi+tc)) * (-tc)/15.0_dp
!    end if

    rh = e/es

    return
  end subroutine calc_rh
!
! dewpoint temperature
!
  subroutine calc_td(q,p,td)
    implicit none
    real(kind=dp),parameter :: a=19.48_dp
    real(kind=dp),parameter :: b=243.5_dp
    real(kind=dp),parameter :: c=440.8_dp
!    real(kind=dp),parameter :: al=6.116441_dp
!    real(kind=dp),parameter :: ml=7.591386_dp
!    real(kind=dp),parameter :: tnl=240.7263_dp
!    real(kind=dp),parameter :: ai=6.114742_dp
!    real(kind=dp),parameter :: mi=9.778707_dp
!    real(kind=dp),parameter :: tni=273.1466_dp

    real(kind=dp),intent(in) :: q,p
    real(kind=dp),intent(out):: td
    real(kind=dp) :: es,lnes!,fact,tc

    if(q.lt.0.0_dp) then
      print *, 'humidity is less than 0, q=',q
      td = 0.0_dp
      return
    end if
    es = q * p * 0.01_dp / (0.378_dp * q + 0.622_dp)
    lnes = log(es)

    td = (b*lnes - c)/(a-lnes)
!    tc = t - t0
!    if(tc.ge.-15.0_dp) then
!      fact = ml / log10(e/al) - 1.0_dp
!      td = tnl / fact
!    else
!      fact = mi / log10(e/ai) - 1.0_dp
!      td = tni / fact
!    end if
    td = td + t0
    return
  end subroutine calc_td
!
! specific humidity from dewpoint temperature
!
  subroutine calc_q(td,p,q)
    implicit none
    real(kind=dp),parameter :: e0=6.112_dp
    real(kind=dp),parameter :: a=17.67_dp
    real(kind=dp),parameter :: b=243.5_dp
!    real(kind=dp),parameter :: al=6.116441_dp
!    real(kind=dp),parameter :: ml=7.591386_dp
!    real(kind=dp),parameter :: tnl=240.7263_dp
!    real(kind=dp),parameter :: ai=6.114742_dp
!    real(kind=dp),parameter :: mi=9.778707_dp
!    real(kind=dp),parameter :: tni=273.1466_dp

    real(kind=dp),intent(in) :: td,p
    real(kind=dp),intent(out):: q
    real(kind=dp) :: e,eps,tc

    tc = td - t0
    e = e0*exp(a*tc/(tc+b)) !Bolton(1980)
    eps = rd / rv
    q = eps*e / (p*0.01_dp - (1.0_dp - eps)*e)

    return
  end subroutine calc_q
!
! specific humidity from relative humidity
!
  subroutine calc_q2(t,rh,p,q)
    implicit none
    real(kind=dp),parameter :: e0=6.112_dp
    real(kind=dp),parameter :: a=17.67_dp
    real(kind=dp),parameter :: b=243.5_dp
!    real(kind=dp),parameter :: al=6.116441_dp
!    real(kind=dp),parameter :: ml=7.591386_dp
!    real(kind=dp),parameter :: tnl=240.7263_dp
!    real(kind=dp),parameter :: ai=6.114742_dp
!    real(kind=dp),parameter :: mi=9.778707_dp
!    real(kind=dp),parameter :: tni=273.1466_dp

    real(kind=dp),intent(in) :: t,rh,p
    real(kind=dp),intent(out):: q
    real(kind=dp) :: e,es,eps,tc

    tc = t - t0
    es = e0*exp(a*tc/(tc+b)) !Bolton(1980)
    e = es * rh
    eps = rd / rv
    q = eps*e / (p*0.01_dp - (1.0_dp - eps)*e)

    return
  end subroutine calc_q2
!
! wind direction
!
  subroutine calc_wd(u,v,wd)
    implicit none
    real(kind=dp), intent(in) :: u,v
    real(kind=dp), intent(out):: wd
    real(kind=dp) :: theta

    if(u.eq.0.0_dp) then
      if(v.ge.0.0_dp) then
        wd = 180.0_dp
      else
        wd = 0.0_dp
      end if
    else if(v.eq.0.0_dp) then
      if(u.ge.0.0_dp) then
        wd = 270.0_dp
      else
        wd = 90.0_dp
      end if
    else
      theta = atan( u/v )
      if(theta.ge.0.0_dp) then
        if(u.lt.0.0_dp) then
          wd = theta * rad2deg
        else
          wd = theta * rad2deg + 180.0_dp
        end if
      else
        if(u.lt.0.0_dp) then
          wd = theta * rad2deg + 180.0_dp
        else
          wd = theta * rad2deg + 360.0_dp
        end if
      end if
    end if

    return
  end subroutine calc_wd
!
! wind components
!
  subroutine calc_uv(ws,wd,u,v)
    implicit none
    real(kind=dp), intent(in) :: ws, wd
    real(kind=dp), intent(out):: u, v
    real(kind=dp) :: theta

    if(wd.eq.0.0_dp) then
      u=0.0_dp
      v=-1.0_dp*ws
    else if(wd.eq.90.0_dp) then
      u=-1.0_dp*ws
      v=0.0_dp
    else if(wd.eq.180.0_dp) then
      u=0.0_dp
      v=ws
    else if(wd.eq.270.0_dp) then
      u=ws
      v=0.0_dp
    else
      theta = wd*deg2rad - pi
      u = ws*sin(theta)
      v = ws*cos(theta)
    end if

    return
  end subroutine calc_uv
!
! pressure adjustment for different height
!
  subroutine prsadj(p,dz,t,q)
    implicit none
    real(kind=dp),intent(inout) :: p !(in)original (out)target
    real(kind=dp),intent(in)    :: dz !target - original level
    real(kind=dp),intent(in)    :: t !original level
    real(kind=dp),intent(in)    :: q !original level
    real(kind=dp) :: tv

    if(dz.ne.0.0) then
      tv = t * (1.0_dp + fvirt*q)
      p = p * ((1.0_dp - lapse/tv*dz)**(grav/rd/lapse))
    end if
    return
  end subroutine prsadj
!
! distance between two points (alon,alat)-(blon,blat)
!
  subroutine distll_1(alon,alat,blon,blat,dist)
    implicit none
    real(kind=dp),intent(in) :: alon, alat, blon, blat
    real(kind=dp),intent(out):: dist
    real(kind=dp) :: lon1,lon2,lat1,lat2,cosd

    lon1 = alon * deg2rad
    lon2 = blon * deg2rad
    lat1 = alat * deg2rad
    lat2 = blat * deg2rad

    cosd = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1)
    cosd = min(1.0_dp,cosd)
    cosd = max(-1.0_dp,cosd)

    dist = acos(cosd) * re

    return
  end subroutine distll_1
!
end module func_module
