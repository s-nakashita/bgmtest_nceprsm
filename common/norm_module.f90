module norm_module
!
! module : norm_module   calculate various norm
!
! history:
! 22-09-13 create
! 
  use kind_module
  use phconst_module
  use read_module
  use write_module

  public :: calc_te, calc_te2, calc_tegrd, calc_teprof
  contains
!======================================================================
! calculate moist total energy
!======================================================================
subroutine calc_te(u,v,t,q,ps,epsq,clat,si,nlon,nlat,kmax,te)
  implicit none
  integer, intent(in) :: nlon, nlat, kmax ! boundaries
  real(kind=dp), intent(in) :: u(:,:,:),v(:,:,:),t(:,:,:),q(:,:,:)
  real(kind=dp), intent(in) :: ps(:,:)
  real(kind=dp), intent(in) :: epsq ! weight for moist term
  real(kind=dp), intent(in) :: clat(:),si(:)
  real(kind=dp), intent(out):: te(4)
  real(kind=dp) :: area,coef
  integer :: igrd1, jgrd1
  integer :: n,i,j,k

  ! calculate energy
  te=0.0d0
  area=0.0d0
  do k=1,kmax
    do j=1,nlat
      coef=(si(k)-si(k+1))*cos(clat(j)*deg2rad)
      do i=1,nlon
        !KE
        te(1)=te(1)+(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k))*coef
        !PE(T)
        te(2)=te(2)+cp/tr*t(i,j,k)*t(i,j,k)*coef
        !LE
        te(3)=te(3)+epsq*lh**2/cp/tr*q(i,j,k)*q(i,j,k)*coef
      end do
    end do
  end do
  do j=1,nlat
    coef=cos(clat(j)*deg2rad)
    do i=1,nlon
      !PE(Ps)
      te(4)=te(4)+rd*tr*ps(i,j)*ps(i,j)/pr/pr*coef
      area=area+coef
    end do
  end do
  do i=1,4
    te(i)=te(i)*0.5d0/area
  end do
  return
end subroutine calc_te

!======================================================================
! calculate moist total energy from 2 different perturbations
!======================================================================
subroutine calc_te2(u1,u2,v1,v2,t1,t2,q1,q2,ps1,ps2,&
                epsq,clat,si,nlon,nlat,kmax,te)
  implicit none
  integer, intent(in) :: nlon, nlat, kmax ! boundaries
  real(kind=dp), intent(in) :: u1(:,:,:),v1(:,:,:),t1(:,:,:),q1(:,:,:)
  real(kind=dp), intent(in) :: u2(:,:,:),v2(:,:,:),t2(:,:,:),q2(:,:,:)
  real(kind=dp), intent(in) :: ps1(:,:),ps2(:,:)
  real(kind=dp), intent(in) :: epsq ! weight for moist term
  real(kind=dp), intent(in) :: clat(:),si(:)
  real(kind=dp), intent(out):: te(4)
  real(kind=dp) :: area,coef
  integer :: igrd1, jgrd1
  integer :: n,i,j,k

  ! calculate energy
  te=0.0d0
  area=0.0d0
  do k=1,kmax
    do j=1,nlat
      coef=(si(k)-si(k+1))*cos(clat(j)*deg2rad)
      do i=1,nlon
        !KE
        te(1)=te(1)+(u1(i,j,k)*u2(i,j,k)+v1(i,j,k)*v2(i,j,k))*coef
        !PE(T)
        te(2)=te(2)+cp/tr*t1(i,j,k)*t2(i,j,k)*coef
        !LE
        te(3)=te(3)+epsq*lh**2/cp/tr*q1(i,j,k)*q2(i,j,k)*coef
      end do
    end do
  end do
  do j=1,nlat
    coef=cos(clat(j)*deg2rad)
    do i=1,nlon
      !PE(Ps)
      te(4)=te(4)+rd*tr*ps1(i,j)*ps2(i,j)/pr/pr*coef
      area=area+coef
    end do
  end do
  do i=1,4
    te(i)=te(i)*0.5d0/area
  end do
  return
end subroutine calc_te2
!======================================================================
! calculate moist total energy for each grids
!======================================================================
subroutine calc_tegrd(u,v,t,q,ps,epsq,clat,si,nlon,nlat,kmax,te)
  implicit none
  integer, intent(in) :: nlon, nlat, kmax! boundaries
  real(kind=dp), intent(in) :: u(:,:,:),v(:,:,:),t(:,:,:),q(:,:,:)
  real(kind=dp), intent(in) :: ps(:,:)
  real(kind=dp), intent(in) :: epsq ! weight for moist term
  real(kind=dp), intent(in) :: clat(:),si(:)
  real(kind=dp), intent(out):: te(:,:,:,:)
  real(kind=dp) :: area,coef
  integer :: igrd1, jgrd1
  integer :: n,i,j,k

  ! calculate energy
  te=0.0d0
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        !KE
        te(i,j,k,1)=0.5d0*(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k))
        !PE(T)
        te(i,j,k,2)=0.5d0*cp/tr*t(i,j,k)*t(i,j,k)
        !LE
        te(i,j,k,3)=0.5d0*epsq*lh**2/cp/tr*q(i,j,k)*q(i,j,k)
      end do
    end do
  end do
  do j=1,nlat
    do i=1,nlon
      !PE(Ps)
      te(i,j,1,4)=0.5d0*rd*tr*ps(i,j)*ps(i,j)/pr/pr
    end do
  end do
  return
end subroutine calc_tegrd
!======================================================================
! calculate moist total energy profile
!======================================================================
subroutine calc_teprof(u,v,t,q,ps,epsq,clat,si,nlon,nlat,kmax,vwgt,te)
  implicit none
  integer, intent(in) :: nlon, nlat, kmax ! boundaries
  real(kind=dp), intent(in) :: u(:,:,:),v(:,:,:),t(:,:,:),q(:,:,:)
  real(kind=dp), intent(in) :: ps(:,:)
  real(kind=dp), intent(in) :: epsq ! weight for moist term
  real(kind=dp), intent(in) :: clat(:),si(:)
  real(kind=dp), intent(out):: vwgt(kmax),te(kmax,4)
  real(kind=dp) :: area,coef
  integer :: igrd1, jgrd1
  integer :: n,i,j,k

  ! calculate energy
  te=0.0d0
  area=0.0d0
  do k=1,kmax
    vwgt(k) = si(k) - si(k+1)
    do j=1,nlat
      coef=cos(clat(j)*deg2rad)
      do i=1,nlon
        !KE
        te(k,1)=te(k,1)+(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k))*coef
        !PE(T)
        te(k,2)=te(k,2)+cp/tr*t(i,j,k)*t(i,j,k)*coef
        !LE
        te(k,3)=te(k,3)+epsq*lh**2/cp/tr*q(i,j,k)*q(i,j,k)*coef
      end do
    end do
  end do
  do j=1,nlat
    coef=cos(clat(j)*deg2rad)
    do i=1,nlon
      !PE(Ps)
      te(1,4)=te(1,4)+rd*tr*ps(i,j)*ps(i,j)/pr/pr*coef
      area=area+coef
    end do
  end do
  do i=1,4
    do k=1,kmax
      te(k,i)=te(k,i)*0.5d0/area
    end do
  end do
  return
end subroutine calc_teprof

end module norm_module
