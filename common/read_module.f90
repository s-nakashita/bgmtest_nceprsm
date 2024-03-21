module read_module
!
! module : read_module   reading sigma, surface, flux files
! history:
! 22-07-28 create
! 
! namelist:
!
  use kind_module
  use phconst_module
  use func_module, only : conv_temp
  implicit none
  private
  integer, parameter, public :: levmax=100, nwext=512-(6+2*levmax)
  integer, parameter, public :: lsoil=2, nfldsfc=26
  integer, parameter, public :: nfldflx=56
  logical, parameter, public :: verbose=.false.

  public :: read_header, read_sig, read_sfc, read_flx
contains
!======================================================================
! read sigma file header
!======================================================================
subroutine read_header(iunit,icld,label,idate,fhour,si,sl,ext,nflds)
! input :
! iunit sigma file (sequential, with header)
! icld  flag for 3D physics
! 
! output: (header components and field number)
! label, idate, fhour, sisl, ext, nflds
!
  implicit none
  integer, intent(in) :: iunit
  integer, intent(in) :: icld !1=include 3D physics, 0=not include
  character(len=8),intent(out) :: label(4)
  integer,intent(out) :: idate(4), nflds
  real(kind=dp),intent(out) :: si(levmax+1), sl(levmax)
  real(kind=sp),intent(out) :: fhour, ext(nwext) 
  integer :: iret
  integer :: iymdh
  real(kind=sp) :: sisl(2*levmax+1)
  ! components of ext
  integer :: iwav1,jwav1,igrd1,jgrd1,levs,nfldx,proj,nonhyd
  real(kind=sp) :: rtruth, rorient, rcenlat, rcenlon, rgrdlft, rgrdbtm, &
 &                delx, dely
  
  print *, 'read and extract header record'
  iret = 0
  rewind(iunit)
! read label
  read(iunit) label
  print '(4a8)', label
! read header
  read(iunit) fhour, idate, sisl, ext
  iymdh = idate(4)*1000000+idate(2)*10000+idate(3)*100+idate(1)
  print *, 'posting date ', iymdh, '+', nint(fhour)
  iwav1 = int(ext(1)); jwav1 = int(ext(2))
  igrd1 = int(ext(3)); jgrd1 = int(ext(4))
  levs  = int(ext(5)); nfldx = int(ext(6))
  proj  = int(ext(7))
  rtruth= ext(8); rorient=ext(9)
  rcenlat=ext(10); rcenlon=ext(11)
  rgrdlft=ext(12); rgrdbtm=ext(13)
  delx  = ext(14); dely = ext(15)
  nonhyd= int(ext(16))
  si(1:levs+1) = sisl(1:levs+1)
  sl(1:levs)   = sisl(levs+2:2*levs+1)
  if(verbose) then
  print *, 'iwav1, jwav1'
  print '(2i4)', iwav1, jwav1
  print *, 'igrd1, jgrd1'
  print '(2i4)', igrd1, jgrd1
  print *, 'levs , nfldx'
  print '(2i4)', levs , nfldx
  print '(a,i4)', 'proj ', proj
  print *, 'rtruth, rorient, rcenlat, rcenlon'
  print '(4f9.3)', rtruth, rorient, rcenlat, rcenlon
  print *, 'rgrdlft, rgrdbtm, delx, dely'
  print '(4f9.3)', rgrdlft, rgrdbtm, delx, dely
  if ( nonhyd.eq.1 ) then
    print *, 'Input is a nonhydrostatic'
  else
    print *, 'Input is a hydrostatic'
  end if
  print *, 'si', si(1:levs+1)
  print *, 'sl', sl(1:levs)
  end if
  !nflds = 8+13*levs
  nflds = 2+9*levs+1 
  !gz,lnps,t(levs),u(levs),v(levs),q(levs),oz(levs),cw(levs),pn(levs),tn(levs),wn(levs+1)
  if((icld.eq.1).and.(fhour.gt.0.0)) then !phys3d
    nflds = nflds + 3*levs
  end if
  print *, 'nflds', nflds
  return
end subroutine read_header
!======================================================================
! read sigma file
!======================================================================
subroutine read_sig(iunit,igrd1,jgrd1,levs,nflds,nonhyd,icld,fhour,sl,&
                  & dfld,mapf,clat,clon,convert)
! input :
! iunit sigma file (sequential, with header)
! 
! output:
! dfld  variables (double precision)
!
  implicit none
  integer, intent(in) :: iunit
  integer, intent(in) :: igrd1, jgrd1, levs, nflds 
  integer, intent(in) :: nonhyd
  integer, intent(in) :: icld !1=include 3D physics, 0=not include
  real(kind=sp), intent(in)  :: fhour
  real(kind=dp), intent(in)  :: sl(levmax)
  real(kind=dp), intent(out) :: dfld(igrd1,jgrd1,nflds)
  real(kind=dp), intent(out) :: mapf(igrd1,jgrd1,3) !map factor
  real(kind=dp), intent(out) :: clat(jgrd1),clon(igrd1)
  logical, intent(in), optional :: convert
  logical :: convert_
  integer :: iret
  integer :: nwf, irec
  integer :: i,j,k,l,m
  integer :: igz, ips, it, iu, iv, iq, ioz, icw, ipn, itn, iwn, &
 &           im2
  integer :: iphys3d(3)
  real(kind=sp), allocatable :: sfld(:)
  real(kind=sp), allocatable :: factor(:,:,:)
  
  convert_=.true.
  if( present(convert) ) then
    convert_=convert
  end if
  iret = 0
  rewind(iunit)
! read label
  read(iunit) 
! read header
  read(iunit)
  if ( nonhyd.eq.1 ) then
    print *, 'Input is a nonhydrostatic'
  else
    print *, 'Input is a hydrostatic'
  end if
  nwf = igrd1*jgrd1
  print *, 'nwf', nwf
  allocate( sfld(nwf) )
  allocate( factor(igrd1,jgrd1,levs) )
  print *, 'start reading sigma data'
  l=1
  ! gz
  igz=1
  read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,igz) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,igz, 'read gz ', dfld(1,1,igz), maxval(dfld(:,:,igz)), minval(dfld(:,:,igz))
  ! ln(ps) [kPa]
  ips=igz+1
  read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ips) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(convert_) then
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,ips) = EXP(dfld(i,j,ips))*1000.0 !kPa=>Pa
      end do
    end do 
  end if
  if(verbose) print *,ips, 'read ps ', dfld(1,1,ips), maxval(dfld(:,:,ips)), minval(dfld(:,:,ips))
  ! Tv
  it=ips+1
  do k=1, levs
    read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,it+k-1) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,it+k-1, 'read Tv at lev=',k, dfld(1,1,it+k-1),&
&  maxval(dfld(:,:,it+k-1)), minval(dfld(:,:,it+k-1))
  end do  
  ! U,V
  iu=it+levs
  iv=iu+levs
  do k=1, levs
    read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,iu+k-1) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
    read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,iv+k-1) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *,iu+k-1, 'read U at lev=',k, dfld(1,1,iu+k-1),&
&  maxval(dfld(:,:,iu+k-1)), minval(dfld(:,:,iu+k-1))
  if(verbose) print *,iv+k-1, 'read V at lev=',k, dfld(1,1,iv+k-1),&
&  maxval(dfld(:,:,iv+k-1)), minval(dfld(:,:,iv+k-1))
  end do
  ! Q
  iq = iv+levs
  do k=1, levs
    read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,iq+k-1) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *,iq+k-1, 'read Q at lev=',k, dfld(1,1,iq+k-1), &
&  maxval(dfld(:,:,iq+k-1)), minval(dfld(:,:,iq+k-1))
  end do
  ! OZ
  ioz=iq+levs
  do k=1,levs
    read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ioz+k-1)=real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *,ioz+k-1, 'read OZ at lev=',k, dfld(1,1,ioz+k-1),&
&  maxval(dfld(:,:,ioz+k-1)), minval(dfld(:,:,ioz+k-1))
  end do
  ! CW
  icw=ioz+levs
  do k=1,levs
    read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,icw+k-1)=real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *,icw+k-1, 'read CW at lev=',k, dfld(1,1,icw+k-1),&
&  maxval(dfld(:,:,icw+k-1)), minval(dfld(:,:,icw+k-1))
  end do
  if(convert_) then
    ! modify the virtual temperature into temperature
    call conv_temp(igrd1,jgrd1,levs,dfld(:,:,it:it+levs-1),dfld(:,:,iq:iq+levs-1),0)
  end if
  ipn=icw+levs
  itn=ipn+levs
  iwn=itn+levs
  if (nonhyd.eq.1) then
    ! pn
    do k=1,levs
      read(iunit) (sfld(i),i=1,nwf)
      do j=1,jgrd1
        do i=1,igrd1
          dfld(i,j,ipn+k-1) = real(sfld(i+(j-1)*igrd1),kind=dp)
!        if (k.eq.1) then
!          dfld(i,j,ips) = dfld(i,j,ipn)/sl(1)
!        end if
        end do
      end do
      if(convert_) then
        do j=1,jgrd1
          do i=1,igrd1
            dfld(i,j,ipn+k-1) = EXP(dfld(i,j,ipn+k-1))*1000.0 !kPa=>Pa
          end do
        end do
      end if
      if(verbose) print *,ipn+k-1, 'read PN at lev=',k, dfld(1,1,ipn+k-1),&
    &  maxval(dfld(:,:,ipn+k-1)), minval(dfld(:,:,ipn+k-1))
    end do
    ! tn
    do k=1,levs
      read(iunit) (sfld(i),i=1,nwf)
      do j=1,jgrd1
        do i=1,igrd1
          dfld(i,j,itn+k-1) = real(sfld(i+(j-1)*igrd1),kind=dp)
          !dfld(i,j,itn+k-1) = dfld(i,j,itn+k-1)/factor(i,j,k) !Tv=>T
        end do
      end do
      if(verbose) print *,itn+k-1, 'read TvN at lev=',k, dfld(1,1,itn+k-1),&
    &  maxval(dfld(:,:,itn+k-1)), minval(dfld(:,:,itn+k-1))
    end do  
    if(convert_) then
      call conv_temp(igrd1,jgrd1,levs,dfld(:,:,itn:itn+levs-1),dfld(:,:,iq:iq+levs-1),0)
    end if
    ! wn
    do k=1,levs+1
      read(iunit) (sfld(i),i=1,nwf)
      do j=1,jgrd1
        do i=1,igrd1
          dfld(i,j,iwn+k-1) = real(sfld(i+(j-1)*igrd1),kind=dp)
        end do
      end do
    if(verbose) print *,iwn+k-1, 'read WN at lev=',k, dfld(1,1,iwn+k-1),&
  &  maxval(dfld(:,:,iwn+k-1)), minval(dfld(:,:,iwn+k-1))
    end do
  else
    ! hydro p
    do k=1,levs
      do j=1,jgrd1
        do i=1,igrd1
          if(convert_) then
            dfld(i,j,ipn+k-1) = dfld(i,j,ips) * sl(k)
          else
            dfld(i,j,ipn+k-1) = dfld(i,j,ips) + log(sl(k))
          end if
        end do
      end do
      if(verbose) print *, 'calc P at lev=',k, maxval(dfld(:,:,ipn+k-1)), minval(dfld(:,:,ipn+k-1))
    end do
    ! hydro t
    do k=1,levs
      do j=1,jgrd1
        do i=1,igrd1
          dfld(i,j,itn+k-1) = dfld(i,j,it+k-1)
        end do
      end do
      if(verbose) print *, 'calc T at lev=',k, maxval(dfld(:,:,itn+k-1)), minval(dfld(:,:,itn+k-1))
    end do  
    ! w=0
    do k=1,levs+1
      do j=1,jgrd1
        do i=1,igrd1
          dfld(i,j,iwn+k-1) = 0.0
        end do
      end do
      if(verbose) print *, 'zero W at lev=',k, maxval(dfld(:,:,iwn+k-1)), minval(dfld(:,:,iwn+k-1))
    end do
  end if
! map factor**2
  im2=iwn+levs+1
  read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      mapf(i,j,1) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *,'read XM2 ', mapf(1,1,1), &
&  maxval(mapf(:,:,1)), minval(mapf(:,:,1))
  if(convert_) then
! modify variables by map factor
 ! U,V
  do k=1,levs
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,iu+k-1)=dfld(i,j,iu+k-1)*sqrt(mapf(i,j,1))
        dfld(i,j,iv+k-1)=dfld(i,j,iv+k-1)*sqrt(mapf(i,j,1))
      end do
    end do
  end do
  end if
! fm2x, fm2y
  read(iunit) (sfld(i),i=1,nwf) ! fm2x
  do j=1,jgrd1
    do i=1,igrd1
      mapf(i,j,2) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *, 'read FM2X ', mapf(1,1,2), &
&  maxval(mapf(:,:,2)), minval(mapf(:,:,2))
  read(iunit) (sfld(i),i=1,nwf) ! fm2y
  do j=1,jgrd1
    do i=1,igrd1
      mapf(i,j,3) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *, 'read FM2Y ', mapf(1,1,3), &
&  maxval(mapf(:,:,3)), minval(mapf(:,:,3))
! latitude and longitude
  read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
    clat(j) = real(sfld(1+(j-1)*igrd1),kind=dp)*rad2deg
    end do
  end do
  !if(verbose) print *, 'latitude ', clat(1), clat(jgrd1)
  print *, 'latitude ', clat(1), clat(jgrd1)
  read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      clon(i) = real(sfld(i+(j-1)*igrd1),kind=dp)*rad2deg
    end do
  end do
  !if(verbose) print *, 'longitude ', clon(1), clon(igrd1)
  print *, 'longitude ', clon(1), clon(igrd1)
! (icld==1 & fhour > 0) 3D physics (f_ice f_rain f_rimef)
  if((icld.eq.1).and.(fhour > 0.0)) then
  iphys3d(1)=iwn+levs+1
  do m=1,3
  do k=1,levs
    read(iunit) (sfld(i),i=1,nwf)
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,iphys3d(m)+k-1)=real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do
  if(verbose) print *,iphys3d(m)+k-1, 'read phys3d at lev=',k, &
&  dfld(1,1,iphys3d(m)+k-1),&
&  maxval(dfld(:,:,iphys3d(m)+k-1)), minval(dfld(:,:,iphys3d(m)+k-1))
  end do
  if(m<3) then
    iphys3d(m+1)=iphys3d(m)+levs
  end if
  end do
  end if
  return
end subroutine
!======================================================================
! read surface file
!======================================================================
subroutine read_sfc(iunit,igrd1,jgrd1,dfld)
! input :
! iunit         surface file (sequential, with header)
! igrd1, jgrd1  grid numbers
! nflds         field numbers
! 
! output:
! dflds         variables
!
  implicit none
  integer, intent(in) :: iunit
  integer,intent(in) :: igrd1, jgrd1
  real(kind=dp), intent(out) :: dfld(igrd1,jgrd1,nfldsfc)
  ! header
  character(len=8) :: label(4)
  integer :: idate(4), iymdh
  real(kind=sp) :: fhour
  integer :: grdtmp(2), version
  !
  integer :: iret
  integer :: nwf, irec
  integer :: i,j,k,l, ifld
  integer :: itsea, ismc, isheleg, istc, itg3, izorl, icv, icvb, &
&            icvt, islmsk, ivfrac, if10m, icanopy, ivtype, istype, iuustar, iffmm, &
&            iffhh, ialvsf, ialvwf, ialnsf, ialnwf, ifacsf, ifacwf
  real(kind=sp), allocatable :: sfld(:), sfldl(:)
  real(kind=sp), allocatable :: tmps2(:,:), tmps4(:,:)

  print *, 'read and extract header record'
  iret = 0
  rewind(iunit)
! read label
  read(iunit) label
  print '(4a8)', label
! read header
  read(iunit) fhour, idate, grdtmp, version
  iymdh = idate(4)*1000000+idate(2)*10000+idate(3)*100+idate(1)
  print *, 'posting date ', iymdh, '+', nint(fhour)
  print *, 'igrd1, jgrd1'
  print '(2i4)', grdtmp(1), grdtmp(2)
  print *, 'version', version
  nwf = igrd1*jgrd1
  print *, 'nwf', nwf
  allocate( sfld(nwf), sfldl(nwf*lsoil) )
  allocate( tmps2(nwf,2), tmps4(nwf,4) )
  print *, 'start reading surface data'
  ifld = 1
! tsea
  itsea = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read tsea ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! smc
  ismc = ifld
  read(iunit) (sfldl(i), i=1,nwf*lsoil)
  l=1
  do k=1,lsoil
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,ifld) = real(sfldl(l),kind=dp)
        l=l+1
      end do
    end do
    if(verbose) print *,ifld, 'read smc ', 'soil_l=',k, dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
    ifld=ifld+1
  end do 
! sheleg
  isheleg = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read sheleg ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! stc
  istc = ifld
  read(iunit) (sfldl(i), i=1,nwf*lsoil)
  l=1
  do k=1,lsoil
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,ifld) = real(sfldl(l),kind=dp)
        l=l+1
      end do
    end do
    if(verbose) print *,ifld, 'read stc ', 'soil_l=',k, dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
    ifld=ifld+1
  end do 
! tg3
  itg3 = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read tg3 ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! zorl
  izorl = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read zorl ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! cv
  icv = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read cv ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! cvb
  icvb = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read cvb ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! cvt
  icvt = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read cvt ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! albedo
  read(iunit) tmps4
  k=1
  ialvsf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(tmps4(l,k),kind=dp)
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'read alvsf ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
  k=k+1
  ialvwf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(tmps4(l,k),kind=dp)
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'read alvwf ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
  k=k+1
  ialnsf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(tmps4(l,k),kind=dp)
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'read alnsf ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
  k=k+1
  ialnwf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(tmps4(l,k),kind=dp)
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'read alnwf ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! slmsk
  islmsk = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read slmsk ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! vfrac
  ivfrac = ifld
  read(iunit) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read vfrac ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! canopy
  icanopy = ifld
  read(iunit,err=5000) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read canopy ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! f10m
  if10m = ifld
  read(iunit,err=5000) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read f10m ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! vtype
  ivtype = ifld
  read(iunit,err=5000) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read vtype ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! stype
  istype = ifld
  read(iunit,err=5000) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read stype ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
! facswf
  read(iunit,err=5000) tmps2
  k=1
  ifacsf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(tmps2(l,k),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read facsf ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld=ifld+1
  k=k+1
  ifacwf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(tmps2(l,k),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read facwf ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
  ifld = ifld + 1
! uustar
  iuustar = ifld
  read(iunit,end=200) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read uustar ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
200 continue
  ifld=ifld+1
! ffmm
  iffmm = ifld
  read(iunit,end=201) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read ffmm ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
201 continue
  ifld=ifld+1
! ffhh
  iffhh = ifld
  read(iunit,end=202) (sfld(i), i=1,nwf)
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,ifld) = real(sfld(l),kind=dp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'read ffhh ', dfld(1,1,ifld), maxval(dfld(:,:,ifld)), minval(dfld(:,:,ifld))
202 continue
! end reading
  print *, 'ifld', ifld, 'nflds ', nfldsfc
  return
5000 print *, ' error reading'
  stop 99
end subroutine read_sfc
!======================================================================
! read flux file
!======================================================================
subroutine read_flx(iunit,igrd1,jgrd1,dfld,ids,iparam,fhour,zhour&
                ,ind_t2m,ind_q2m,ind_u10m,ind_v10m)
!
! read flux file
!
! input :
! iunit  flux file (sequential, with header)
! 
! output:
! dfld   variables
! ids    GRIB scaling
! iparam    GRIB index
!
  implicit none
  integer, intent(in) :: iunit
  integer, intent(in) :: igrd1, jgrd1
  real(kind=dp), intent(out) :: dfld(igrd1,jgrd1,nfldflx)
  integer, intent(out)      :: ids(255)
  integer, intent(out)      :: iparam(nfldflx)
  real(kind=sp), intent(out) :: fhour, zhour
  integer, intent(out), optional :: ind_t2m,ind_q2m,ind_u10m,ind_v10m
  logical :: lt2m, lq2m, lu10m, lv10m
  integer, parameter :: iprs=1, itemp=11, iznlw=33, imerw=34, isphum=51, ipwat=54, &
  &                     ipcpr=59, isnowd=65, icldf=71, iccldf=72, islmsk=81, izorl=83, &
  &                     ialbdo=84, isoilm=144, icemsk=91, ilhflx=121, ishflx=122, izws=124, &
  &                     imws=125, ighflx=155, iuswfc=160, idswfc=161, iulwfc=162, idlwfc=163, &
  &                     inswfc=164, inlwfc=165, idswvb=166, idswvd=167, idswnb=168, idswnd=169, &
  &                     itmx=15, itmn=16, irnof=90, iep=145, iqmx=118, iqmn=119, icldwk=14, &
  &                     izgw=147, imgw=148, ihpbl=221, idswf=204, idlwf=205, iuswf=211, iulwf=212, &
  &                     icpcpr=214
  integer, parameter :: isfc=1, itoa=8, ielev=105, isglev=107, idbls=111, i2dbls=112, icolmn=200, &
  &                     ilcbl=212, ilctl=213, ilclyr=214, imcbl=222, imctl=223, imclyr=224, &
  &                     ihcbl=232, ihctl=233, ihclyr=234
  integer, parameter :: nfld=16
  integer, parameter :: lflux=27
  integer, dimension(nfld) :: ipur, itlr
  data ipur/iulwf, iuswf, iuswf, idswf, icldf, iprs, iprs, itemp, icldf, iprs, iprs, itemp,&
          & icldf, iprs, iprs, itemp/
  data itlr/itoa, itoa, isfc, isfc, ihclyr, ihctl, ihcbl, ihctl, imclyr, imctl, imcbl, imctl,&
          & ilclyr, ilctl, ilcbl, ilctl/ 
  character(len=8) :: label(4)
  integer :: idate(4), iymdh
  integer :: iyr, imo, ida, ihr
  integer :: iens(5), idstmp
  integer :: iret
  logical, allocatable :: lbm(:)
  real, allocatable :: fluxf(:,:)
  ! components of flux file
  integer :: maxbit, iptv, icen, igen, ibm, il1k, il2k, ip1, ip2, ina, inm, &
  &          icen2, inst, iavg, iacc, ifhour, ifday, ifhr, ithr, ilpds, idrt,&
  &          igrd2, jgrd2
  real(kind=sp) :: colat, rlat1, rlon1, rlat2, rlon2, delx, dely, ortru, proj
  integer :: nflds, nwf, irec
  integer :: itype, ilev, itime
  integer :: i,j,k,l,k4
  real(kind=sp), allocatable :: sfld(:)
  real(kind=sp), parameter :: rd=2.8705e2, rv=4.6150e2, fvirt=rv/rd-1.0
  real(kind=sp), parameter :: pi=3.141592, rad2deg=180.0/pi

  lt2m=.false.
  lq2m=.false.
  lu10m=.false.
  lv10m=.false.
  if(present(ind_t2m)) then
    lt2m=.true.
  end if
  if(present(ind_q2m)) then
    lq2m=.true.
  end if
  if(present(ind_u10m)) then
    lu10m=.true.
  end if
  if(present(ind_v10m)) then
    lv10m=.true.
  end if
  if(verbose) print *, lt2m, lq2m, lu10m, lv10m
  
  nwf = igrd1*jgrd1
  print *, 'nwf', nwf
  allocate( sfld(nwf) )
  allocate( lbm(nwf) )
  print *, 'start reading flux data'
  rewind(iunit)
  ! dusfc
  l=1
  if(verbose) print *, 'itype ',izws,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  fhour=real(ithr)
  zhour=real(ifhr)
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read dusfc ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! dvsfc
  l=l+1
  if(verbose) print *, 'itype ',imws,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read dvsfc ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! dtsfc
  l=l+1
  if(verbose) print *, 'itype ',ishflx,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read dtsfc ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! dqsfc
  l=l+1
  if(verbose) print *, 'itype ',ilhflx,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read dqsfc ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! tsea
  l=l+1
  if(verbose) print *, 'itype ',itemp,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read tsea ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! smc(1)
  l=l+1
  if(verbose) print *, 'itype ',isoilm,' ilev ',i2dbls
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read smc(:,1) ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! smc(2)
  l=l+1
  if(verbose) print *, 'itype ',isoilm,' ilev ',i2dbls
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read smc(:,2) ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! stc(1)
  l=l+1
  if(verbose) print *, 'itype ',itemp,' ilev ',i2dbls
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read stc(:,1) ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! stc(2)
  l=l+1
  if(verbose) print *, 'itype ',itemp,' ilev ',i2dbls
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read stc(:,2) ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! sheleg
  l=l+1
  if(verbose) print *, 'itype ',isnowd,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read sheleg ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! dlwsfc
  l=l+1
  if(verbose) print *, 'itype ',idlwf,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read dlwsfc ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! ulwsfc
  l=l+1
  if(verbose) print *, 'itype ',iulwf,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read ulwsfc ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! raw fluxes
  do k=1,4
    l=l+1
    if(verbose) print *, 'itype ',ipur(k),' ilev ',itlr(k)
    read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    ids(itype)=idstmp
    iparam(l)=itype
    if(verbose) then
    print *,'lbm ',lbm(1:min(10,nwf))
    print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
    print *,'maxbit ',maxbit,' colat ',colat
    print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
    print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
    print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
    print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
    print *,'delx ',delx,' dely ',dely
    print *,'ortru ',ortru,' proj ',proj
    print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
    print *,'il1k ',il1k,' il2k ',il2k
    print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
    print *,'ids ',ids(itype),' iens ',iens
    end if
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
      end do
    end do 
    if(verbose) print *,l,'read raw flux ',k, dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  end do
  ! fixed fluxes for approx diurnal cycle
  do k=5,7
    l=l+1
    k4=4+(k-5)*4
    if(verbose) print *, 'itype ',ipur(k4+1),' ilev ',itlr(k4+1)
    read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    ids(itype)=idstmp
    iparam(l)=itype
    if(verbose) then
    print *,'lbm ',lbm(1:min(10,nwf))
    print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
    print *,'maxbit ',maxbit,' colat ',colat
    print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
    print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
    print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
    print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
    print *,'delx ',delx,' dely ',dely
    print *,'ortru ',ortru,' proj ',proj
    print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
    print *,'il1k ',il1k,' il2k ',il2k
    print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
    print *,'ids ',ids(itype),' iens ',iens
    end if
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
      end do
    end do 
    if(verbose) print *,l,'read fixed flux ',k4+1, dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
    l=l+1
    if(verbose) print *, 'itype ',ipur(k4+2),' ilev ',itlr(k4+2)
    read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    ids(itype)=idstmp
    iparam(l)=itype
    if(verbose) then
    print *,'lbm ',lbm(1:min(10,nwf))
    print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
    print *,'maxbit ',maxbit,' colat ',colat
    print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
    print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
    print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
    print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
    print *,'delx ',delx,' dely ',dely
    print *,'ortru ',ortru,' proj ',proj
    print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
    print *,'il1k ',il1k,' il2k ',il2k
    print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
    print *,'ids ',ids(itype),' iens ',iens
    end if
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
      end do
    end do 
    if(verbose) print *,l,'read fixed flux ',k4+2, dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
    l=l+1
    if(verbose) print *, 'itype ',ipur(k4+3),' ilev ',itlr(k4+3)
    read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    ids(itype)=idstmp
    iparam(l)=itype
    if(verbose) then
    print *,'lbm ',lbm(1:min(10,nwf))
    print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
    print *,'maxbit ',maxbit,' colat ',colat
    print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
    print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
    print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
    print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
    print *,'delx ',delx,' dely ',dely
    print *,'ortru ',ortru,' proj ',proj
    print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
    print *,'il1k ',il1k,' il2k ',il2k
    print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
    print *,'ids ',ids(itype),' iens ',iens
    end if
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
      end do
    end do 
    if(verbose) print *,l,'read fixed flux ',k4+3, dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
    l=l+1
    if(verbose) print *, 'itype ',ipur(k4+4),' ilev ',itlr(k4+4)
    read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    ids(itype)=idstmp
    iparam(l)=itype
    if(verbose) then
    print *,'lbm ',lbm(1:min(10,nwf))
    print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
    print *,'maxbit ',maxbit,' colat ',colat
    print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
    print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
    print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
    print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
    print *,'delx ',delx,' dely ',dely
    print *,'ortru ',ortru,' proj ',proj
    print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
    print *,'il1k ',il1k,' il2k ',il2k
    print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
    print *,'ids ',ids(itype),' iens ',iens
    end if
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
      end do
    end do 
    if(verbose) print *,l,'read fixed flux ',k4+4, dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  end do 
  ! geshem
  l=l+1
  if(verbose) print *, 'itype ',ipcpr,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read geshem ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! bengsh
  l=l+1
  if(verbose) print *, 'itype ',icpcpr,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read bengsh ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! gflux
  l=l+1
  if(verbose) print *, 'itype ',ighflx,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read gflux ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! slmsk
  l=l+1
  if(verbose) print *, 'itype ',islmsk,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read slmsk ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! cemsk
  l=l+1
  if(verbose) print *, 'itype ',icemsk,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read cemsk ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! u10
  l=l+1
  if(lu10m) then
    ind_u10m=l
  end if
  if(verbose) print *, 'itype ',iznlw,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read u10 ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! v10
  l=l+1
  if(lv10m) then
    ind_v10m=l
  end if
  if(verbose) print *, 'itype ',imerw,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read v10 ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! t2
  l=l+1
  if(lt2m) then
    ind_t2m=l
  end if
  if(verbose) print *, 'itype ',itemp,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read t2 ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! q2
  l=l+1
  if(lq2m) then
    ind_q2m=l
  end if
  if(verbose) print *, 'itype ',isphum,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read q2 ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! psurf
  l=l+1
  if(verbose) print *, 'itype ',iprs,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read psurf ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! t2max
  l=l+1
  if(verbose) print *, 'itype ',itmx,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read t2max ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! q2max
  l=l+1
  if(verbose) print *, 'itype ',iqmx,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read q2max ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! t2min
  l=l+1
  if(verbose) print *, 'itype ',itmn,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read t2min ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! q2min
  l=l+1
  if(verbose) print *, 'itype ',iqmn,' ilev ',ielev
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read q2min ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! runoff
  l=l+1
  if(verbose) print *, 'itype ',irnof,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read runoff ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! ep
  l=l+1
  if(verbose) print *, 'itype ',iep,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read ep ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! cldwrk
  l=l+1
  if(verbose) print *, 'itype ',icldwk,' ilev ',icolmn
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read cldwrk ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! dugwd
  l=l+1
  if(verbose) print *, 'itype ',izgw,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read dugwd ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! dvgwd
  l=l+1
  if(verbose) print *, 'itype ',imgw,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read dvgwd ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! hpbl
  l=l+1
  if(verbose) print *, 'itype ',ihpbl,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read hpbl ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! pwat
  l=l+1
  if(verbose) print *, 'itype ',ipwat,' ilev ',icolmn
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read pwat ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! albedo(percent)
  l=l+1
  if(verbose) print *, 'itype ',ialbdo,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read albedo ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! cldf
  l=l+1
  if(verbose) print *, 'itype ',icldf,' ilev ',icolmn
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read cldf ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! wvuflx
  l=l+1
  if(verbose) print *, 'itype ',242,' ilev ',icolmn
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
!  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',idstmp,' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read wvuflx ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! wvvflx
  l=l+1
  if(verbose) print *, 'itype ',243,' ilev ',icolmn
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
!  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose)then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',idstmp,' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read wvvflx ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! srunoff
  l=l+1
  if(verbose) print *, 'itype ',235,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
!  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',idstmp,' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read srunoff ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! soilm
  l=l+1
  if(verbose) print *, 'itype ',86,' ilev ',i2dbls
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read soilm ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  ! snwdph
  l=l+1
  if(verbose) print *, 'itype ',66,' ilev ',isfc
  read(iunit) sfld,lbm,idrt,igrd2,jgrd2,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  ids(itype)=idstmp
  iparam(l)=itype
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd2,' jgrd ',jgrd2
  print *,'maxbit ',maxbit,' colat ',colat
  print *,'ilpds ',ilpds,' iptv ',iptv,' icen ',icen,' igen ',igen
  print *,'iyr ',iyr,' imo ',imo,' ida ',ida,' ihr ',ihr
  print *,'ifhour ',ifhour,' ifhr ',ifhr,' ithr ',ithr
  print *,'rlat1 ',rlat1,' rlon1 ',rlon1,' rlat2 ',rlat2,' rlon2 ',rlon2
  print *,'delx ',delx,' dely ',dely
  print *,'ortru ',ortru,' proj ',proj
  print *,'ibm ',ibm,' itype ',itype,' ilev ',ilev
  print *,'il1k ',il1k,' il2k ',il2k
  print *,'itime ',itime,' ina ',ina,' inm ',inm,' icen2 ',icen2
  print *,'ids ',ids(itype),' iens ',iens
  end if
  do j=1,jgrd1
    do i=1,igrd1
      dfld(i,j,l) = real(sfld(i+(j-1)*igrd1),kind=dp)
    end do
  end do 
  if(verbose) print *,l,'read snwdph ', dfld(1,1,l), maxval(dfld(:,:,l)), minval(dfld(:,:,l))
  print *,'end reading flux file'
  return
end subroutine read_flx
!
end module read_module
