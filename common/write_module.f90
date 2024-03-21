module write_module
!
! module : write_module   writing sigma, surface, flux files
! history:
! 22-07-28 create
! 
! namelist:
!
  use kind_module
  use phconst_module, only: pi,deg2rad,fvirt
  use func_module, only: conv_temp
  use read_module, only: levmax, nwext, lsoil, nfldsfc, nfldflx, verbose
  
  implicit none
  private

  public :: write_sig, write_sfc, write_flx
contains
!======================================================================
! write sigma file
!======================================================================
subroutine write_sig(ounit,label,idate,fhour,si,sl,ext,&
&                    igrd1,jgrd1,levs,nflds,nonhyd,icld,dfld,mapf,clat,clon,&
&                    convert)
! input :
! ounit sigma file (sequential, with header)
! dfld  variables (double precision)
!
  implicit none
  integer, intent(in) :: ounit
  character(len=8),intent(in) :: label(4)
  integer,intent(in) :: idate(4)
  real(kind=sp),intent(in) :: fhour, ext(nwext) 
  real(kind=dp),intent(in) :: si(levmax+1), sl(levmax)
  integer, intent(in) :: igrd1, jgrd1, levs, nflds 
  integer, intent(in) :: nonhyd
  integer, intent(in) :: icld !1=include 3D physics, 0=not include
  real(kind=dp), intent(inout) :: dfld(igrd1,jgrd1,nflds)
  real(kind=dp), intent(in) :: mapf(igrd1,jgrd1,3)
  real(kind=dp), intent(in) :: clat(jgrd1), clon(igrd1)
  logical, intent(in), optional :: convert
  logical :: convert_
  real(kind=sp) :: sisl(2*levmax+1)
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
! write label
  write(ounit) label
! write header
  print *, 'posting date = ', idate(4), idate(2), idate(3), idate(1), '+', nint(fhour)
  sisl=0.0
  sisl(1:levs+1) = si(1:levs+1)
  sisl(levs+2:2*levs+1) = sl(1:levs)
  write(ounit) fhour, idate, sisl, ext
  if ( nonhyd.eq.1 ) then
    print *, 'Input is a nonhydrostatic'
  else
    print *, 'Input is a hydrostatic'
  end if
  nwf = igrd1*jgrd1
  print *, 'nwf', nwf
  allocate( sfld(nwf) )
  allocate( factor(igrd1,jgrd1,levs) )
  print *, 'start writing sigma data'
  l=1
  ! gz
  igz=1
  do j=1,jgrd1
    do i=1,igrd1
       sfld(i+(j-1)*igrd1) = real(dfld(i,j,igz),kind=sp)
    end do
  end do 
  write(ounit) (sfld(i),i=1,nwf)
  if(verbose) print *,igz, 'write gz ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! ln(ps)[kPa]
  ips=igz+1
  do j=1,jgrd1
    do i=1,igrd1
       sfld(i+(j-1)*igrd1) = real(dfld(i,j,ips),kind=sp)
       if(convert_) then
         sfld(i+(j-1)*igrd1) = LOG(sfld(i+(j-1)*igrd1)/1000.0) !Pa=>kPa
       end if
    end do
  end do 
  write(ounit) (sfld(i),i=1,nwf)
  if(verbose) print *,ips, 'write lnps ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! T
  it=ips+1
  iq=it+3*levs
  if(convert_) then
  ! temperature => virtual temperature
  call conv_temp(igrd1,jgrd1,levs,dfld(:,:,it:it+levs-1),dfld(:,:,iq:iq+levs-1),1)
  end if
  do k=1, levs
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,it+k-1),kind=sp)
      end do
    end do 
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,it+k-1, 'write Tv at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
  end do  
  ! U,V
  iu=it+levs
  iv=iu+levs
  do k=1, levs
    if(convert_) then
! modify variables by map factor
    do j=1,jgrd1
      do i=1,igrd1
        dfld(i,j,iu+k-1)=dfld(i,j,iu+k-1)/sqrt(mapf(i,j,1))
        dfld(i,j,iv+k-1)=dfld(i,j,iv+k-1)/sqrt(mapf(i,j,1))
      end do
    end do
    end if
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,iu+k-1),kind=sp) 
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,iu+k-1, 'write U at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,iv+k-1),kind=sp) 
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,iv+k-1, 'write V at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
  end do
  ! Q
  iq = iv+levs
  do k=1, levs
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,iq+k-1),kind=sp)
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,iq+k-1, 'write Q at lev=',k, sfld(1), &
&    maxval(sfld(:)), minval(sfld(:))
  end do
  ! OZ
  ioz=iq+levs
  do k=1,levs
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,ioz+k-1),kind=sp)
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,ioz+k-1, 'write OZ at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
  end do
  ! CW
  icw=ioz+levs
  do k=1,levs
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1)= real(dfld(i,j,icw+k-1),kind=sp)
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,icw+k-1, 'write CW at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
  end do
  ipn=icw+levs
  itn=ipn+levs
  iwn=itn+levs
  if (nonhyd.eq.1) then
  ! pn
  do k=1,levs
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,ipn+k-1),kind=sp)
        if(convert_) then
          sfld(i+(j-1)*igrd1) = LOG(sfld(i+(j-1)*igrd1)/1000.0) !Pa=>kPa
        end if
!        if (k.eq.1) then
!          dfld(i,j,ips) = dfld(i,j,ipn)/sl(1)
!        end if
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,ipn+k-1, 'write PN at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
  end do
  ! tn
  if (convert_) then
  ! temperature => virtual temperature
  call conv_temp(igrd1,jgrd1,levs,dfld(:,:,itn:itn+levs-1),dfld(:,:,iq:iq+levs-1),1)
  end if
  do k=1,levs
    do j=1,jgrd1
      do i=1,igrd1
         sfld(i+(j-1)*igrd1) = real(dfld(i,j,itn+k-1),kind=sp)
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,itn+k-1, 'write TN at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
  end do  
  ! wn
  do k=1,levs+1
    do j=1,jgrd1
      do i=1,igrd1
         sfld(i+(j-1)*igrd1) = real(dfld(i,j,iwn+k-1),kind=sp)
      end do
    end do
    write(ounit) (sfld(i),i=1,nwf)
    if(verbose) print *,iwn+k-1, 'write WN at lev=',k, sfld(1),&
&    maxval(sfld(:)), minval(sfld(:))
  end do
  end if
! map factor**2
  im2=iwn+levs+1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(mapf(i,j,1),kind=sp)
    end do
  end do
  write(ounit) (sfld(i),i=1,nwf)
  if(verbose) print *, 'write XM2 ', sfld(1), &
&  maxval(sfld(:)), minval(sfld(:))
! fm2x, fm2y
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(mapf(i,j,2),kind=sp)
    end do
  end do
  write(ounit) (sfld(i),i=1,nwf) !fm2x
  if(verbose) print *, 'write FM2X ', sfld(1), &
&  maxval(sfld(:)), minval(sfld(:))
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(mapf(i,j,3),kind=sp)
    end do
  end do
  write(ounit) (sfld(i),i=1,nwf) !fm2y
  if(verbose) print *, 'write FM2Y ', sfld(1), &
&  maxval(sfld(:)), minval(sfld(:))
! latitude and longitude
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(clat(j)*deg2rad,kind=sp)
    end do
  end do
  write(ounit) (sfld(i),i=1,nwf)
  if(verbose) print *, 'latitude ', sfld(1), sfld(nwf)
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(clon(i)*deg2rad,kind=sp)
    end do
  end do
  write(ounit) (sfld(i),i=1,nwf)
  if(verbose) print *, 'longitude ', sfld(1), sfld(nwf)
! (icld==1 & fhour > 0) 3D physics (f_ice f_rain f_rimef)
  if((icld.eq.1).and.(fhour > 0.0)) then
  iphys3d(1)=iwn+levs+1
  iphys3d(2)=iphys3d(1)+levs
  iphys3d(3)=iphys3d(3)+levs
  do m=1,3
    do k=1,levs
      do j=1,jgrd1
        do i=1,igrd1
          sfld(i+(j-1)*igrd1)=real(dfld(i,j,iphys3d(m)+k-1),kind=sp)
        end do
      end do
      write(ounit) (sfld(i),i=1,nwf)
      if(verbose) print *,iphys3d(m)+k-1, 'write phys3d at lev=',k, &
&      sfld(1),&
&      maxval(sfld(:)), minval(sfld(:))
    end do
  end do
  end if
  return
end subroutine
!======================================================================
! write surface file
!======================================================================
subroutine write_sfc(ounit,igrd1,jgrd1,dfld,label,idate,fhour)
! input :
! ounit         surface file (sequential, with header)
! igrd1, jgrd1  grid numbers
! dflds         variables
!
  implicit none
  integer, intent(in) :: ounit
  integer,intent(in) :: igrd1, jgrd1
  real(kind=dp), intent(in) :: dfld(igrd1,jgrd1,nfldsfc)
  ! header
  character(len=8),intent(in) :: label(4)
  integer,intent(in) :: idate(4)
  real(kind=sp),intent(in) :: fhour
  integer :: version
  integer :: iymdh
  !
  integer :: iret
  integer :: nwf, irec
  integer :: i,j,k,l, ifld
  integer :: itsea, ismc, isheleg, istc, itg3, izorl, icv, icvb, &
&            icvt, islmsk, ivfrac, if10m, icanopy, ivtype, istype, iuustar, iffmm, &
&            iffhh, ialvsf, ialvwf, ialnsf, ialnwf, ifacsf, ifacwf
  real(kind=sp), allocatable :: sfld(:), sfldl(:)
  real(kind=sp), allocatable :: tmps2(:,:), tmps4(:,:)
  
  print *, 'write header record'
  iret = 0
! write label
  write(ounit) label
! write header
  print *, 'posting date = ', idate(4), idate(2), idate(3), idate(1), '+', nint(fhour)
  version=199802
  write(ounit) fhour, idate, igrd1, jgrd1, version
  nwf = igrd1*jgrd1
  print *, 'nwf', nwf
  allocate( sfld(nwf), sfldl(nwf*lsoil) )
  allocate( tmps2(nwf,2), tmps4(nwf,4) )
  print *, 'start writing surface data'
  ifld = 1
! tsea
  itsea = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write tsea ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! smc
  ismc = ifld
  l=1
  do k=1,lsoil
    do j=1,jgrd1
      do i=1,igrd1
        sfldl(l) = real(dfld(i,j,ifld),kind=sp)
        l=l+1
      end do
    end do
    if(verbose) print *,ifld, 'write smc ', 'soil_l=',k, sfldl(1+(k-1)*nwf), &
    &  maxval(sfldl(1+(k-1)*nwf:nwf*k)), minval(sfldl(1+(k-1)*nwf:nwf*k))
    ifld=ifld+1
  end do 
  write(ounit) (sfldl(i), i=1,nwf*lsoil)
! sheleg
  isheleg = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp) 
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write sheleg ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! stc
  istc = ifld
  l=1
  do k=1,lsoil
    do j=1,jgrd1
      do i=1,igrd1
        sfldl(l) = real(dfld(i,j,ifld),kind=sp)
        l=l+1
      end do
    end do
    if(verbose) print *,ifld, 'write stc ', 'soil_l=',k, sfldl(1+(k-1)*nwf), &
    &  maxval(sfldl(1+(k-1)*nwf:nwf*k)), minval(sfldl(1+(k-1)*nwf:nwf*k))
    ifld=ifld+1
  end do 
  write(ounit) (sfldl(i), i=1,nwf*lsoil)
! tg3
  itg3 = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp) 
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write tg3 ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! zorl
  izorl = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp) 
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write zorl ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! cv
  icv = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp) 
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write cv ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! cvb
  icvb = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write cvb ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! cvt
  icvt = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp) 
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write cvt ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! albedo
  k=1
  ialvsf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      tmps4(l,k) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'write alvsf ', tmps4(1,k), maxval(tmps4(:,k)), minval(tmps4(:,k))
  ifld=ifld+1
  k=k+1
  ialvwf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      tmps4(l,k) = real(dfld(i,j,ifld),kind=sp) 
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'write alvwf ', tmps4(1,k), maxval(tmps4(:,k)), minval(tmps4(:,k))
  ifld=ifld+1
  k=k+1
  ialnsf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      tmps4(l,k) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'write alnsf ', tmps4(1,k), maxval(tmps4(:,k)), minval(tmps4(:,k))
  ifld=ifld+1
  k=k+1
  ialnwf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      tmps4(l,k) = real(dfld(i,j,ifld),kind=sp) 
      l=l+1
    end do
  end do
  if(verbose) print *,ifld, 'write alnwf ', tmps4(1,k), maxval(tmps4(:,k)), minval(tmps4(:,k))
  write(ounit) tmps4
  ifld=ifld+1
! slmsk
  islmsk = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write slmsk ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! vfrac
  ivfrac = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write vfrac ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! canopy
  icanopy = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write canopy ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! f10m
  if10m = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write f10m ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! vtype
  ivtype = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write vtype ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! stype
  istype = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write stype ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! facswf
  k=1
  ifacsf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      tmps2(l,k) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'write facsf ', tmps2(1,k), maxval(tmps2(:,k)), minval(tmps2(:,k))
  ifld=ifld+1
  k=k+1
  ifacwf = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      tmps2(l,k) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  if(verbose) print *,ifld, 'write facwf ', tmps2(1,k), maxval(tmps2(:,k)), minval(tmps2(:,k))
  write(ounit) tmps2
  ifld = ifld + 1
! uustar
  iuustar = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write uustar ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! ffmm
  iffmm = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write ffmm ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ifld=ifld+1
! ffhh
  iffhh = ifld
  l=1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(l) = real(dfld(i,j,ifld),kind=sp)
      l=l+1
    end do
  end do 
  write(ounit) (sfld(i), i=1,nwf)
  if(verbose) print *,ifld, 'write ffhh ', sfld(1), maxval(sfld(:)), minval(sfld(:))
! end writing
  print *, 'ifld', ifld, 'nflds ', nfldsfc
  return
end subroutine write_sfc
!======================================================================
! write flux file
!======================================================================
subroutine write_flx(ounit,igrd1,jgrd1,dfld,ids,iparam,slmsk,&
                   & idate,fhour,zhour,&
                   & rdelx,rdely,rlat,rlon,rtruth,rorient,rproj)
!
! write flux file
!
! input :
! ounit  flux file (sequential, with header)
! dfld   variables
! ids    GRIB scaling
! iparam    GRIB index
!
  implicit none
  integer, intent(in) :: ounit
  integer, intent(in) :: igrd1, jgrd1
  real(kind=dp), intent(inout) :: dfld(igrd1,jgrd1,nfldflx)
  integer, intent(inout)      :: ids(255)
  integer, intent(inout)      :: iparam(nfldflx)
  real(kind=dp), intent(in)    :: slmsk(igrd1,jgrd1) ! read from surface file
  integer, intent(in) :: idate(4)
  real(kind=sp), intent(in) :: fhour, zhour
  real(kind=dp), intent(in) :: rdelx, rdely, rtruth, rorient, rproj
  real(kind=dp), intent(in) :: rlat(jgrd1), rlon(igrd1)
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
  integer :: iymdh
  integer :: iyr, imo, ida, ihr
  integer :: iens(5), idstmp
  integer :: iret
  logical, allocatable :: lbm(:)
  real, allocatable :: fluxf(:,:)
  real, allocatable :: slmask(:)
  ! components of flux file
  integer :: maxbit, iptv, icen, igen, ibm, ibm0, ibm1, il1k, il2k, ip1, ip2, ina, inm, &
  &          icen2, inst, iavg, iacc, ifhour, ifday, ifhr, ithr, ilpds, idrt
  real(kind=sp) :: colat, rlat1, rlon1, rlat2, rlon2, delx, dely, ortru, proj
  integer :: nflds, nwf, irec
  integer :: itype, ilev, itime
  integer :: i,j,k,l,k4
  real(kind=sp), allocatable :: sfld(:)
  
  nwf = igrd1*jgrd1
  print *, 'nwf', nwf
  allocate( sfld(nwf) )
  allocate( lbm(nwf) )
  allocate( fluxf(nwf,4) )
  allocate( slmask(nwf) )
! constants
  maxbit=16
  colat=0.
  iptv=2
  icen=7
  igen=99
  ibm0=0
  ibm1=1
  il1k=0
  il2k=0
  ip1=0
  ip2=0
  ina=0
  inm=0
  icen2=0
  inst=10
  iavg=3
  iacc=4
  ifhour=1
  ifday=2
!
  ilpds=28
  iens(1)=1; iens(2)=0; iens(3)=0; iens(4)=1; iens(5)=255
  iyr=idate(4); imo=idate(2); ida=idate(3); ihr=idate(1)
  ifhr=nint(zhour); ithr=nint(fhour)
  delx=rdelx; dely=rdely
  rlat1=rlat(1)*deg2rad; rlat2=rlat(1)*deg2rad
  rlon1=rlon(1)*deg2rad; rlon2=rlon(igrd1)*deg2rad
  if( rproj .eq. 0.0 ) then !mercater
    idrt=1
    ortru=rtruth
    proj=0.0
  else ! polar
    idrt=5
    ortru=rorient
    proj=rproj
  end if
  do j=1,nwf
    lbm(j)=.TRUE.
  end do
  l=1
  do j=1,jgrd1
    do i=1, igrd1
      slmask(l) = real(slmsk(i,j),kind=sp)
      l=l+1
    end do
  end do
!  
  print *, 'start writing flux data'
  ! dusfc
  l=1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) then
  print *,'lbm ',lbm(1:min(10,nwf))
  print *,'idrt ',idrt,' igrd ',igrd1,' jgrd ',jgrd1
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
  print *,l,'write dusfc ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  end if
  ! dvsfc
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write dvsfc ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! dtsfc
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write dtsfc ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! dqsfc
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write dqsfc ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! tsea
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write tsea ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! smc(1)
  l=l+1
  itype=iparam(l)
  ilev=i2dbls
  idstmp=ids(itype)
  itime=inst
  ibm=ibm1
  il2k=10
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write smc(:,1) ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! smc(2)
  l=l+1
  itype=iparam(l)
  ilev=i2dbls
  idstmp=ids(itype)
  itime=inst
  ibm=ibm1
  il1k=10
  il2k=200
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write smc(:,2) ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! stc(1)
  l=l+1
  itype=iparam(l)
  ilev=i2dbls
  idstmp=ids(itype)
  itime=inst
  ibm=ibm1
  il1k=0
  il2k=10
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write stc(:,1) ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! stc(2)
  l=l+1
  itype=iparam(l)
  ilev=i2dbls
  idstmp=ids(itype)
  itime=inst
  ibm=ibm1
  il1k=10
  il2k=200
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write stc(:,2) ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! sheleg
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il1k=0
  il2k=0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write sheleg ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! dlwsfc
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write dlwsfc ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! ulwsfc
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write ulwsfc ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! raw fluxes
  do k=1,4
    l=l+1
  itype=iparam(l)
  ilev=itlr(k)
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
    write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    if(verbose) print *,l,'write raw flux ',k, sfld(1), maxval(sfld(:)), minval(sfld(:))
  end do
  ! fixed fluxes for approx diurnal cycle
  do k=5,7
    l=l+1
    k4=4+(k-5)*4
    itype=iparam(l)
    ilev=itlr(k4+1)
    idstmp=ids(itype)
    itime=iavg
    ibm=ibm0
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
      end do
    end do 
    do j=1,nwf
      lbm(j)=sfld(j).gt.0.5
    end do
    write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    if(verbose) print *,l,'write fixed flux ',k4+1, sfld(1), maxval(sfld(:)), minval(sfld(:))
    l=l+1
    itype=iparam(l)
    ilev=itlr(k4+2)
    idstmp=ids(itype)
    itime=iavg
    ibm=ibm1
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
      end do
    end do 
    write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    if(verbose) print *,l,'write fixed flux ',k4+2, sfld(1), maxval(sfld(:)), minval(sfld(:))
    l=l+1
    itype=iparam(l)
    ilev=itlr(k4+3)
    idstmp=ids(itype)
    itime=iavg
    ibm=ibm1
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
      end do
    end do 
    write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    if(verbose) print *,l,'write fixed flux ',k4+3, sfld(1), maxval(sfld(:)), minval(sfld(:))
    l=l+1
    itype=iparam(l)
    ilev=itlr(k4+4)
    idstmp=ids(itype)
    itime=iavg
    ibm=ibm1
    do j=1,jgrd1
      do i=1,igrd1
        sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
      end do
    end do 
    write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
    &           ilpds,iptv,icen,igen,&
    &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
    &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
    &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
    if(verbose) print *,l,'write fixed flux ',k4+4, sfld(1), maxval(sfld(:)), minval(sfld(:))
  end do 
  ! geshem
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write geshem ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! bengsh
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write bengsh ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! gflux
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  do j=1,nwf
    lbm(j)=slmask(j).ne.0
  end do
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write gflux ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! slmsk
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write slmsk ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! cemsk
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write cemsk ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! u10
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=10
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write u10 ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! v10
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=10
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write v10 ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! t2
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=2
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write t2 ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! q2
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=2
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write q2 ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! psurf
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write psurf ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! t2max
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=2
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write t2max ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! q2max
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=2
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write q2max ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! t2min
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=2
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write t2min ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! q2min
  l=l+1
  itype=iparam(l)
  ilev=ielev
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  il2k=2
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write q2min ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! runoff
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iacc
  ibm=ibm1
  il2k=0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  do j=1,nwf
    lbm(j)=slmask(j).ne.0.
  end do
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write runoff ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! ep
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write ep ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! cldwrk
  l=l+1
  itype=iparam(l)
  ilev=icolmn
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write cldwrk ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! dugwd
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write dugwd ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! dvgwd
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write dvgwd ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! hpbl
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write hpbl ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! pwat
  l=l+1
  itype=iparam(l)
  ilev=icolmn
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write pwat ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! albedo(percent)
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write albedo ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! cldf
  l=l+1
  itype=iparam(l)
  ilev=icolmn
  idstmp=ids(itype)
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  do j=1,nwf
    lbm(j)=sfld(j).gt.0
  end do
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write cldf ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! wvuflx
  l=l+1
  itype=iparam(l)
  ilev=icolmn
  idstmp=0
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,131,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write wvuflx ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! wvvflx
  l=l+1
  itype=iparam(l)
  ilev=icolmn
  idstmp=0
  itime=iavg
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,131,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write wvvflx ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! srunoff
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(irnof)
  itime=iacc
  ibm=ibm1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  do j=1,nwf
    lbm(j)=slmask(j).ne.0
  end do
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ifhr,ithr,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write srunoff ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! soilm
  l=l+1
  itype=iparam(l)
  ilev=i2dbls
  idstmp=ids(itype)
  itime=inst
  ibm=ibm1
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  do j=1,nwf
    lbm(j)=slmask(j).ne.0
  end do
  il1k=0
  il2k=200
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,0,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write soilm ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  ! snwdph
  l=l+1
  itype=iparam(l)
  ilev=isfc
  idstmp=ids(itype)
  itime=inst
  ibm=ibm0
  do j=1,jgrd1
    do i=1,igrd1
      sfld(i+(j-1)*igrd1) = real(dfld(i,j,l),kind=sp)
    end do
  end do 
  do j=1,nwf
    lbm(j)=.true.
  end do
  il1k=0
  il2k=0
  write(ounit) sfld,lbm,idrt,igrd1,jgrd1,maxbit,colat, &
  &           ilpds,iptv,icen,igen,&
  &           ibm,itype,ilev,il1k,il2k,iyr,imo,ida,ihr,&
  &           ifhour,ithr,ip2,itime,ina,inm,icen2,idstmp,iens,&
  &           rlat1,rlon1,rlat2,rlon2,delx,dely,ortru,proj
  if(verbose) print *,l,'write snwdph ', sfld(1), maxval(sfld(:)), minval(sfld(:))
  print *,'end writing flux file'
  return
end subroutine write_flx
!
end module write_module
