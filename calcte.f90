program calcte
!
! calculate full or perturbation moist total energy
!
  use kind_module
  use phconst_module
  use read_module
  use write_module
  use norm_module, only: calc_te
  use date_module, only: ndate
  implicit none
  logical :: lprtb=.true. ! False=>calculate for full field
  real(kind=dp) :: epsq=1.0d0 ! weight for moist term
  real(kind=dp) :: lonw=-999.9d0, lone=-999.9d0 ! calculation region
  real(kind=dp) :: lats=-999.9d0, latn=-999.9d0 ! calculation region
  integer :: kmax=21 ! upper limit for energy calculation
  namelist /namlst_prtb/ lprtb, epsq, lonw, lone, lats, latn, kmax
  integer :: ilonw, ilone, jlats, jlatn ! calculation region indexes
  integer :: nlon, nlat
  ! for energy calculation
  real(kind=dp), parameter :: p0=1000.0d2, ptheta=rd/cp ! potential temperature
  real(kind=dp), allocatable :: u(:,:,:),v(:,:,:),t(:,:,:),q(:,:,:)
  real(kind=dp), allocatable :: theta(:,:,:),fact(:,:,:)
  real(kind=dp), allocatable :: ps(:,:)
  real(kind=dp) :: area,te(4),coef
  integer :: ips,it,iu,iv,iq
  character(len=6) :: ofile='te.dat'
  ! input files' units (initial, plus 1 for 12h forecast)
  integer, parameter :: nisig=11, nisfc=21, niflx=31
  integer            :: nsig,     nsfc,     nflx
  real(kind=dp), allocatable :: dfld(:,:,:)
  real(kind=dp), allocatable :: mapf(:,:,:), clat(:), clon(:), slmsk(:,:)
  character(len=8) :: label(4)
  integer :: idate(4), nfldsig
  integer :: idate2(4), iymdh, iymdh2
  ! for ndate
  integer :: date1(5),date2(5),dtmin
  real(kind=sp) :: fhour, ext(nwext) 
  real(kind=sp) :: fhour2
  real(kind=dp) :: si(levmax+1), sl(levmax)
  real(kind=dp) :: rdelx, rdely, rtruth, rorient, rproj
  integer :: ids(255), iparam(nfldflx)
  integer :: igrd1, jgrd1, levs, nonhyd, icld
  integer :: n,i,j,k

  read(5,namlst_prtb)
  write(6,namlst_prtb)
  
!!! sigma files (r_sig.fNN)
  ! headers are assumed to be identical for all initial time
  icld=1
  call read_header(nisig,icld,label,idate,fhour,si,sl,ext,nfldsig)
  print*, label
  print*, idate
  print*, fhour
  dtmin=nint(fhour)*60
  date1(1)=idate(4)
  date1(2)=idate(2)
  date1(3)=idate(3)
  date1(4)=idate(1)
  date1(5)=0
  call ndate(date1,dtmin,date2)
  idate(4)=date2(1)
  idate(2)=date2(2)
  idate(3)=date2(3)
  idate(1)=date2(4)
  iymdh = idate(4)*1000000+idate(2)*10000+idate(3)*100+idate(1)
  !print*, ext(1:16)
!  print*, nfldsig
  igrd1 = int(ext(3))
  jgrd1 = int(ext(4))
  print*, igrd1, jgrd1
  levs = int(ext(5))
  print*, levs
  rproj = ext(7)
  rtruth=ext(8); rorient=ext(9)
  rdelx=ext(14); rdely=ext(15)
  nonhyd=int(ext(16))
  print*, nonhyd
!  print*, si(1:levs+1)
!  print*, sl(1:levs)
  allocate( dfld(igrd1,jgrd1,nfldsig) )
  allocate( mapf(igrd1,jgrd1,3) )
  allocate( clat(jgrd1), clon(igrd1) )
  
  ips=2
  if(nonhyd.eq.1) then 
    it=2+7*levs
  else
    it=3
  end if
  iu=3+levs
  iv=iu+levs
  iq=iv+levs
  ! 0h forecast (analysis)
  nsig=nisig
!  call read_sig(nsig,igrd1,jgrd1,levs,nfldsig,nonhyd,fhour,sl,dfld,mapf,clat,clon)
  call read_sig(nsig,igrd1,jgrd1,levs,nfldsig,nonhyd,icld,0.0,sl,dfld,mapf,clat,clon)
  !! setting boundaries
  if ((lonw.gt.-999.9d0).and.(lone.gt.-999.9d0)) then
    do i=1,igrd1
      if(clon(i).ge.lonw) then
        ilonw=i
        exit
      end if
    end do 
    do i=1,igrd1
      if(clon(i).ge.lone) then
        ilone=i
        exit
      end if
    end do 
  else
    ilonw=1
    ilone=igrd1
  end if
  if ((lats.gt.-999.9d0).and.(latn.gt.-999.9d0)) then
    do j=1,jgrd1
      if(clat(j).ge.lats) then
        jlats=j
        exit
      end if
    end do 
    do j=1,jgrd1
      if(clat(j).ge.latn) then
        jlatn=j
        exit
      end if
    end do 
  else
    jlats=1
    jlatn=jgrd1
  end if
  print *, 'boundary ',ilonw,'-',ilone,' lon ',clon(ilonw),'-',clon(ilone)
  print *, 'boundary ',jlats,'-',jlatn,' lat ',clat(jlats),'-',clat(jlatn)
  nlon = ilone - ilonw + 1
  nlat = jlatn - jlats + 1
  print *, 'nlon ',nlon,' nlat ',nlat
  allocate( u(nlon,nlat,kmax),v(nlon,nlat,kmax) )
  allocate( t(nlon,nlat,kmax),q(nlon,nlat,kmax) )
  allocate( theta(nlon,nlat,kmax),fact(nlon,nlat,kmax) )
  allocate( ps(nlon,nlat) )
  u=0.0d0
  v=0.0d0
  t=0.0d0
  q=0.0d0
  ps=0.0d0

  do j=1,nlat
    do i=1,nlon
      ps(i,j)=dfld(i+ilonw-1,j+jlats-1,ips)
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        t(i,j,k) = dfld(i+ilonw-1,j+jlats-1,it+k-1)
        theta(i,j,k) = t(i,j,k)*(p0/dfld(i+ilonw-1,j+jlats-1,ips)/sl(k))**ptheta
      end do
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        u(i,j,k) = dfld(i+ilonw-1,j+jlats-1,iu+k-1)
      end do
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        v(i,j,k) = dfld(i+ilonw-1,j+jlats-1,iv+k-1)
      end do
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        q(i,j,k) = dfld(i+ilonw-1,j+jlats-1,iq+k-1)
      end do
    end do
  end do
  do k=1,kmax
    print*, k
    print *, 'u(full,max)', maxval(u(:,:,k)),maxloc(u(:,:,k))
    print *, 'u(full,min)', minval(u(:,:,k)),minloc(u(:,:,k))
    print *, 'v(full,max)', maxval(v(:,:,k)),maxloc(v(:,:,k))
    print *, 'v(full,min)', minval(v(:,:,k)),minloc(v(:,:,k))
    print *, 't(full,max)', maxval(t(:,:,k)),maxloc(t(:,:,k))
    print *, 't(full,min)', minval(t(:,:,k)),minloc(t(:,:,k))
!    print *, 'th(full,max)', maxval(theta(:,:,k)),maxloc(theta(:,:,k))
!    print *, 'th(full,min)', minval(theta(:,:,k)),minloc(theta(:,:,k))
    print *, 'q(full,max)', maxval(q(:,:,k)),maxloc(q(:,:,k))
    print *, 'q(full,min)', minval(q(:,:,k)),minloc(q(:,:,k))
  end do
  print *, 'ps(full,max)', maxval(ps(:,:)),maxloc(ps(:,:))
  print *, 'ps(full,min)', minval(ps(:,:)),minloc(ps(:,:))
  if(lprtb) then
  ! 12h forecast
  nsig=nisig+1
  !! consistency check
  call read_header(nsig,icld,label,idate2,fhour2,si,sl,ext,nfldsig)
  dtmin=nint(fhour2)*60
  date1(1)=idate2(4)
  date1(2)=idate2(2)
  date1(3)=idate2(3)
  date1(4)=idate2(1)
  date1(5)=0
  call ndate(date1,dtmin,date2)
  idate2(4)=date2(1)
  idate2(2)=date2(2)
  idate2(3)=date2(3)
  idate2(1)=date2(4)
  iymdh2 = idate2(4)*1000000+idate2(2)*10000+idate2(3)*100+idate2(1)
  if (iymdh.ne.iymdh2) then
    print *, 'valid dates are different, ',iymdh,' ',iymdh2
    stop 99
  end if
!  call read_sig(nsig,igrd1,jgrd1,levs,nfldsig,nonhyd,fhour,sl,dfld,mapf,clat,clon)
  call read_sig(nsig,igrd1,jgrd1,levs,nfldsig,nonhyd,icld,0.0,sl,dfld,mapf,clat,clon)
  do j=1,nlat
    do i=1,nlon
     ps(i,j)=ps(i,j)-dfld(i+ilonw-1,j+jlats-1,ips)
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        t(i,j,k) = t(i,j,k)-dfld(i+ilonw-1,j+jlats-1,it+k-1)
!        t(i,j,k) = dfld(i+ilonw-1,j+jlats-1,it+k-1)
!        theta(i,j,k) = theta(i,j,k) - t(i,j,k)*(p0/dfld(i+ilonw-1,j+jlats-1,ips)/sl(k))**ptheta
      end do
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        u(i,j,k) = u(i,j,k)-dfld(i+ilonw-1,j+jlats-1,iu+k-1)
      end do
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        v(i,j,k) = v(i,j,k)-dfld(i+ilonw-1,j+jlats-1,iv+k-1)
      end do
    end do
  end do
  do k=1,kmax
    do j=1,nlat
      do i=1,nlon
        q(i,j,k) = q(i,j,k)-dfld(i+ilonw-1,j+jlats-1,iq+k-1)
      end do
    end do
  end do
  do k=1,kmax
    print*, k
    print *, 'u(prtb,max)', maxval(u(:,:,k)),maxloc(u(:,:,k))
    print *, 'u(prtb,min)', minval(u(:,:,k)),minloc(u(:,:,k))
    print *, 'v(prtb,max)', maxval(v(:,:,k)),maxloc(v(:,:,k))
    print *, 'v(prtb,min)', minval(v(:,:,k)),minloc(v(:,:,k))
    print *, 't(prtb,max)', maxval(t(:,:,k)),maxloc(t(:,:,k))
    print *, 't(prtb,min)', minval(t(:,:,k)),minloc(t(:,:,k))
!    print *, 'th(prtb,max)', maxval(theta(:,:,k)),maxloc(theta(:,:,k))
!    print *, 'th(prtb,min)', minval(theta(:,:,k)),minloc(theta(:,:,k))
    print *, 'q(prtb,max)', maxval(q(:,:,k)),maxloc(q(:,:,k))
    print *, 'q(prtb,min)', minval(q(:,:,k)),minloc(q(:,:,k))
  end do
  print *, 'ps(prtb,max)', maxval(ps(:,:)),maxloc(ps(:,:))
  print *, 'ps(prtb,min)', minval(ps(:,:)),minloc(ps(:,:))
  end if
!  do k=1,10
!    print *, u(:,67,k)
!  end do
  ! calculate energy
  !call calc_te(u,v,theta,q,ps,epsq,clat(jlats:jlatn),si,nlon,nlat,kmax,te)
  call calc_te(u,v,t,q,ps,epsq,clat(jlats:jlatn),si,nlon,nlat,kmax,te)
  open(55,file=ofile)
  write(55,'(A1,A12,4A13)') '#','ke','pe(t)','lh','pe(ps)','sum'
  write(55,'(5F13.5)') te(1),te(2),te(3),te(4),sum(te)
  close(55)
  deallocate( dfld,u,v,t,q,ps,theta,fact ) 
end program
