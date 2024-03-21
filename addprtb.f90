program addprtb
!
! add rescaled perturbations based on error moist total energy
!
  use kind_module
  use phconst_module
  use rsmcom_module
  use read_module
  use write_module, only: write_sig, write_sfc
  use norm_module, only: calc_te, calc_te2
  use func_module, only: calc_rh, calc_q2 !=> calc_q
  use date_module, only: ndate
  implicit none
  ! for energy calculation
  !real(kind=dp), parameter :: teref=5.6d0 !rescaled total energy [J/kg/m2]
  !! Reference values below are following Saito et al. (2011, Tellus A).
  real(kind=dp), parameter :: uscl=1.0d0,vscl=1.0d0&
                           &,thetascl=0.4d0,rhscl=5.0d-2,psscl=3.5d1
  real(kind=dp), parameter :: p0=1000.0d2, ptheta=rd/cp ! potential temperature
  real(kind=dp)            :: tscl,qscl,pscl,tbase,qbase

  integer       :: member=10     !number of perturbations
  real(kind=dp) :: teref=0.0d0 !rescaled total energy [J/kg/m2]
  real(kind=dp) :: epsq=1.0d0   !weight for moist term (0.0=dry)
  logical       :: setnorm=.FALSE. !whether rescaling norm magnitude is given from namelist or not
  real(kind=dp) :: alpha=0.0d0 !rescaled factor
  real(kind=dp) :: lonw=-999.9d0, lone=-999.9d0 !calculation region
  real(kind=dp) :: lats=-999.9d0, latn=-999.9d0 !calculation region
  integer       :: ilonw,ilone,jlats,jlatn !calculation region
  integer       :: nlonl,nlatl
  integer       :: kmax=21
  logical       :: adjust_q=.false. !whether super saturation and super dry are removed or not
  logical       :: orth=.false. !orthogonalization
  logical       :: lsub=.false. !subtract perturbations
  namelist /namlst_prtb/ member,&
          setnorm,teref,epsq,&
          lonw,lone,lats,latn,kmax,&
          adjust_q,orth,lsub
  real(kind=dp), allocatable :: u(:,:,:),v(:,:,:),t(:,:,:),q(:,:,:)
  real(kind=dp), allocatable :: ps(:,:)
!  real(kind=dp), allocatable :: fact(:,:,:),theta(:,:,:)
  real(kind=dp), allocatable :: u2(:,:,:),v2(:,:,:),t2(:,:,:),q2(:,:,:)
  real(kind=dp), allocatable :: ps2(:,:)
  real(kind=dp) :: tecmp(4),tecross
  real(kind=dp), allocatable :: te(:)
  real(kind=dp) :: area,coef
  real(kind=dp) :: t1,p1,q1,rh1,cw1,qlim,tulim,tllim !for q adjustment
  integer :: ips,it,iu,iv,iq,icw
  ! input filenames
  character(len=15) :: file_basename='r_.@@@@.LEV.grd'
  character(len=15) :: filename
  ! input files' units (base, prtb)
  integer :: nisigb, nisigp1, nisigp2
!  integer :: nisigib=15 !intermediate file
  integer :: nisfc
  logical :: lexist
  ! output file's unit
  integer :: nosig, nosigm !plus and minus
  integer :: nosfc 
  real(kind=dp), allocatable :: dfld(:,:,:)
  real(kind=dp), allocatable :: dfldb(:,:,:),dfldm(:,:,:),dfldp(:,:,:,:)
  !real(kind=dp), allocatable :: mapf(:,:,:), rlat(:), rlon(:), slmsk(:,:)
  real(kind=dp), allocatable :: dummapf(:,:,:), dumlat(:), dumlon(:)
!  character(len=8) :: label(4)
!  integer :: idate(4)
  real(kind=sp) :: ext(nwext) 
!  real(kind=sp) :: fhour, zhour
!  real(kind=dp) :: si(levmax+1), sl(levmax)
!  real(kind=dp) :: rdelx, rdely, rtruth, rorient, rproj
!  integer :: igrd1, jgrd1, levs, nonhyd, icld
  integer :: nfldsig, levs
  ! for ndate
  integer :: date1(5),date2(5),dtmin
  !
  integer :: n,i,j,k,l,im,km

  read(5,namlst_prtb)
  write(6,namlst_prtb)

!!! sigma files (r_sig.fNN)
  filename=file_basename
  write(filename(1:2),'(a2)') 'rb'
  write(filename(4:7),'(i4.4)') 0
  write(filename(9:11),'(a3)') 'sig'
  call set_rsmparm(filename(1:7))
!  icld=1
  levs=nlev
  nfldsig=nflds
  ext(1) = real(iwav1,kind=sp)
  ext(2) = real(jwav1,kind=sp)
  ext(3) = real(igrd1,kind=sp)
  ext(4) = real(jgrd1,kind=sp)
  ext(5) = real(nlev,kind=sp)
  ext(6) = real(2+nlev*4+5,kind=sp)
  ext(7) = real(rproj,kind=sp)
  ext(8) = real(rtruth,kind=sp)
  ext(9) = real(rorient,kind=sp)
  ext(10)= real(rcenlat,kind=sp)
  ext(11)= real(rcenlon,kind=sp)
  ext(12)= real(rlftgrd,kind=sp)
  ext(13)= real(rbtmgrd,kind=sp)
  ext(14)= real(rdelx,kind=sp)
  ext(15)= real(rdely,kind=sp)
  ext(16)= real(nonhyd,kind=sp)
  allocate( dfldb(igrd1,jgrd1,nfldsig) ) !base(control)
  allocate( dfldp(igrd1,jgrd1,nfldsig,member) ) !perturbation
  allocate( dfldm(igrd1,jgrd1,nfldsig) ) !perturbation mean
  allocate( dummapf(igrd1,jgrd1,3) )
  allocate( dumlat(jgrd1), dumlon(igrd1) )
  ips=2
  if(nonhyd.eq.1) then
    it=2+7*levs
  else
    it=3
  end if
  iu=3+levs
  iv=iu+levs
  iq=iv+levs
  icw=iq+2*levs
  !! set boundaries
  if((lonw.gt.-999.9d0).and.(lone.gt.-999.9d0)) then
    do i=1,igrd1
      ilonw=i
      if(rlon(i).ge.lonw) exit
    end do
    do i=1,igrd1
      ilone=i
      if(rlon(i).ge.lone) exit
    end do
  else
    ilonw=1
    ilone=igrd1
  end if
  if((lats.gt.-999.9d0).and.(latn.gt.-999.9d0)) then
    do j=1,jgrd1
      jlats=j
      if(rlat(j).ge.lats) exit
    end do
    do j=1,jgrd1
      jlatn=j
      if(rlat(j).ge.latn) exit
    end do
  else
    jlats=1
    jlatn=jgrd1
  end if
  print *, "boundary ",ilonw,"-",ilone," lon ",rlon(ilonw),"-",rlon(ilone)
  print *, "boundary ",jlats,"-",jlatn," lat ",rlat(jlats),"-",rlat(jlatn)
  nlonl = ilone - ilonw + 1
  nlatl = jlatn - jlats + 1
  print *, 'nlon ',nlonl,' nlat ',nlatl
  allocate( u(nlonl,nlatl,kmax),v(nlonl,nlatl,kmax) )
  allocate( t(nlonl,nlatl,kmax),q(nlonl,nlatl,kmax) )
!  allocate( theta(nlonl,nlatl,kmax),fact(nlonl,nlatl,kmax) )
  allocate( ps(nlonl,nlatl) )
  if(orth) then
    allocate( u2(nlonl,nlatl,kmax),v2(nlonl,nlatl,kmax) )
    allocate( t2(nlonl,nlatl,kmax),q2(nlonl,nlatl,kmax) )
    allocate( ps2(nlonl,nlatl) )
  end if
  if(.not.setnorm) then
    ! calculate reference energy
    teref=0.0d0
    area=0.0d0
    do k=1,kmax
      pscl=p0*sig(k)
      tscl=thetascl*sig(k)**(1.0d0/ptheta)
      tbase=tr
      call calc_q2(tbase,rhscl,pscl,qbase)
      tbase=tr+tscl
      call calc_q2(tbase,rhscl,pscl,qscl)
      qscl=qscl-qbase
      print '(5(a,es11.4))', 'uscl ',uscl,' vscl ',vscl,' tscl ',tscl,' qscl ',qscl,' psscl ',psscl 
      do j=1,nlatl
        coef=(sigh(k)-sigh(k+1))*cos(rlat(j+jlats-1)*deg2rad)
        do i=1,nlonl
          !KE
          teref=teref+(uscl*uscl+vscl*vscl)*coef
          !PE(T)
          teref=teref+cp/tr*thetascl*thetascl*coef
          !LE
          teref=teref+epsq*lh**2/cp/tr*qscl*qscl*coef
        end do
      end do
    end do
    do j=1,nlatl
      coef=cos(rlat(j+jlats-1)*deg2rad)
      do i=1,nlonl
        !PE(Ps)
        teref=teref+rd*tr*psscl*psscl/pr/pr*coef
        area=area+coef
      end do
    end do
    teref=teref*0.5d0/area
    print*, 'reference total energy = ', teref
  end if
  
  ! base field
  nisigb=11
  filename=file_basename
  write(filename(1:2),'(a2)') 'rb'
  write(filename(4:7),'(i4.4)') 0
  write(filename(9:11),'(a3)') 'sig'
  write(6,'(2a)') 'input base= ',filename
  open(nisigb,file=filename,form='unformatted',access='sequential',action='read')
  call read_sig(nisigb,igrd1,jgrd1,levs,nfldsig,nonhyd,icld,fhour,sig,&
    dfldb,dummapf,dumlat,dumlon)

  allocate( dfld(igrd1,jgrd1,nfldsig) )
  dfldm = 0.0d0
  nisigp1=13
  nisigp2=14
  do im=1,member
    ! perturbation (nisigp1 - nisigp2)
    filename=file_basename
    write(filename(1:2),'(a2)') 'ri'
    write(filename(4:7),'(i4.4)') 2*im-1
    write(filename(9:11),'(a3)') 'sig'
    write(6,'(2a)') 'input prtb1= ',filename
    open(nisigp1,file=filename,form='unformatted',access='sequential',action='read')
    dfldp(:,:,:,im)=0.0
    call read_sig(nisigp1,igrd1,jgrd1,levs,nfldsig,nonhyd,icld,fhour,sig,&
    dfldp(:,:,:,im),dummapf,dumlat,dumlon)
    close(nisigp1)
    ! perturbation
    filename=file_basename
    write(filename(1:2),'(a2)') 'ri'
    write(filename(4:7),'(i4.4)') 2*im
    write(filename(9:11),'(a3)') 'sig'
    write(6,'(2a)') 'input prtb2= ',filename
    open(nisigp2,file=filename,form='unformatted',access='sequential',action='read')
    call read_sig(nisigp2,igrd1,jgrd1,levs,nfldsig,nonhyd,icld,fhour,sig,&
    dfld,dummapf,dumlat,dumlon)
    close(nisigp2)
    dfldp(:,:,:,im) = dfldp(:,:,:,im) - dfld
    dfldm = dfldm + dfldp(:,:,:,im)
    nisigp1=nisigp1+2
    nisigp2=nisigp2+2
  end do !im=1,member   
  deallocate( dfld )

  dfldm = dfldm / real(member,kind=dp)

  ! perturbation energy
  allocate( te(member) )
  nosig=51
  nosfc=52
  do im=1,member
    u=0.0d0
    v=0.0d0
    t=0.0d0
    q=0.0d0
    ps=0.0d0
    do j=1,nlatl
      do i=1,nlonl
        ps(i,j)=dfldp(i+ilonw-1,j+jlats-1,ips,im)
      end do
    end do
    do k=1,kmax
      do j=1,nlatl
        do i=1,nlonl
          t(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,it+k-1,im)
          !theta(i,j,k) = t(i,j,k) * (p0/dfld(i+ilonw-1,j+jlats-1,ips)/sig(k))**ptheta
        end do
      end do
    end do
    do k=1,kmax
      do j=1,nlatl
        do i=1,nlonl
          u(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iu+k-1,im)
        end do
      end do
    end do
    do k=1,kmax
      do j=1,nlatl
        do i=1,nlonl
          v(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iv+k-1,im)
        end do
      end do
    end do
    do k=1,kmax
      do j=1,nlatl
        do i=1,nlonl
          q(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iq+k-1,im)
        end do
      end do
    end do
    if(orth.and.im.gt.1) then
      do km=1,im-1
        u2=0.0d0
        v2=0.0d0
        t2=0.0d0
        q2=0.0d0
        ps2=0.0d0
        do j=1,nlatl
          do i=1,nlonl
            ps2(i,j)=dfldp(i+ilonw-1,j+jlats-1,ips,km)
          end do
        end do
        do k=1,kmax
          do j=1,nlatl
            do i=1,nlonl
              t2(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,it+k-1,km)
            end do
          end do
        end do
        do k=1,kmax
          do j=1,nlatl
            do i=1,nlonl
              u2(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iu+k-1,km)
            end do
          end do
        end do
        do k=1,kmax
          do j=1,nlatl
            do i=1,nlonl
              v2(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iv+k-1,km)
            end do
          end do
        end do
        do k=1,kmax
          do j=1,nlatl
            do i=1,nlonl
              q2(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iq+k-1,km)
            end do
          end do
        end do
        call calc_te2(u,u2,v,v2,t,t2,q,q2,ps,ps2,&
              epsq,rlat(jlats:jlatn),sigh,nlonl,nlatl,kmax,tecmp)
        tecross=sum(tecmp)
        print*, 'cross perturbation energy = ',tecross
        dfldp(:,:,2:,im) = dfldp(:,:,2:,im) - dfldp(:,:,2:,km)*tecross/te(km)
      end do
      u=0.0d0
      v=0.0d0
      t=0.0d0
      q=0.0d0
      ps=0.0d0
      do j=1,nlatl
        do i=1,nlonl
          ps(i,j)=dfldp(i+ilonw-1,j+jlats-1,ips,im)
        end do
      end do
      do k=1,kmax
        do j=1,nlatl
          do i=1,nlonl
            t(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,it+k-1,im)
            !theta(i,j,k) = t(i,j,k) * (p0/dfld(i+ilonw-1,j+jlats-1,ips)/sig(k))**ptheta
          end do
        end do
      end do
      do k=1,kmax
        do j=1,nlatl
          do i=1,nlonl
            u(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iu+k-1,im)
          end do
        end do
      end do
      do k=1,kmax
        do j=1,nlatl
          do i=1,nlonl
            v(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iv+k-1,im)
          end do
        end do
      end do
      do k=1,kmax
        do j=1,nlatl
          do i=1,nlonl
            q(i,j,k) = dfldp(i+ilonw-1,j+jlats-1,iq+k-1,im)
          end do
        end do
      end do
    end if
    print*, 'u(prtb)', maxval(u),minval(u)
    print*, 'v(prtb)', maxval(v),minval(v)
    print*, 't(prtb)', maxval(t),minval(t)
!    print*, 'theta(prtb)', maxval(theta),minval(theta)
    print*, 'q(prtb)', maxval(q),minval(q)
    print*, 'ps(prtb)', maxval(ps),minval(ps)
    
    !if(alpha.eq.0.0d0) then
    ! determining rescale factor based on perturbation energy
    print*, 'normalized total energy = ', teref
    ! calculate energy
    !call calc_te(u,v,theta,q,ps,epsq,rlat(jlats:jlatn),sigh,nlonl,nlatl,tecmp)
    call calc_te(u,v,t,q,ps,epsq,rlat(jlats:jlatn),sigh,nlonl,nlatl,kmax,tecmp)
    te(im)=sum(tecmp)
    print *, tecmp
    print*, im,'th perturbation total energy = ', te(im)
    ! rescaling
    alpha = sqrt(teref / te(im))
    if(lsub) then
      alpha = -1.0d0 * alpha
    end if
    !end if ! alpha.eq.0.0d0
    print*, 'rescaling factor = ', alpha
    ! write output
    allocate( dfld(igrd1,jgrd1,nfldsig) )
    !! reset forecast hour
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
    fhour=0.0
    print *, idate(4),idate(2),idate(3),idate(1),'+',nint(fhour)
    !! add perturbations
    dfld(:,:,1) = dfldb(:,:,1) !surface elevation
    dfld(:,:,2:) = dfldb(:,:,2:) + dfldp(:,:,2:,im) * alpha
    if(adjust_q) then
      ! super saturation(dry) adjustment
      tllim = t0 - 30.0_dp
      tulim = t0 + 35.0_dp
      do k=1,levs
        do j=1,jgrd1
          do i=1,igrd1
            t1 = dfld(i,j,it+k-1)
            q1 = dfld(i,j,iq+k-1)
            cw1 = dfld(i,j,icw+k-1)
            p1 = dfld(i,j,ips)*sig(k)
            if(q1.lt.0.0_dp) then !super dry
              write(0,'(a,f10.2,a,es10.2,a)') &
                  'super dry adjustment: p=',p1,' q=',q1,' < 0.0'
              dfld(i,j,iq+k-1)=0.0_dp
            else if(p1.gt.20000.0_dp.and.(t1.gt.tllim.and.t1.lt.tulim)) then
            !saturation water vapor accuracy is acceptable for p > 200mb, -30 celsius < T < 35 celsius
              rh1=1.2_dp
              call calc_q2(t1,rh1,p1,qlim)
              if(q1.gt.qlim) then !super saturation
                write(0,'(a,f10.2,x,f10.2,a,f10.2,a,f10.2,x,a,es10.2,a,es10.2)') &
                  'super saturation adjustment: p=',p1,&
                  tllim,' < t=',t1,' < ',tulim,&
                  ' q=',q1,' > ',qlim
                dfld(i,j,iq+k-1)=qlim
              end if
            end if
            if(cw1.lt.0.0d0) then !negative cloud water
              write(0,'(a,f10.2,a,es10.2,a)') &
                 'super dry adjustment: p=',p1,' cw=',cw1,'<0.0'
              dfld(i,j,icw+k-1)=0.0_dp
            end if
          end do
        end do
      end do
    end if
    ! write output
    filename=file_basename
    write(filename(1:2),'(a2)') 'ro'
    write(filename(4:7),'(i4.4)') im
    write(filename(9:11),'(a3)') 'sig'
    write(6,'(2a)') 'output= ',filename
    open(nosig,file=filename,form='unformatted',access='sequential',action='write')
    call write_sig(nosig,label,idate,fhour,sigh,sig,ext,&
&                    igrd1,jgrd1,levs,nfldsig,nonhyd,icld,dfld,mapf,rlat,rlon)
    close(nosig)
    ! read surface and change forecast date and hour, then write out
    deallocate( dfld )
    allocate( dfld(igrd1,jgrd1,nfldsfc) )
    nisfc=12
    filename=file_basename
    write(filename(1:2),'(a2)') 'rb'
    write(filename(4:7),'(i4.4)') 0
    write(filename(9:11),'(a3)') 'sfc'
    write(6,'(2a)') 'input= ',filename
    open(nisfc,file=filename,form='unformatted',access='sequential',action='read')
    call read_sfc(nisfc,igrd1,jgrd1,dfld)
    close(nisfc)
    filename=file_basename
    write(filename(1:2),'(a2)') 'ro'
    write(filename(4:7),'(i4.4)') im
    write(filename(9:11),'(a3)') 'sfc'
    write(6,'(2a)') 'output= ',filename
    open(nosfc,file=filename,form='unformatted',access='sequential',action='write')
    call write_sfc(nosfc,igrd1,jgrd1,dfld,label,idate,fhour)
    close(nosfc)
    deallocate( dfld )
    nosig=nosig+2
    nosfc=nosfc+2
  end do
  deallocate( dfldm,dfldb,dfldp,u,v,t,q,ps ) 
end program
