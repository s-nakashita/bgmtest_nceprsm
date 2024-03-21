module rsmcom_module
!
! RSM common information and routines
!
! history:
! 22-08-04 SN create
!
  use kind_module
  use phconst_module, only : rad2deg
  use read_module, only: levmax,nwext,read_header,read_sig,nfldsfc,read_sfc &
          ,nfldflx,read_flx
  use write_module, only: write_sig,write_sfc
  implicit none
  public
! 
! parameters read from header
!
  !!! horizontal grid numbers
  integer,save :: igrd1, jgrd1
  integer,save :: nlon, nlat
  integer,save :: lngrd
  !!! horizontal wave numbers
  integer,save :: iwav1, jwav1
  integer,save :: lnwav
  !!! vertical levels
  integer,save :: nlev
  !!! projection parameters
  integer,save :: rproj !0=mercator, 1.or.-1=polar
  real(kind=dp),save :: rtruth,rorient,rdelx,rdely,rcenlat,rcenlon
  real(kind=dp),save :: rlftgrd, rbtmgrd
  !!! other header information (saved for write restart file)
  real(kind=sp),save :: fhour
  character(len=8),save :: label(4)
  integer,save :: idate(4)
  !!! model version
  integer,save :: nonhyd !0=hydrostatic, 1=nonhydrostatic
  integer,parameter   :: icld=0 !1=output includes 3D physics, 0=not includes
!
! grid information
!
  real(kind=dp),allocatable,save :: rlon(:),rlat(:)
  real(kind=dp),allocatable,save :: sig(:),sigh(:)
!
! map factors
!
  real(kind=dp),allocatable,save :: mapf(:,:) !1:square of map factor
                                                     !2:d(fm2)/dx
                                                     !3:d(fm2)/dy
! 
! parameters read from base header
!
  !!! horizontal grid numbers
  integer,save :: cigrd1, cjgrd1
  integer,save :: clngrd
  !!! horizontal wave numbers
  integer,save :: ciwav1, cjwav1
  integer,save :: clnwav
  !!! vertical levels
  integer,save :: cnlev
  !!! projection parameters
  integer,save :: cproj !0=mercator, 1.or.-1=polar
  real(kind=dp),save :: ctruth,corient,cdelx,cdely,ccenlat,ccenlon
  real(kind=dp),save :: clftgrd, cbtmgrd
  !!! model version
  integer,save :: cnonhyd !0=hydrostatic, 1=nonhydrostatic
  real(kind=dp),allocatable,save :: clon(:),clat(:)
  real(kind=dp),allocatable,save :: csig(:),csigh(:)
!
! variables' indexes
!
  !!! in r_sig.fNN
  integer,save :: nflds
  !!! 2D variables
  integer,save      :: nv2d
  integer,save      :: nv2d_sig    !gz,ps,(wbtm)
  integer,parameter :: iv2d_gz=1 ! terrain height
  integer,parameter :: iv2d_ps=2 ! surface pressure
  integer,parameter :: iv2d_wb=3 ! vertical velocity at bottom level
  integer,parameter :: nv2d_sfc=nfldsfc !r_sfc
  integer,parameter :: iv2d_tsfc=1
  integer,parameter :: iv2d_smc_1=2
  integer,parameter :: iv2d_smc_2=3
  integer,parameter :: iv2d_sheleg=4
  integer,parameter :: iv2d_stc_1=5
  integer,parameter :: iv2d_stc_2=6
  integer,parameter :: iv2d_tg3=7
  integer,parameter :: iv2d_zorl=8
  integer,parameter :: iv2d_cv=9
  integer,parameter :: iv2d_cvb=10
  integer,parameter :: iv2d_cvt=11
  integer,parameter :: iv2d_alvsf=12
  integer,parameter :: iv2d_alvwf=13
  integer,parameter :: iv2d_alnsf=14
  integer,parameter :: iv2d_alnwf=15
  integer,parameter :: iv2d_slmsk=16
  integer,parameter :: iv2d_vfrac=17
  integer,parameter :: iv2d_canopy=18
  integer,parameter :: iv2d_f10m=19
  integer,parameter :: iv2d_vtype=20
  integer,parameter :: iv2d_stype=21
  integer,parameter :: iv2d_facsf=22
  integer,parameter :: iv2d_facwf=23
  integer,parameter :: iv2d_uustar=24
  integer,parameter :: iv2d_ffmm=25
  integer,parameter :: iv2d_ffhh=26
  character(len=6),parameter :: varnames_sfc(nfldsfc) = (/&
  '  TSFC',' SMC_1',' SMC_2','SHELEG',' STC_1',' STC_2','   TG3',&
  '  ZORL','    CV','   CVB','   CVT',' ALVSF',' ALVWF',' ALNSF',' ALNWF',&
  ' SLMSK',' VFRAC','CANOPY','  F10m',&
  ' VTYPE',' STYPE',' FACSF',' FACWF','UUSTAR','  FFMM','  FFHH'/)
  integer,parameter :: nv2d_flx=4 !r_flx(diagnostic,no need to write)
  integer,parameter :: iv2d_t2m=1
  integer,parameter :: iv2d_q2m=2
  integer,parameter :: iv2d_u10m=3
  integer,parameter :: iv2d_v10m=4
  character(len=6),parameter :: varnames_flx(nv2d_flx) = (/&
          '    Ts','    Qs','    Us','    Vs'/)
  integer,save      :: ind_t2m,ind_q2m,ind_u10m,ind_v10m
  !!! 3D variables
  !!! pn,tn,ww are only for nonhydrostatic
  integer,save      :: nv3d
  integer,parameter :: nv3d_hyd=6    !th,u,v,q,oz,cw 
  integer,parameter :: nv3d_nonhyd=3 !pn,tn,ww(halflevel)
  integer,parameter :: iv3d_th=1 ! hydrostatic temperature
  integer,parameter :: iv3d_u=2  ! x-direction wind
  integer,parameter :: iv3d_v=3  ! y-direction wind
  integer,parameter :: iv3d_q=4  ! specific humidity
  integer,parameter :: iv3d_oz=5 ! ozone mixing ratio
  integer,parameter :: iv3d_cw=6 ! cloud water
  integer,parameter :: iv3d_pn=7 ! full pressure
  integer,parameter :: iv3d_tn=8 ! non-hydrostatic temperature
  integer,parameter :: iv3d_wn=9 ! vertical velocity at half levels (except bottom)
  integer,save      :: iv3d_t,iv3d_tt !update and not update temperature indexes
  
  integer,save :: nlevall
  character(len=6),allocatable :: varnames(:)
!
! IO
!
  character(len=4) :: filesuffix='.grd'
!
  contains
!
! set grid information
!
  subroutine set_rsmparm(cfile)
    implicit none
!    integer,intent(in) :: nsig ! input sigma file unit
    character(len=*), intent(in) :: cfile !input sigma file
    integer :: nsig
    character(len=filelenmax) :: filename
    integer :: nskip
    real(kind=sp) :: ext(nwext)
    real(kind=dp) :: si(levmax+1),sl(levmax)
    real(kind=sp),allocatable :: sfld(:)
    integer :: iwav, jwav
    integer :: i,j,k
    
    write(6,'(A)') 'start set_rsmparm'
    
    !nsig=70
    call search_fileunit(nsig)
    filename=trim(cfile)//'.sig'//filesuffix
    write(6,'(3a,i3)') 'open file ',trim(filename),' unit=',nsig
    open(nsig,file=filename,access='sequential',form='unformatted',action='read')
    call read_header(nsig,icld,label,idate,fhour,si,sl,ext,nflds)
    iwav1  = int(ext(1))
    jwav1  = int(ext(2))
    igrd1  = int(ext(3))
    jgrd1  = int(ext(4))
    nlev   = int(ext(5))
    !nfldx  = int(ext(6))
    rproj  = int(ext(7))
    rtruth =real(ext(8),kind=dp)
    rorient=real(ext(9),kind=dp)
    rcenlat=real(ext(10),kind=dp) 
    rcenlon=real(ext(11),kind=dp)
    rlftgrd=real(ext(12),kind=dp)
    rbtmgrd=real(ext(13),kind=dp)
    rdelx  =real(ext(14),kind=dp)
    rdely  =real(ext(15),kind=dp)
    nonhyd =int(ext(16))
    lngrd=igrd1*jgrd1
    write(6,'(A,2i6)') 'igrd1 jgrd1=',igrd1,jgrd1
    write(6,'(A,I8)') 'lngrd=',lngrd
    lnwav = iwav1*jwav1
    write(6,'(A,2i6)') 'iwav1 jwav1=',iwav1,jwav1
    write(6,'(A,I8)') 'lnwav=',lnwav
    nlon = igrd1
    nlat = jgrd1
    allocate( rlon(nlon), rlat(nlat) )
    allocate( sig(nlev), sigh(nlev+1) )
    allocate( mapf(lngrd,3) )
    do i=1,nlev
      sig(i) = sl(i)
      sigh(i) = si(i)
    end do
    sigh(nlev+1) = si(nlev+1)
    write(6,'(a,2f7.4)') 'sig=',minval(sig),maxval(sig)
    write(6,'(a,2f7.4)') 'sigh=',minval(sigh),maxval(sigh)

    if(nonhyd.eq.1) then
      nv2d_sig=3
      nv2d=nv2d_sig+nv2d_sfc+nv2d_flx
      nv3d=nv3d_hyd+nv3d_nonhyd
      write(6,*) 'model version is nonhydrostatic'
      allocate( varnames(nv3d+nv2d) )
      varnames(1:nv3d+nv2d_sig) = (/&
              '    Th','     U','     V','     Q','    OZ','    CW',&
              '    Pn','    Tn','    Wn','    GZ','    Ps','    Wb'/)
    else
      nv2d_sig=2
      nv2d=nv2d_sig+nv2d_sfc+nv2d_flx
      nv3d=nv3d_hyd
      write(6,*) 'model version is hydrostatic'
      allocate( varnames(nv3d+nv2d) )
      varnames(1:nv3d+nv2d_sig) = (/&
              '     T','     U','     V','     Q','    OZ','    CW',&
              '    GZ','    Ps'/)
    end if
    varnames(nv3d+nv2d_sig+1:nv3d+nv2d_sig+nv2d_sfc) = varnames_sfc
    varnames(nv3d+nv2d_sig+nv2d_sfc+1:) = varnames_flx
    nlevall=nv3d*nlev+nv2d
!    nskip=2+(nlevall-nv2d_sfc)!+3
    nskip=2+(nlevall-nv2d_sfc-nv2d_flx)!+3
    allocate( sfld(lngrd) )
    rewind(nsig)
    do i=1,nskip
      read(nsig)
    end do
    !fm2
    read(nsig) (sfld(i),i=1,lngrd)
    do i=1,lngrd
      mapf(i,1) = real(sfld(i),kind=dp)
    end do
    !fm2x
    read(nsig) (sfld(i),i=1,lngrd)
    do i=1,lngrd
      mapf(i,2) = real(sfld(i),kind=dp)
    end do
    !fm2y
    read(nsig) (sfld(i),i=1,lngrd)
    do i=1,lngrd
      mapf(i,3) = real(sfld(i),kind=dp)
    end do
    !rlat(S->N)
    read(nsig) (sfld(i),i=1,lngrd)
    k=1
    do j=1,nlat
      rlat(j) = real(sfld(k),kind=dp)*rad2deg
      k=k+nlon
    end do
    write(6,'(a,2f9.2)') 'rlat=',minval(rlat),maxval(rlat)
    !rlon(W->E)
    read(nsig) (sfld(i),i=1,lngrd)
    do i=1,nlon
      rlon(i) = real(sfld(i),kind=dp)*rad2deg
    end do
    write(6,'(a,2f9.2)') 'rlon=',minval(rlon),maxval(rlon)
    deallocate(sfld)
    close(nsig)
    return

  end subroutine set_rsmparm
!
! set base grid information
!
  subroutine set_rsmparm_base(cfile)
    implicit none
!    integer,intent(in) :: nsig ! input sigma file unit
    character(len=*), intent(in) :: cfile !input sigma file
    integer :: nsig
    character(len=filelenmax) :: filename
    integer :: nskip
    real(kind=sp) :: ext(nwext)
    real(kind=dp) :: si(levmax+1),sl(levmax)
    real(kind=sp),allocatable :: sfld(:)
    integer :: iwav, jwav
    integer :: nv3dc, nv2dc
    !!! dummy header information
    real(kind=sp) :: dumfhour
    character(len=8) :: dumlabel(4)
    integer :: dumidate(4)
    integer :: dumnflds
    integer :: i,j,k
    
    write(6,'(A)') 'start set_rsmparm_base'
    
    !nsig=70
    call search_fileunit(nsig)
    filename=trim(cfile)//'.sig'//filesuffix
    write(6,'(3a,i3)') 'open file ',trim(filename),' unit=',nsig
    open(nsig,file=filename,access='sequential',form='unformatted',action='read')
    call read_header(nsig,icld,dumlabel,dumidate,dumfhour,si,sl,ext,dumnflds)
    ciwav1  = int(ext(1))
    cjwav1  = int(ext(2))
    cigrd1  = int(ext(3))
    cjgrd1  = int(ext(4))
    cnlev   = int(ext(5))
    !nfldx  = int(ext(6))
    cproj  = int(ext(7))
    ctruth =real(ext(8),kind=dp)
    corient=real(ext(9),kind=dp)
    ccenlat=real(ext(10),kind=dp) 
    ccenlon=real(ext(11),kind=dp)
    clftgrd=real(ext(12),kind=dp)
    cbtmgrd=real(ext(13),kind=dp)
    cdelx  =real(ext(14),kind=dp)
    cdely  =real(ext(15),kind=dp)
    cnonhyd =int(ext(16))
    clngrd=cigrd1*cjgrd1
    write(6,'(A,2i6)') 'cigrd1 cjgrd1=',cigrd1,cjgrd1
    write(6,'(A,I8)') 'clngrd=',clngrd
    clnwav = ciwav1*cjwav1
    write(6,'(A,2i6)') 'ciwav1 cjwav1=',ciwav1,cjwav1
    write(6,'(A,I8)') 'clnwav=',clnwav
    allocate( clon(cigrd1), clat(cjgrd1) )
    allocate( csig(cnlev), csigh(cnlev+1) )
    do i=1,cnlev
      csig(i) = sl(i)
      csigh(i) = si(i)
    end do
    csigh(cnlev+1) = si(cnlev+1)
    write(6,'(a,2f7.4)') 'csig=',minval(csig),maxval(csig)
    write(6,'(a,2f7.4)') 'csigh=',minval(csigh),maxval(csigh)

    if(cnonhyd.eq.1) then
      nv2dc=3+nv2d_sfc
      nv3dc=nv3d_hyd+nv3d_nonhyd
      write(6,*) 'base model version is nonhydrostatic'
    else
      nv2dc=2+nv2d_sfc
      nv3dc=nv3d_hyd
      write(6,*) 'base model version is hydrostatic'
    end if
    nskip=2+(nv3dc*cnlev+nv2dc-nv2d_sfc)+3
    allocate( sfld(clngrd) )
    rewind(nsig)
    do i=1,nskip
      read(nsig)
    end do
    !clat(S->N)
    read(nsig) (sfld(i),i=1,clngrd)
    k=1
    do j=1,nlat
      clat(j) = real(sfld(k),kind=dp)*rad2deg
      k=k+nlon
    end do
    write(6,'(a,2f9.2)') 'clat=',minval(clat),maxval(clat)
    !clon(W->E)
    read(nsig) (sfld(i),i=1,clngrd)
    do i=1,nlon
      clon(i) = real(sfld(i),kind=dp)*rad2deg
    end do
    write(6,'(a,2f9.2)') 'clon=',minval(clon),maxval(clon)
    deallocate(sfld)
    close(nsig)
    return

  end subroutine set_rsmparm_base
!
! clean up
!
  subroutine clean_rsmparm
    implicit none

    deallocate(rlon,rlat,sig,sigh,mapf)
  end subroutine clean_rsmparm
!
! ensemble mean
!
  subroutine ensmean_grd(mem,ni,nj,v3d,v2d,v3dm,v2dm)
    implicit none
    integer, intent(in) :: mem, ni,nj
    real(kind=dp),intent(in) :: v3d(ni,nj,nlev,mem,nv3d)
    real(kind=dp),intent(in) :: v2d(ni,nj,     mem,nv2d)
    real(kind=dp),intent(out):: v3dm(ni,nj,nlev,nv3d)
    real(kind=dp),intent(out):: v2dm(ni,nj,     nv2d)
    integer :: i,j,k,m,n

    do n=1,nv3d
!$OMP PARALLEL DO PRIVATE(i,j,k,m)
      do k=1,nlev
        do j=1,nj
          do i=1,ni
            v3dm(i,j,k,n)=v3d(i,j,k,1,n)
            do m=2,mem
              v3dm(i,j,k,n)=v3dm(i,j,k,n)+v3d(i,j,k,m,n)
            end do
            v3dm(i,j,k,n)=v3dm(i,j,k,n)/real(mem,kind=dp)
          end do
        end do
      end do
!$OMP END PARALLEL DO
    end do

    do n=1,nv2d
!$OMP PARALLEL DO PRIVATE(i,j,m)
      do i=1,ni
        do j=1,nj
          v2dm(i,j,n)=v2d(i,j,1,n)
          do m=2,mem
            v2dm(i,j,n)=v2dm(i,j,n)+v2d(i,j,m,n)
          end do
          v2dm(i,j,n)=v2dm(i,j,n)/real(mem,kind=dp)
        end do
      end do
!$OMP END PARALLEL DO
    end do

    return
  end subroutine ensmean_grd
!
! read restart file
!
  subroutine read_restart(cfile,v3dg,v2dg,convert)
    implicit none
    character(len=*), intent(in) :: cfile
    !integer, intent(in) :: nsig
    real(kind=dp), intent(out) :: v3dg(nlon,nlat,nlev,nv3d)
    real(kind=dp), intent(out) :: v2dg(nlon,nlat,nv2d)
    logical, intent(in), optional :: convert

    integer :: nsig,nsfc,nflx
    character(len=filelenmax) :: filename
    character(len=3) :: clev
    real(kind=dp), allocatable :: dfld(:,:,:)
    real(kind=dp), allocatable :: dummp(:,:,:),dumlat(:),dumlon(:) !dummy
    integer :: k,kk
    ! read_flx
    integer      :: ids(255)
    integer      :: iparam(nfldflx)
    real(kind=sp) :: fhr, zhr

    allocate( dfld(igrd1,jgrd1,nflds) )
    allocate( dummp(igrd1,jgrd1,3) )
    allocate( dumlat(jgrd1), dumlon(igrd1) )

    !nsig=70
    call search_fileunit(nsig)
    clev='sig'
    filename=trim(cfile)//'.'//clev//filesuffix
    write(6,'(3a,i3)') 'open file ',trim(filename),' unit=',nsig
    open(nsig,file=filename,access='sequential',form='unformatted',action='read')
    if( present(convert) ) then
    call read_sig( nsig,igrd1,jgrd1,nlev,nflds,nonhyd,icld,fhour,sig,&
      &  dfld,dummp,dumlat,dumlon,convert=convert )
    else
    call read_sig( nsig,igrd1,jgrd1,nlev,nflds,nonhyd,icld,fhour,sig,&
      &  dfld,dummp,dumlat,dumlon )
    end if
    kk=1
    v2dg(:,:,iv2d_gz)=dfld(:,:,kk)
    kk=kk+1
    v2dg(:,:,iv2d_ps)=dfld(:,:,kk)
    kk=kk+1
    do k=1,nlev
      v3dg(:,:,k,iv3d_th) = dfld(:,:,kk) 
      kk=kk+1
    end do
    do k=1,nlev
      v3dg(:,:,k,iv3d_u) = dfld(:,:,kk)
      kk=kk+1
    end do
    do k=1,nlev
      v3dg(:,:,k,iv3d_v) = dfld(:,:,kk)
      kk=kk+1
    end do
    do k=1,nlev
      v3dg(:,:,k,iv3d_q) = dfld(:,:,kk)
      kk=kk+1
    end do
    do k=1,nlev
      v3dg(:,:,k,iv3d_oz) = dfld(:,:,kk)
      kk=kk+1
    end do
    do k=1,nlev
      v3dg(:,:,k,iv3d_cw) = dfld(:,:,kk)
      kk=kk+1
    end do
    if(nonhyd.eq.1) then
      do k=1,nlev
        v3dg(:,:,k,iv3d_pn) = dfld(:,:,kk)
        kk=kk+1
      end do
      do k=1,nlev
        v3dg(:,:,k,iv3d_tn) = dfld(:,:,kk)
        kk=kk+1
      end do
      v2dg(:,:,iv2d_wb) = dfld(:,:,kk)
      kk=kk+1
      do k=1,nlev
        v3dg(:,:,k,iv3d_wn) = dfld(:,:,kk)
        kk=kk+1
      end do
    end if
    close(nsig)

    deallocate( dfld )
    allocate( dfld(igrd1,jgrd1,nfldsfc) )
    call search_fileunit(nsfc)
    clev='sfc'
    filename=trim(cfile)//'.'//clev//filesuffix
    write(6,'(3a,i3)') 'open file ',trim(filename),' unit=',nsfc
    open(nsfc,file=filename,access='sequential',form='unformatted',action='read')
    call read_sfc(nsfc,igrd1,jgrd1,dfld)
    do k=1,nfldsfc
      v2dg(:,:,nv2d_sig+k) = dfld(:,:,k)
    end do
    close(nsfc)

    deallocate( dfld )
    allocate( dfld(igrd1,jgrd1,nfldflx) )
    call search_fileunit(nflx)
    clev='flx'
    filename=trim(cfile)//'.'//clev//filesuffix
    write(6,'(3a,i3)') 'open file ',trim(filename),' unit=',nflx
    open(nflx,file=filename,access='sequential',form='unformatted',action='read')
    call read_flx(nflx,igrd1,jgrd1,dfld,ids,iparam,fhr,zhr&
                ,ind_t2m,ind_q2m,ind_u10m,ind_v10m)
    v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_t2m) = dfld(:,:,ind_t2m)
    v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_q2m) = dfld(:,:,ind_q2m)
    v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_u10m) = dfld(:,:,ind_u10m)
    v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_v10m) = dfld(:,:,ind_v10m)
    !print *, 't2m ',maxval(v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_t2m))&
    !               ,minval(v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_t2m))
    !print *, 'q2m ',maxval(v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_q2m))&
    !               ,minval(v2dg(:,:,nv2d_sig+nv2d_sfc+iv2d_q2m))
    close(nflx)
    deallocate(dfld)
    return
  end subroutine read_restart
!
! write restart file
!
  subroutine write_restart(cfile,v3dg,v2dg,convert)
    implicit none
    character(len=*), intent(in) :: cfile
    !integer, intent(in) :: nsig
    real(kind=dp), intent(in) :: v3dg(nlon,nlat,nlev,nv3d)
    real(kind=dp), intent(in) :: v2dg(nlon,nlat,nv2d)
    logical, intent(in), optional :: convert

    integer :: nsig,nsfc
    character(len=filelenmax) :: filename
    character(len=3) :: clev
    real(kind=dp), allocatable :: dfld(:,:,:)
    real(kind=dp) :: si(levmax+1),sl(levmax)
    real(kind=sp) :: ext(nwext)
    integer :: k,kk

    allocate( dfld(igrd1,jgrd1,nflds) )
    si=0.0
    sl=0.0
    si(1:nlev+1) = sigh
    sl(1:nlev) = sig
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
    kk=1
    dfld(:,:,kk)=v2dg(:,:,iv2d_gz)
    kk=kk+1
    dfld(:,:,kk)=v2dg(:,:,iv2d_ps)
    kk=kk+1
    do k=1,nlev
      dfld(:,:,kk)=v3dg(:,:,k,iv3d_th)
      kk=kk+1
    end do
    do k=1,nlev
      dfld(:,:,kk)=v3dg(:,:,k,iv3d_u)
      kk=kk+1
    end do
    do k=1,nlev
      dfld(:,:,kk)=v3dg(:,:,k,iv3d_v)
      kk=kk+1
    end do
    do k=1,nlev
      dfld(:,:,kk)=v3dg(:,:,k,iv3d_q)
      kk=kk+1
    end do
    do k=1,nlev
      dfld(:,:,kk)=v3dg(:,:,k,iv3d_oz)
      kk=kk+1
    end do
    do k=1,nlev
      dfld(:,:,kk)=v3dg(:,:,k,iv3d_cw)
      kk=kk+1
    end do
    if(nonhyd.eq.1) then
      do k=1,nlev
        dfld(:,:,kk)=v3dg(:,:,k,iv3d_pn)
        kk=kk+1
      end do
      do k=1,nlev
        dfld(:,:,kk)=v3dg(:,:,k,iv3d_tn)
        kk=kk+1
      end do
      dfld(:,:,kk)=v2dg(:,:,iv2d_wb)
      kk=kk+1
      do k=1,nlev
        dfld(:,:,kk)=v3dg(:,:,k,iv3d_wn)
        kk=kk+1
      end do
    end if
    !nsig=80
    call search_fileunit(nsig)
    clev='sig'
    filename=trim(cfile)//'.'//clev//filesuffix
    write(6,'(3a,i3)') 'open file ',trim(filename),' unit=',nsig
    open(nsig,file=filename,form='unformatted',access='sequential')
    if( present(convert) )then
    call write_sig(nsig,label,idate,fhour,si,sl,ext,&
     &  igrd1,jgrd1,nlev,nflds,nonhyd,icld,dfld,mapf,rlat,rlon,convert=convert)
    else
    call write_sig(nsig,label,idate,fhour,si,sl,ext,&
     &  igrd1,jgrd1,nlev,nflds,nonhyd,icld,dfld,mapf,rlat,rlon)
    end if
    close(nsig)

    deallocate( dfld )
    allocate( dfld(igrd1,jgrd1,nfldsfc) )
    do k=1,nfldsfc
      dfld(:,:,k) = v2dg(:,:,nv2d_sig+k)
    end do
    call search_fileunit(nsfc)
    clev='sfc'
    filename=trim(cfile)//'.'//clev//filesuffix
    write(6,'(3a,i3)') 'open file ',trim(filename),' unit=',nsfc
    open(nsfc,file=filename,access='sequential',form='unformatted')
    call write_sfc(nsfc,igrd1,jgrd1,dfld,label,idate,fhour)
    close(nsfc)

    return
  end subroutine write_restart
!
! search empty file unit (referring to JMA GSM)
!
  subroutine search_fileunit(iunit)
    implicit none
    integer, intent(out) :: iunit
    integer, parameter :: iunit_s=11, iunit_e=99
    integer :: iu
    logical :: lopened, lexist

    iunit=-999

    do iu=iunit_s,iunit_e
      inquire(unit=iu,opened=lopened,exist=lexist)
      if(lexist.and.(.not.lopened)) then
        iunit=iu
        exit
      end if
    end do

    if(iunit.eq.-999) then
      write(6,'(a,i3,a,i3,a)') 'search_fileunit : error : unit number from ', &
       & iunit_s, ' to ', iunit_e, ' are opened.'
      stop 999
    end if
    return
  end subroutine search_fileunit 
!
end module rsmcom_module
