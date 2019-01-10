program circuit
  !###################################################################!
  !###################################################################!
  !####                                                           ####!
  !####                PROGRAM CIRCUIT                            ####!
  !####                                                           ####!
  !####  Given a center (xc,yc,zc), a radius, and a               ####!
  !####  spacing between parcels (original_distance), creates a   ####!
  !####  circuit of parcels in the horizontal plane.  Then com-   ####!
  !####  putes trajectories of parcels.  Every time the parcels   ####!
  !####  get farther away from each other than the specified min- ####!
  !####  imum distance (minddist), a new parcel is added halfway  ####!
  !####  between them.  Circulation and related quantities are    ####!
  !####  computed.                                                ####!
  !####                                                           ####!
  !####  Based on code written by Paul Markowski.                 ####!
  !####===========================================================####!
  !#### Ryan Hastings, 15 Sep 2015                                ####!
  !####-----------------------------------------------------------####!
  !#### v1.0, updated comments and eliminated diabatic stuff,     ####!
  !#### Ryan Hastings, 16 Dec 2016                                ####!
  !###################################################################!
  !###################################################################!
  use globals
  use nc_inout
  use traj_utils
  implicit none
  !===================================================================!
  !###################################################################!
  !####################### VARIABLES #################################!
  !###################################################################!
  !===================================================================!
  !-------------------------------------------------------------------!
  ! configuration variables
  real :: begintime, endtime, xc, yc, zc, radius, original_distance, mindist
  integer :: maxparcels

  character(len=100),dimension(:),allocatable :: infiles, dynfiles
  character(len=100) :: cm1list, outfile, dynlist

  !-------------------------------------------------------------------!
  ! dimensional variables
  real :: tt, xx, yy, zz, dtt
  real,dimension(:),allocatable :: zzh
  real,dimension(:),allocatable :: times, tp

  integer ti, tstep

  integer nfiles, nparcels

  real :: circumference, incr, dl

  !-----------------------------------------------------------------------------!
  ! netcdf variables
  integer rcode, ncid_out, varid
  integer :: ncid
  integer tpid, xpid, ypid, zpid, upid, vpid, wpid
  integer qvpid, qcpid, qspid, qgpid, qhlpid, qipid, qrpid, thpid, thrpid, buoypid
  integer ncpid, nspid, ngpid, nipid, nrpid, rhodid, pipid
  integer pgfzpid, pgfxpid, pgfypid
  integer pgfxbpid, pgfxdlpid, pgfxdnpid
  integer pgfybpid, pgfydlpid, pgfydnpid
  integer pgfzbpid, pgfzdlpid, pgfzdnpid
  integer diabid, tkepid
  integer xvortpid, yvortpid, zvortpid, svortpid, cvortpid, psipid
  integer zvort_stretchid, zvort_tiltid
  integer svort_tiltid, svort_stretchid, svort_bclnid, svort_xchgid
  integer dpsidtid
  integer cvort_tiltid, cvort_stretchid, cvort_bclnid, cvort_xchgid
  integer circid, intBdzid, vdotl1id, vdotl2id

  !-----------------------------------------------------------------------------!
  ! field variables
  real,dimension(:,:,:),allocatable :: u, v, w, qv, qc, qi, qs, qr, qg, qhl, th, thr, rhod
  real,dimension(:,:,:),allocatable :: ncc, ncr, nci, ncs, ncg, ppi
  real,dimension(:,:,:),allocatable :: utem, vtem, wtem
  real,dimension(:,:,:),allocatable :: buoy, pgfz, pgfz_b, pgfz_dl, pgfz_dn
  real,dimension(:,:,:),allocatable :: pgfx, pgfx_b, pgfx_dl, pgfx_dn
  real,dimension(:,:,:),allocatable :: pgfy, pgfy_b, pgfy_dl, pgfy_dn
  real,dimension(:),allocatable :: pi0, th0, thr0, qv0
  real,dimension(:,:,:),allocatable :: tke

  real,dimension(:,:,:),allocatable :: xvort, yvort, zvort
  real,dimension(:,:,:),allocatable :: dpsi_dx, dpsi_dy, dpsi_dz
  real,dimension(:,:,:),allocatable :: dvhdx, dvhdy, dvhdz
  real,dimension(:,:,:),allocatable :: dbdx, dbdy
  real,dimension(:,:,:),allocatable :: dwdx, dwdy, dwdz


  !-----------------------------------------------------------------------------!
  ! parcel variables
  real,dimension(:,:),allocatable :: xpc, ypc, zpc
  real,dimension(:,:),allocatable :: up, vp, wp
  real,dimension(:,:),allocatable :: qvp, qcp, qrp, qsp, qip, qgp, qhlp
  real,dimension(:,:),allocatable :: ncp, nrp, nsp, nip, ngp
  real,dimension(:,:),allocatable :: thp, thrp, rhodp, pip
  real,dimension(:,:),allocatable :: buoyp, pgfzp, pgfxp, pgfyp
  real,dimension(:,:),allocatable :: pgfxbp, pgfxdlp, pgfxdnp
  real,dimension(:,:),allocatable :: pgfybp, pgfydlp, pgfydnp
  real,dimension(:,:),allocatable :: pgfzbp, pgfzdlp, pgfzdnp
  real,dimension(:,:),allocatable :: tkep
  real,dimension(:,:),allocatable :: xvortp, yvortp, zvortp, svortp, cvortp
  real,dimension(:,:),allocatable :: psip, zvort_stretchp, zvort_tiltp
  real,dimension(:,:),allocatable :: svort_tiltp, svort_stretchp, svort_bclnp
  real,dimension(:,:),allocatable :: svort_xchgp, dpsidt
  real,dimension(:,:),allocatable :: cvort_tiltp, cvort_stretchp, cvort_bclnp
  real,dimension(:,:),allocatable :: cvort_xchgp
  real,dimension(:),allocatable   :: circ, intBdz
  real,dimension(:,:),allocatable :: vdotl1, vdotl2

  !-----------------------------------------------------------------------------!
  ! looping and dummy variables
  integer i, t, np, p
  integer,dimension(1) :: tmp


  !-----------------------------------------------------------------------------!
  ! namelist

  namelist /filenames/  cm1list, dynlist, outfile
  namelist /parameters/ begintime, endtime, dt, xc, yc, zc, radius, original_distance, mindist, maxparcels, ptype
  namelist /outputs/    output_dyn, output_vor, output_z_only

  !=============================================================================!
  !#############################################################################!
  !########################### MAIN BODY #######################################!
  !#############################################################################!

  write(*,*)
  write(*,*) '                     CIRCUIT'
  write(*,*)

  !#############################################################################!
  !#                   INITIALIZATION                                          #!
  !#############################################################################!

  write(*,*) 'opening circuit.input'
  open(8,file='circuit.input')

  read(8,nml=filenames)
  write(*,*) 'cm1list =',cm1list ! file containing list of CM1 filenames
  write(*,*) 'dynlist =',dynlist ! file containing list of dyn filenames
  write(*,*) 'outfile =',outfile ! name of outfile

  read(8,nml=parameters)
  write(*,*) 'begintime         =',begintime ! time to regress parcels to
  write(*,*) 'endtime           =',endtime   ! time of parcels
  write(*,*) 'dt                =',dt ! time step
  write(*,*) 'xc                =',xc ! center of circuit, x-coordinate
  write(*,*) 'yc                =',yc ! center of circuit, y-coordinate
  write(*,*) 'zc                =',zc ! center of circuit, z-coordinate
  write(*,*) 'radius            =',radius ! radius of circuit
  write(*,*) 'original_distance =',original_distance ! starting distance of parcels
  write(*,*) 'minimum distance  =',mindist ! minimum distance for parcels
  write(*,*) 'maximum parcels   =',maxparcels ! maximum number of parcels
  write(*,*) 'ptype             =',ptype ! precipitation type

  read(8,nml=outputs)
  write(*,*) 'output_dyn   =',output_dyn ! output dynamics?
  write(*,*) 'output_vor   =',output_vor ! output vorticity dynamics?
  write(*,*) 'output_z_only=',output_z_only ! output z stuff only?

  close(8)

  !----------------------------------------------------------------------------!
  ! set up arrays for times                                                    !
  !----------------------------------------------------------------------------!
  ! set up array for times of parcels
  write(*,*)
  write(*,*) 'setting up arrays for parcel times'
  ntimes = int(1+(endtime-begintime)/dt)

  allocate( tp(ntimes) )
  do t=1,ntimes
    tp(t) = begintime + (t-1)*dt
    write(*,*) 't,tp(t)=',t,tp(t)
  enddo

  !----------------------------------------------------------------------------!
  ! set up arrays for times of files

  ! read infilelist
  write(*,*) 'reading list of cm1 files:'
  open(8,file=cm1list)

  read(8,*) nfiles
  allocate( infiles(nfiles) )
  do i=1,nfiles
    read(8,*) infiles(i)
    write(*,*) infiles(i)
  enddo
  close(8)

  write(*,*)
  ! read times
  allocate( times(nfiles) )
  write(*,*) 'reading times:'
  do i=1,nfiles
    call cm1_nc_gettime_traj( infiles(i), times(i) )
    write(*,*) infiles(i),times(i)
  enddo

  write(*,*)
  ! read dynlist
  if(output_dyn.eq.1)then
    write(*,*) 'reading list of dynamics files:'
    open(8,file=dynlist)

    allocate( dynfiles(nfiles) )
    do i=1,nfiles
      read(8,*) dynfiles(i)
      write(*,*) dynfiles(i)
    enddo
    close(8)
  endif

  !------------------------------------------------------------------------------!
  ! get dimensions and set up arrays                                             !
  !------------------------------------------------------------------------------!
  ! get dimensions

  write(*,*) 'reading dimensions..'
  ncid = ncopn( infiles(1), NCNOWRIT, rcode )

  call cm1_nc_getdims( ncid )

  !------------------------------------------------------------------------------!
  ! set up spacial arrays

  write(*,*) 'reading space..'
  allocate( xh(ni) )
  allocate( xf(nip1) )
  allocate( yh(nj) )
  allocate( yf(njp1) )
  allocate( zh(nk) )
  allocate( zf(nkp1) )

  ! we will be assuming free-slip conditions, so the horizontal winds will remain
  ! constant below the lowest grid point.  zzh(0) represents the surface for
  ! mass grid points, then.  
  allocate( zzh(0:nk) )

  allocate( mzh(1:nk) )
  allocate( mzf(1:nk+1) )

  write(*,*) 'reading space-files..'
  call cm1_nc_getspace( ncid )

  dz = xh(2)-xh(1)
  rdz = 1.0/dz

  write(*,*) 'setting up stretching arrays'
  do i=2,nk
    mzh(i) = dz/(zf(i+1)-zf(i))
  enddo
  mzh(1)=1.0
  do i=2,nk
    mzf(i) = dz/(zh(i)-zh(i-1))
  enddo
  mzf(1)=1.0
  mzf(nk+1) = 1.0

  zzh(1:nk) = zh(1:nk)

  zzh(0)    = -zzh(1)
  write(*,*) zzh(nk)

  !-----------------------------------------------------------------------!
  ! Read in variables

  write(*,*) 'reading in base state'
  allocate( th0(1:nk) )
  allocate( qv0(1:nk) )
  allocate( pi0(1:nk) )
  allocate( thr0(1:nk) )

  write(*,*) '..qv0'
  varid =  ncvid( ncid, 'qv0', rcode )
  call ncvgt( ncid, varid, (/1,1,1,1/), (/1,1,nk,1/), qv0, rcode )

  write(*,*) '..th0'
  varid =  ncvid( ncid, 'th0', rcode )
  call ncvgt( ncid, varid, (/1,1,1,1/), (/1,1,nk,1/), th0, rcode )

  thr0=th0*(1+qv0*reps)/(1+qv0)

  write(*,*) '..pi0'
  varid = ncvid( ncid, 'pi0', rcode )
  call ncvgt( ncid ,varid, (/1,1,1,1/), (/1,1,nk,1/), pi0, rcode )

  call ncclos( ncid, rcode )

  !---------------------------------------------------------------------------!
  ! set up other arrays, now that we have dimensions                          !
  !---------------------------------------------------------------------------!
  write(*,*) 'setting up some wind arrays'
  allocate( u(1:ni+1,1:nj,0:nk) )
  allocate( v(1:ni,1:nj+1,0:nk) )
  allocate( w(1:ni,1:nj,1:nk+1) )
  allocate( utem(1:ni+1,1:nj,0:nk) )
  allocate( vtem(1:ni,1:nj+1,0:nk) )
  allocate( wtem(1:ni,1:nj,1:nk+1) )

  !------------------------------------------!
  ! variable for all scalars

  allocate( qv(ni,nj,0:nk) )
  allocate( qc(ni,nj,0:nk) )
  allocate( ncc(ni,nj,0:nk) )
  allocate( qr(ni,nj,0:nk) )
  allocate( ncr(ni,nj,0:nk) )
  allocate( qi(ni,nj,0:nk) )
  allocate( nci(ni,nj,0:nk) )
  allocate( qs(ni,nj,0:nk) )
  allocate( ncs(ni,nj,0:nk) )
  allocate( qg(ni,nj,0:nk) )
  allocate( qhl(ni,nj,0:nk) )
  allocate( ncg(ni,nj,0:nk) )
  allocate( th(ni,nj,0:nk) )
  allocate( ppi(ni,nj,0:nk) )

  allocate( thr(ni,nj,0:nk) )
  allocate( rhod(ni,nj,0:nk) )

  allocate( buoy(ni,nj,0:nk) )
  allocate( pgfx(ni+1,nj,0:nk) )
  allocate( pgfy(ni,nj+1,0:nk) )
  allocate( pgfz(ni,nj,1:nk+1) )

  allocate( pgfx_b(1:ni+1,1:nj,0:nk) )
  allocate( pgfy_b(1:ni,1:nj+1,0:nk) )
  allocate( pgfz_b(1:ni,1:nj,1:nk+1) )

  allocate( pgfx_dl(1:ni+1,1:nj,0:nk) )
  allocate( pgfy_dl(1:ni,1:nj+1,0:nk) )
  allocate( pgfz_dl(1:ni,1:nj,1:nk+1) )

  allocate( pgfx_dn(1:ni+1,1:nj,1:nk) )
  allocate( pgfy_dn(1:ni,1:nj+1,1:nk) )
  allocate( pgfz_dn(1:ni,1:nj,1:nk+1) )

  allocate( tke(1:ni,1:nj,1:nk+1) )

  allocate( xvort(1:ni,1:nj+1,1:nk+1) )
  allocate( yvort(1:ni+1,1:nj,1:nk+1) )
  allocate( zvort(1:ni+1,1:nj+1,0:nk) )

  allocate( dpsi_dx(1:ni+1,1:nj,0:nk) )
  allocate( dpsi_dy(1:ni,1:nj+1,0:nk) )
  allocate( dpsi_dz(1:ni,1:nj,1:nk+1) )

  allocate( dvhdx(1:ni+1,1:nj,0:nk) )
  allocate( dvhdy(1:ni,1:nj+1,0:nk) )
  allocate( dvhdz(1:ni,1:nj,1:nk+1) )

  allocate( dbdx(1:ni+1,1:nj,0:nk) )
  allocate( dbdy(1:ni,1:nj+1,0:nk) )

  allocate( dwdx(1:ni+1,1:nj,1:nk+1) )
  allocate( dwdy(1:ni,1:nj+1,1:nk+1) )
  allocate( dwdz(1:ni,1:nj,0:nk) )


  !-----------------------------------------------------------!
  ! set up trajectory arrays
  write(*,*) 'setting up trajectory arrays'
  allocate( xpc(ntimes,maxparcels) )
  allocate( ypc(ntimes,maxparcels) )
  allocate( zpc(ntimes,maxparcels) )

  xpc = 9.9e9 ! assign all variables to a missing value
  ypc = 9.9e9
  zpc = 9.9e9

  allocate( up(ntimes,maxparcels) )
  allocate( vp(ntimes,maxparcels) )
  allocate( wp(ntimes,maxparcels) )

  up = 9.9e9
  vp = 9.9e9
  wp = 9.9e9

  allocate( qvp(ntimes,maxparcels) )
  allocate( qcp(ntimes,maxparcels) )
  allocate( ncp(ntimes,maxparcels) )
  allocate( qrp(ntimes,maxparcels) )
  allocate( nrp(ntimes,maxparcels) )
  allocate( qip(ntimes,maxparcels) )
  allocate( nip(ntimes,maxparcels) )
  allocate( qsp(ntimes,maxparcels) )
  allocate( nsp(ntimes,maxparcels) )
  allocate( qgp(ntimes,maxparcels) )
  allocate( qhlp(ntimes,maxparcels) )
  allocate( ngp(ntimes,maxparcels) )
  allocate( thp(ntimes,maxparcels) )
  allocate( pip(ntimes,maxparcels) )

  qvp = 9.9e9
  qcp = 9.9e9
  ncp = 9.9e9
  qrp = 9.9e9
  nrp = 9.9e9
  qip = 9.9e9
  nip = 9.9e9
  qsp = 9.9e9
  nsp = 9.9e9
  qgp = 9.9e9
  qhlp = 9.9e9
  ngp = 9.9e9
  thp = 9.9e9
  pip = 9.9e9

  allocate( thrp(ntimes,maxparcels) )
  allocate( rhodp(ntimes,maxparcels) )

  thrp = 9.9e9
  rhodp = 9.9e9

  allocate( buoyp(ntimes,maxparcels) )

  buoyp = 9.9e9

  allocate( pgfxp(ntimes,maxparcels) )
  allocate( pgfyp(ntimes,maxparcels) )
  allocate( pgfzp(ntimes,maxparcels) )

  pgfxp = 9.9e9
  pgfyp = 9.9e9
  pgfzp = 9.9e9

  allocate( pgfxbp(ntimes,maxparcels) )
  allocate( pgfybp(ntimes,maxparcels) )
  allocate( pgfzbp(ntimes,maxparcels) )

  pgfxbp = 9.9e9
  pgfybp = 9.9e9
  pgfzbp = 9.9e9

  allocate( pgfxdlp(ntimes,maxparcels) )
  allocate( pgfydlp(ntimes,maxparcels) )
  allocate( pgfzdlp(ntimes,maxparcels) )

  pgfxdlp = 9.9e9
  pgfydlp = 9.9e9
  pgfzdlp = 9.9e9

  allocate( pgfxdnp(ntimes,maxparcels) )
  allocate( pgfydnp(ntimes,maxparcels) )
  allocate( pgfzdnp(ntimes,maxparcels) )

  pgfxdnp = 9.9e9
  pgfydnp = 9.9e9
  pgfzdnp = 9.9e9

  allocate( tkep( ntimes, maxparcels) )

  tkep = 9.9e9

  allocate( xvortp(ntimes,maxparcels) )
  allocate( yvortp(ntimes,maxparcels) )
  allocate( zvortp(ntimes,maxparcels) )
  allocate( svortp(ntimes,maxparcels) )
  allocate( cvortp(ntimes,maxparcels) )

  xvortp = 9.9e9
  yvortp = 9.9e9
  zvortp = 9.9e9

  allocate( psip(ntimes,maxparcels) )

  psip = 9.9e9

  allocate( zvort_stretchp(ntimes,maxparcels) )
  allocate( zvort_tiltp(ntimes,maxparcels) )

  zvort_stretchp = 9.9e9
  zvort_tiltp = 9.9e9

  allocate( svort_tiltp(ntimes,maxparcels) )
  allocate( svort_stretchp(ntimes,maxparcels) )
  allocate( svort_xchgp(ntimes,maxparcels) )
  allocate( svort_bclnp(ntimes,maxparcels) )

  svort_tiltp = 9.9e9
  svort_stretchp = 9.9e9
  svort_xchgp = 9.9e9
  svort_bclnp = 9.9e9

  allocate( dpsidt(ntimes,maxparcels) )

  dpsidt = 9.9e9

  allocate( cvort_tiltp(ntimes,maxparcels) )
  allocate( cvort_stretchp(ntimes,maxparcels) )
  allocate( cvort_xchgp(ntimes,maxparcels) )
  allocate( cvort_bclnp(ntimes,maxparcels) )

  cvort_tiltp = 9.9e9
  cvort_stretchp = 9.9e9
  cvort_xchgp = 9.9e9
  cvort_bclnp = 9.9e9

  allocate( circ(ntimes) )
  allocate( intBdz(ntimes) )
  allocate( vdotl1(ntimes,maxparcels) )
  allocate( vdotl2(ntimes,maxparcels) )
  vdotl1=9.9e9
  vdotl2=9.9e9

  ! set up parcels
  circumference = 2*pi*radius
  nparcels = int( circumference/original_distance )
  incr = 360./float(nparcels)

  do np=1,nparcels
    xpc(ntimes,np) = xc + radius*sin(incr*float(np)*pi/180)
    ypc(ntimes,np) = yc + radius*cos(incr*float(np)*pi/180)
    zpc(ntimes,np) = zc
  enddo

  !================================================================!
  !################################################################!
  !############ COMPUTE TRAJECTORIES ##############################!
  !################################################################!
  !================================================================!

  write(*,*) '####################################################'
  write(*,*) 'starting backward trajectory calculations'

  dt=-dt

  !------------- set up initial times -----------------------------!
  tmp = minloc( abs( endtime - times ) )
  ti = tmp(1)

  write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

  call traject_nc_read_cm1( infiles(ti), u, v, w, qv, qc, qr, qi, qs, qg, qhl, th, thr0, pi0, &
    pgfx, pgfy, pgfz, thr, buoy,&
    tke, ncc, ncr, nci, ncs, ncg, rhod, ppi )

  if (output_dyn.eq.1) then
    call traject_nc_read_dyn( dynfiles(ti), pgfx_b, pgfy_b, pgfz_b, pgfx_dl, pgfy_dl, pgfz_dl, &
      pgfx_dn, pgfy_dn, pgfz_dn, thr, buoy )
  endif



  do tstep=ntimes,2,-1

    tt=tp(tstep)
    write(*,*) '--------------------------------------------------------------'
    write(*,*) 'time step ',tstep,' or ',tt,'s'

    if(output_vor.eq.1)then
      call compute_vorts( buoy, u, v, w, xvort, yvort, zvort, &
        dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
        dwdy, dwdx, dwdz )
    endif

    do np=1,nparcels-1

      ! compute distance between parcels
      dl = sqrt( (xpc(tstep,np+1)-xpc(tstep,np))**2 + &
                 (ypc(tstep,np+1)-ypc(tstep,np))**2 + &
                 (zpc(tstep,np+1)-zpc(tstep,np))**2 )

      if (dl.gt.mindist) then ! parcel distance exceeds mindist
        if(nparcels.eq.maxparcels)then
          write(*,*)
          write(*,*) 'Parcels needed to be added to the material circuit, but MAXPARCELS was exceeded'
          write(*,*) 'Program stopped'
          stop
        else
          nparcels = nparcels + 1
          write(*,*) 'adding parcel'
        endif

        do p=nparcels,np+2,-1 ! shift all parcels over
          xpc(1:tstep,p) = xpc(1:tstep,p-1)
          ypc(1:tstep,p) = ypc(1:tstep,p-1)
          zpc(1:tstep,p) = zpc(1:tstep,p-1)
          up(1:tstep,p)  = up(1:tstep,p-1)
          vp(1:tstep,p)  = vp(1:tstep,p-1)
          wp(1:tstep,p)  = wp(1:tstep,p-1)
          qvp(1:tstep,p) = qvp(1:tstep,p-1)
          qcp(1:tstep,p) = qcp(1:tstep,p-1)
          ncp(1:tstep,p) = ncp(1:tstep,p-1)
          qrp(1:tstep,p) = qrp(1:tstep,p-1)
          nrp(1:tstep,p) = nrp(1:tstep,p-1)
          qip(1:tstep,p) = qip(1:tstep,p-1)
          nip(1:tstep,p) = nip(1:tstep,p-1)
          qsp(1:tstep,p) = qsp(1:tstep,p-1)
          nsp(1:tstep,p) = nsp(1:tstep,p-1)
          qgp(1:tstep,p) = qgp(1:tstep,p-1)
          qhlp(1:tstep,p) = qhlp(1:tstep,p-1)
          ngp(1:tstep,p) = ngp(1:tstep,p-1)
          thp(1:tstep,p) = thp(1:tstep,p-1)
          pip(1:tstep,p) = pip(1:tstep,p-1)
          thrp(1:tstep,p) = thrp(1:tstep,p-1)
          rhodp(1:tstep,p) = rhodp(1:tstep,p-1)
          buoyp(1:tstep,p) = buoyp(1:tstep,p-1)
          pgfxp(1:tstep,p) = pgfxp(1:tstep,p-1)
          pgfyp(1:tstep,p) = pgfyp(1:tstep,p-1)
          pgfzp(1:tstep,p) = pgfzp(1:tstep,p-1)
          pgfxbp(1:tstep,p) = pgfxbp(1:tstep,p-1)
          pgfybp(1:tstep,p) = pgfybp(1:tstep,p-1)
          pgfzbp(1:tstep,p) = pgfzbp(1:tstep,p-1)
          pgfxdlp(1:tstep,p) = pgfxdlp(1:tstep,p-1)
          pgfydlp(1:tstep,p) = pgfydlp(1:tstep,p-1)
          pgfzdlp(1:tstep,p) = pgfzdlp(1:tstep,p-1)
          pgfxdnp(1:tstep,p) = pgfxdnp(1:tstep,p-1)
          pgfydnp(1:tstep,p) = pgfydnp(1:tstep,p-1)
          pgfzdnp(1:tstep,p) = pgfzdnp(1:tstep,p-1)
          tkep(1:tstep,p) = tkep(1:tstep,p-1)
          xvortp(1:tstep,p) = xvortp(1:tstep,p-1)
          yvortp(1:tstep,p) = yvortp(1:tstep,p-1)
          zvortp(1:tstep,p) = zvortp(1:tstep,p-1)
          svortp(1:tstep,p) = svortp(1:tstep,p-1)
          cvortp(1:tstep,p) = cvortp(1:tstep,p-1)
          psip(1:tstep,p) = psip(1:tstep,p-1)
          zvort_stretchp(1:tstep,p) = zvort_stretchp(1:tstep,p-1)
          zvort_tiltp(1:tstep,p) = zvort_tiltp(1:tstep,p-1)
          svort_tiltp(1:tstep,p) = svort_tiltp(1:tstep,p-1)
          svort_stretchp(1:tstep,p) = svort_stretchp(1:tstep,p-1)
          svort_xchgp(1:tstep,p) = svort_xchgp(1:tstep,p-1)
          svort_bclnp(1:tstep,p) = svort_bclnp(1:tstep,p-1)
          cvort_tiltp(1:tstep,p) = cvort_tiltp(1:tstep,p-1)
          cvort_stretchp(1:tstep,p) = cvort_stretchp(1:tstep,p-1)
          cvort_xchgp(1:tstep,p) = cvort_xchgp(1:tstep,p-1)
          cvort_bclnp(1:tstep,p) = cvort_bclnp(1:tstep,p-1)
        enddo

        ! define coordinates of new parcel
        xpc(tstep,np) = 0.5*(xpc(tstep,np)+xpc(tstep,np+2))
        ypc(tstep,np) = 0.5*(ypc(tstep,np)+ypc(tstep,np+2))
        zpc(tstep,np) = 0.5*(zpc(tstep,np)+zpc(tstep,np+2))

      endif

    enddo

    dl = sqrt( (xpc(tstep,1)-xpc(tstep,nparcels))**2 + &
               (ypc(tstep,1)-ypc(tstep,nparcels))**2 + &
               (zpc(tstep,1)-zpc(tstep,nparcels))**2 )

    if (dl.gt.mindist) then
      if(nparcels.eq.maxparcels)then
        write(*,*)
        write(*,*) 'Parcels needed to be added to the material circuit, but MAXPARCELS was exceeded'
        write(*,*) 'Program stopped'
        stop
      else
        nparcels = nparcels + 1
        write(*,*) 'adding parcel'
      endif

      xpc(tstep,nparcels) = 0.5*(xpc(tstep,nparcels)+xpc(tstep,1))
      ypc(tstep,nparcels) = 0.5*(ypc(tstep,nparcels)+ypc(tstep,1))
      zpc(tstep,nparcels) = 0.5*(zpc(tstep,nparcels)+zpc(tstep,1))
    endif

    do np=1,nparcels

      call interp_all( xpc(tstep,np), ypc(tstep,np), zpc(tstep,np), zzh, u, up(tstep,np), v, vp(tstep,np), &
        w, wp(tstep,np), buoy, buoyp(tstep,np), pgfx, pgfxp(tstep,np), pgfy, pgfyp(tstep,np), &
        pgfz, pgfzp(tstep,np), &
        tke, tkep(tstep,np), thr0, thr, thrp(tstep,np), qv, qvp(tstep,np), qc, qcp(tstep,np), &
        qr, qrp(tstep,np), qi, qip(tstep,np), qs, qsp(tstep,np), qg, qgp(tstep,np), qhl, qhlp(tstep,np), &
        ncc, ncp(tstep,np), ncr, nrp(tstep,np), nci, nip(tstep,np), &
        ncs, nsp(tstep,np), ncg, ngp(tstep,np), rhod, rhodp(tstep,np), &
        th, thp(tstep,np), ppi, pip(tstep,np) )      

      if(output_dyn.eq.1)then
        call interp_dynz( xpc(tstep,np), ypc(tstep,np), zpc(tstep,np), &
          pgfx_b, pgfxbp(tstep,np), pgfy_b, pgfybp(tstep,np), pgfz_b, pgfzbp(tstep,np), &
          pgfx_dl, pgfxdlp(tstep,np), pgfy_dl, pgfydlp(tstep,np), pgfz_dl, pgfzdlp(tstep,np), &
          pgfx_dn, pgfxdnp(tstep,np), pgfz_dn, pgfzdnp(tstep,np), pgfz_dn, pgfzdnp(tstep,np) )
      endif

      if(output_vor.eq.1)then
        call interp_vort3( xpc(tstep,np), ypc(tstep,np), &
          zpc(tstep,np), up(tstep,np), vp(tstep,np), wp(tstep,np), &
          xvort, yvort, zvort, dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, &
          dpsi_dz, dvhdz, dwdy, dwdx, dwdz, &
          xvortp(tstep,np), yvortp(tstep,np), zvortp(tstep,np), &
          svortp(tstep,np), cvortp(tstep,np), psip(tstep,np), &
          zvort_tiltp(tstep,np), zvort_stretchp(tstep,np), &
          svort_tiltp(tstep,np), svort_stretchp(tstep,np), &
          svort_bclnp(tstep,np), cvort_stretchp(tstep,np), &
          cvort_tiltp(tstep,np), cvort_bclnp(tstep,np) )
      endif

    enddo

    call compute_circ( nparcels, xpc(tstep,:), ypc(tstep,:), zpc(tstep,:), up(tstep,:), vp(tstep,:), wp(tstep,:), &
      circ(tstep), vdotl1(tstep,:), vdotl2(tstep,:) )
    call compute_intBdz( nparcels, zpc(tstep,:), buoyp(tstep,:), intBdz(tstep) )

    ti=ti-1
    write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

    call traject_nc_read_cm1( infiles(ti), utem, vtem, wtem, qv, qc, qr, qi, qs, qg, qhl, th, &
      thr0, pi0, pgfx, pgfy, pgfz, thr, buoy, tke, ncc, ncr, nci, ncs, ncg, &
      rhod, ppi )

    do np=1,nparcels
      call rk4( xpc(tstep-1:tstep,np), ypc(tstep-1:tstep,np), zpc(tstep-1:tstep,np), &
        up(tstep,np), vp(tstep,np), wp(tstep,np), u, v, w, utem, vtem, wtem, dt )
    enddo

    u = utem
    v = vtem
    w = wtem

    if (output_dyn.eq.1) then
      call traject_nc_read_dyn( dynfiles(ti), pgfx_b, pgfy_b, pgfz_b, pgfx_dl, pgfy_dl, pgfz_dl, &
        pgfx_dn, pgfy_dn, pgfz_dn, thr, buoy )
    endif

  enddo ! do tstep=start_time,2,-1

  if(output_vor.eq.1)then
    call compute_vorts( buoy, u, v, w, xvort, yvort, zvort, &
      dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
      dwdy, dwdx, dwdz )
  endif

  do np=1,nparcels

    call interp_all( xpc(1,np), ypc(1,np), zpc(1,np), zzh, u, up(1,np), v, vp(1,np), &
      w, wp(1,np), buoy, buoyp(1,np), pgfx, pgfxp(1,np), pgfy, pgfyp(1,np), pgfz, pgfzp(1,np), &
      tke, tkep(1,np), thr0, thr, thrp(1,np), qv, qvp(1,np), qc, qcp(1,np), &
      qr, qrp(1,np), qi, qip(1,np), qs, qsp(1,np), qg, qgp(1,np), qhl, qhlp(1,np), &
      ncc, ncp(1,np), ncr, nrp(1,np), nci, nip(1,np), &
      ncs, nsp(1,np), ncg, ngp(1,np), rhod, rhodp(1,np), th, thp(1,np), ppi, pip(1,np) )

    if(output_dyn.eq.1)then
      call interp_dynz( xpc(1,np), ypc(1,np), zpc(1,np), &
        pgfx_b, pgfxbp(1,np), pgfy_b, pgfybp(1,np), pgfz_b, pgfzbp(1,np), &
        pgfx_dl, pgfxdlp(1,np), pgfy_dl, pgfydlp(1,np), pgfz_dl, pgfzdlp(1,np), &
        pgfx_dn, pgfxdnp(1,np), pgfy_dn, pgfydnp(1,np), pgfz_dn, pgfzdnp(1,np) )
    endif

    if(output_vor.eq.1)then
      call interp_vort3( xpc(1,np), ypc(1,np), zpc(1,np), up(1,np), &
        vp(1,np), wp(1,np), &
        xvort, yvort, zvort, dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, &
        dpsi_dz, dvhdz, dwdy, dwdx, dwdz, &
        xvortp(1,np), yvortp(1,np), zvortp(1,np), &
        svortp(1,np), cvortp(1,np), psip(1,np), zvort_tiltp(1,np), &
        zvort_stretchp(1,np), svort_tiltp(1,np), svort_stretchp(1,np), &
        svort_bclnp(1,np), cvort_stretchp(1,np), cvort_tiltp(1,np), &
        cvort_bclnp(1,np) )
    endif

  enddo

  call compute_circ( nparcels, xpc(1,:), ypc(1,:), zpc(1,:), up(1,:), vp(1,:), wp(1,:), &
    circ(1), vdotl1(1,:), vdotl2(1,:) )
  call compute_intBdz( nparcels, zpc(1,:), buoyp(1,:), intBdz(1) )


  if(output_vor.eq.1)then
    !---------------------------------------------------------------!
    ! compute exchange terms for horizontal vorticity               !
    !---------------------------------------------------------------!
      dpsidt(:,:) = 0.0
      do np=1,nparcels
        do tstep=2,ntimes-1
          dpsidt(tstep,np) = compute_dpsi2( psip(tstep-1,np), psip(tstep,np), &
            psip(tstep+1,np), 2*dt )
        enddo
     enddo

     svort_xchgp =  cvortp*dpsidt
     cvort_xchgp = -svortp*dpsidt

  endif

  !##############################################################################!
  !#                                WRITE OUTPUT                                #!
  !##############################################################################!

  write(*,*) '###################################################################'
  write(*,*) '# write out netcdf file'
  write(*,*)

  !----------------------------------------------------------!
  ! create netcdf file
  call nc_circuit_vardef(  outfile, ncid_out, ntimes, nparcels, &
      tpid, xpid, ypid, zpid, upid, vpid, wpid, qvpid, qcpid, qrpid, qipid, &
      qspid, qgpid, qhlpid, ncpid, nrpid, nipid, nspid, ngpid, &
      thpid, thrpid, buoypid, pgfxpid, pgfypid, pgfzpid, &
      pgfxbpid, pgfybpid, pgfzbpid, pgfxdlpid, pgfydlpid, pgfzdlpid, &
      pgfxdnpid, pgfydnpid, pgfzdnpid, &
      tkepid, rhodid, pipid, xvortpid, yvortpid, zvortpid, &
      svortpid, cvortpid, psipid, zvort_stretchid, zvort_tiltid,&
      svort_tiltid, svort_stretchid, svort_bclnid, svort_xchgid, &
      dpsidtid, cvort_tiltid, cvort_stretchid, cvort_bclnid, &
      cvort_xchgid, circid, vdotl1id, vdotl2id, intBdzid )!

  call nc_circuit_writeout( ntimes, nparcels, ncid_out, tpid, tp, xpid, xpc, ypid, ypc, zpid, zpc, &
      upid, up, vpid, vp, wpid, wp, buoypid, buoyp, &
      pgfxpid, pgfxp, pgfypid, pgfyp, pgfzpid, pgfzp, &
      pgfxbpid, pgfxbp, pgfybpid, pgfybp, pgfzbpid, pgfzbp, &
      pgfxdlpid, pgfxdlp, pgfydlpid, pgfydlp, pgfzdlpid, pgfzdlp, &
      pgfxdnpid, pgfxdnp, pgfydnpid, pgfydnp, pgfzdnpid, pgfzdnp, &
      tkepid, tkep, thrpid, thrp, &
      qvpid, qvp, qcpid, qcp, qrpid, qrp, qipid, qip, qspid, qsp, qgpid, qgp, qhlpid, qhlp, &
      ncpid, ncp, nrpid, nrp, nipid, nip, nspid, nsp, ngpid, ngp, rhodid, rhodp, &
      thpid, thp, pipid, pip, xvortpid, xvortp, yvortpid, yvortp, zvortpid, zvortp, &
      svortpid, svortp, cvortpid, cvortp, psipid, psip, zvort_tiltid, &
      zvort_tiltp, zvort_stretchid, zvort_stretchp, svort_tiltid, &
      svort_tiltp, svort_stretchid, svort_stretchp, svort_bclnid, &
      svort_bclnp, svort_xchgid, svort_xchgp, dpsidtid, dpsidt, &
      cvort_tiltid, cvort_tiltp, cvort_stretchid, cvort_stretchp, &
      cvort_bclnid, cvort_bclnp, cvort_xchgid, cvort_xchgp, &
      circid, circ, vdotl1id, vdotl1, vdotl2id, vdotl2, intBdzid, intbdz )




  write(*,*) 'done writing out netcdf file'


  write(*,*) 'program terminated normally'

contains
  !##################################################################!
  !#                    SUBROUTINE COMPUTE_CIRC                     #!
  !#                                                                #!
  !#  Compute circulation around circuit.                           #!
  !#  Passed variables:                                             #!
  !#    integer nparcels : number of parcels                        #!
  !#    real xpc, ypc, zpc : parcel positions                       #!
  !#    real up, vp, wp : winds at parcel positions                 #!
  !#  Returned variables:                                           #!
  !#    real circ : circulation                                     #!
  !#    real vdotl1 : dot product between wind and distance         #!
  !#    real vdotl2 : vdotl1 normalized by distance                 #!
  !#================================================================#!
  !#  v1.0, Ryan Hastings, 16 Dec 2016                              #!
  !##################################################################!
  subroutine compute_circ( nparcels, xpc, ypc, zpc, up, vp, wp, circ, vdotl1, vdotl2 )

    implicit none
    integer nparcels
    real,dimension(:) :: xpc, ypc, zpc, up, vp, wp, vdotl1, vdotl2
    real :: circ

    integer np
    real,dimension(nparcels) :: uavg, vavg, wavg, deltax, deltay, deltaz

    circ = 0.0

    do np=1,nparcels
      if ( (up(np).eq.9.9e9).or.(vp(np).eq.9.9e9).or.(wp(np).eq.9.9e9).or.(xpc(np).eq.9.9e9).or. &
           (ypc(np).eq.9.9e9).or.(zpc(np).eq.9.9e9) )then
        circ = 9.9e9
        vdotl1(np)=9.9e9
        vdotl2(np)=9.9e9
      endif
    enddo

    do np=1,nparcels-1
      uavg(np) = 0.5*(up(np)+up(np+1))
      vavg(np) = 0.5*(vp(np)+vp(np+1))
      wavg(np) = 0.5*(wp(np)+wp(np+1))
      deltax(np) = xpc(np+1)-xpc(np)
      deltay(np) = ypc(np+1)-ypc(np)
      deltaz(np) = zpc(np+1)-zpc(np)
    enddo
    uavg(nparcels) = 0.5*(up(1)+up(nparcels))
    vavg(nparcels) = 0.5*(vp(1)+vp(nparcels))
    wavg(nparcels) = 0.5*(wp(1)+wp(nparcels))
    deltax(nparcels) = xpc(1)-xpc(nparcels)
    deltay(nparcels) = ypc(1)-ypc(nparcels)
    deltaz(nparcels) = zpc(1)-zpc(nparcels)

    do np=1,nparcels
      circ = circ - uavg(np)*deltax(np) - vavg(np)*deltay(np) - wavg(np)*deltaz(np)
      vdotl1(np) = -uavg(np)*deltax(np) - vavg(np)*deltay(np) - wavg(np)*deltaz(np)
      vdotl2(np) = vdotl1(np)/sqrt( deltax(np)**2 + deltay(np)**2 + deltaz(np)**2 )
    enddo

end subroutine

!##################################################################!
!#              SUBROUTINE COMPUTE_INTBDZ                         #!
!#                                                                #!
!#  Compute circuit integral of buoyancy with height              #!
!#  Passed variables:                                             #!
!#    integer nparcels : number of parcels                        #!
!#    real zpc : height of parcels                                #!
!#    real buoyp : parcel buoyancy                                #!
!#  Returned variable:                                            #!
!#    real intBdz : integral of buoyancy with height              #!
!#================================================================#!
!#  v1.0, Ryan Hastings, 16 Dec 2016                              #!
!##################################################################!
subroutine compute_intBdz( nparcels, zpc, buoyp, intBdz )
  implicit none
  integer nparcels
  real,dimension(:) :: zpc, buoyp
  real intBdz

  real,dimension(nparcels) :: bavg, deltaz
  integer np

  do np=1,nparcels
    if ( (buoyp(np).eq.9.9e9).or.(zpc(np).eq.9.9e9) ) then
      intBdz = 9.9e9
    endif
  enddo
  do np=1,nparcels-1
    bavg(np) = 0.5*(buoyp(np)+buoyp(np+1))
    deltaz(np) = zpc(np+1)-zpc(np)
  enddo
  bavg(nparcels)=0.5*(buoyp(nparcels)+buoyp(1))
  deltaz(nparcels)=zpc(1)-zpc(nparcels)

  intBdz = 0.0
  do np=1,nparcels
    intBdz = intBdz - deltaz(np)*bavg(np)
  enddo

end subroutine

end program
