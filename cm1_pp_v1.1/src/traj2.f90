program traj
  !#############################################################################!
  !#############################################################################!
  !####                                                                     ####!
  !####                         PROGRAM TRAJ                                ####!
  !####                                                                     ####!
  !####  This program computes trajectories of parcels from specified       ####!
  !####  starting positions.  The CM1 output must be netCDF with each his-  ####!
  !####  tory dump in a different file.  The time-step variable dt must     ####!
  !####  equal the difference in times between the history dumps.  That is, ####!
  !####  no fourth-dimensional interpolated is used to find the parcel      ####!
  !####  values between history files.                                      ####!
  !####                                                                     ####!
  !####  Note that this does compute vorticity dynamics, but use at your    ####!
  !####  own risk, especially the horizontal vorticity.                     ####!
  !####=====================================================================####!
  !#### Ryan Hastings, 8 Sep 2011                                           ####!
  !####---------------------------------------------------------------------####!
  !#### v1.0 completed on 10 December 2016                                  ####!
  !#############################################################################!
  !#############################################################################!
  use globals
  use nc_inout
  use traj_utils
  implicit none
  !=============================================================================!
  !#############################################################################!
  !################################# VARIABLES #################################!
  !#############################################################################!
  !=============================================================================!
  !-----------------------------------------------------------------------------!
  ! configuration variables
  real                                        :: begintime ! time to integrate back to
  real                                        :: parcelt ! initialization time
  real                                        :: endtime ! time to integrate forward to

  character(len=100),dimension(:),allocatable :: infiles ! array to hold list of files
  character(len=100),dimension(:),allocatable :: dynfiles ! array to hold list of dyn files
  character(len=100)                          :: cm1list ! name of file containing list of
                                                         ! CM1 output files
  character(len=100)                          :: plist ! file containing list of parcel
                                                       ! initialization positions
  character(len=100)                          :: dynlist ! name of file containing list of
  character(len=100)                          :: outfile ! output filename

  integer :: integration_method ! 0=euler, 1=runge-kutta

  real :: file_timestep ! time step of the files...dt will now be the time step of the parcel advance

  !-----------------------------------------------------------------------------!
  ! dimensional variables
  integer                       :: nftimes ! number of forward integration steps
  integer                       :: nbtimes ! number of backward integration steps
  integer                       :: start_time ! index of parcel initialization time
  integer                       :: nfiles ! number of files in infiles 
  integer                       :: nparcels ! number of parcels
  integer                       :: np ! looping variable

  real                          :: tt ! time step
  real,dimension(:),allocatable :: zzh ! scalar height array (see below for explanation)

  real,dimension(:),allocatable :: times ! times of files in infiles
  real,dimension(:),allocatable :: tp ! integration times
  integer                       :: ti, tstep ! looping variables

  !-----------------------------------------------------------------------------!
  ! netcdf variables
  integer rcode, ncid_out, varid, ncid
  integer tpid, xpid, ypid, zpid, upid, vpid, wpid
  integer qvpid, qcpid, qspid, qgpid, qhlpid, qipid, qrpid, thpid, thrpid, buoypid
  integer ncpid, nspid, ngpid, nipid, nrpid, rhodid, pipid
  integer pgfzpid, pgfxpid, pgfypid
  integer pgfxbpid, pgfxdlpid, pgfxdnpid
  integer pgfybpid, pgfydlpid, pgfydnpid
  integer pgfzbpid, pgfzdlpid, pgfzdnpid
  integer tkepid
  integer xvortpid, yvortpid, zvortpid, svortpid, cvortpid, psipid
  integer zvort_stretchid, zvort_tiltid
  integer svort_tiltid, svort_stretchid, svort_bclnid, svort_xchgid
  integer dpsidtid
  integer cvort_tiltid, cvort_stretchid, cvort_bclnid, cvort_xchgid

  !-----------------------------------------------------------------------------!
  ! field variables
  real,dimension(:,:,:),allocatable :: u, v, w ! wind (m/s)
  real,dimension(:,:,:),allocatable :: u1, v1, w1 ! wind (m/s)
  real,dimension(:,:,:),allocatable :: u2, v2, w2 ! wind (m/s)

  real,dimension(:,:,:),allocatable :: qv, qc, qi, qs, qr, qg, qhl ! mixing ratios (kg/kg)
  real,dimension(:,:,:),allocatable :: qv1, qc1, qi1, qs1, qr1, qg1, qhl1 ! mixing ratios (kg/kg)
  real,dimension(:,:,:),allocatable :: qv2, qc2, qi2, qs2, qr2, qg2, qhl2 ! mixing ratios (kg/kg)


  real,dimension(:,:,:),allocatable :: th, thr ! potential and density potential temp (K)
  real,dimension(:,:,:),allocatable :: th1, thr1 ! potential and density potential temp (K)
  real,dimension(:,:,:),allocatable :: th2, thr2 ! potential and density potential temp (K)

  real,dimension(:,:,:),allocatable :: rhod, rhod1, rhod2 ! dry air density (kg/m^3)
  real,dimension(:,:,:),allocatable :: ncc, ncc1, ncc2 ! number concentration of cloud drops (m^-3)
  real,dimension(:,:,:),allocatable :: ncr, ncr1, ncr2 ! number concentration of rain drops
  real,dimension(:,:,:),allocatable :: nci, nci1, nci2 ! number concertration of ice particles
  real,dimension(:,:,:),allocatable :: ncs, ncs1, ncs2 ! number concertration of snowflakes
  real,dimension(:,:,:),allocatable :: ncg, ncg1, ncg2 ! number concentration of graupel
  real,dimension(:,:,:),allocatable :: ppi, ppi1, ppi2 ! nondimensionalized pressure
  real,dimension(:,:,:),allocatable :: utem, vtem, wtem ! temporary wind arrays
  real,dimension(:,:,:),allocatable :: buoy, buoy1, buoy2 ! buoyancy (m/s^2)

  real,dimension(:,:,:),allocatable :: pgfx, pgfy, pgfz ! total PGF (m/s^2)
  real,dimension(:,:,:),allocatable :: pgfx1, pgfy1, pgfz1 ! total PGF (m/s^2)
  real,dimension(:,:,:),allocatable :: pgfx2, pgfy2, pgfz2 ! total PGF (m/s^2)

  real,dimension(:,:,:),allocatable :: pgfx_b, pgfy_b, pgfz_b ! buoyancy PPGF
  real,dimension(:,:,:),allocatable :: pgfx_b1, pgfy_b1, pgfz_b1 ! buoyancy PPGF
  real,dimension(:,:,:),allocatable :: pgfx_b2, pgfy_b2, pgfz_b2 ! buoyancy PPGF

  real,dimension(:,:,:),allocatable :: pgfx_dl, pgfy_dl, pgfz_dl ! linear dynamic PPGF
  real,dimension(:,:,:),allocatable :: pgfx_dl1, pgfy_dl1, pgfz_dl1 ! linear dynamic PPGF
  real,dimension(:,:,:),allocatable :: pgfx_dl2, pgfy_dl2, pgfz_dl2 ! linear dynamic PPGF

  real,dimension(:,:,:),allocatable :: pgfx_dn, pgfy_dn, pgfz_dn ! nonlinear dynamic PPGF
  real,dimension(:,:,:),allocatable :: pgfx_dn1, pgfy_dn1, pgfz_dn1 ! nonlinear dynamic PPGF
  real,dimension(:,:,:),allocatable :: pgfx_dn2, pgfy_dn2, pgfz_dn2 ! nonlinear dynamic PPGF

  real,dimension(:),allocatable     :: pi0, th0, thr0, qv0 ! base state variables
  real,dimension(:,:,:),allocatable :: tke, tke1, tke2 ! tubulent kinetic energy

  real,dimension(:,:,:),allocatable :: xvort, yvort, zvort ! vorticity
  real,dimension(:,:,:),allocatable :: dpsi_dx, dpsi_dy, dpsi_dz ! for exchange term
  real,dimension(:,:,:),allocatable :: dvhdx, dvhdy, dvhdz ! horizontal wind speed
  real,dimension(:,:,:),allocatable :: dbdx, dbdy ! horizontal buoyancy gradients
  real,dimension(:,:,:),allocatable :: dwdx, dwdy, dwdz ! w-gradients


  !-----------------------------------------------------------------------------!
  ! parcel variables
  real,dimension(:,:),allocatable :: xpc, ypc, zpc ! parcel position (m)
  real,dimension(:,:),allocatable :: up, vp, wp ! parcel winds
  real,dimension(:,:),allocatable :: qvp, qcp, qrp, qsp, qip, qgp, qhlp ! parcel mixing ratios
  real,dimension(:,:),allocatable :: ncp, nrp, nsp, nip, ngp ! parcel number concentrations
  real,dimension(:,:),allocatable :: thp, thrp, rhodp, pip
  real,dimension(:,:),allocatable :: buoyp, pgfzp, pgfxp, pgfyp
  real,dimension(:,:),allocatable :: pgfxbp, pgfxdlp, pgfxdnp
  real,dimension(:,:),allocatable :: pgfybp, pgfydlp, pgfydnp
  real,dimension(:,:),allocatable :: pgfzbp, pgfzdlp, pgfzdnp
  real,dimension(:,:),allocatable :: tkep
  real,dimension(:,:),allocatable :: xvortp, yvortp, zvortp ! parcel vorticities
  real,dimension(:,:),allocatable :: svortp, cvortp ! streamwise and crosswise vorticity
  real,dimension(:,:),allocatable :: psip ! angle between vorticity and wind vectors
  real,dimension(:,:),allocatable :: zvort_stretchp ! stretching of vertical vorticity
  real,dimension(:,:),allocatable :: zvort_tiltp ! tilting into vertical vorticity
  real,dimension(:,:),allocatable :: svort_tiltp ! tilting into streamwise vorticity
  real,dimension(:,:),allocatable :: svort_stretchp ! stretching of streamwise vorticity
  real,dimension(:,:),allocatable :: svort_bclnp ! baroclinic production of streamwise vorticity
  real,dimension(:,:),allocatable :: svort_xchgp ! streamwise exchange term
  real,dimension(:,:),allocatable :: dpsidt ! change in angle between vorticity and wind vectors
  real,dimension(:,:),allocatable :: cvort_tiltp, cvort_stretchp, cvort_bclnp ! crosswise terms
  real,dimension(:,:),allocatable :: cvort_xchgp

  !-----------------------------------------------------------------------------!
  ! looping and dummy variables
  integer i, t
  integer,dimension(1) :: tmp


  !-----------------------------------------------------------------------------!
  ! namelist

  namelist /filenames/  cm1list, dynlist, plist, outfile
  namelist /parameters/ begintime, parcelt, endtime, dz, dt, integration_method, ptype
  namelist /outputs/    output_dyn, output_vor, output_z_only, noncq

  !=============================================================================!
  !#############################################################################!
  !################################# MAIN BODY #################################!
  !#############################################################################!
  !=============================================================================!

  write(*,*)
  write(*,*) '                  TRAJECT_CM1'
  write(*,*)

  !##############################################################################!
  !#                              INITIALIZATION                                #!
  !##############################################################################!

  !------------------------------------------------------------------------------!
  ! read namelist.input                                                          !
  !------------------------------------------------------------------------------!

  write(*,*) 'opening traject.input'
  open(8,file='traject.input')

  read(8,nml=filenames)
  write(*,*) 'cm1list    =',cm1list     ! list of cm1 files
  write(*,*) 'dynlist    =',dynlist     ! list of dynamics files
  write(*,*) 'plist      =',plist       ! file with list of parcel starting positions
  write(*,*) 'outfile    =',outfile     ! outfile

  read(8,nml=parameters)
  write(*,*) 'begintime  =',begintime   ! beginning time for trajectory run
  write(*,*) 'parcelt=   ',parcelt      ! parcel beginning time
  write(*,*) 'endtime    =',endtime     ! ending time for trajectory run
  write(*,*) 'dt         =',dt          ! time step

  read(8,nml=outputs)
  write(*,*) 'output_dyn =',output_dyn  ! output dynamics?
  write(*,*) 'output_vor =',output_vor  ! output vorticity dynamics?
  write(*,*) 'output_z_oly = ',output_z_only

  close(8)

  !------------------------------------------------------------------------------!
  ! set up arrays for times                                                      !
  !------------------------------------------------------------------------------!
  ! set up array for times of parcels
  write(*,*)
  write(*,*) 'setting up arrays for parcel times'
  ntimes     = int( 1 + (endtime-begintime)/dt ) ! total number of times
  nftimes    = int( 1 + (endtime-parcelt)/dt )   ! number of forward integrations
  nbtimes    = int( (parcelt-begintime)/dt )     ! number of backward integrations
  start_time = nbtimes + 1                       ! start_time is an integer index to the tp array

  allocate( tp(ntimes) )          ! parcel times
  do t=1,ntimes
    tp(t) = begintime + (t-1)*dt
    write(*,*) 't,tp(t)=',t,tp(t)
  enddo

  !------------------------------------------------------------------------------!
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

!  stop
  !-----------------------------------------------------------------------------!
  ! read in base state

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
  allocate( u1(1:ni+1,1:nj,0:nk) )
  allocate( v1(1:ni,1:nj+1,0:nk) )
  allocate( w1(1:ni,1:nj,1:nk+1) )
  allocate( u2(1:ni+1,1:nj,0:nk) )
  allocate( v2(1:ni,1:nj+1,0:nk) )
  allocate( w2(1:ni,1:nj,1:nk+1) )


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

  allocate( qv1(ni,nj,0:nk) )
  allocate( qc1(ni,nj,0:nk) )
  allocate( ncc1(ni,nj,0:nk) )
  allocate( qr1(ni,nj,0:nk) )
  allocate( ncr1(ni,nj,0:nk) )
  allocate( qi1(ni,nj,0:nk) )
  allocate( nci1(ni,nj,0:nk) )
  allocate( qs1(ni,nj,0:nk) )
  allocate( ncs1(ni,nj,0:nk) )
  allocate( qg1(ni,nj,0:nk) )
  allocate( qhl1(ni,nj,0:nk) )
  allocate( ncg1(ni,nj,0:nk) )
  allocate( th1(ni,nj,0:nk) )
  allocate( ppi1(ni,nj,0:nk) )

  allocate( qv2(ni,nj,0:nk) )
  allocate( qc2(ni,nj,0:nk) )
  allocate( ncc2(ni,nj,0:nk) )
  allocate( qr2(ni,nj,0:nk) )
  allocate( ncr2(ni,nj,0:nk) )
  allocate( qi2(ni,nj,0:nk) )
  allocate( nci2(ni,nj,0:nk) )
  allocate( qs2(ni,nj,0:nk) )
  allocate( ncs2(ni,nj,0:nk) )
  allocate( qg2(ni,nj,0:nk) )
  allocate( qhl2(ni,nj,0:nk) )
  allocate( ncg2(ni,nj,0:nk) )
  allocate( th2(ni,nj,0:nk) )
  allocate( ppi2(ni,nj,0:nk) )


  allocate( thr(ni,nj,0:nk) )
  allocate( rhod(ni,nj,0:nk) )

  allocate( thr1(ni,nj,0:nk) )
  allocate( rhod1(ni,nj,0:nk) )

  allocate( thr2(ni,nj,0:nk) )
  allocate( rhod2(ni,nj,0:nk) )


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

  allocate( buoy1(ni,nj,0:nk) )
  allocate( pgfx1(ni+1,nj,0:nk) )
  allocate( pgfy1(ni,nj+1,0:nk) )
  allocate( pgfz1(ni,nj,1:nk+1) )

  allocate( pgfx_b1(1:ni+1,1:nj,0:nk) )
  allocate( pgfy_b1(1:ni,1:nj+1,0:nk) )
  allocate( pgfz_b1(1:ni,1:nj,1:nk+1) )

  allocate( pgfx_dl1(1:ni+1,1:nj,0:nk) )
  allocate( pgfy_dl1(1:ni,1:nj+1,0:nk) )
  allocate( pgfz_dl1(1:ni,1:nj,1:nk+1) )

  allocate( pgfx_dn1(1:ni+1,1:nj,1:nk) )
  allocate( pgfy_dn1(1:ni,1:nj+1,1:nk) )
  allocate( pgfz_dn1(1:ni,1:nj,1:nk+1) )

  allocate( tke1(1:ni,1:nj,1:nk+1) )

  allocate( buoy2(ni,nj,0:nk) )
  allocate( pgfx2(ni+1,nj,0:nk) )
  allocate( pgfy2(ni,nj+1,0:nk) )
  allocate( pgfz2(ni,nj,1:nk+1) )

  allocate( pgfx_b2(1:ni+1,1:nj,0:nk) )
  allocate( pgfy_b2(1:ni,1:nj+1,0:nk) )
  allocate( pgfz_b2(1:ni,1:nj,1:nk+1) )

  allocate( pgfx_dl2(1:ni+1,1:nj,0:nk) )
  allocate( pgfy_dl2(1:ni,1:nj+1,0:nk) )
  allocate( pgfz_dl2(1:ni,1:nj,1:nk+1) )

  allocate( pgfx_dn2(1:ni+1,1:nj,1:nk) )
  allocate( pgfy_dn2(1:ni,1:nj+1,1:nk) )
  allocate( pgfz_dn2(1:ni,1:nj,1:nk+1) )

  allocate( tke2(1:ni,1:nj,1:nk+1) )

! not doing any of vort stuff yet 8 mar 2018
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

  !----------------------------------------------------------!
  ! read in parcels
  write(*,*)
  write(*,*) 'reading parcels'
  open(8,file=plist)

  read(8,*) nparcels
  write(*,*) nparcels,' parcels'

  !----------------------------------------------------------!
  ! create netcdf file
  call nc_traject_vardef(  outfile, ncid_out, ntimes, nparcels, &
      tpid, xpid, ypid, zpid, upid, vpid, wpid, qvpid, qcpid, qrpid, qipid, &
      qspid, qgpid, qhlpid, ncpid, nrpid, nipid, nspid, ngpid, &
      thpid, thrpid, buoypid, pgfxpid, pgfypid, pgfzpid, &
      pgfxbpid, pgfybpid, pgfzbpid, pgfxdlpid, pgfydlpid, pgfzdlpid, &
      pgfxdnpid, pgfydnpid, pgfzdnpid, &
      tkepid, rhodid, pipid, xvortpid, yvortpid, zvortpid, &
      svortpid, cvortpid, psipid, zvort_stretchid, zvort_tiltid,&
      svort_tiltid, svort_stretchid, svort_bclnid, svort_xchgid, &
      dpsidtid, cvort_tiltid, cvort_stretchid, cvort_bclnid, &
      cvort_xchgid )

  !-----------------------------------------------------------!
  ! set up trajectory arrays
  write(*,*) 'setting up trajectory arrays'
  allocate( xpc(ntimes,nparcels) )
  allocate( ypc(ntimes,nparcels) )
  allocate( zpc(ntimes,nparcels) )

  allocate( up(ntimes,nparcels) )
  allocate( vp(ntimes,nparcels) )
  allocate( wp(ntimes,nparcels) )

  allocate( qvp(ntimes,nparcels) )
  allocate( qcp(ntimes,nparcels) )
  allocate( ncp(ntimes,nparcels) )
  allocate( qrp(ntimes,nparcels) )
  allocate( nrp(ntimes,nparcels) )
  allocate( qip(ntimes,nparcels) )
  allocate( nip(ntimes,nparcels) )
  allocate( qsp(ntimes,nparcels) )
  allocate( nsp(ntimes,nparcels) )
  allocate( qgp(ntimes,nparcels) )
  allocate( qhlp(ntimes,nparcels) )
  allocate( ngp(ntimes,nparcels) )
  allocate( thp(ntimes,nparcels) )
  allocate( pip(ntimes,nparcels) )

  allocate( thrp(ntimes,nparcels) )
  allocate( rhodp(ntimes,nparcels) )

  allocate( buoyp(ntimes,nparcels) )

  allocate( pgfxp(ntimes,nparcels) )
  allocate( pgfyp(ntimes,nparcels) )
  allocate( pgfzp(ntimes,nparcels) )

  allocate( pgfxbp(ntimes,nparcels) )
  allocate( pgfybp(ntimes,nparcels) )
  allocate( pgfzbp(ntimes,nparcels) )

  allocate( pgfxdlp(ntimes,nparcels) )
  allocate( pgfydlp(ntimes,nparcels) )
  allocate( pgfzdlp(ntimes,nparcels) )

  allocate( pgfxdnp(ntimes,nparcels) )
  allocate( pgfydnp(ntimes,nparcels) )
  allocate( pgfzdnp(ntimes,nparcels) )

  allocate( tkep( ntimes, nparcels) )

  allocate( xvortp(ntimes,nparcels) )
  allocate( yvortp(ntimes,nparcels) )
  allocate( zvortp(ntimes,nparcels) )
  allocate( svortp(ntimes,nparcels) )
  allocate( cvortp(ntimes,nparcels) )

  allocate( psip(ntimes,nparcels) )

  allocate( zvort_stretchp(ntimes,nparcels) )
  allocate( zvort_tiltp(ntimes,nparcels) )

  allocate( svort_tiltp(ntimes,nparcels) )
  allocate( svort_stretchp(ntimes,nparcels) )
  allocate( svort_xchgp(ntimes,nparcels) )
  allocate( svort_bclnp(ntimes,nparcels) )

  allocate( dpsidt(ntimes,nparcels) )

  allocate( cvort_tiltp(ntimes,nparcels) )
  allocate( cvort_stretchp(ntimes,nparcels) )
  allocate( cvort_xchgp(ntimes,nparcels) )
  allocate( cvort_bclnp(ntimes,nparcels) )

  write(*,*) 'number of backward integration steps = ',nbtimes
  write(*,*) 'number of forward integration steps = ',nftimes

  !----------------------------------------------------------------------------!
  ! read in parcels
  do t=1,nparcels
    read(8,*) xpc(start_time,t), ypc(start_time,t), zpc(start_time,t)
    write(*,*) 'parcel ',t,' starts at x,y,z=',xpc(start_time,t),'km,',ypc(start_time,t),'km,',zpc(start_time,t),'km'
  enddo

  write(*,*)
  close(8)

  !##############################################################################!
  !#                           TRAJECTORY CALCULATIONS                          #!
  !##############################################################################!

  !------------------------------------------------------------------------------!
  ! do forward trajectory calculations

  if(endtime.gt.parcelt)then

    write(*,*) '################################################################'
    write(*,*) 'starting forward trajectory calculations'
    write(*,*) 

    !----------------- set up initial times ------------------------!
    tmp = minloc( abs( parcelt - times ) )
    ti = tmp(1)

    write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

    call traject_nc_read_cm1( infiles, ti, u1, u2, v1, v2, w1, w2, qv1, qv2, &
      qc1, qc2, qr1, qr2, qi1, qi2, qs1, qs2, qg1, qg2, qhl1, qhl2, th1, th2, &
      thr0, pi0, &
      pgfx1, pgfx2, pgfy1, pgfy2, pgfz1, pgfz2, thr1, thr2, buoy1, buoy2,&
      tke1, tke2, ncc1, ncc2, ncr1, vcr2, nci1, nci2, ncs1, ncs2, ncg1, ng2, &
      rhod1, rhod2, ppi1, ppi2,dt )

    if (output_dyn.eq.1) then
      call traject_nc_read_dyn( dynfiles, ti, pgfx_b1, pgfx_b2, &
        pgfy_b1, pgfy_b2, pgfz_b1, pgfz_b2, pgfx_dl1, pgfx_dl2, &
        pgfy_dl1, pgfy_dl2, pgfz_dl1, pgfz_dl2, &
        pgfx_dn1, pgfx_dn2, pgfy_dn1, pgfy_dn2, pgfz_dn1, pgfz_dn1, thr1, 
        thr2, buoy1, buoy2,dt )
    endif

    !-------------------- loop through forward times ----------------!
    do tstep=start_time,ntimes-1

      tt=tp(tstep)
      write(*,*) '--------------------------------------------------------------'
      write(*,*) 'time step ',tstep,' or ',tt,'s'

      if(output_vor.eq.1)then
        call compute_vorts( buoy, u, v, w, xvort, yvort, zvort, &
          dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
          dwdy, dwdx, dwdz ) 
      endif

      call interpintime_winds( u1, u2, u, v1, v2, v, w1, w2, w, times(ti), times(ti+1), tt )
      call interpintime_scalars( buoy1, buoy2, buoy, tke1, tke2, tke, thr1, thr2, thr, &
        qv1, qv2, qv, qc1, qc2, qc, qr1, qr2, qr, qi1, qi2, qi, qs1, qs2, qs, &
        qg1, qg2, qg, qi1, qi2, qi, qs1, qs2, qs, qg1, qg2, qg, qhl1, qhl2, qhl, &
        ncc1, ncc2, ncc, ncr1, ncr2, ncr, nci1, nci2, nci, ncs1, ncs2, ncs, &
        ngi1, ngi1, ngi, rhod1, rhod2, rhod, th1, th2, th, ppi1, ppi1, ppi, times(ti), times(ti+1), tt )
      call interpintime_winds( pgfx1, pgfx2, pgfx, pgfy1, pgfy2, pgfz1, pgfz2, pgfz, times(ti), times(ti+1), tt )

      if(output_dyn.eq.1)then
        call interpintime_winds( pgfx_b1, pgfx_b2, pgfx_b, pgfy_b1, pgfy_b2, pgfy_b, &
          pgfz_b1, pgfz_b2, pgfz_b, times(ti), times(ti+1), tt )
        call interpintime_winds( pgfx_dl1, pgfx_dl2, pgfx_dl, pgfy_dl1, pgfy_dl2, pgfy_dl, &
          pgfz_dl1, pgfz_dl2, pgfz_dl, times(ti), times(ti+1), tt )
        call interpintime_winds( pgfx_dn1, pgfx_dn2, pgfx_dn, pgfy_dn1, pgfy_dn2, pgfy_dn, &
          pgfz_dn1, pgfz_dn2, pgfz_dn, times(ti), times(ti+1), tt )
      endif

      do np=1,nparcels
 
        call interp_all( xpc(tstep,np), ypc(tstep,np), zpc(tstep,np), zzh, u, up(tstep,np),&
          v, vp(tstep,np), &
          w, wp(tstep,np), buoy, buoyp(tstep,np), pgfx, pgfxp(tstep,np), &
          pgfy, pgfyp(tstep,np), &
          pgfz, pgfzp(tstep,np), &
          tke, tkep(tstep,np), thr0, thr, thrp(tstep,np), qv, qvp(tstep,np), &
          qc, qcp(tstep,np), &
          qr, qrp(tstep,np), qi, qip(tstep,np), qs, qsp(tstep,np), &
          qg, qgp(tstep,np), qhl, qhlp(tstep,np), &
          ncc, ncp(tstep,np), ncr, nrp(tstep,np), nci, nip(tstep,np), &
          ncs, nsp(tstep,np), ncg, ngp(tstep,np), rhod, rhodp(tstep,np),&
          th, thp(tstep,np),&
          ppi, pip(tstep,np) )

        if(output_dyn.eq.1)then
          call interp_dynz( xpc(tstep,np), ypc(tstep,np), zpc(tstep,np), &
            pgfx_b, pgfxbp(tstep,np), pgfy_b, pgfybp(tstep,np), pgfz_b, pgfzbp(tstep,np), &
            pgfx_dl, pgfxdlp(tstep,np), pgfy_dl, pgfydlp(tstep,np), pgfz_dl, pgfzdlp(tstep,np), &
            pgfx_dn, pgfxdnp(tstep,np), pgfy_dn, pgfydnp(tstep,np), pgfz_dn, pgfzdnp(tstep,np) )
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

        IF(INTEGRATION_METHOD.EQ.0)THEN
          xpc(tstep+1,np)=xpc(tstep,np)+dt*up(tstep,np)
          ypc(tstep+1,np)=ypc(tstep,np)+dt*vp(tstep,np)
          zpc(tstep+1,np)=zpc(tstep,np)+dt*wp(tstep,np)
        ENDIF

      enddo

!      ti=ti+1
!      write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

        IF(INTEGRATION_METHOD.EQ.0)THEN
      if(tp(tstep).gt.infile(ti+1))then
        ti=ti+1
        call traject_nc_read_cm1( infiles, ti, u1, u2, v1, v2, w1, w2, qv1, qv2, &
    qc1, qc2, qr1, qr2, qi1, qi2, qs1, qs2, qg1, qg2, qhl1, qhl2, &
    th1, th2, thr0, pi0, pgfx1, pgfx2, pgfy1, pgfy2, pgfz1, pgfz2, thr1, thr2,&
    buoy1, buoy2, tke1, tke2,&
    ncc1, ncc2, ncr1, ncr2, nci1, nci2, ncs1, ncs2, ncg1, ncg2, rhod1, rhod2, &
    ppi1, ppi2, dt )
        if(output_dyn.eq.1)then
      call traject_nc_read_dyn( dynfiles, ti, pgfx_b1, pgfx_b2, &
        pgfy_b1, pgfy_b2, pgfz_b1, pgfz_b2, pgfx_dl1, pgfx_dl2, &
        pgfy_dl1, pgfy_dl2, pgfz_dl1, pgfz_dl2, &
        pgfx_dn1, pgfx_dn2, pgfy_dn1, pgfy_dn2, pgfz_dn1, pgfz_dn1, thr1,
        thr2, buoy1, buoy2, dt )
        endif
      endif


        ELSEIF(INTEGRATION_METHOD.EQ.1)THEN

!      call traject_nc_read_cm1( infiles(ti), utem, vtem, wtem, qv, qc, qr, qi, qs, qg, qhl, th, &
!           thr0, pi0, pgfx, pgfy, pgfz, thr, buoy, tke, ncc, ncr, nci, ncs, ncg, &
!           rhod, ppi )

      if(tp(tstep).gt.infiles(ti+1))then
        ti=ti+1
        call traject_nc_read_cm1( infiles, ti, u1, u2, v1, v2, w1, w2, qv1, qv2, &
    qc1, qc2, qr1, qr2, qi1, qi2, qs1, qs2, qg1, qg2, qhl1, qhl2, &
    th1, th2, thr0, pi0, pgfx1, pgfx2, pgfy1, pgfy2, pgfz1, pgfz2, thr1, thr2,&
    buoy1, buoy2, tke1, tke2,&
    ncc1, ncc2, ncr1, ncr2, nci1, nci2, ncs1, ncs2, ncg1, ncg2, rhod1, rhod2, &
    ppi1, ppi2, dt )
        if(output_dyn.eq.1)then
      call traject_nc_read_dyn( dynfiles, ti, pgfx_b1, pgfx_b2, &
        pgfy_b1, pgfy_b2, pgfz_b1, pgfz_b2, pgfx_dl1, pgfx_dl2, &
        pgfy_dl1, pgfy_dl2, pgfz_dl1, pgfz_dl2, &
        pgfx_dn1, pgfx_dn2, pgfy_dn1, pgfy_dn2, pgfz_dn1, pgfz_dn1, thr1,
        thr2, buoy1, buoy2, dt )
        endif
      endif

      call interpintime_winds( u1, u2, utem, v1, v2, vtem, w1, w2, wtem, times(ti), times(ti+1), tp(tstep+1) )
 
      do np=1,nparcels
        call rk4( xpc(tstep:tstep+1,np), ypc(tstep:tstep+1,np), zpc(tstep:tstep+1,np), &
          up(tstep,np), vp(tstep,np), wp(tstep,np), u, v, w, utem, vtem, wtem, dt )
      enddo

      u = utem
      v = vtem
      w = wtem

        ENDIF 



!      if (output_dyn.eq.1) then
!        call traject_nc_read_dyn( dynfiles(ti), pgfx_b, pgfy_b, pgfz_b, pgfx_dl, pgfy_dl, pgfz_dl, &
!          pgfx_dn, pgfy_dn, pgfz_dn, thr,buoy )
!      endif

    enddo ! loop over times

    !------------------ computations for final time step -------------------!
    if(output_vor.eq.1)then
      call compute_vorts( buoy, u, v, w, xvort, yvort, zvort, &
        dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
        dwdy, dwdx, dwdz ) 
    endif
      tt=tp(ntimes)
      call interpintime_winds( u1, u2, u, v1, v2, v, w1, w2, w, times(nfiles-1), times(nfiles), tt )
      call interpintime_scalars( buoy1, buoy2, buoy, tke1, tke2, tke, thr1, thr2, thr, &
        qv1, qv2, qv, qc1, qc2, qc, qr1, qr2, qr, qi1, qi2, qi, qs1, qs2, qs, &
        qg1, qg2, qg, qi1, qi2, qi, qs1, qs2, qs, qg1, qg2, qg, qhl1, qhl2, qhl, &
        ncc1, ncc2, ncc, ncr1, ncr2, ncr, nci1, nci2, nci, ncs1, ncs2, ncs, &
        ngi1, ngi1, ngi, rhod1, rhod2, rhod, th1, th2, th, ppi1, ppi1, ppi, times(nfiles-1), times(nfiles), tt )
      call interpintime_winds( pgfx1, pgfx2, pgfx, pgfy1, pgfy2, pgfz1, pgfz2, pgfz, times(nfiles-1), times(nfiles), tt )

      if(output_dyn.eq.1)then
        call interpintime_winds( pgfx_b1, pgfx_b2, pgfx_b, pgfy_b1, pgfy_b2, pgfy_b, &
          pgfz_b1, pgfz_b2, pgfz_b, times(nfiles-1), times(nfiles), tt )
        call interpintime_winds( pgfx_dl1, pgfx_dl2, pgfx_dl, pgfy_dl1, pgfy_dl2, pgfy_dl, &
          pgfz_dl1, pgfz_dl2, pgfz_dl, times(nfiles-1), times(nfiles), tt )
        call interpintime_winds( pgfx_dn1, pgfx_dn2, pgfx_dn, pgfy_dn1, pgfy_dn2, pgfy_dn, &
          pgfz_dn1, pgfz_dn2, pgfz_dn, times(nfiles-1), times(nfiles), tt )
      endif


    do np=1,nparcels

      call interp_all( xpc(ntimes,np), ypc(ntimes,np), zpc(ntimes,np), zzh, u, up(ntimes,np), v, vp(ntimes,np), &
        w, wp(ntimes,np), buoy, buoyp(ntimes,np), pgfx, pgfxp(ntimes,np), pgfy, pgfyp(ntimes,np), &
        pgfz, pgfzp(ntimes,np), &
        tke, tkep(ntimes,np), thr0, thr, thrp(ntimes,np), qv, qvp(ntimes,np), &
        qc, qcp(ntimes,np), qr, qrp(ntimes,np), qi, qip(ntimes,np), qs, qsp(ntimes,np), &
        qg, qgp(ntimes,np), qhl, qhlp(ntimes,np),  &
          ncc, ncp(ntimes,np), ncr, nrp(ntimes,np), nci, nip(ntimes,np), &
          ncs, nsp(ntimes,np), ncg, ngp(ntimes,np), rhod, rhodp(ntimes,np), &
          th, thp(ntimes,np), ppi, pip(ntimes,np) )

      if(output_dyn.eq.1)then
        call interp_dynz( xpc(ntimes,np), ypc(ntimes,np), zpc(ntimes,np), &
          pgfx_b, pgfxbp(ntimes,np), pgfy_b, pgfybp(ntimes,np), pgfz_b, pgfzbp(ntimes,np), &
          pgfx_dl, pgfxdlp(ntimes,np), pgfy_dl, pgfydlp(ntimes,np), pgfz_dl, pgfzdlp(ntimes,np), &
          pgfx_dn, pgfxdnp(ntimes,np), pgfy_dn, pgfydnp(ntimes,np), pgfz_dn, pgfzdnp(ntimes,np) )
      endif

      if(output_vor.eq.1)then
        call interp_vort3( xpc(ntimes,np), ypc(ntimes,np), &
          zpc(ntimes,np), up(ntimes,np), vp(ntimes,np), wp(ntimes,np), &
          xvort, yvort, zvort, dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, &
          dpsi_dz, dvhdz, dwdy, dwdx, dwdz, &
          xvortp(ntimes,np), yvortp(ntimes,np), zvortp(ntimes,np), &
          svortp(ntimes,np), cvortp(ntimes,np), psip(ntimes,np), &
          zvort_tiltp(ntimes,np), zvort_stretchp(ntimes,np), &
          svort_tiltp(ntimes,np), svort_stretchp(ntimes,np), &
          svort_bclnp(ntimes,np), cvort_stretchp(ntimes,np), &
          cvort_tiltp(ntimes,np), cvort_bclnp(ntimes,np) )
      endif


    enddo

  endif

  !------------------------------------------------------------------------------!
  ! do backward trajectory calculations

  if(begintime.lt.parcelt)then

    write(*,*) '################################################################'
    write(*,*) 'starting backward trajectory calculations'
    write(*,*) 

    dt=-dt

    !----------------- set up initial times ------------------------!
    tmp = minloc( abs( parcelt - times ) )
    ti = tmp(1)

    write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

!    call traject_nc_read_cm1( infiles(ti), u, v, w, qv, qc, qr, qi, qs, qg, qhl, th, thr0, pi0, &
!      pgfx, pgfy, pgfz, thr, buoy,&
!      tke, ncc, ncr, nci, ncs, ncg, rhod, ppi )
!
!    if (output_dyn.eq.1) then
!      call traject_nc_read_dyn( dynfiles(ti), pgfx_b, pgfy_b, pgfz_b, pgfx_dl, pgfy_dl, pgfz_dl, &
!        pgfx_dn, pgfy_dn, pgfz_dn, thr, buoy )
!    endif

    call traject_nc_read_cm1( infiles, ti, u1, u2, v1, v2, w1, w2, qv1, qv2, &
      qc1, qc2, qr1, qr2, qi1, qi2, qs1, qs2, qg1, qg2, qhl1, qhl2, th1, th2, &
      thr0, pi0, &
      pgfx1, pgfx2, pgfy1, pgfy2, pgfz1, pgfz2, thr1, thr2, buoy1, buoy2,&
      tke1, tke2, ncc1, ncc2, ncr1, vcr2, nci1, nci2, ncs1, ncs2, ncg1, ng2, &
      rhod1, rhod2, ppi1, ppi2,dt )

    if (output_dyn.eq.1) then
      call traject_nc_read_dyn( dynfiles, ti, pgfx_b1, pgfx_b2, &
        pgfy_b1, pgfy_b2, pgfz_b1, pgfz_b2, pgfx_dl1, pgfx_dl2, &
        pgfy_dl1, pgfy_dl2, pgfz_dl1, pgfz_dl2, &
        pgfx_dn1, pgfx_dn2, pgfy_dn1, pgfy_dn2, pgfz_dn1, pgfz_dn1, thr1,
        thr2, buoy1, buoy2,dt )
    endif



    !----------------- loop through backwards times ----------------!
    do tstep=start_time,2,-1

      tt=tp(tstep)
      write(*,*) '--------------------------------------------------------------'
      write(*,*) 'time step ',tstep,' or ',tt,'s'

      if(output_vor.eq.1)then
        call compute_vorts( buoy, u, v, w, xvort, yvort, zvort, &
          dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
          dwdy, dwdx, dwdz ) 
      endif

      call interpintime_winds( u1, u2, u, v1, v2, v, w1, w2, w, times(ti-1), times(ti), tt )
      call interpintime_scalars( buoy1, buoy2, buoy, tke1, tke2, tke, thr1, thr2, thr, &
        qv1, qv2, qv, qc1, qc2, qc, qr1, qr2, qr, qi1, qi2, qi, qs1, qs2, qs, &
        qg1, qg2, qg, qi1, qi2, qi, qs1, qs2, qs, qg1, qg2, qg, qhl1, qhl2, qhl, &
        ncc1, ncc2, ncc, ncr1, ncr2, ncr, nci1, nci2, nci, ncs1, ncs2, ncs, &
        ngi1, ngi1, ngi, rhod1, rhod2, rhod, th1, th2, th, ppi1, ppi1, ppi, times(ti-1), times(ti), tt )
      call interpintime_winds( pgfx1, pgfx2, pgfx, pgfy1, pgfy2, pgfz1, pgfz2, pgfz, times(ti-1), times(ti), tt )

      if(output_dyn.eq.1)then
        call interpintime_winds( pgfx_b1, pgfx_b2, pgfx_b, pgfy_b1, pgfy_b2, pgfy_b, &
          pgfz_b1, pgfz_b2, pgfz_b, times(ti-1), times(ti), tt )
        call interpintime_winds( pgfx_dl1, pgfx_dl2, pgfx_dl, pgfy_dl1, pgfy_dl2, pgfy_dl, &
          pgfz_dl1, pgfz_dl2, pgfz_dl, times(ti-1), times(ti), tt )
        call interpintime_winds( pgfx_dn1, pgfx_dn2, pgfx_dn, pgfy_dn1, pgfy_dn2, pgfy_dn, &
          pgfz_dn1, pgfz_dn2, pgfz_dn, times(ti-1), times(ti), tt )
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
            pgfx_dn, pgfxdnp(tstep,np), pgfy_dn, pgfydnp(tstep,np), pgfz_dn, pgfzdnp(tstep,np) )
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

           IF(INTEGRATION_METHOD.EQ.0)THEN

        xpc(tstep-1,np) = xpc(tstep,np) + dt*up(tstep,np)
        ypc(tstep-1,np) = ypc(tstep,np) + dt*vp(tstep,np)
        zpc(tstep-1,np) = zpc(tstep,np) + dt*wp(tstep,np)

           ENDIF

      enddo ! do np=1,nparcels

!      ti=ti-1
      write(*,*) 'opening ',infiles(ti),'...time ',times(ti)

            IF(INTEGRATION_METHOD.EQ.0)THEN

!      call traject_nc_read_cm1( infiles(ti), u, v, w, qv, qc, qr, qi, qs, qg, qhl, th, &
!           thr0, pi0, pgfx, pgfy, pgfz, thr, buoy, tke, ncc, ncr, nci, ncs, ncg,&
!           rhod, ppi )

      if(tp(tstep).lt.infile(ti-1))then
        ti=ti-1
        call traject_nc_read_cm1( infiles, ti, u1, u2, v1, v2, w1, w2, qv1, qv2, &
    qc1, qc2, qr1, qr2, qi1, qi2, qs1, qs2, qg1, qg2, qhl1, qhl2, &
    th1, th2, thr0, pi0, pgfx1, pgfx2, pgfy1, pgfy2, pgfz1, pgfz2, thr1, thr2,&
    buoy1, buoy2, tke1, tke2,&
    ncc1, ncc2, ncr1, ncr2, nci1, nci2, ncs1, ncs2, ncg1, ncg2, rhod1, rhod2, &
    ppi1, ppi2, dt )
        if(output_dyn.eq.1)then
      call traject_nc_read_dyn( dynfiles, ti, pgfx_b1, pgfx_b2, &
        pgfy_b1, pgfy_b2, pgfz_b1, pgfz_b2, pgfx_dl1, pgfx_dl2, &
        pgfy_dl1, pgfy_dl2, pgfz_dl1, pgfz_dl2, &
        pgfx_dn1, pgfx_dn2, pgfy_dn1, pgfy_dn2, pgfz_dn1, pgfz_dn1, thr1,
        thr2, buoy1, buoy2, dt )
        endif
      endif



            ELSEIF(INTEGRATION_METHOD.EQ.1)THEN

!      call traject_nc_read_cm1( infiles(ti), utem, vtem, wtem, qv, qc, qr, qi, qs, qg, qhl, th, &
!           thr0, pi0, pgfx, pgfy, pgfz, thr, buoy, tke, ncc, ncr, nci, ncs, ncg, &
!           rhod, ppi )
! 

      if(tp(tstep).gt.infiles(ti-1))then
        ti=ti-1
        call traject_nc_read_cm1( infiles, ti, u1, u2, v1, v2, w1, w2, qv1, qv2, &
    qc1, qc2, qr1, qr2, qi1, qi2, qs1, qs2, qg1, qg2, qhl1, qhl2, &
    th1, th2, thr0, pi0, pgfx1, pgfx2, pgfy1, pgfy2, pgfz1, pgfz2, thr1, thr2,&
    buoy1, buoy2, tke1, tke2,&
    ncc1, ncc2, ncr1, ncr2, nci1, nci2, ncs1, ncs2, ncg1, ncg2, rhod1, rhod2, &
    ppi1, ppi2, dt )
        if(output_dyn.eq.1)then
      call traject_nc_read_dyn( dynfiles, ti, pgfx_b1, pgfx_b2, &
        pgfy_b1, pgfy_b2, pgfz_b1, pgfz_b2, pgfx_dl1, pgfx_dl2, &
        pgfy_dl1, pgfy_dl2, pgfz_dl1, pgfz_dl2, &
        pgfx_dn1, pgfx_dn2, pgfy_dn1, pgfy_dn2, pgfz_dn1, pgfz_dn1, thr1,
        thr2, buoy1, buoy2, dt )
        endif
      endif

      call interpintime_winds( u1, u2, utem, v1, v2, vtem, w1, w2, wtem, times(ti-1, times(ti), tp(tstep-1) )

      do np=1,nparcels
        call rk4( xpc(tstep:tstep+1,np), ypc(tstep:tstep+1,np), zpc(tstep:tstep+1,np), &
          up(tstep,np), vp(tstep,np), wp(tstep,np), u, v, w, utem, vtem, wtem, dt )
      enddo




      u = utem
      v = vtem
      w = wtem

        ENDIF 



      if (output_dyn.eq.1) then
        call traject_nc_read_dyn( dynfiles(ti), pgfx_b, pgfy_b, pgfz_b, pgfx_dl, pgfy_dl, pgfz_dl, &
          pgfx_dn, pgfy_dn, pgfz_dn, thr, buoy )
      endif

    enddo ! do tstep=start_time,2,-1
    tt=tp(1)
    !--------------- do first time ---------------------------!
    if(output_vor.eq.1)then
      call compute_vorts( buoy, u, v, w, xvort, yvort, zvort, &
        dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
        dwdy, dwdx, dwdz ) 
    endif

      call interpintime_winds( u1, u2, u, v1, v2, v, w1, w2, w, times(1), times(2), tt )
      call interpintime_scalars( buoy1, buoy2, buoy, tke1, tke2, tke, thr1, thr2, thr, &
        qv1, qv2, qv, qc1, qc2, qc, qr1, qr2, qr, qi1, qi2, qi, qs1, qs2, qs, &
        qg1, qg2, qg, qi1, qi2, qi, qs1, qs2, qs, qg1, qg2, qg, qhl1, qhl2, qhl, &
        ncc1, ncc2, ncc, ncr1, ncr2, ncr, nci1, nci2, nci, ncs1, ncs2, ncs, &
        ngi1, ngi1, ngi, rhod1, rhod2, rhod, th1, th2, th, ppi1, ppi1, ppi, times(1), times(2), tt )
      call interpintime_winds( pgfx1, pgfx2, pgfx, pgfy1, pgfy2, pgfz1, pgfz2, pgfz, times(1), times(2), tt )

      if(output_dyn.eq.1)then
        call interpintime_winds( pgfx_b1, pgfx_b2, pgfx_b, pgfy_b1, pgfy_b2, pgfy_b, &
          pgfz_b1, pgfz_b2, pgfz_b, times(1), times(2) tt )
        call interpintime_winds( pgfx_dl1, pgfx_dl2, pgfx_dl, pgfy_dl1, pgfy_dl2, pgfy_dl, &
          pgfz_dl1, pgfz_dl2, pgfz_dl, times(1), times(2), tt )
        call interpintime_winds( pgfx_dn1, pgfx_dn2, pgfx_dn, pgfy_dn1, pgfy_dn2, pgfy_dn, &
          pgfz_dn1, pgfz_dn2, pgfz_dn, times(1), times(2), tt )
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

  endif

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


  call nc_traject_writeout( ntimes, nparcels, ncid_out, tpid, tp, xpid, xpc, ypid, ypc, zpid, zpc, &
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
      cvort_bclnid, cvort_bclnp, cvort_xchgid, cvort_xchgp )




  write(*,*) 'done writing out netcdf file'


  write(*,*) 'program terminated normally'

contains


  real function compute_dpsi( u1, v1, u2, v2, dx )

    real :: u1, v1, u2, v2, dx
    real :: dpsi, udotv_uv
    real :: u1mag, u2mag

    u1mag=sqrt(u1**2+v1**2)
    u2mag=sqrt(u2**2+v2**2)
    if((u1mag.ne.0).and.(u2mag.ne.0))then
      udotv_uv = ( u1*u2+v1*v2 )/(u1mag*u2mag)
    else
      udotv_uv = 0
    endif

    if( (v1.ge.0).and.(v2.ge.0) )then
      if(u1.gt.u2)then
        dpsi = -acos( udotv_uv )
      elseif(u1.lt.u2)then
        dpsi =  acos( udotv_uv )
      else
        if(v1.gt.v2)then
          dpsi = -acos( udotv_uv )
        else
          dpsi =  acos( udotv_uv )
        endif
      endif
    elseif( (v1.lt.0).and.(v2.lt.0) )then
      if(u1.gt.u2)then
        dpsi =  acos( udotv_uv )
      elseif(u1.lt.u2)then
        dpsi = -acos( udotv_uv )
      else
        if(v1.gt.v2)then
          dpsi =  acos( udotv_uv )
        else
          dpsi = -acos( udotv_uv )
        endif
      endif
    elseif( (v1.le.0).and.(v2.gt.0))then
      dpsi =  acos( udotv_uv )
    elseif( (v1.gt.0).and.(v2.le.0))then
      dpsi = -acos( udotv_uv )
    endif

    compute_dpsi = dpsi/dx

  end function


end program
