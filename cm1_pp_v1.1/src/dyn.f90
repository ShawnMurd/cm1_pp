program dyn
  !#############################################################################!
  !#############################################################################!
  !####                                                                     ####!
  !####                            PROGRAM DYN                              ####!
  !####                                                                     ####!
  !####  Program DYN post-processes netCDF output from CM1, creating new    ####!
  !####  files containing data of dynamical interest.  This version         ####!
  !####  specifically works with microphysics that have six (cloud and rain ####!
  !####  water, cloud ice, snow, and graupel) or seven (hail, ZVD scheme),  ####!
  !####  however can be run for dry simulations (imoist=0).                 ####!
  !####  It is otherwise compatible with output going back to at least r14  ####!
  !####  It is written for single time-single file outputs only.            ####!
  !####  CM1 files are expected to include base state and total variables,  ####!
  !####  and perturbations are calculated from those.  E.g., both th and    ####!
  !####  th0 must be present; thpert will be calculated from those.         ####!
  !####                                                                     ####!
  !####    In the output files are 3-D grids of nondimensionalized pressure ####!
  !####  fields, and associated gradient forces, as well as buoyancy.  Also ####!
  !####  included are density and equivalent potential temperatures.  The   ####!
  !####  Poisson equation for buoyancy and linear dynamic pressures are     ####!
  !####  solved for Neuman boundary conditions, and nonlinear dynamic pres- ####!
  !####  sure is found by subtracting those fields from the total pertur-   ####!
  !####  bation pressure.  The nonuniqueness, while not a serious issue for ####!
  !####  dynamical calculations, is aesthetically distasteful.  The linear  ####!
  !####  dynamic pressure is adjusted so that the southeast bottom corner   ####!
  !####  is 0.0; the buoyancy pressure adjusted so the sum is zero.         ####!
  !####                                                                     ####!
  !####  Note that this version requires the Fastest Fourier Transform in   ####!
  !####  the West, available at fftw.org.                                   ####!
  !####=====================================================================####!
  !#### v0.4 written by Ryan Hastings, completed 25 Jan 2012.               ####!
  !####---------------------------------------------------------------------####!
  !#### v0.5 completed 13 March 2012.  Switched from red-and-black Jacobi   ####!
  !#### method to DCT.                                                      ####!
  !####---------------------------------------------------------------------####!
  !#### v1.0 completed 10 Dec 2016.  Cleaned up code and comments.          ####!
  !#############################################################################!
  !#############################################################################!

  use globals
  use nc_inout
  implicit none
  !=============================================================================!
  !#############################################################################!
  !################################# VARIABLES #################################!
  !#############################################################################!
  !=============================================================================!

  !-----------------------------------------------------------------------------!
  !  I/O variables
  character(len=100)                :: infile, outfile
  integer ncid, rcode ! netcdf file identifier and error code

  !-----------------------------------------------------------------------------!
  ! Poisson solver variables
  real,dimension(:),allocatable      :: cfa, cfc, rf0, ad1, ad2
  real,dimension(:,:,:),allocatable  :: buoyf
  real,dimension(:,:,:),allocatable  :: cfb
  real,dimension(:,:,:),allocatable  :: ff

  !-----------------------------------------------------------------------------!
  ! kinematic variables
  real,dimension(:),allocatable     :: u0, v0 ! base state (m/s)
  real,dimension(:,:,:),allocatable :: u, v, w
  real,dimension(:,:,:),allocatable :: ui, vi, wi ! interpolated to scalar points
  real,dimension(:,:,:),allocatable :: upert, vpert

  !-----------------------------------------------------------------------------!
  ! thermodynamic variables
  real,dimension(:,:,:),allocatable :: th, thpert ! potential temperature (K)
  real,dimension(:,:,:),allocatable :: thr, thrpert ! density potential temperature
  real,dimension(:,:,:),allocatable :: qv ! water vapor mixing ratio (kg/kg)
  real,dimension(:,:,:),allocatable :: qc ! cloud water mixing ratio
  real,dimension(:,:,:),allocatable :: qr ! rain water mixing ratio
  real,dimension(:,:,:),allocatable :: qi ! cloud ice mixing ratio
  real,dimension(:,:,:),allocatable :: qs ! snow mixing ratio
  real,dimension(:,:,:),allocatable :: qg ! graupel mixing ratio
  real,dimension(:,:,:),allocatable :: qhl ! hail mixing ratio (ptype=27 only)
  real,dimension(:,:,:),allocatable :: qt ! total precipitation mixing ratio
  real,dimension(:,:,:),allocatable :: ppi, pipert ! nondimensionalized pressure
  real,dimension(:,:,:),allocatable :: prs, prspert ! pressure (Pa)
  real,dimension(:,:,:),allocatable :: buoy ! buoyancy (m/s^2)
  real,dimension(:),allocatable     :: th0, thr0, pi0, prs0, rho0, qv0 ! base state

  !-----------------------------------------------------------------------------!
  ! perturbation pressure variables

  real,dimension(:,:,:),allocatable :: pi_dl, p_dl ! linear dynamic pressure
  real,dimension(:,:,:),allocatable :: pi_dn, p_dn ! nonlinear dynamic pressure
  real,dimension(:,:,:),allocatable :: pi_b, p_b   ! buoyancy pressure

  real,dimension(:,:,:,:),allocatable :: pgfb, pgfdn, pgfdl ! PGF vectors (m/s^2)

  real,dimension(:),allocatable     :: Sx, Sy ! vertical shear of horizontal base state winds
  real,dimension(:,:,:),allocatable :: wx, wy ! horizontal gradient of vertical velocity


  !-----------------------------------------------------------------------------!
  ! other variables
  integer i, j, k ! looping variables
  real sumpb ! adjustment to pressure fields for uniqueness
  integer imoist ! imoist=0, dry; imoist=1, moist

  !-----------------------------------------------------------------------------!
  ! namelists

  namelist /inputs/ infile, outfile, dz, imoist, ptype
  namelist /outputs/ output_dl, output_dn, output_b

  !=============================================================================!
  !#############################################################################!
  !################################# MAIN BODY #################################!
  !#############################################################################!
  !=============================================================================!

  write(*,*)
  write(*,*) '                CM1 DYNAMICS'
  write(*,*)

  !##############################################################################!
  !#                              INITIALIZATION
  !#!
  !##############################################################################!


  !------------------------------------------------------------------------------!
  ! read namelist.input                                                          !
  !------------------------------------------------------------------------------!

  write(*,*) 'reading in namelist'

  open(8,file='dyn.input')
  write(*,*) '..reading inputs'
  read(8,nml=inputs)
  write(*,*) '..reading outputs'
  read(8,nml=outputs)
  close(8)

  write(*,*)
  write(*,*) 'inputs:'
  write(*,*) '   infile      = ',infile
  write(*,*) '   outfile     = ',outfile
  write(*,*) '   dz          = ',dz
  write(*,*) '   imoist      = ',imoist
  write(*,*) '   ptype       = ',ptype
  write(*,*)

  write(*,*) 'outputs:'
  write(*,*) '   output_dl   = ',output_dl
  write(*,*) '   output_b    = ',output_b
  write(*,*) '   output_dn   = ',output_dn

  rdz=1.0/dz

  !------------------------------------------------------------------------------!
  ! read cm1 dimensions
  ! !
  !------------------------------------------------------------------------------!

  write(*,*)
  write(*,*) 'opening ',infile
  ncid=ncopn(infile,NCNOWRIT,rcode)

  call cm1_nc_getdims( ncid )

  write(*,*) 'ni,nj,nk=',ni,nj,nk

  !------------------------------------------------------------------------------!
  ! get spatial dimesions                                                        !
  !------------------------------------------------------------------------------!

  allocate( xh(1:ni) )
  allocate( xf(1:nip1) )
  allocate( yh(1:nj) )
  allocate( yf(njp1) )
  allocate( zh(1:nk) )
  allocate( zf(1:nk+1) )
  allocate( mzh(1:nk) )   ! transformation between computational space and stretched
  allocate( mzf(1:nk+1) ) ! coordinate grid

  call cm1_nc_getspace( ncid )
  call cm1_nc_gettime( ncid )

  !------------------------------------------------------------------------------!
  ! get kinematics variables                                                     !
  !------------------------------------------------------------------------------!

  allocate( u0(1:nk) )
  allocate( v0(1:nk) )

  allocate( u(ni+1,nj,nk) )
  allocate( v(ni,nj+1,nk) )
  allocate( w(ni,nj,nk+1) )

  allocate( ui(ni,nj,nk) )
  allocate( vi(ni,nj,nk) )
  allocate( wi(ni,nj,nk) )

  allocate( upert(ni+1,nj,nk) )
  allocate( vpert(ni,nj+1,nk) )

  call dyn_nc_getwinds( ncid, u0, v0, u, v, w, ui, vi, wi, upert, vpert )

  !------------------------------------------------------------------------------!
  ! get thermodynamic variables

  allocate( rho0(nk) )

  allocate( th(ni,nj,nk) )
  allocate( th0(nk) )
  allocate( thpert(ni,nj,nk) )

  allocate( thr(ni,nj,nk) )
  allocate( thr0(0:nk+1) )
  allocate( thrpert(ni,nj,nk) )

  allocate( ppi(ni,nj,nk) )
  allocate( pipert(ni,nj,nk) )
  allocate( pi0(nk) )

  allocate( prs(ni,nj,nk) )
  allocate( prs0(nk) )
  allocate( prspert(ni,nj,nk) )

  allocate( qv(ni,nj,nk) )
  allocate( qv0(nk) )
  allocate( qc(ni,nj,nk) )
  allocate( qr(ni,nj,nk) )
  allocate( qi(ni,nj,nk) )
  allocate( qg(ni,nj,nk) )
  allocate( qhl(ni,nj,nk) )
  allocate( qs(ni,nj,nk) )
  allocate( qt(ni,nj,nk) )

  allocate( buoy(ni,nj,0:nk+1) )

  call cm1_nc_getthermo( ncid, th, th0, thpert, ppi, pi0, pipert, &
    prs, prs0, prspert, qv, qv0, qc, qr, qs, qi, qg, qhl, qt, imoist )

  do k=1,nk
    !thr0(k) = th0(k)*(1+reps*qv0(k))/(1+qv0(k)) !hastings version
    thr0(k) = th0(k)*(1+reps*qv0(k))  !markowski version

    rho0(k) = prs0(k)/(rd*thr0(k)*pi0(k))
    do i=1,ni
    do j=1,nj
      
     ! thr(i,j,k) = th(i,j,k)*(1+reps*qv(i,j,k))/(1+qv(i,j,k)+qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qg(i,j,k)) !hastings version
      thr(i,j,k) = th(i,j,k)*(1.+reps*qv(i,j,k)-(qc(i,j,k)+qr(i,j,k)+qi(i,j,k)+qs(i,j,k)+qg(i,j,k))) ! markowski version
      buoy(i,j,k) = g*(thr(i,j,k)/thr0(k) - 1.0 )
      thrpert(i,j,k) = thr(i,j,k)-thr0(k);

    enddo
    enddo
  enddo
  do i=1,ni
  do j=1,nj
   buoy(i,j,0)=-buoy(i,j,1) ! we are assuming the buoyancy is constant below the
                           ! the lowest scalar grid point.  this is because,
                           ! when computing p_b, we get physically unrealistic
                           ! results near the surface, because B is forced to be
                           ! zero.  not only does this result in an
                           ! overestimation of dB/dz and thus p_b, but, because
                           ! we are computing p_dn from p'-p_dl-p_b, it leads to
                           ! anomalously negative values in p_dn
    buoy(i,j,nk+1)=-buoy(i,j,nk)
  enddo
  enddo

  thr0(0)=thr0(1)
  thr0(nk+1)=thr0(nk)

  ! close netcdf file
  write(*,*) 'closing ',infile
  call ncclos( ncid, rcode )

  ! set up arrays for elliptical solver
  allocate( cfa(nk) )
  allocate( cfc(nk) )
  allocate( cfb(ni,nj,nk) )
  allocate( ad1(nk) )
  allocate( ad2(nk) )
  allocate( rf0(0:nk+1) )

  rf0(1:nk+1)=0.0
  do k=2,nk
    rf0(k) = 0.5*(rho0(k)+rho0(k-1))
  enddo
!  rf0(1)=rho0(1)     ! note that this produces incorrect results
!  rf0(nk+1)=rho0(nk)
  rf0(1)=0.0
  rf0(nk+1)=0.0

  do k=1,nk
    cfa(k) = mzh(k)*mzf(k)*rf0(k)*0.5*(thr0(k-1)+thr0(k))/ &
           & (dz*dz*rho0(k)*thr0(k))
    cfc(k) = mzh(k)*mzf(k+1)*rf0(k+1)*0.5*(thr0(k)+thr0(k+1)) / &
           & (dz*dz*rho0(k)*thr0(k))
    ad1(k) = 1.0/(cp*rho0(k)*thr0(k))
    ad2(k) = 1.0
  enddo

  do k=1,nk
  do j=1,nj
  do i=1,ni
    cfb(i,j,k) = 2.0*( cos(pi*(i-1)/ni) + cos(pi*(j-1)/nj) - &
               & 2.0 )*rdz*rdz - cfa(k) - cfc(k)
  enddo
  enddo
  enddo


  !##############################################################################!
  !#                               COMPUTE DYNAMICS                             #!
  !##############################################################################!

  !------------------------------------------------------------------------------!
  ! initialize arrays

  ! forcing and perturbation arrays
  allocate( ff(ni,nj,nk) ) ! forcing function

  ff(:,:,:) = 0.0

  allocate( pi_dl(ni,nj,nk) )
  allocate( p_dl(ni,nj,nk) )

  p_dl(:,:,:) = 0.0

  allocate( pi_dn(ni,nj,nk) )
  allocate( p_dn(ni,nj,nk) )

  p_dn(:,:,:) = 0.0

  allocate( pi_b(ni,nj,nk) )
  allocate( p_b(ni,nj,nk) )

  p_b(:,:,:) = 0.0

  allocate( pgfb(ni,nj,nk+1,3) )
  allocate( pgfdn(ni,nj,nk+1,3) )
  allocate( pgfdl(ni,nj,nk+1,3) )

  write(*,*)
  write(*,*) 'and now for some dynamical calculations:'
  write(*,*)

 
  !------------------------------------------------------------------------------!
  ! compute linear dynamic pressure                                              !
  !------------------------------------------------------------------------------!

        IF(OUTPUT_DL.EQ.1)THEN

  !-----------------------------------------------------------------------------!
  ! compute forcing function for linear dynamic perturbation pressure

  allocate( Sx(nk+1) )
  allocate( Sy(nk+1) )
  allocate( wx(ni,nj,nk+1) )
  allocate( wy(ni,nj,nk+1) )

  ! vertical shear of horizontal wind
  do k=2,nk
    Sx(k) = (u0(k)-u0(k-1))*rdz*mzf(k)
    Sy(k) = (v0(k)-v0(k-1))*rdz*mzf(k)
  enddo
  Sx(1) = 0.0
  Sy(1) = 0.0
  Sx(nk+1) = 0.0
  Sy(nk+1) = 0.0

  ! horizontal shear of vertical wind
  wx(:,:,:)=0.0
  wy(:,:,:)=0.0

  do j=1,nj
  do k=1,nk+1
    wx(1,j,k)  = Sx(k)*( w(2,j,k)  - w(   1,j,k) )*rdz
    wx(ni,j,k) = Sx(k)*( w(ni,j,k) - w(ni-1,j,k) )*rdz
    do i=2,ni-1
      wx(i,j,k) = Sx(k)*0.5*( w(i+1,j,k) - w(i-1,j,k) )*rdz
    enddo
  enddo
  enddo
 
  do k=1,nk+1
  do i=1,ni
    wy(i,1,k) = Sy(k)*( w(i,2,k) - w(i,1,k) )*rdz
    wy(i,nj,k) = Sy(k)*( w(i,nj,k) - w(i,nj-1,k) )*rdz
    do j=2,nj-1
      wy(i,j,k) = Sy(k)*0.5*( w(i,j+1,k) - w(i,j-1,k) )*rdz
    enddo
  enddo
  enddo

  ! set up forcing function
  do i=1,ni
  do j=1,nj
  do k=1,nk
    ff(i,j,k) = -rho0(k)*( (wx(i,j,k)+wx(i,j,k+1)) + &
                         & (wy(i,j,k)+wy(i,j,k+1)) )
  enddo
  enddo
  enddo


  !-----------------------------------------------------------------------------!
  ! solve for linear dynamic perturbation pressure

  write(*,*) 'solving for linear dynamic perturbation pressure'
  call poiss( ff, cfb, cfa, cfc, ad1, ad2, pi_dl )

  ! adjust so that the sum over the domain is zero
  sumpb=0.0
  do i=1,ni
  do j=1,nj
  do k=1,nk
    sumpb=sumpb+pi_dl(i,j,k)
  enddo
  enddo
  enddo
  pi_dl = pi_dl - sumpb/(ni*nj*nk)

  write(*,*)

  deallocate( Sx )
  deallocate( Sy )
  deallocate( wx )
  deallocate( wy )

  ! convert exner function to Pascals
  do i=1,ni
  do j=1,nj
  do k=1,nk
    p_dl(i,j,k) = prs(i,j,k) - p00*((ppi(i,j,k)-pi_dl(i,j,k))**rkappa)
  enddo
  enddo
  enddo

  write(*,*) 'done with linear dynamic perturbation pressure'

  ! compute pgf
  write(*,*) 'compute pgf'
  call compute_pgfxyz( thr, pi_dl, pgfdl )

          ENDIF

  !-----------------------------------------------------------------------------!
  ! buoyancy perturbation pressure                                              !
  !-----------------------------------------------------------------------------!


          IF(OUTPUT_B.EQ.1)THEN

  !-----------------------------------------------------------------------------!
  ! set up forcing function and boundary conditions

  ff(:,:,:) = 0.0

  allocate( buoyf(ni,nj,nk+1) )
  do i=1,ni
  do j=1,nj
  do k=1,nk+1
    buoyf(i,j,k) = 0.5*( buoy(i,j,k-1)+buoy(i,j,k) )
  enddo
  enddo
  enddo

  do i=1,ni
  do j=1,nj
  do k=1,nk
    ff(i,j,k) = ( rf0(k+1)*buoyf(i,j,k+1) - rf0(k)*buoyf(i,j,k) )*rdz*mzh(k)
  enddo
  enddo
  enddo

  !-----------------------------------------------------------------------------!
  ! solve poisson equation

  write(*,*) 'solving for buoyancy pressure perturbation'

  call poiss( ff, cfb, cfa, cfc, ad1, ad2, pi_b )

  write(*,*) maxval(pi_b)

  ! adjust so sum is zero
  sumpb = 0.0
  do i=1,ni
  do j=1,nj
  do k=1,nk
    sumpb = sumpb+pi_b(i,j,k)
  enddo
  enddo
  enddo
  pi_b = pi_b - sumpb/(ni*nj*nk)

  ! convert to Pascals
  do i=1,ni
  do j=1,nj
  do k=1,nk
    p_b(i,j,k) = prs(i,j,k)-p00*((ppi(i,j,k)-pi_b(i,j,k))**rkappa)
  enddo
  enddo
  enddo

  write(*,*) 'solved buoyancy pressure.'

!stop

  ! compute PGF
  write(*,*) 'computing pgf..'
  call compute_pgfxyz( thr, pi_b, pgfb )
  

          ENDIF

  !-----------------------------------------------------------------------------!
  ! compute nonlinear dynamic perturbation pressure                             !
  !-----------------------------------------------------------------------------!


          IF(OUTPUT_DN.EQ.1)THEN

  write(*,*) 'nonlinear dynamic pressure'


      if((output_dl.eq.1).and.(output_b.eq.1))then

  do i=1,ni
  do j=1,nj
  do k=1,nk
    pi_dn(i,j,k) = ppi(i,j,k) - pi0(k) - pi_b(i,j,k) - pi_dl(i,j,k)
  enddo
  enddo
  enddo

  do i=1,ni
  do j=1,nj
  do k=1,nk
    p_dn(i,j,k) = prs(i,j,k) - p00*((ppi(i,j,k)-pi_dn(i,j,k))**rkappa)
  enddo
  enddo
  enddo

  ! compute pgf
  call compute_pgfxyz( thr, pi_dn, pgfdn )

      else
        write(*,*) 'error:  cannot compute nonlinear dynamic pressure without linear and buoyancy pressures'
        stop
      endif


          ENDIF


  !#############################################################################!
  !#                     WRITE OUT RESULTS IN NETCDF FILE                      #!
  !#############################################################################!



  write(*,*) '-----------------------------------------------------'
  write(*,*) 'writing out ',outfile

  call nc_dyn_writeout( outfile, thr0, thrpert, &
        pi_dl, p_dl, pi_b, p_b, pi_dn, p_dn, &
        pgfb, pgfdl, pgfdn )


  ! deallocation block

  deallocate( xh )
  deallocate( xf )
  deallocate( yh )
  deallocate( yf )
  deallocate( zh )
  deallocate( zf )
  deallocate( mzh )
  deallocate( mzf )

  deallocate( u0 )
  deallocate( v0 )

  deallocate( u )
  deallocate( v )
  deallocate( w )

  deallocate( ui )
  deallocate( vi )
  deallocate( wi )

  deallocate( rho0 )

  deallocate( th )
  deallocate( th0 )
  deallocate( thpert )
  deallocate( thr )
  deallocate( thr0 )
  deallocate( thrpert )

  deallocate( ppi )
  deallocate( pipert )
  deallocate( pi0 )
  deallocate( prs )
  deallocate( prs0 )
  deallocate( prspert )

  deallocate( buoy )

  deallocate( cfa )
  deallocate( cfc )
  deallocate( rf0 )
  deallocate( buoyf )

  !deallocate( prs_h )
  !deallocate( prs_nh )

  deallocate( ff )

  deallocate( p_dl )

  deallocate( p_dn )

  deallocate( p_b )

  write(*,*)
  write(*,*) 'program terminated successfully!  w00t!'
  write(*,*)

contains
  !=============================================================================!
  !#############################################################################!
  !################################ SUBROUTINES ################################!
  !#############################################################################!
  !=============================================================================!
  !#############################################################################!
  !#                          SUBROUTINE POISS                                 #!
  !#                                                                           #!
  !#  Solves pressure equations.  Based on George Bryan's poiss.F code, except #!
  !#  using a DCT rather than DFT.                                             #!
  !#---------------------------------------------------------------------------#!
  !# Ryan Hastings, 13 Mar 2012                                                #!
  !#############################################################################!
  subroutine poiss( f, cfb, cfa, cfc, d1, d2, ppi )
    implicit none
    include 'fftw3.f'
    !###########################################################################!
    !################################ VARIABLES ################################!
    !---------------------------- input variables ------------------------------!
    real,dimension(ni,nj,nk) :: f, cfb
    real,dimension(nk)       :: cfa, cfc, d1, d2

    !--------------------------- output variables ------------------------------! 
    real,dimension(ni,nj,nk) :: ppi

    !----------------------------- local variables -----------------------------!
    real,dimension(0:nk+1) :: r1

    real,dimension(ni,nj,nk) :: pdt, deft
    real*8,dimension(ni,nj) :: rhs, trans
    real,dimension(nk) :: lgbth, lgbph

    integer i, j, k, ii, jj, kk, wi, wj, pdti, pdtj

    integer*8 :: plan


    ! DCT the total forcing

    write(*,*) 'DCT the forcing'
    do k=1,nk
      do i=1,ni
      do j=1,nj
          rhs(i,j) = dble(f(i,j,k)*d1(k))
      enddo
      enddo


      !DCT
      call dfftw_plan_r2r_2d(plan,ni,nj,rhs,trans,FFTW_REDFT10,FFTW_REDFT10,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)

      do j=1,nj
      do i=1,ni
        deft(i,j,k) = trans(i,j) 
      enddo
      enddo

    enddo

    ! solve tri-diagonal matrix
    write(*,*) 'solve the matrix'

    do j=1,nj
    do i=1,ni

      if((i.eq.1).and.(j.eq.1))then
        r1(nk+1)=0.0
        r1(nk)=0.0
        do k=nk,2,-1
          r1(k-1)=(deft(i,j,k)-cfc(k)*r1(k+1)-cfb(i,j,k)*r1(k))/cfa(k)
        enddo
        do k=1,nk
           pdt(i,j,k)=r1(k)
        enddo
      else
        lgbth(1)=-cfc(1)/cfb(i,j,1)
        lgbph(1)= deft(i,j,1)/cfb(i,j,1)
        do k=2,nk
          lgbth(k)=-cfc(k)/(cfa(k)*lgbth(k-1)+cfb(i,j,k))
          lgbph(k)=(deft(i,j,k)-cfa(k)*lgbph(k-1))/ &
                  &(cfa(k)*lgbth(k-1)+cfb(i,j,k))
        enddo
        pdt(i,j,nk)=lgbph(nk)
        do k=nk-1,1,-1
          pdt(i,j,k)=lgbth(k)*pdt(i,j,k+1)+lgbph(k)
        enddo
      endif

    enddo
    enddo

    ! reverse transform
    write(*,*) 'invert DCT'
    do k=1,nk

      do i=1,ni
      do j=1,nj
        rhs(i,j) = dble(pdt(i,j,k))
      enddo
      enddo

      !DCT
      call dfftw_plan_r2r_2d(plan,ni,nj,rhs,trans,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)

      do j=1,nj
      do i=1,ni
        ppi(i,j,k)=real(trans(i,j))*d2(k)/(4*ni*nj)
      enddo
      enddo

    enddo

  end subroutine

  !#############################################################################!
  !#                      SUBROUTINE COMPUTE_PGFXYZ                            #!
  !#                                                                           #!
  !#  Computes PGF.  Output is a four-dimensional variable, with each of the   #!
  !#  three components given in the last dimension.                            #!
  !#---------------------------------------------------------------------------#!
  !# Ryan Hastings, 13 Mar 2012                                                #!
  !#############################################################################!
  subroutine compute_pgfxyz( thr, ppi, pgf )

    implicit none
    real,intent(in),dimension(ni,nj,nk) :: thr
    real,intent(in),dimension(ni,nj,nk) :: ppi
    real,intent(out),dimension(ni,nj,nk+1,3) :: pgf

    integer i, j, k
    real tem

    pgf(:,:,:,:) = 0.0

    tem = -cp*rdz*0.5


    do i=2,ni
    do j=1,nj
    do k=1,nk
      pgf(i,j,k,1) = tem*(thr(i,j,k)+thr(i-1,j,k))*(ppi(i,j,k)-ppi(i-1,j,k))
    enddo
    enddo
    enddo

    write(*,*) 'pgfx'

    do i=1,ni
    do j=2,nj
    do k=1,nk
      pgf(i,j,k,2) = tem*(thr(i,j,k)+thr(i,j-1,k))*(ppi(i,j,k)-ppi(i,j-1,k))
    enddo
    enddo
    enddo

    write(*,*) 'pgfy'

    do i=1,ni
    do j=1,nj
      do k=2,nk
        pgf(i,j,k,3) = tem*(thr(i,j,k)+thr(i,j,k-1))*(ppi(i,j,k)-ppi(i,j,k-1))*mzf(k)
      enddo
    enddo
    enddo

    write(*,*) 'pgfz'

  end subroutine


end program
