program mass_flux
  !########################################################################!
  !########################################################################!
  !##                                                                    ##!
  !##                         MASS_FLUX                                  ##!
  !##                                                                    ##!
  !## Computes updraft, downdraft, and total mass flux at user-specified ##!
  !## grid height for a number of files.                                 ##!
  !##====================================================================##!
  !## v1.0, Ryan Hastings, 19 Jan 2017                                   ##!
  !########################################################################!
  use globals
  use nc_inout
  implicit none

  !##############################################!
  !##### SPECIFICATION PART #####################!
  !##############################################!

  !----------------------------------------------!
  ! configuration variables
  character(len=100) :: outfile ! name of output file
  integer t1, t2 ! time index for first and last files
  integer tstep  ! step for time indices
  integer gpx, gpy ! number of grid points in x- and y-direction for wholedomain=0
  integer kud, kdd ! vertical grid points at which updraft (kud) and downdraft (kdd)
                   ! flux are to be taken
  real dx, dy      ! horizontal grid spacing
  real toffset     ! time offset (for cm1r16 and earlier)
  integer wholedomain ! compute over whole domain, or just part?

  namelist /configuration/ outfile, t1, tstep, t2, gpx, gpy, dx, dy, kud, kdd, toffset, wholedomain, ptype

  !----------------------------------------------!
  ! netcdf variables
  integer ncid, rcode
  character(len=100) :: filename
  character(len=6) :: timetag
  integer ntid, tid, mfluid, mfldid, mfltid

  !----------------------------------------------!
  ! dimensions
  integer nt
  integer,dimension(2) :: wmax_loc
  integer imax, i1, i2, jmax, j1, j2, tt, t

  !----------------------------------------------!
  ! field variables
  real,dimension(:,:,:),allocatable :: qv, qr, qi, qs, qc, qg, qhl
  real,dimension(:,:,:),allocatable :: rho, rhot, mflux_all, w

  !----------------------------------------------!
  ! output variables
  real,dimension(:),allocatable :: times, mflux_up, mflux_dn, mflux_to

  !##############################################!
  !########### EXECUTION PART ###################!
  !##############################################!

  open(8,file='mass_flux.input')
  read(8,nml=configuration)
  close(8)

  nt=t2-t1+1
  allocate( times(nt) )
  allocate( mflux_up(nt) )
  allocate( mflux_dn(nt) )
  allocate( mflux_to(nt) )

  write(timetag,'(I6.6)') t1
  filename='cm1out_'//timetag//'.nc'
  ncid=ncopn(filename,NCNOWRIT,rcode)
  call cm1_nc_getdims( ncid )
  call ncclos(ncid,rcode)

  allocate( mflux_all(ni,nj,nk) )
  allocate( rho(ni,nj,nk) )
  allocate( qv(ni,nj,nk) )
  allocate( qc(ni,nj,nk) )
  allocate( qr(ni,nj,nk) )
  allocate( qi(ni,nj,nk) )
  allocate( qs(ni,nj,nk) )
  allocate( qg(ni,nj,nk) )
  allocate( qhl(ni,nj,nk) )
  allocate( rhot(ni,nj,nk) )
  allocate( w(ni,nj,nk) )

  do t=1,nt!t1,t2,tstep

    tt=t1+tstep*(t-1)
    write(timetag,'(I6.6)') tt
    filename='cm1out_'//timetag//'.nc'

    write(*,*) 'opening ',filename
    ncid=ncopn(filename,NCNOWRIT,rcode)

    call cm1_nc_gettime( ncid )
    times(t) = time+toffset
    call read_cm1_var( ncid, 'rho     ', rho )
    call read_cm1_var( ncid, 'qv      ', qv )
    call read_cm1_var( ncid, 'qc      ', qc )
    call read_cm1_var( ncid, 'qr      ', qr )
    call read_cm1_var( ncid, 'qi      ', qi )
    call read_cm1_var( ncid, 'qg      ', qg )
    if(ptype.eq.27)then
      call read_cm1_var( ncid, 'qhl     ', qhl )
    endif
    call read_cm1_var( ncid, 'qs      ', qs )
    call read_cm1_var( ncid, 'winterp ', w  )

    call ncclos(ncid,rcode)
    !---------------- mass flux -------------------!
    ! compute total density
    if(ptype.eq.27)then
      rhot = rho*(1+qv+qc+qr+qi+qg+qs+qhl)
    else
      rhot = rho*(1+qv+qc+qr+qi+qg+qs)
    endif
    ! compute total mass flux
    mflux_all=rhot*w*dx*dy

      IF(WHOLEDOMAIN.EQ.0)THEN
    ! maximum w
    wmax_loc = maxloc( w(:,:,kud) )
    imax=wmax_loc(1)
    i1=max(imax-gpx,1)
    i2=min(imax+gpx,ni)
    jmax=wmax_loc(2)
    j1=max(1,jmax-gpy)
    j2=min(nj,jmax+gpy)

    ! mass fluxes
    mflux_up(t) = sum( mflux_all(i1:i2,j1:j2,kud), mflux_all(i1:i2,j1:j2,kud).gt.0.0 )
    mflux_dn(t) = sum( mflux_all(i1:i2,j1:j2,kdd), mflux_all(i1:i2,j1:j2,kdd).lt.0.0 )
    mflux_to(t) = sum( mflux_all(i1:i2,j1:j2,kud) )
      ELSEIF(WHOLEDOMAIN.EQ.1)THEN
    mflux_up(t) = sum( mflux_all(:,:,kud), mflux_all(:,:,kud).gt.0.0 )
    mflux_dn(t) = sum( mflux_all(:,:,kdd), mflux_all(:,:,kdd).lt.0.0 )
    mflux_to(t) = sum( mflux_all(:,:,kud) )
      ENDIF
    

  enddo

  !--! WRITE OUT RESULTS
  write(*,*) 'writing out ',outfile
  ncid=nccre(outfile,NCNOCLOB,rcode)

  ntid=ncddef(ncid,'ntimes',nt,rcode)

  tid=ncvdef(ncid,'time',nf_real,1,ntid,rcode)
  mfluid=ncvdef(ncid,'updraft_mass_flux',nf_real,1,ntid,rcode)
  mfldid=ncvdef(ncid,'downdraft_mass_flux',nf_real,1,ntid,rcode)
  mfltid=ncvdef(ncid,'total_mass_flux',nf_real,1,ntid,rcode)
  call ncendf(ncid,rcode)

  call ncvpt(ncid,tid,1,nt,times,rcode)
  call ncvpt(ncid,mfluid,1,nt,mflux_up,rcode)
  call ncvpt(ncid,mfldid,1,nt,mflux_dn,rcode)
  call ncvpt(ncid,mfltid,1,nt,mflux_to,rcode)

  call ncclos(ncid,rcode)

contains

subroutine read_cm1_var( ncid, varname, var )
  implicit none
  integer,intent(in) :: ncid
  character(len=8),intent(in) :: varname
  real,dimension(ni,nj,nk),intent(out) :: var

  real,dimension(ni,nj,nk,1) :: dummyvar
  integer rcode, varid, i, j, k

  varid = ncvid( ncid, varname, rcode )
  call ncvgt( ncid, varid, (/1,1,1,1/), (/ni,nj,nk,1/), dummyvar, rcode )
  do i=1,ni
  do j=1,nj
  do k=1,nk
    var(i,j,k) = dummyvar(i,j,k,1)
  enddo
  enddo
  enddo

end subroutine

subroutine read_cm1_z( ncid )
  implicit none
  integer,intent(in) :: ncid
  real,dimension(nk+1) :: zf
  integer zfid, rcode

  zfid=ncvid(ncid,'zf      ',rcode)
  call ncvgt(ncid,zfid,1,nk+1,zf,rcode)

end subroutine

end program
