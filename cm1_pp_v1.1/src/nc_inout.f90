module nc_inout
  !#############################################################################!
  !############################################################################!
  !####                                                                    ####!
  !####                          MODULE NC_INOUT                           ####!
  !####                                                                    ####!
  !####  Module containing subroutines for input and output from netCDF    ####!
  !####  files useful for CM1_PP.                                          ####!
  !####====================================================================####!
  !####  Ryan Hastings, 7 September 2011                                   ####!
  !####                                                                    ####!
  !############################################################################!
  !####                                                                    ####!
  !####  SUBROUTINE INDEX:                                                 ####!
  !####                                                                    ####!
  !####      cm1_nc_getdims             : get dimensions                   ####!
  !####                                                                    ####!
  !####      cm1_nc_getspace            : get spatial variables            ####!
  !####                                                                    ####!
  !####      cm1_nc_gettime             : get time for DYN                 ####!
  !####                                                                    ####!
  !####      cm1_nc_gettime_traj        : get time for TRAJ                ####!
  !####                                                                    ####!
  !####      cm1_nc_getwinds            : get winds                        ####!
  !####                                                                    ####!
  !####      cm1_nc_getthermo           : get thermodynamic variables      ####!
  !####                                                                    ####!
  !####      cm1_nc_getvar3             : get generic 3d variable          ####!
  !####                                                                    ####!
  !####      cm1_nc_getvar4             : get generic 4d var, store as 3d  ####!
  !####                                                                    ####!
  !####      nc_dyn_writeout            : write variables for dyn          ####!
  !####                                                                    ####!
  !####      nc_handle_error            : handle error for netcdf          ####!
  !####                                                                    ####!
  !####      dyn_nc_getwinds            : get winds for dyn                ####!
  !####                                                                    ####!
  !####      nc_traject_vardef          : define variables for traj output ####!
  !####                                                                    ####!
  !####      nc_circuit_vardef          : define variables for circuit     ####!
  !####                                                                    ####!
  !####      traject_nc_read_cm1        : read CM1 file for traj           ####!
  !####                                                                    ####!
  !####      traject_nc_read_dyn        : read dyn output file for traj    ####!
  !####                                                                    ####!
  !####      nc_readvar4                : see description below            ####!
  !####                                                                    ####!
  !####      nc_traject_writeout        : write out traj output file       ####!
  !####                                                                    ####!
  !####--------------------------------------------------------------------####!
  !#### Modified by Ryan Hastings on 26 September 2012                     ####!
  !####      added nc_error and modified existing subroutines to use       ####!
  !####      F77 netcdf libraries (i.e., nf_inq_dimlen rather than ncdinq, ####!
  !####      which doesn't seem to be working)                             ####!
  !####--------------------------------------------------------------------####!
  !#### v0.8 28 Mar 2013                                                   ####!
  !####--------------------------------------------------------------------####!
  !#### v1.0 13 Dec 2016                                                   ####!
  !############################################################################!
  !############################################################################!
  use globals
  implicit none
  include 'netcdf.inc'

  contains

    !###########################################################################!
    !#                    SUBROUTINE CM1_NC_GETTIME_TRAJ                       #!
    !#                                                                         #!
    !#  Gets the time from a cm1 output file.                                  #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 9 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_gettime_traj( filename, tt )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ------------------------------!
      character(len=100) :: filename
      !------------------------- output variables ------------------------------!
      real tt

      !-------------------------- local variables ------------------------------!
      integer ncid, rcode, timedid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      ncid = ncopn( filename, NCNOWRIT, rcode )
      timedid = ncvid( ncid, 'time', rcode )
      call ncvgt( ncid, timedid, 1, 1, tt, rcode )
      call ncclos( ncid, rcode )

    end subroutine

    !#########################################################################!
    !#                        SUBROUTINE NC_HANDLE_ERROR                     #!
    !#                                                                       #!
    !#  Handle error.                                                        #!
    !#-----------------------------------------------------------------------#!
    !# Ryan Hastings, 26 September, 2012                                     #!
    !#########################################################################!
      subroutine nc_handle_error( retval )
      implicit none

      integer retval

      if (retval.ne.nf_noerr) then
        write(*,*) 'error: ',nf_strerror(retval)
        stop 2
      endif

    end subroutine



    !###########################################################################!
    !#                        SUBROUTINE CM1_NC_GETDIMS                        #!
    !#                                                                         #!
    !#  Get dimensions from CM1 output file.                                   #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_getdims( ncid,vcm1 )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!

      !--------------------------- input variables -----------------------------!
      integer ncid
      real vcm1

      !--------------------------- local variables -----------------------------!
      integer rcode, nidid, njdid, nkdid, nip1did, njp1did, nkp1did
      character(len=2) char2
      character(len=4) char4

      !#########################################################################!
      !############################## MAIN BODY ################################!

      write(*,*)
      write(*,*) 'reading dimensions'
      write(*,*)

      if (vcm1.lt.20) then
        nidid   = ncdid(ncid,'ni',rcode)
        nip1did = ncdid(ncid,'nip1',rcode)
        njdid   = ncdid(ncid,'nj',rcode)
        njp1did = ncdid(ncid,'njp1',rcode)
        nkdid   = ncdid(ncid,'nk',rcode)
        nkp1did = ncdid(ncid,'nkp1',rcode)
      else
        nidid   = ncdid(ncid,'xh',rcode)
        nip1did = ncdid(ncid,'xf',rcode)
        njdid   = ncdid(ncid,'yh',rcode)
        njp1did = ncdid(ncid,'yf',rcode)
        nkdid   = ncdid(ncid,'zh',rcode)
        nkp1did = ncdid(ncid,'zf',rcode)
      endif

      call nc_handle_error( nf_inq_dimlen(ncid,nidid,ni) )
      call nc_handle_error( nf_inq_dimlen(ncid,nip1did,nip1) )
      call nc_handle_error( nf_inq_dimlen(ncid,njdid,nj) )
      call nc_handle_error( nf_inq_dimlen(ncid,njp1did,njp1) )
      call nc_handle_error( nf_inq_dimlen(ncid,nkdid,nk) )
      call nc_handle_error( nf_inq_dimlen(ncid,nkp1did,nkp1) )

      write(*,*)
      write(*,*) 'ni,nj,nk=',ni,nj,nk
      write(*,*) 'nip1,njp1,nkp1=',nip1,njp1,nkp1
      write(*,*)

      return

    end subroutine
    !===========================================================================!
    !###########################################################################!
    !#                        SUBROUTINE CM1_NC_GETTIME                        #!
    !#                                                                         #!
    !#  Gets the time from a cm1 output file.                                  #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 9 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_gettime( ncid )
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ------------------------------!
      integer ncid

      !-------------------------- local variables ------------------------------!
      integer rcode, timedid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      timedid = ncvid( ncid, 'time', rcode )
      call ncvgt( ncid, timedid, 1, 1, time, rcode )

    end subroutine
   !===========================================================================!
    !===========================================================================!
    !###########################################################################!
    !#                        SUBROUTINE CM1_NC_GETSPACE                       #!
    !#                                                                         #!
    !#  Reads in spacial variables from cm1 file:  xh, yh, zh, zf.  Calculates #!
    !# mzh and mfh, too.  These are the Jacobian functions. 'h' indicates pos- #!
    !# ition on mass grid, 'f' on momentum grid.  mzh(k) is the Jacobian for   #!
    !# the kth point on zh, mzf for zf.                                        #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 9 Sep 2011                                               #!
    !###########################################################################!
    subroutine cm1_nc_getspace( ncid,vcm1 )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ------------------------------!
      integer ncid
      real vcm1

      !-------------------------- local variables ------------------------------! 
      integer i, j, k
      integer rcode
      integer xhid, xfid, yhid, yfid, zhid, zfid, timeid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      write(*,*)
      write(*,*) 'reading in spatial variables'
      write(*,*)

      !-------------------------------------------------------------------------!
      ! get xh, yh, zh, zf

            write(*,*) '....xh....'
      if (vcm1.gt.20) then
        xhid = ncvid(ncid, 'xh', rcode)
      else
        xhid = ncvid(ncid, 'ni', rcode)
      endif
      call ncvgt( ncid, xhid, 1,   ni, xh(1:ni), rcode )

      write(*,*) '....xf....'
      if (vcm1.gt.20) then
        xfid = ncvid(ncid, 'xf', rcode)
      else
        xfid = ncvid(ncid, 'nip1', rcode)
      endif
      call ncvgt( ncid, xfid, 1, nip1, xf(1:nip1), rcode )

      if(nj.eq.1)then
        yh(1)=0.5
        yh(2)=1.5
        yf(1)=0.0
        yf(2)=1.0
        yf(3)=2.0
      else

        if (vcm1.gt.20) then

          write(*,*) '....yh....'
          yhid = ncvid(ncid, 'yh', rcode)
          call ncvgt( ncid, yhid, 1, nj,   yh(1:nj),   rcode )

          write(*,*) '....yf....'
          yfid = ncvid(ncid, 'yf', rcode)
          call ncvgt( ncid, yfid, 1, njp1, yf(1:njp1),   rcode )

          write(*,*) '....zh....'
          zhid = ncvid(ncid, 'zh', rcode)
          call ncvgt( ncid, zhid, 1, nk,   zh(1:nk),   rcode )

          write(*,*) '....zf....'
          zfid = ncvid(ncid, 'zf', rcode)
          call ncvgt( ncid, zfid, 1, nk+1, zf(1:nkp1), rcode )

        else

          write(*,*) '....yh....'
          yhid = ncvid(ncid, 'nj', rcode)
          call ncvgt( ncid, yhid, 1, nj,   yh(1:nj),   rcode )

          write(*,*) '....yf....'
          yfid = ncvid(ncid, 'njp1', rcode)
          call ncvgt( ncid, yfid, 1, njp1, yf(1:njp1),   rcode )

          write(*,*) '....zh....'
          zhid = ncvid(ncid, 'nk', rcode)
          call ncvgt( ncid, zhid, 1, nk,   zh(1:nk),   rcode )

          write(*,*) '....zf....'
          zfid = ncvid(ncid, 'nkp1', rcode)
          call ncvgt( ncid, zfid, 1, nk+1, zf(1:nkp1), rcode )

        endif

        ! convert km to m
        xh=xh*1000
        xf=xf*1000
        yh=yh*1000
        yf=yf*1000
        zh=zh*1000
        zf=zf*1000

      endif

      ! calculate mzh and mzf
      write(*,*) '....mzh....'
      do k=1,nk
        mzh(k) = dz/(zf(k+1)-zf(k))
      enddo

      do k=2,nk
        mzf(k) = dz/(zh(k)-zh(k-1))
      enddo
      mzf(1)=1.0
      mzf(nk+1)=1.0

      ! write out variables
      write(*,*)
      do i=1,ni
        write(*,*) 'i,xh(i),xf(i)=',i,xh(i),xf(i)
      enddo
      write(*,*) 'i,xf(i)=',nip1,xf(nip1)
      write(*,*)
      do j=1,nj
        write(*,*) 'j,yh(j),yf(j)=',j,yh(j),yf(j)
      enddo
      write(*,*) 'j,yf(j)=',njp1,yf(njp1)
      write(*,*)
      do k=1,nk
        write(*,*) 'k,zh(k),zf(k)=',k,zh(k),zf(k)
      enddo
      write(*,*) 'k,zh(k),zf(k)=',nk+1,zf(nk+1)
      write(*,*) '---------------------------------------------------------'

      write(*,*)
      write(*,*) '.......done reading in space'
      write(*,*)

    end subroutine

    !===========================================================================!
    !###########################################################################!
    !#                           SUBROUTINE NC_GETVAR3                         #!
    !#                                                                         #!
    !#  Get a generic 3-dimensional variable from an netCDF file and store it  #!
    !#  in a 3-dimensional variable.                                           #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !###########################################################################!
    subroutine nc_getvar3( ncid, varname, strt, cnt, var )
      !#########################################################################!
      !############################## VARIABLES ################################!

      !-------------------------- input variables ------------------------------!
      integer                 ncid
      character(len=8)        varname
      integer,dimension(3) :: strt, cnt

      !------------------------- output variables ------------------------------!
      real,dimension(ni,nj,nk) :: var

      !-------------------------- local variables ------------------------------!
      integer rcode, varid

      !#########################################################################!
      !############################## MAIN BODY ################################!

      varid = ncvid( ncid, varname, rcode )
      call ncvgt( ncid, varid, strt, cnt, var, rcode )

    end subroutine

    !===========================================================================!
    !###########################################################################!
    !#                          SUBROUTINE NC_GETVAR4                          #!
    !#                                                                         #!
    !#  Get a generic 4-dimensional variable and store it in a 3-dimensional   #!
    !# variable.                                                               #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !###########################################################################!
    subroutine nc_getvar4( ncid, varname, strt, cnt, var )
      !#########################################################################!
      !############################## VARIABLES ################################!

      !-------------------------- input variables ------------------------------!
      integer                  ncid
      character(len=8)         varname
      integer,dimension(4)  :: strt, cnt

      !------------------------- output variables ------------------------------!
      real,dimension(ni,nj,nk) :: var

      !-------------------------- local variables ------------------------------!
      real,dimension(ni,nj,nk,1) :: dummyvar
      integer                       rcode, varid, i, j, k

      !#########################################################################!
      !############################## MAIN BODY ################################!

      varid = ncvid( ncid, varname, rcode )
      call ncvgt( ncid, varid, strt, cnt, dummyvar, rcode )
      do i=1,ni
      do j=1,nj
      do k=1,nk
        var(i,j,k) = dummyvar(i,j,k,1)
      enddo
      enddo
      enddo

    end subroutine

    !===========================================================================!
    !===========================================================================!
    !###########################################################################!
    !#                      SUBROUTINE CM1_NC_GETTHERMO                        #!
    !#                                                                         #!
    !#  Get thermodynamic variables for DYN.                                   #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !#    Modified 22 January 2012                                             #!
    !###########################################################################!
    subroutine cm1_nc_getthermo( ncid, th, th0, thpert, ppi, pi0, pipert, &
      prs, prs0, prspert, qv, qv0, qc, qr, qs, qi, qi2, qg, qhl, qt, imoist, &
      ptype )

      !#########################################################################!
      !############################## VARIABLES ################################!

      !-------------------------- input variables ------------------------------!
      integer ncid, imoist, ptype

      !------------------------- output variables ------------------------------!
      real,dimension(nk)       :: th0, pi0, prs0, thr0, qv0
      real,dimension(ni,nj,nk) :: th, thpert, thr, ppi, pipert, prs, prspert
      real,dimension(ni,nj,nk) :: qv, qc, qr, qi, qi2, qs, qg, qhl, qt
!      real,dimension(ni,nj,nk) :: tke, khh, khv, kmh, kmv
 
      !-------------------------- local variables ------------------------------!
      integer,dimension(4)     :: strt3, cnt3
      integer,dimension(4)     :: strt,  cnt
      integer                  :: i, j, k
!      integer                  :: output_turb

      !#########################################################################!
      !############################## MAIN BODY ################################!

      write(*,*) 'reading in thermodynamic data'
      write(*,*)

      !-------------------------------------------------------------------------!
      ! set up counting arrays
      strt3 = (/ 1, 1,  1, 1/)
      cnt3  = (/ 1,  1, nk, 1/)

      strt  = (/1, 1,  1,  1/)
      cnt   = (/ni, nj, nk,  1/)

      !-------------------------------------------------------------------------!
      ! get variables

      ! potential temperature
      write(*,*) '...potential temperature'
      call nc_getvar4( ncid, 'th      ',  strt,  cnt,  th )
      write(*,*) '...base state potential temperature'
      call nc_getvar3( ncid, 'th0     ', strt3, cnt3, th0 )

      ! pressure
      write(*,*) '...nondimensional pressure...'
      call nc_getvar4( ncid, 'pi      ',  strt,  cnt, ppi )
      write(*,*) '...base state nondimensional pressure'
      call nc_getvar3( ncid, 'pi0     ', strt3, cnt3, pi0 )
      write(*,*) '...pressure'
      call nc_getvar4( ncid, 'prs     ',  strt,  cnt, prs )
      write(*,*) '...base state pressure'
      call nc_getvar3( ncid, 'prs0    ', strt3, cnt3, prs0 )

      qv0(:)=0.0
      qv(:,:,:)=0.0
      qc(:,:,:)=0.0
      qr(:,:,:)=0.0
      qi(:,:,:)=0.0
      qi2(:,:,:)=0.0
      qs(:,:,:)=0.0
      qg(:,:,:)=0.0
      qhl(:,:,:)=0.0
      qt(:,:,:)=0.0

      if (imoist.eq.1) then
        ! mixing ratios

        write(*,*) '...base state water vapor'
        call nc_getvar3( ncid, 'qv0     ', strt3, cnt3, qv0  )
        write(*,*) '...water vapor'
        call nc_getvar4( ncid, 'qv      ',  strt,  cnt,  qv  )
        write(*,*) '...cloud water'
        call nc_getvar4( ncid, 'qc      ',  strt,  cnt,  qc  )
        write(*,*) '...rain water'
        call nc_getvar4( ncid, 'qr      ',  strt,  cnt,  qr  )
        write(*,*) '...cloud ice'
        call nc_getvar4( ncid, 'qi      ',  strt,  cnt,  qi  )

        if (ptype.eq.27) then
          write(*,*) '...graupel'
          call nc_getvar4( ncid, 'qg      ',  strt,  cnt,  qg  )
          write(*,*) '...snow'
          call nc_getvar4( ncid, 'qs      ',  strt,  cnt,  qs  )
          write(*,*) '...hail'
          call nc_getvar4( ncid, 'qhl     ', strt, cnt, qhl )
        elseif (ptype.eq.5) then
          write(*,*) '...graupel'
          call nc_getvar4( ncid, 'qg      ',  strt,  cnt,  qg  )
          write(*,*) '...snow'
          call nc_getvar4( ncid, 'qs      ',  strt,  cnt,  qs  )
        elseif (ptype.eq.52) then
          write(*,*) '...2nd free ice category'
          call nc_getvar4( ncid, 'qi2     ',  strt,  cnt,  qi2 )
        endif
       
        qt = qc + qr + qi + qi2 + qg + qs + qhl

      endif

      do i=1,ni
      do j=1,nj
      do k=1,nk
        thpert(i,j,k) = th(i,j,k)-th0(k)
        pipert(i,j,k) = ppi(i,j,k)-pi0(k)
        prspert(i,j,k) = prs(i,j,k)-prs0(k)
      enddo
      enddo
      enddo

      write(*,*) 'thermodynamic data complete'
      write(*,*) '--------------------------------------------------'

    end subroutine

    !===========================================================================!

    !===========================================================================!
    !###########################################################################!
    !#                   SUBROUTINE DYN_NC_GETWINDS                            #!
    !#                                                                         #!
    !#  Returns base state winds, total wind field, winds interpolated to      #!
    !#  scalar grid, and perturbation winds.                                   #!
    !#  Passed variable:                                                       #!
    !#    integer ncid : netcdf file identifier                                #!
    !#  Returned variables:                                                    #!
    !#    real,dimension(:) u0, v0 : base state winds                          #!
    !#    real,dimension(:,:,:) u, v, w : total wind field on momentum grids   #!
    !#    real,dimension(:,:,:) ui, vi, wi : wind interpolated to scalar grid  #!
    !#    real,dimension(:,:,:) upert, vpert: perturbation winds               #!
    !#=========================================================================#!
    !# Ryan Hastings, 6 Sep 2011                                               #!
    !#-------------------------------------------------------------------------#!
    !# v1.0, updated commentary, Ryan Hastings 13 Dec 2016                     #!
    !###########################################################################!
    subroutine dyn_nc_getwinds( ncid, u0, v0, u, v, w, ui, vi, wi, upert, vpert )

      !#########################################################################!
      !############################## VARIABLES ################################!

      !-------------------------- input variables ------------------------------!
      integer ncid

      !------------------------- output variables ------------------------------!
      real,dimension(1:nk)     :: u0, v0
      real,dimension(1:ni,1:nj,1:nk) :: ui, vi, wi
      real,dimension(ni+1,nj,nk) :: u, upert
      real,dimension(ni,nj+1,nk) :: v, vpert
      real,dimension(ni,nj,nk+1) :: w

      !-------------------------- local variables ------------------------------!

      integer,dimension(4)             :: strt4

      integer                             i, j, k, varid, rcode

      !#########################################################################!
      !############################## MAIN BODY ################################!

      write(*,*) 'reading in winds..'
      write(*,*)

      !-------------------------------------------------------------------------!
      ! set up start arrays

      strt4 = (/1,1,1,1/)

      !-------------------------------------------------------------------------!
      ! get base state winds on momentum grid
      write(*,*) 'base state zonal wind',nk
      call nc_getvar3( ncid, 'u0      ', strt4, (/1,1,nk,1/), u0 )
      write(*,*) 'base state meridional wind'
      call nc_getvar3( ncid, 'v0      ', strt4, (/1,1,nk,1/), v0 )

      !-------------------------------------------------------------------------!
      ! get winds on mass grid
      write(*,*) 'zonal wind on mass grid'
      call nc_getvar4( ncid, 'uinterp ', strt4, (/ni,nj,nk,1/), ui )
      write(*,*) 'meridional wind on mass grid'
      call nc_getvar4( ncid, 'vinterp ', strt4, (/ni,nj,nk,1/), vi )
      write(*,*) 'vertical wind on mass grid'
      call nc_getvar4( ncid, 'winterp ', strt4, (/ni,nj,nk,1/), wi )

      !-------------------------------------------------------------------------!
      ! get winds on momentum grid
      write(*,*) 'zonal wind on momentum grid'
      varid = ncvid( ncid, 'u       ', rcode )
      call ncvgt( ncid, varid, strt4, (/ni+1,nj,nk,1/), u, rcode )
      write(*,*) 'meridional wind on momentum grid'
      varid = ncvid( ncid, 'v       ', rcode )
      call ncvgt( ncid, varid, strt4, (/ni,nj+1,nk,1/), v, rcode )
      write(*,*) 'vertical wind on momentum grid'
      varid = ncvid( ncid, 'w       ', rcode )
      call ncvgt( ncid, varid, strt4, (/ni,nj,nk+1,1/), w, rcode )

      !-------------------------------------------------------------------------!
      ! get perturbation winds on momentum grid
      do k=1,nk
        do i=1,ni+1
        do j=1,nj
          upert(i,j,k) = u(i,j,k)-u0(k)
        enddo
        enddo
        do j=1,nj+1
        do i=1,ni
          vpert(i,j,k) = v(i,j,k)-v0(k)
        enddo
        enddo
      enddo

      write(*,*) 'winds read'
      write(*,*) '-------------------------------------------------'

    end subroutine

    !##########################################################################!
    !#                      SUBROUTINE NC_TRAJECT_VARDEF                      #!
    !#                                                                        #!
    !#  Sets up netCDF output file for traj.  Defines all variables.          #!
    !#  Passed variable:                                                      #!
    !#    character outfile : name of output file                             #!
    !#  Returned variables:                                                   #!
    !#    integers that are netCDF variable identifiers                       #!
    !#========================================================================#!
    !# v1.0, Ryan Hastings, 14 Dec 2016                                       #!
    !##########################################################################!
    subroutine nc_traject_vardef( outfile, ncid, ntimes, nparcels, &
      tpid, xpid, ypid, zpid, upid, vpid, wpid, qvpid, qcpid, qrpid, qipid, &
      qspid, qgpid, qhlpid, ncpid, nrpid, nipid, nspid, ngpid, &
      thpid, thrpid, buoypid, pgfxpid, pgfypid, pgfzpid, &
      pgfxbpid, pgfybpid, pgfzbpid, pgfxdlpid, pgfydlpid, pgfzdlpid, &
      pgfxdnpid, pgfydnpid, pgfzdnpid, tkepid, rhodid, pipid, &
      xvortpid, yvortpid, zvortpid, svortpid, cvortpid, psipid, &
      zvort_stretchid, zvort_tiltid, svort_tiltid, svort_stretchid, svort_bclnid,&
      svort_xchgid, dpsidtid, cvort_tiltid, cvort_stretchid, cvort_bclnid, &
      cvort_xchgid )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ----------------------------!
      character(len=100)             :: outfile
      integer                        :: nparcels, ntimes
  
      !---------------------------- output variables ---------------------------!
      integer                            :: ncid
      integer tpid, xpid, ypid, zpid, upid, vpid, wpid
      integer qvpid, qcpid, qrpid, qipid, qspid, qgpid, qhlpid, thpid, thrpid, buoypid
      integer ncpid, nrpid, nipid, nspid, ngpid, pipid
      integer pgfxpid, pgfypid, pgfzpid  
      integer pgfxbpid, pgfxdlpid, pgfxdnpid
      integer pgfybpid, pgfydlpid, pgfydnpid
      integer pgfzbpid, pgfzdlpid, pgfzdnpid
      integer tkepid, rhodid
      integer xvortpid, yvortpid, zvortpid, svortpid, cvortpid, psipid
      integer zvort_stretchid, zvort_tiltid
      integer svort_tiltid, svort_stretchid, svort_bclnid, svort_xchgid
      integer dpsidtid
      integer cvort_tiltid, cvort_stretchid, cvort_bclnid, cvort_xchgid

      !---------------------------- local variables -------------------------!
      integer rcode, i
      integer tid, npid
      integer,dimension(2) :: dimid2


     !#########################################################################!
     !############################## MAIN BODY ################################!
     !----------------------------------------------------------!
     ! create netcdf file
     write(*,*) 'creating ',outfile
     ncid = nccre( outfile, NCNOCLOB, rcode )

     !----------------------------------------------------------!
     ! define dimensions
       write(*,*) 'ntimes=',ntimes
       tid  = ncddef( ncid, 'time',    ntimes,   rcode )
       write(*,*) 'nparcels=',nparcels
       npid = ncddef( ncid, 'parcel', ncunlim, rcode )

    dimid2 = (/tid,npid/)

    !----------------------------------------------------------!
    ! define variables
    tpid = ncvdef( ncid, 'tp', nf_real, 1,    tid, rcode )

    ! spacial variables
    xpid = ncvdef( ncid, 'xp', nf_real, 2, dimid2, rcode )
    ypid = ncvdef( ncid, 'yp', nf_real, 2, dimid2, rcode )
    zpid = ncvdef( ncid, 'zp', nf_real, 2, dimid2, rcode )

    ! wind
    upid = ncvdef( ncid, 'up', nf_real, 2, dimid2, rcode )
    vpid = ncvdef( ncid, 'vp', nf_real, 2, dimid2, rcode )
    wpid = ncvdef( ncid, 'wp', nf_real, 2, dimid2, rcode )
 
    ! mixing ratios and number concentrations
    qvpid = ncvdef( ncid, 'qv', nf_real, 2, dimid2, rcode )
    qcpid = ncvdef( ncid, 'qc', nf_real, 2, dimid2, rcode )
    ncpid = ncvdef( ncid, 'ncc', nf_real, 2, dimid2, rcode )
    qrpid = ncvdef( ncid, 'qr', nf_real, 2, dimid2, rcode )
    nrpid = ncvdef( ncid, 'ncr', nf_real, 2, dimid2, rcode )
    qipid = ncvdef( ncid, 'qi', nf_real, 2, dimid2, rcode )
    nipid = ncvdef( ncid, 'nci', nf_real, 2, dimid2, rcode )
    qspid = ncvdef( ncid, 'qs', nf_real, 2, dimid2, rcode )
    nspid = ncvdef( ncid, 'ncs', nf_real, 2, dimid2, rcode )
    qgpid = ncvdef( ncid, 'qg', nf_real, 2, dimid2, rcode )
    ngpid = ncvdef( ncid, 'ncg', nf_real, 2, dimid2, rcode )
    if(ptype.eq.27)then
      qhlpid = ncvdef( ncid, 'qhl', nf_real, 2, dimid2, rcode )
    endif


    ! tempertatures
    thpid = ncvdef( ncid, 'th', nf_real, 2, dimid2, rcode )
    thrpid = ncvdef( ncid, 'thrpert', nf_real, 2, dimid2, rcode )
    rhodid = ncvdef( ncid, 'rhod', nf_real, 2, dimid2, rcode )
    pipid = ncvdef( ncid, 'pi', nf_real, 2, dimid2, rcode )

    if(output_dyn.eq.1)then
    ! forces
      buoypid = ncvdef( ncid, 'buoy', nf_real, 2, dimid2, rcode )

      pgfxpid = ncvdef( ncid, 'pgfx', nf_real, 2, dimid2, rcode )
      pgfypid = ncvdef( ncid, 'pgfy', nf_real, 2, dimid2, rcode )
      pgfzpid = ncvdef( ncid, 'pgfz', nf_real, 2, dimid2, rcode )

      pgfxbpid = ncvdef( ncid, 'pgfx_b', nf_real, 2, dimid2, rcode )
      pgfybpid = ncvdef( ncid, 'pgfy_b', nf_real, 2, dimid2, rcode )
      pgfzbpid = ncvdef( ncid, 'pgfz_b', nf_real, 2, dimid2, rcode )

      pgfxdlpid = ncvdef( ncid, 'pgfx_dl', nf_real, 2, dimid2, rcode )
      pgfydlpid = ncvdef( ncid, 'pgfy_dl', nf_real, 2, dimid2, rcode )
      pgfzdlpid = ncvdef( ncid, 'pgfz_dl', nf_real, 2, dimid2, rcode )

      pgfxdnpid = ncvdef( ncid, 'pgfx_dn', nf_real, 2, dimid2, rcode )
      pgfydnpid = ncvdef( ncid, 'pgfy_dn', nf_real, 2, dimid2, rcode )
      pgfzdnpid = ncvdef( ncid, 'pgfz_dn', nf_real, 2, dimid2, rcode )

    else

      buoypid=0
      pgfxpid=0
      pgfypid=0
      pgfzpid=0
      pgfxbpid=0
      pgfybpid=0
      pgfzbpid=0
      pgfxdlpid=0
      pgfydlpid=0
      pgfzdlpid=0
      pgfxdnpid=0
      pgfydnpid=0
      pgfzdnpid=0

    endif 

    ! extra stuff
    tkepid = ncvdef( ncid, 'tke', nf_real, 2, dimid2, rcode )

    if(output_vor.eq.1)then
    ! vorticity
      xvortpid = ncvdef( ncid, 'xvort', nf_real, 2, dimid2, rcode )
      yvortpid = ncvdef( ncid, 'yvort', nf_real, 2, dimid2, rcode )
      zvortpid = ncvdef( ncid, 'zvort', nf_real, 2, dimid2, rcode )
      svortpid = ncvdef( ncid, 'svort', nf_real, 2, dimid2, rcode )
      cvortpid = ncvdef( ncid, 'cvort', nf_real, 2, dimid2, rcode )

      psipid = ncvdef( ncid, 'psi', nf_real, 2, dimid2, rcode )

      zvort_stretchid = ncvdef( ncid, 'zvort_stretch', nf_real, 2, dimid2, &
        rcode )
      zvort_tiltid = ncvdef( ncid, 'zvort_tilt', nf_real, 2, dimid2, rcode )

      svort_stretchid = ncvdef( ncid, 'svort_stretch', nf_real, 2, dimid2, &
        rcode )
      svort_tiltid = ncvdef( ncid, 'svort_tilt', nf_real, 2, dimid2, rcode )
      svort_xchgid = ncvdef( ncid, 'svort_xchg', nf_real, 2, dimid2, rcode )
      svort_bclnid = ncvdef( ncid, 'svort_bcln', nf_real, 2, dimid2, rcode )

      dpsidtid = ncvdef( ncid, 'dpsidt', nf_real, 2, dimid2, rcode )

      cvort_stretchid = ncvdef( ncid, 'cvort_stretch', nf_real, 2, dimid2, &
        rcode )
      cvort_tiltid = ncvdef( ncid, 'cvort_tilt', nf_real, 2, dimid2, rcode )
      cvort_xchgid = ncvdef( ncid, 'cvort_xchg', nf_real, 2, dimid2, rcode )
      cvort_bclnid = ncvdef( ncid, 'cvort_bcln', nf_real, 2, dimid2, rcode )

    else

      xvortpid=0
      yvortpid=0
      zvortpid=0
      svortpid=0
      cvortpid=0

      psipid=0

      zvort_stretchid = 0
      zvort_tiltid = 0

    endif

    ! end definitions
    call ncendf( ncid, rcode )

  end subroutine
  !############################################################################!
  !#                     SUBROUTINE NC_CIRCUIT_VARDEF                         #!
  !#                                                                          #!
  !#  Creates netCDF output file for circuit, and defines variables.          #!
  !#  Passed variables:                                                       #!
  !#    character outfile : name of output file                               #!
  !#    integer nparcels : number of parcels                                  #!
  !#    integer ntimes : number of times                                      #!
  !#  Returned variables:                                                     #!
  !#    integers that are variable identifiers                                #!
  !#==========================================================================#!
  !# v1.0, Ryan Hastings, 14 Dec 2016                                         #!
  !############################################################################!
  subroutine nc_circuit_vardef( outfile, ncid, ntimes, nparcels, &
      tpid, xpid, ypid, zpid, upid, vpid, wpid, qvpid, qcpid, qrpid, qipid, &
      qspid, qgpid, qhlpid, ncpid, nrpid, nipid, nspid, ngpid, &
      thpid, thrpid, buoypid, pgfxpid, pgfypid, pgfzpid, &
      pgfxbpid, pgfybpid, pgfzbpid, pgfxdlpid, pgfydlpid, pgfzdlpid, &
      pgfxdnpid, pgfydnpid, pgfzdnpid, tkepid, rhodid, pipid, &
      xvortpid, yvortpid, zvortpid, svortpid, cvortpid, psipid, &
      zvort_stretchid, zvort_tiltid, svort_tiltid, svort_stretchid, svort_bclnid,&
      svort_xchgid, dpsidtid, cvort_tiltid, cvort_stretchid, cvort_bclnid, &
      cvort_xchgid, circid, vdotl1id, vdotl2id, intbdzid )
      implicit none
      !#########################################################################!
      !############################## VARIABLES ################################!
      !-------------------------- input variables ----------------------------!
      character(len=100)             :: outfile
      integer                        :: nparcels, ntimes

      !---------------------------- output variables ---------------------------!
      integer                            :: ncid
      integer tpid, xpid, ypid, zpid, upid, vpid, wpid
      integer qvpid, qcpid, qrpid, qipid, qspid, qgpid, qhlpid, thpid, thrpid, buoypid
      integer ncpid, nrpid, nipid, nspid, ngpid, pipid
      integer pgfxpid, pgfypid, pgfzpid
      integer pgfxbpid, pgfxdlpid, pgfxdnpid
      integer pgfybpid, pgfydlpid, pgfydnpid
      integer pgfzbpid, pgfzdlpid, pgfzdnpid
      integer tkepid, diabid, rhodid
      integer xvortpid, yvortpid, zvortpid, svortpid, cvortpid, psipid
      integer zvort_stretchid, zvort_tiltid
      integer svort_tiltid, svort_stretchid, svort_bclnid, svort_xchgid
      integer dpsidtid
      integer cvort_tiltid, cvort_stretchid, cvort_bclnid, cvort_xchgid
      integer circid, vdotl1id, vdotl2id, intbdzid

      !---------------------------- local variables -------------------------!
      integer rcode, i
      integer tid, npid
      integer,dimension(2) :: dimid2


     !#########################################################################!
     !############################## MAIN BODY ################################!
     !----------------------------------------------------------!
     ! create netcdf file
     write(*,*) 'creating ',outfile
     ncid = nccre( outfile, NCNOCLOB, rcode )

     !----------------------------------------------------------!
     ! define dimensions
       write(*,*) 'ntimes=',ntimes
       tid  = ncddef( ncid, 'time',    ntimes,   rcode )
       write(*,*) 'nparcels=',nparcels
       npid = ncddef( ncid, 'parcel', ncunlim, rcode )

    dimid2 = (/tid,npid/)

    !----------------------------------------------------------!
    ! define variables
    tpid = ncvdef( ncid, 'tp', nf_real, 1,    tid, rcode )

    ! spacial variables
    xpid = ncvdef( ncid, 'xp', nf_real, 2, dimid2, rcode )
    ypid = ncvdef( ncid, 'yp', nf_real, 2, dimid2, rcode )
    zpid = ncvdef( ncid, 'zp', nf_real, 2, dimid2, rcode )

    ! wind
    upid = ncvdef( ncid, 'up', nf_real, 2, dimid2, rcode )
    vpid = ncvdef( ncid, 'vp', nf_real, 2, dimid2, rcode )
    wpid = ncvdef( ncid, 'wp', nf_real, 2, dimid2, rcode )

    ! mixing ratios and number concentrations
    qvpid = ncvdef( ncid, 'qv', nf_real, 2, dimid2, rcode )
    qcpid = ncvdef( ncid, 'qc', nf_real, 2, dimid2, rcode )
    ncpid = ncvdef( ncid, 'ncc', nf_real, 2, dimid2, rcode )
    qrpid = ncvdef( ncid, 'qr', nf_real, 2, dimid2, rcode )
    nrpid = ncvdef( ncid, 'ncr', nf_real, 2, dimid2, rcode )
    qipid = ncvdef( ncid, 'qi', nf_real, 2, dimid2, rcode )
    nipid = ncvdef( ncid, 'nci', nf_real, 2, dimid2, rcode )
    qspid = ncvdef( ncid, 'qs', nf_real, 2, dimid2, rcode )
    nspid = ncvdef( ncid, 'ncs', nf_real, 2, dimid2, rcode )
    qgpid = ncvdef( ncid, 'qg', nf_real, 2, dimid2, rcode )

    if(ptype.eq.27)then
      qhlpid = ncvdef( ncid, 'qhl', nf_real, 2, dimid2, rcode )
    endif
    ngpid = ncvdef( ncid, 'ncg', nf_real, 2, dimid2, rcode )

    ! tempertatures
    thpid = ncvdef( ncid, 'th', nf_real, 2, dimid2, rcode )
    thrpid = ncvdef( ncid, 'thrpert', nf_real, 2, dimid2, rcode )
    rhodid = ncvdef( ncid, 'rhod', nf_real, 2, dimid2, rcode )
    pipid = ncvdef( ncid, 'pi', nf_real, 2, dimid2, rcode )

    circid = ncvdef( ncid, 'circ', nf_real, 1, tid, rcode )
    vdotl1id = ncvdef( ncid, 'vdotl1', nf_real, 2, dimid2, rcode )
    vdotl2id = ncvdef( ncid, 'vdotl2', nf_real, 2, dimid2, rcode )
    intBdzid = ncvdef( ncid, 'intBdz', nf_real, 1, tid, rcode )

    if(output_dyn.eq.1)then
    ! forces
      buoypid = ncvdef( ncid, 'buoy', nf_real, 2, dimid2, rcode )

      pgfxpid = ncvdef( ncid, 'pgfx', nf_real, 2, dimid2, rcode )
      pgfypid = ncvdef( ncid, 'pgfy', nf_real, 2, dimid2, rcode )
      pgfzpid = ncvdef( ncid, 'pgfz', nf_real, 2, dimid2, rcode )

      pgfxbpid = ncvdef( ncid, 'pgfx_b', nf_real, 2, dimid2, rcode )
      pgfybpid = ncvdef( ncid, 'pgfy_b', nf_real, 2, dimid2, rcode )
      pgfzbpid = ncvdef( ncid, 'pgfz_b', nf_real, 2, dimid2, rcode )

      pgfxdlpid = ncvdef( ncid, 'pgfx_dl', nf_real, 2, dimid2, rcode )
      pgfydlpid = ncvdef( ncid, 'pgfy_dl', nf_real, 2, dimid2, rcode )
      pgfzdlpid = ncvdef( ncid, 'pgfz_dl', nf_real, 2, dimid2, rcode )

      pgfxdnpid = ncvdef( ncid, 'pgfx_dn', nf_real, 2, dimid2, rcode )
      pgfydnpid = ncvdef( ncid, 'pgfy_dn', nf_real, 2, dimid2, rcode )
      pgfzdnpid = ncvdef( ncid, 'pgfz_dn', nf_real, 2, dimid2, rcode )

    else

      buoypid=0
      pgfxpid=0
      pgfypid=0
      pgfzpid=0
      pgfxbpid=0
      pgfybpid=0
      pgfzbpid=0
      pgfxdlpid=0
      pgfydlpid=0
      pgfzdlpid=0
      pgfxdnpid=0
      pgfydnpid=0
      pgfzdnpid=0

    endif

    ! extra stuff
    tkepid = ncvdef( ncid, 'tke', nf_real, 2, dimid2, rcode )

    if(output_vor.eq.1)then
    ! vorticity
      xvortpid = ncvdef( ncid, 'xvort', nf_real, 2, dimid2, rcode )
      yvortpid = ncvdef( ncid, 'yvort', nf_real, 2, dimid2, rcode )
      zvortpid = ncvdef( ncid, 'zvort', nf_real, 2, dimid2, rcode )
      svortpid = ncvdef( ncid, 'svort', nf_real, 2, dimid2, rcode )
      cvortpid = ncvdef( ncid, 'cvort', nf_real, 2, dimid2, rcode )

      psipid = ncvdef( ncid, 'psi', nf_real, 2, dimid2, rcode )

      zvort_stretchid = ncvdef( ncid, 'zvort_stretch', nf_real, 2, dimid2, &
        rcode )
      zvort_tiltid = ncvdef( ncid, 'zvort_tilt', nf_real, 2, dimid2, rcode )

      svort_stretchid = ncvdef( ncid, 'svort_stretch', nf_real, 2, dimid2, &
        rcode )
      svort_tiltid = ncvdef( ncid, 'svort_tilt', nf_real, 2, dimid2, rcode )
      svort_xchgid = ncvdef( ncid, 'svort_xchg', nf_real, 2, dimid2, rcode )
      svort_bclnid = ncvdef( ncid, 'svort_bcln', nf_real, 2, dimid2, rcode )

      dpsidtid = ncvdef( ncid, 'dpsidt', nf_real, 2, dimid2, rcode )

      cvort_stretchid = ncvdef( ncid, 'cvort_stretch', nf_real, 2, dimid2, &
        rcode )
      cvort_tiltid = ncvdef( ncid, 'cvort_tilt', nf_real, 2, dimid2, rcode )
      cvort_xchgid = ncvdef( ncid, 'cvort_xchg', nf_real, 2, dimid2, rcode )
      cvort_bclnid = ncvdef( ncid, 'cvort_bcln', nf_real, 2, dimid2, rcode )

    else

      xvortpid=0
      yvortpid=0
      zvortpid=0
      svortpid=0
      cvortpid=0

      psipid=0

      zvort_stretchid = 0
      zvort_tiltid = 0

    endif

    ! end definitions
    call ncendf( ncid, rcode )

  end subroutine

  !############################################################################!
  !#                  SUBROUTINE TRAJECT_NC_READ_CM1                          #!
  !#                                                                          #!
  !#  Reads variables from CM1 file for traj.                                 #!
  !#  Passed variables:                                                       #!
  !#    character filename:  name of CM1 file                                 #!
  !#    real,dimension(:) pi0, thr0 : base state variables                    #!
  !#  Returned variables:                                                     #!
  !#    real,dimension(:,:,:) u, v, w:  total winds                           #!
  !#    real,dimension(:,:,:) qv, qc, qr, qi, qs, qg, qhl : mixing ratios     #!
  !#    real,dimension(:,:,:) th : potential temperature                      #!
  !#    real,dimension(:,:,:) thr : density potential temperature             #!
  !#    real,dimension(:,:,:) pgfx, pgfy, pgfz:  PGFs                         #!
  !#    real,dimension(:,:,:) buoy : buoyancy                                 #!
  !#    real,dimension(:,:,:) ncc, ncr, nci, ncs, ncg : number concentrations #!
  !#    real,dimension(:,:,:) rhod : dry air density                          #!
  !#    real,dimension(:,:,:) ppi : nondimensionalized pressure               #!
  !############################################################################!
  subroutine traject_nc_read_cm1( infiles, ti, u1, u2, v1, v2, w1, w2, qv1, qv2, &
    qc1, qc2, qr1, qr2, qi1, qi2, qs1, qs2, qg1, qg2, qhl1, qhl2, &
    th1, th2, thr0, pi0, pgfx1, pgfx2, pgfy1, pgfy2, pgfz1, pgfz2, thr1, thr2,&
    buoy1, buoy2, tke1, tke2,&
    ncc1, ncc2, ncr1, ncr2, nci1, nci2, ncs1, ncs2, ncg1, ncg2, rhod1, rhod2, &
    ppi1, ppi2,dt )

    implicit none

    character(len=100),intent(in),dimension(:) :: infiles
    integer,intent(in) :: ti
    real,intent(in) :: dt

    real,intent(inout),dimension(1:ni+1,1:nj,0:nk) :: u1, u2, pgfx1, pgfx2
    real,intent(inout),dimension(1:ni,1:nj+1,0:nk) :: v1, v2, pgfy1, pgfy2
    real,intent(inout),dimension(1:ni,1:nj,1:nk+1) :: w1, w2

    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: qv1, qv2, qc1, qc2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: qr1, qr2, qi1, qi2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: qs1, qs2, qg1, qg2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: qhl1, qhl2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: ncc1, ncc2, ncr1, ncr2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: nci1, nci2, ncs1, ncs2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: ncg1, ncg2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: th1, th2, thr1, thr2
    real,intent(inout),dimension(1:ni,1:nj,0:nk) :: buoy1, buoy2
    real,intent(inout),dimension(1:ni,1:nj,1:nk+1) :: pgfz1, pgfz2, tke1, tke2

    real,intent(in),dimension(1:nk) :: pi0, thr0

    integer ncid, rcode, varid 
    integer,dimension(4) :: st
    integer i, j, k

    real,dimension(ni,nj,0:nk) :: ppi1, ppi2, pipert1, pipert2, rhod1, rhod2

!    if(dt.gt.0)then
      ncid = ncopn( infiles(ti), NCNOWRIT, rcode )
!    else
!      ncid = ncopn( infile(ti-1),NCNOWRIT, rcode )
!    endif

    st=(/1,1,1,1/)

    varid = ncvid( ncid, 'u', rcode )
    call ncvgt( ncid, varid, st, (/nip1,nj,nk,1/), u1(:,:,1:nk), rcode )
    varid = ncvid( ncid, 'v', rcode )
    call ncvgt( ncid, varid, st, (/ni,njp1,nk,1/), v1(:,:,1:nk), rcode )
    varid = ncvid( ncid, 'w', rcode )
    call ncvgt( ncid, varid, st, (/ni,nj,nkp1,1/), w1(:,:,:), rcode )

    u1(:,:,0)=u1(:,:,1)
    v1(:,:,0)=v1(:,:,1)

    call nc_readvar4( ncid, 'qv      ', qv1 )
    call nc_readvar4( ncid, 'qc      ', qc1 )
   ! call nc_readvar4( ncid, 'ncc     ', ncc )
    ncc1(:,:,:)=0.0 ! oops.  cloud droplets have a fixed number concentration..
                   ! i didn't want to have to go through and eliminate this everywhere,
                   ! so as a kluge, just setting it to zero

    if(noncq.eq.1)then
      ncr1(:,:,:)=0.0
      nci1(:,:,:)=0.0
      ncs1(:,:,:)=0.0
      ncg1(:,:,:)=0.0
    else
      call nc_readvar4( ncid, 'ncr     ', ncr1 )
      call nc_readvar4( ncid, 'nci     ', nci1 )
      call nc_readvar4( ncid ,'ncs     ', ncs1 )
      call nc_readvar4( ncid, 'ncg     ', ncg1 )
      call nc_readvar4( ncid, 'rho     ', rhod1 )
    endif

    call nc_readvar4( ncid, 'qr      ', qr1 )
    call nc_readvar4( ncid, 'qi      ', qi1 )
    call nc_readvar4( ncid, 'qs      ', qs1 )
    call nc_readvar4( ncid, 'qg      ', qg1 )
    if(ptype.eq.27)then
      call nc_readvar4( ncid, 'qhl     ', qhl1 )
    endif
    call nc_readvar4( ncid, 'th      ', th1 )
    call nc_readvar4( ncid, 'pi      ', ppi1 )

    do i=1,ni
    do j=1,nj
    do k=1,nk
      pipert1(i,j,k) = ppi1(i,j,k) - pi0(k)
    enddo
    enddo
    enddo

    if(ptype.eq.27)then
      thr1 = th1*(1+reps*qv1)/(1+qv1+qc1+qr1+qi1+qs1+qg1+qhl1)
    else
      thr1 = th1*(1+reps*qv1)/(1+qv1+qc1+qr1+qi1+qs1+qg1)
    endif

    do i=1,ni
    do j=1,nj
      do k=1,nk
        buoy1(i,j,k) = g*( thr1(i,j,k)/thr0(k) - 1.0 )
      enddo
      buoy1(i,j,0) = buoy1(i,j,1)
    enddo
    enddo

    pgfx1(:,:,:) = 0.0
    pgfy1(:,:,:) = 0.0
    pgfz1(:,:,:) = 0.0

    do i=1,ni
    do j=1,nj
      do k=2,nk
        pgfz1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i,j,k-1))*(pipert1(i,j,k)-pipert1(i,j,k-1))*mzf(k)
      enddo
      pgfz1(i,j,1)    = -buoy1(i,j,1)
      pgfz1(i,j,nk+1) = -buoy1(i,j,nk)
    enddo
    enddo

    if(output_z_only.eq.0)then

      do i=2,ni
      do j=1,nj
        do k=1,nk
          pgfx1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i-1,j,k))*&
            (pipert1(i,j,k)-pipert1(i-1,j,k))
        enddo
        pgfx1(i,j,0) = pgfy1(1,j,1)
      enddo
      enddo

      do i=1,ni
      do j=2,nj
        do k=1,nk
          pgfy1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i,j-1,k))*&
            (pipert1(i,j,k)-pipert1(i,j-1,k))
        enddo
        pgfy1(i,j,0)=pgfy1(i,j,1)
      enddo
      enddo

    endif
      
    call nc_readvar4( ncid, 'tke     ', tke1 )

    call ncclos(ncid,rcode)

!    if(dt.gt.0)then
      ncid = ncopn( infiles(ti+1), NCNOWRIT, rcode )
!    else
!      ncid = ncopn( infile(ti), NCNOWRIT, rcode )
!    endif


    st=(/1,1,1,1/)

    varid = ncvid( ncid, 'u', rcode )
    call ncvgt( ncid, varid, st, (/nip1,nj,nk,1/), u2(:,:,1:nk), rcode )
    varid = ncvid( ncid, 'v', rcode )
    call ncvgt( ncid, varid, st, (/ni,njp1,nk,1/), v2(:,:,1:nk), rcode )
    varid = ncvid( ncid, 'w', rcode )
    call ncvgt( ncid, varid, st, (/ni,nj,nkp1,1/), w2(:,:,:), rcode )

    u2(:,:,0)=u2(:,:,1)
    v2(:,:,0)=v2(:,:,1)

    call nc_readvar4( ncid, 'qv      ', qv2 )
    call nc_readvar4( ncid, 'qc      ', qc2 )
   ! call nc_readvar4( ncid, 'ncc     ', ncc )
    ncc2(:,:,:)=0.0 ! oops.  cloud droplets have a fixed number concentration..
                   ! i didn't want to have to go through and eliminate this everywhere,
                   ! so as a kluge, just setting it to zero

    if(noncq.eq.1)then
      ncr2(:,:,:)=0.0
      nci2(:,:,:)=0.0
      ncs2(:,:,:)=0.0
      ncg2(:,:,:)=0.0
    else
      call nc_readvar4( ncid, 'ncr     ', ncr2 )
      call nc_readvar4( ncid, 'nci     ', nci2 )
      call nc_readvar4( ncid ,'ncs     ', ncs2 )
      call nc_readvar4( ncid, 'ncg     ', ncg2 )
      call nc_readvar4( ncid, 'rho     ', rhod2 )
    endif

    call nc_readvar4( ncid, 'qr      ', qr2 )
    call nc_readvar4( ncid, 'qi      ', qi2 )
    call nc_readvar4( ncid, 'qs      ', qs2 )
    call nc_readvar4( ncid, 'qg      ', qg2 )
    if(ptype.eq.27)then
      call nc_readvar4( ncid, 'qhl     ', qhl2 )
    endif
    call nc_readvar4( ncid, 'th      ', th2 )
    call nc_readvar4( ncid, 'pi      ', ppi2 )

    do i=1,ni
    do j=1,nj
    do k=1,nk
      pipert2(i,j,k) = ppi2(i,j,k) - pi0(k)
    enddo
    enddo
    enddo

    if(ptype.eq.27)then
      thr2 = th2*(1+reps*qv2)/(1+qv2+qc2+qr2+qi2+qs2+qg2+qhl2)
    else
      thr2 = th2*(1+reps*qv2)/(1+qv2+qc2+qr2+qi2+qs2+qg2)
    endif

    do i=1,ni
    do j=1,nj
      do k=1,nk
        buoy2(i,j,k) = g*( thr2(i,j,k)/thr0(k) - 1.0 )
      enddo
      buoy2(i,j,0) = buoy2(i,j,1)
    enddo
    enddo

    pgfx2(:,:,:) = 0.0
    pgfy2(:,:,:) = 0.0
    pgfz2(:,:,:) = 0.0

    do i=1,ni
    do j=1,nj
      do k=2,nk
        pgfz2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i,j,k-1))*(pipert2(i,j,k)-pipert2(i,j,k-1))*mzf(k)
      enddo
      pgfz2(i,j,1)    = -buoy2(i,j,1)
      pgfz2(i,j,nk+1) = -buoy2(i,j,nk)
    enddo
    enddo

    if(output_z_only.eq.0)then

      do i=2,ni
      do j=1,nj
        do k=1,nk
          pgfx2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i-1,j,k))*&
            (pipert2(i,j,k)-pipert2(i-1,j,k))
        enddo
        pgfx2(i,j,0) = pgfy2(1,j,1)
      enddo
      enddo

      do i=1,ni
      do j=2,nj
        do k=1,nk
          pgfy2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i,j-1,k))*&
            (pipert2(i,j,k)-pipert2(i,j-1,k))
        enddo
        pgfy2(i,j,0)=pgfy2(i,j,1)
      enddo
      enddo

    endif


    call ncclos(ncid,rcode)

  end subroutine

  !############################################################################!
  !#               SUBROUTINE TRAJECT_NC_READ_DYN                             #!
  !#                                                                          #!
  !#  Reads in dyn output files for traj.                                     #!
  !#  Passed variables:                                                       #!
  !#    character dynfile : name of dyn file to read                          #!
  !#    real,dimension(:,:,:) thr : density potential temperature             #!
  !#    real,dimension(:,:,:) buoy : buoyancy                                 #!
  !#  Returned variables:                                                     #!
  !#    real,dimension(:,:,:) pgfx_b, pgfy_b, pgfz_b : buoyancy PPGF          #!
  !#    real,dimension(:,:,:) pgfx_dl, pgfy_dl, pgfz_dl : linear dynamic PPGF #!
  !#    real,dimension(:,:,:) pgfx_dn, pgfy_dn, pgfz_dn : nonlinear dynamic   #!
  !#      PPGF                                                                #!
  !#==========================================================================#!
  !# v1.0, Ryan Hastings, 14 Dec 2016                                         #!
  !############################################################################!
  subroutine traject_nc_read_dyn( dynfiles, ti, pgfx_b1, pgfx_b2, pgfy_b1, pgfy_b2,&
    pgfz_b1, pgfz_b2, &
    pgfx_dl1, pgfx_dl2, pgfy_dl1, pgfy_dl2, pgfz_dl1, pgfz_dl2, pgfx_dn1, pgfx_dn2,&
    pgfy_dn1, pgfy_dn2, pgfz_dn1, pgfz_dn2, thr1, thr2, buoy1, buoy2, dt )

    implicit none
    !######################### VARIABLES ######################################!
    ! Passed variables
    character(len=100),intent(in),dimension(:) :: dynfiles
    integer ti
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: thr1,thr2,buoy1,buoy2
    real,intent(in) :: dt

    ! Returned variables
    real,intent(out),dimension(1:ni+1,1:nj,0:nk) :: pgfx_b1, pgfx_b2
    real,intent(out),dimension(1:ni+1,1:nj,0:nk) :: pgfx_dl1, pgfx_dl2
    real,intent(out),dimension(1:ni+1,1:nj,0:nk) :: pgfx_dn1, pgfx_dn2
    real,intent(out),dimension(1:ni,1:nj+1,0:nk) :: pgfy_b1, pgfy_b2
    real,intent(out),dimension(1:ni,1:nj+1,0:nk) :: pgfy_dl1, pgfy_dl2
    real,intent(out),dimension(1:ni,1:nj+1,0:nk) :: pgfy_dn1, pgfy_dn2
    real,intent(out),dimension(1:ni,1:nj,1:nk+1) :: pgfz_b1, pgfz_b2
    real,intent(out),dimension(1:ni,1:nj,1:nk+1) :: pgfz_dl1, pgfz_dl2
    real,intent(out),dimension(1:ni,1:nj,1:nk+1) :: pgfz_dn1, pgfz_dn2

    ! Local variables
    real,dimension(1:ni,1:nj,0:nk) :: pi_b1, pi_b2, pi_dl1, pi_dl2
    real,dimension(1:ni,1:nj,0:nk) :: pi_dn1, pi_dn2

    integer ncid, rcode
    integer i, j, k

    !################################ MAIN BODY ###############################!

    ! open netcdf file
!    if(dt.gt.0)then
      ncid = ncopn( dynfiles(ti), NCNOWRIT, rcode )
!    else
!      ncid = ncopn( dynfiles(ti-1), NCNOWRIT, rcode )
!    endif


    ! read in pressure fields
    call nc_readvar4( ncid, 'pi_b     ', pi_b1 )
    call nc_readvar4( ncid, 'pi_dl    ', pi_dl1 )
    call nc_readvar4( ncid, 'pi_dn    ', pi_dn1 )

    ! compute PGFs
    pgfz_b1(:,:,:) = 0.0
    pgfz_dl1(:,:,:) = 0.0
    pgfz_dn1(:,:,:) = 0.0

    do i=1,ni
    do j=1,nj
      do k=2,nk
        pgfz_b1(i,j,k) = -0.5*rdz*cp*(thr1(i,j,k)+thr1(i,j,k-1))* &
          (pi_b1(i,j,k)-pi_b1(i,j,k-1))*mzf(k)
        pgfz_dl1(i,j,k) = -0.5*rdz*cp*(thr1(i,j,k)+thr1(i,j,k-1))* &
          (pi_dl1(i,j,k)-pi_dl1(i,j,k-1))*mzf(k)
        pgfz_dn1(i,j,k) = -0.5*rdz*cp*(thr1(i,j,k)+thr1(i,j,k-1))* &
          (pi_dn1(i,j,k)-pi_dn1(i,j,k-1))*mzf(k)
      enddo
      pgfz_b1(i,j,1) = -buoy1(i,j,1)
      pgfz_b1(i,j,nk+1)=-buoy1(i,j,nk)
    enddo
    enddo

    if(output_z_only.eq.0)then
      
      do i=2,ni
      do j=1,nj
        do k=1,nk
          pgfx_b1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i-1,j,k))*&
            (pi_b1(i,j,k)-pi_b1(i-1,j,k))
          pgfx_dl1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i-1,j,k))*&
            (pi_dl1(i,j,k)-pi_dl1(i-1,j,k))
          pgfx_dn1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i-1,j,k))*&
            (pi_dn1(i,j,k)-pi_dn1(i-1,j,k))

        enddo
        pgfx_b1(i,j,0) = pgfx_b1(i,j,1)
        pgfx_dl1(i,j,0) = pgfx_dl1(i,j,1)
        pgfx_dn1(i,j,0) = pgfx_dn1(i,j,1)
      enddo
      enddo
    
      do i=1,ni
      do j=2,nj
        do k=1,nk
          pgfy_b1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i,j-1,k))*&
            (pi_b1(i,j,k)-pi_b1(i,j-1,k))
          pgfy_dl1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i,j-1,k))*&
            (pi_dl1(i,j,k)-pi_dl1(i,j-1,k))
          pgfy_dn1(i,j,k) = -0.5*cp*rdz*(thr1(i,j,k)+thr1(i,j-1,k))*&
            (pi_dn1(i,j,k)-pi_dn1(i,j-1,k))
        enddo
        pgfy_b1(i,j,0)=pgfy_b1(i,j,1)
        pgfy_dl1(i,j,0)=pgfy_dl1(i,j,1)
        pgfy_dn1(i,j,0)=pgfy_dn1(i,j,1)
      enddo
      enddo
   
    endif

    ! close netcdf file
    call ncclos(ncid,rcode)

    ! open netcdf file
!    if(dt.gt.0)then
      ncid = ncopn( dynfiles(ti+1), NCNOWRIT, rcode )
!    else
!      ncid = ncopn( dynfiles(ti),NCNOWRIT,rcode)
!    endif

    ! read in pressure fields
    call nc_readvar4( ncid, 'pi_b     ', pi_b2 )
    call nc_readvar4( ncid, 'pi_dl    ', pi_dl2 )
    call nc_readvar4( ncid, 'pi_dn    ', pi_dn2 )

    ! compute PGFs
    pgfz_b2(:,:,:) = 0.0
    pgfz_dl2(:,:,:) = 0.0
    pgfz_dn2(:,:,:) = 0.0

    do i=1,ni
    do j=1,nj
      do k=2,nk
        pgfz_b2(i,j,k) = -0.5*rdz*cp*(thr2(i,j,k)+thr2(i,j,k-1))* &
          (pi_b2(i,j,k)-pi_b2(i,j,k-1))*mzf(k)
        pgfz_dl2(i,j,k) = -0.5*rdz*cp*(thr2(i,j,k)+thr2(i,j,k-1))* &
          (pi_dl2(i,j,k)-pi_dl2(i,j,k-1))*mzf(k)
        pgfz_dn2(i,j,k) = -0.5*rdz*cp*(thr2(i,j,k)+thr2(i,j,k-1))* &
          (pi_dn2(i,j,k)-pi_dn2(i,j,k-1))*mzf(k)
      enddo
      pgfz_b2(i,j,1) = -buoy2(i,j,1)
      pgfz_b2(i,j,nk+1)=-buoy2(i,j,nk)
    enddo
    enddo

    if(output_z_only.eq.0)then

      do i=2,ni
      do j=1,nj
        do k=1,nk
          pgfx_b2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i-1,j,k))*&
            (pi_b2(i,j,k)-pi_b2(i-1,j,k))
          pgfx_dl2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i-1,j,k))*&
            (pi_dl2(i,j,k)-pi_dl2(i-1,j,k))
          pgfx_dn2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i-1,j,k))*&
            (pi_dn2(i,j,k)-pi_dn2(i-1,j,k))

        enddo
        pgfx_b2(i,j,0) = pgfx_b2(i,j,1)
        pgfx_dl2(i,j,0) = pgfx_dl2(i,j,1)
        pgfx_dn2(i,j,0) = pgfx_dn2(i,j,1)
      enddo
      enddo


      do i=1,ni
      do j=2,nj
        do k=1,nk
          pgfy_b2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i,j-1,k))*&
            (pi_b2(i,j,k)-pi_b2(i,j-1,k))
          pgfy_dl2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i,j-1,k))*&
            (pi_dl2(i,j,k)-pi_dl2(i,j-1,k))
          pgfy_dn2(i,j,k) = -0.5*cp*rdz*(thr2(i,j,k)+thr2(i,j-1,k))*&
            (pi_dn2(i,j,k)-pi_dn2(i,j-1,k))
        enddo
        pgfy_b2(i,j,0)=pgfy_b2(i,j,1)
        pgfy_dl2(i,j,0)=pgfy_dl2(i,j,1)
        pgfy_dn2(i,j,0)=pgfy_dn2(i,j,1)
      enddo
      enddo

    endif

    ! close netcdf file
    call ncclos(ncid,rcode)



  end subroutine

  !############################################################################!
  !#                   SUBROUTINE NC_READVAR4                                 #!
  !#                                                                          #!
  !#  Reads 4-D netcdf variable on scalar grid.  This differs from            #!
  !#  CM1_NC_GETVAR4 in that it reads to an extended scalar vertical grid.    #!
  !#  I.e., the lowest grid point is below the ground.  Assumes values remain #!
  !#  constant below lowest level.                                            #!
  !############################################################################!
  subroutine nc_readvar4( ncid, varname, var )

    implicit none

    integer,intent(in) :: ncid
    character(len=8),intent(in) :: varname
    real,intent(out),dimension(ni,nj,0:nk) :: var

    integer rcode, varid
    integer i, j

    varid = ncvid( ncid, varname, rcode )
    call ncvgt( ncid, varid, (/1,1,1,1/), (/ni,nj,nk,1/), var(:,:,1:nk), rcode )

    do i=1,ni
    do j=1,nj
      var(i,j,0) = var(i,j,1)
    enddo
    enddo
 
  end subroutine

  !############################################################################!
  !#                     SUBROUTINE NC_TRAJECT_WRITEOUT                       #!
  !#                                                                          #!
  !#  Writes out variables for traj.                                          #!
  !############################################################################!
  subroutine nc_traject_writeout( ntimes, nparcels, ncid, tpid, tp, &
    xpid, xpc, ypid, ypc, zpid, zpc, upid, up, vpid, vp, wpid, wp, &
    buoypid, buoyp, pgfxpid, pgfxp, pgfypid, pgfyp, pgfzpid, pgfzp, &
    pgfxbpid, pgfxbp, pgfybpid, pgfybp, pgfzbpid, pgfzbp, &
    pgfxdlpid, pgfxdlp, pgfydlpid, pgfydlp, pgfzdlpid, pgfzdlp, &
    pgfxdnpid, pgfxdnp, pgfydnpid, pgfydnp, pgfzdnpid, pgfzdnp, &
    tkepid, tkep, &
      thrpid, thrp, qvpid, qvp, qcpid, qcp, qrpid, qrp, qipid, qip, &
      qspid, qsp, qgpid, qgp, qhlpid, qhlp, &
      ncpid, ncp, nrpid, nrp, nipid, nip, nspid, nsp, ngpid, ngp, rhodid, rhodp, &
      thpid, thp, pipid, pip, xvortpid, xvortp, yvortpid, yvortp, zvortpid, &
      zvortp, svortpid, svortp, cvortpid, cvortp, psipid, psip, &
      zvort_tiltid, zvort_tiltp, zvort_stretchid, zvort_stretchp, &
      svort_tiltid, svort_tiltp, svort_stretchid, svort_stretchp, &
      svort_bclnid, svort_bclnp, svort_xchgid, svort_xchgp, dpsidtid, dpsidt,&
      cvort_tiltid, cvort_tiltp, cvort_stretchid, cvort_stretchp, &
      cvort_bclnid, cvort_bclnp, cvort_xchgid, cvort_xchgp ) 
    implicit none

    integer,intent(in) :: ntimes, nparcels
    integer,intent(in) :: ncid
    integer,intent(in) :: tpid, xpid, ypid, zpid, upid, vpid, wpid, buoypid, pgfzpid
    integer,intent(in) :: pgfxpid, pgfypid
    integer,intent(in) :: pgfxbpid, pgfxdlpid, pgfxdnpid
    integer,intent(in) :: pgfybpid, pgfydlpid, pgfydnpid
    integer,intent(in) :: pgfzbpid, pgfzdlpid, pgfzdnpid
    integer,intent(in) :: tkepid, thrpid, rhodid, thpid, pipid
    integer,intent(in) :: qvpid, qcpid, qrpid, qipid, qspid, qgpid, qhlpid
    integer,intent(in) :: ncpid, nrpid, nipid, nspid, ngpid
    integer,intent(in) :: xvortpid, yvortpid, zvortpid, svortpid, cvortpid
    integer,intent(in) :: psipid
    integer,intent(in) :: zvort_tiltid, zvort_stretchid
    integer,intent(in) :: svort_stretchid, svort_tiltid, svort_bclnid
    integer,intent(in) :: svort_xchgid, dpsidtid
    integer,intent(in) :: cvort_stretchid, cvort_tiltid, cvort_bclnid
    integer,intent(in) :: cvort_xchgid 

    real,intent(in),dimension(ntimes) :: tp
    real,intent(in),dimension(ntimes,nparcels) :: xpc, ypc, zpc, up, vp, wp, buoyp
    real,intent(in),dimension(ntimes,nparcels) :: pgfxp, pgfyp, pgfzp
    real,intent(in),dimension(ntimes,nparcels) :: pgfxbp, pgfxdlp, pgfxdnp
    real,intent(in),dimension(ntimes,nparcels) :: pgfybp, pgfydlp, pgfydnp
    real,intent(in),dimension(ntimes,nparcels) :: pgfzbp, pgfzdlp, pgfzdnp
    real,intent(in),dimension(ntimes,nparcels) :: tkep, thrp, rhodp, thp, pip
    real,intent(in),dimension(ntimes,nparcels) :: qvp, qcp, qrp, qip, qsp, qgp, qhlp
    real,intent(in),dimension(ntimes,nparcels) :: ncp, nrp, nip, nsp, ngp
    real,intent(in),dimension(ntimes,nparcels) :: xvortp, yvortp, zvortp
    real,intent(in),dimension(ntimes,nparcels) :: svortp, cvortp, psip
    real,intent(in),dimension(ntimes,nparcels) :: zvort_tiltp, zvort_stretchp
    real,intent(in),dimension(ntimes,nparcels) :: svort_tiltp, svort_stretchp
    real,intent(in),dimension(ntimes,nparcels) :: svort_bclnp, svort_xchgp
    real,intent(in),dimension(ntimes,nparcels) :: dpsidt
    real,intent(in),dimension(ntimes,nparcels) :: cvort_tiltp, cvort_stretchp
    real,intent(in),dimension(ntimes,nparcels) :: cvort_bclnp, cvort_xchgp

    integer,dimension(2) :: s2, c2
    integer rcode


    s2 = (/1,1/)
    c2 = (/ntimes,nparcels/)

    call ncvpt( ncid, tpid, 1, ntimes, tp, rcode )
    call ncvpt( ncid, xpid, s2, c2, xpc, rcode )
    call ncvpt( ncid, ypid, s2, c2, ypc, rcode )
    call ncvpt( ncid, zpid, s2, c2, zpc, rcode )

    call ncvpt( ncid, upid, s2, c2, up, rcode )
    call ncvpt( ncid, vpid, s2, c2, vp, rcode )
    call ncvpt( ncid, wpid, s2, c2, wp, rcode )

    call ncvpt( ncid, qvpid, s2, c2, qvp, rcode )
    call ncvpt( ncid, qcpid, s2, c2, qcp, rcode )
    call ncvpt( ncid, qrpid, s2, c2, qrp, rcode )
    call ncvpt( ncid, qipid, s2, c2, qip, rcode )
    call ncvpt( ncid, qspid, s2, c2, qsp, rcode )
    call ncvpt( ncid, qgpid, s2, c2, qgp, rcode )
    if(ptype.eq.27)then
      call ncvpt( ncid, qhlpid, s2, c2, qhlp, rcode )
    endif

    call ncvpt( ncid, ncpid, s2, c2, ncp, rcode )
    call ncvpt( ncid, nrpid, s2, c2, nrp, rcode )
    call ncvpt( ncid, nipid, s2, c2, nip, rcode )
    call ncvpt( ncid, nspid, s2, c2, nsp, rcode )
    call ncvpt( ncid, ngpid, s2, c2, ngp, rcode )

    if(output_dyn.eq.1)then
      call ncvpt( ncid, buoypid, s2, c2, buoyp, rcode )

      call ncvpt( ncid, pgfxpid, s2, c2, pgfxp, rcode )
      call ncvpt( ncid, pgfypid, s2, c2, pgfyp, rcode )
      call ncvpt( ncid, pgfzpid, s2, c2, pgfzp, rcode )

      call ncvpt( ncid, pgfxbpid, s2, c2, pgfxbp, rcode )
      call ncvpt( ncid, pgfybpid, s2, c2, pgfybp, rcode )
      call ncvpt( ncid, pgfzbpid, s2, c2, pgfzbp, rcode )

      call ncvpt( ncid, pgfxdlpid, s2, c2, pgfxdlp, rcode )
      call ncvpt( ncid, pgfydlpid, s2, c2, pgfydlp, rcode )
      call ncvpt( ncid, pgfzdlpid, s2, c2, pgfzdlp, rcode )

      call ncvpt( ncid, pgfxdnpid, s2, c2, pgfxdnp, rcode )
      call ncvpt( ncid, pgfydnpid, s2, c2, pgfydnp, rcode )
      call ncvpt( ncid, pgfzdnpid, s2, c2, pgfzdnp, rcode )
    endif

    call ncvpt( ncid, tkepid, s2, c2, tkep,  rcode )

    call ncvpt( ncid, thrpid, s2, c2, thrp, rcode )
    call ncvpt( ncid, rhodid, s2, c2, rhodp, rcode )
    call ncvpt( ncid, thpid,  s2, c2, thp, rcode )
    call ncvpt( ncid, pipid,  s2, c2, pip, rcode )

    if(output_vor.eq.1)then
      call ncvpt( ncid, xvortpid, s2, c2, xvortp, rcode )
      call ncvpt( ncid, yvortpid, s2, c2, yvortp, rcode )
      call ncvpt( ncid, zvortpid, s2, c2, zvortp, rcode )
      call ncvpt( ncid, svortpid, s2, c2, svortp, rcode )
      call ncvpt( ncid, cvortpid, s2, c2, cvortp, rcode )

      call ncvpt( ncid, psipid, s2, c2, psip, rcode )

      call ncvpt( ncid, zvort_tiltid, s2, c2, zvort_tiltp, rcode )
      call ncvpt( ncid, zvort_stretchid, s2, c2, zvort_stretchp, rcode )

      call ncvpt( ncid, svort_tiltid, s2, c2, svort_tiltp, rcode )
      call ncvpt( ncid, svort_stretchid, s2, c2, svort_stretchp, rcode )
      call ncvpt( ncid, svort_bclnid, s2, c2, svort_bclnp, rcode )
      call ncvpt( ncid, svort_xchgid, s2, c2, svort_xchgp, rcode )

      call ncvpt( ncid, dpsidtid, s2, c2, dpsidt, rcode )

      call ncvpt( ncid, cvort_tiltid, s2, c2, cvort_tiltp, rcode )
      call ncvpt( ncid, cvort_stretchid, s2, c2, cvort_stretchp, rcode )
      call ncvpt( ncid, cvort_bclnid, s2, c2, cvort_bclnp, rcode )
      call ncvpt( ncid, cvort_xchgid, s2, c2, cvort_xchgp, rcode )

    endif

    call ncclos( ncid, rcode )

  end subroutine
  !############################################################################!
  !#                SUBROUTINE NC_CIRCUIT_WRITEOUT                            #!
  !#                                                                          #!
  !#  Writes out variables for circuit.                                       #!
  !#==========================================================================#!
  !# v1.0, Ryan Hastings, 14 Dec 2016                                         #!
  !############################################################################!
  subroutine nc_circuit_writeout( ntimes, nparcels, ncid, tpid, tp, &
    xpid, xpc, ypid, ypc, zpid, zpc, upid, up, vpid, vp, wpid, wp, &
    buoypid, buoyp, pgfxpid, pgfxp, pgfypid, pgfyp, pgfzpid, pgfzp, &
    pgfxbpid, pgfxbp, pgfybpid, pgfybp, pgfzbpid, pgfzbp, &
    pgfxdlpid, pgfxdlp, pgfydlpid, pgfydlp, pgfzdlpid, pgfzdlp, &
    pgfxdnpid, pgfxdnp, pgfydnpid, pgfydnp, pgfzdnpid, pgfzdnp, &
    tkepid, tkep, &
      thrpid, thrp, qvpid, qvp, qcpid, qcp, qrpid, qrp, qipid, qip, &
      qspid, qsp, qgpid, qgp, qhlpid, qhlp, &
      ncpid, ncp, nrpid, nrp, nipid, nip, nspid, nsp, ngpid, ngp, rhodid, rhodp, &
      thpid, thp, pipid, pip, xvortpid, xvortp, yvortpid, yvortp, zvortpid, &
      zvortp, svortpid, svortp, cvortpid, cvortp, psipid, psip, &
      zvort_tiltid, zvort_tiltp, zvort_stretchid, zvort_stretchp, &
      svort_tiltid, svort_tiltp, svort_stretchid, svort_stretchp, &
      svort_bclnid, svort_bclnp, svort_xchgid, svort_xchgp, dpsidtid, dpsidt,&
      cvort_tiltid, cvort_tiltp, cvort_stretchid, cvort_stretchp, &
      cvort_bclnid, cvort_bclnp, cvort_xchgid, cvort_xchgp, &
      circid, circ, vdotl1id, vdotl1, vdotl2id, vdotl2, intBdzid, intBdz )
    implicit none
      
    integer,intent(in) :: ntimes, nparcels
    integer,intent(in) :: ncid
    integer,intent(in) :: tpid, xpid, ypid, zpid, upid, vpid, wpid, buoypid, pgfzpid
    integer,intent(in) :: pgfxpid, pgfypid
    integer,intent(in) :: pgfxbpid, pgfxdlpid, pgfxdnpid
    integer,intent(in) :: pgfybpid, pgfydlpid, pgfydnpid
    integer,intent(in) :: pgfzbpid, pgfzdlpid, pgfzdnpid
    integer,intent(in) :: tkepid, thrpid, rhodid, thpid, pipid
!removed diabid
    integer,intent(in) :: qvpid, qcpid, qrpid, qipid, qspid, qgpid, qhlpid
    integer,intent(in) :: ncpid, nrpid, nipid, nspid, ngpid
    integer,intent(in) :: xvortpid, yvortpid, zvortpid, svortpid, cvortpid
    integer,intent(in) :: psipid
    integer,intent(in) :: zvort_tiltid, zvort_stretchid
    integer,intent(in) :: svort_stretchid, svort_tiltid, svort_bclnid
    integer,intent(in) :: svort_xchgid, dpsidtid
    integer,intent(in) :: cvort_stretchid, cvort_tiltid, cvort_bclnid
    integer,intent(in) :: cvort_xchgid
    integer,intent(in) :: circid, vdotl1id, vdotl2id, intBdzid

    real,intent(in),dimension(ntimes) :: tp
    real,intent(in),dimension(ntimes,nparcels) :: xpc, ypc, zpc, up, vp, wp, buoyp
    real,intent(in),dimension(ntimes,nparcels) :: pgfxp, pgfyp, pgfzp
    real,intent(in),dimension(ntimes,nparcels) :: pgfxbp, pgfxdlp, pgfxdnp
    real,intent(in),dimension(ntimes,nparcels) :: pgfybp, pgfydlp, pgfydnp
    real,intent(in),dimension(ntimes,nparcels) :: pgfzbp, pgfzdlp, pgfzdnp
    real,intent(in),dimension(ntimes,nparcels) :: tkep, thrp, rhodp, thp, pip
!removed diabp
    real,intent(in),dimension(ntimes,nparcels) :: qvp, qcp, qrp, qip, qsp, qgp, qhlp
    real,intent(in),dimension(ntimes,nparcels) :: ncp, nrp, nip, nsp, ngp
    real,intent(in),dimension(ntimes,nparcels) :: xvortp, yvortp, zvortp
    real,intent(in),dimension(ntimes,nparcels) :: svortp, cvortp, psip
    real,intent(in),dimension(ntimes,nparcels) :: zvort_tiltp, zvort_stretchp
    real,intent(in),dimension(ntimes,nparcels) :: svort_tiltp, svort_stretchp
    real,intent(in),dimension(ntimes,nparcels) :: svort_bclnp, svort_xchgp
    real,intent(in),dimension(ntimes,nparcels) :: dpsidt
    real,intent(in),dimension(ntimes,nparcels) :: cvort_tiltp, cvort_stretchp
    real,intent(in),dimension(ntimes,nparcels) :: cvort_bclnp, cvort_xchgp
    real,intent(in),dimension(ntimes)          :: circ, intBdz
    real,intent(in),dimension(ntimes,nparcels) :: vdotl1, vdotl2
    
    integer,dimension(2) :: s2, c2 
    integer rcode 
      
    
    s2 = (/1,1/)
    c2 = (/ntimes,nparcels/)
    
    call ncvpt( ncid, tpid, 1, ntimes, tp, rcode ) 
    call ncvpt( ncid, xpid, s2, c2, xpc, rcode ) 
    call ncvpt( ncid, ypid, s2, c2, ypc, rcode ) 
    call ncvpt( ncid, zpid, s2, c2, zpc, rcode )
    
    call ncvpt( ncid, upid, s2, c2, up, rcode )
    call ncvpt( ncid, vpid, s2, c2, vp, rcode )
    call ncvpt( ncid, wpid, s2, c2, wp, rcode )
      
    call ncvpt( ncid, qvpid, s2, c2, qvp, rcode ) 
    call ncvpt( ncid, qcpid, s2, c2, qcp, rcode )
    call ncvpt( ncid, qrpid, s2, c2, qrp, rcode )
    call ncvpt( ncid, qipid, s2, c2, qip, rcode )
    call ncvpt( ncid, qspid, s2, c2, qsp, rcode )
    call ncvpt( ncid, qgpid, s2, c2, qgp, rcode )
    if(ptype.eq.27)then
      call ncvpt( ncid, qhlpid, s2, c2, qhlp, rcode )
    endif
    
    call ncvpt( ncid, ncpid, s2, c2, ncp, rcode )
    call ncvpt( ncid, nrpid, s2, c2, nrp, rcode )
    call ncvpt( ncid, nipid, s2, c2, nip, rcode )
    call ncvpt( ncid, nspid, s2, c2, nsp, rcode )
    call ncvpt( ncid, ngpid, s2, c2, ngp, rcode )

    if(output_dyn.eq.1)then
      call ncvpt( ncid, buoypid, s2, c2, buoyp, rcode )

      call ncvpt( ncid, pgfxpid, s2, c2, pgfxp, rcode )
      call ncvpt( ncid, pgfypid, s2, c2, pgfyp, rcode )
      call ncvpt( ncid, pgfzpid, s2, c2, pgfzp, rcode )

      call ncvpt( ncid, pgfxbpid, s2, c2, pgfxbp, rcode )
      call ncvpt( ncid, pgfybpid, s2, c2, pgfybp, rcode )
      call ncvpt( ncid, pgfzbpid, s2, c2, pgfzbp, rcode )

      call ncvpt( ncid, pgfxdlpid, s2, c2, pgfxdlp, rcode )
      call ncvpt( ncid, pgfydlpid, s2, c2, pgfydlp, rcode )
      call ncvpt( ncid, pgfzdlpid, s2, c2, pgfzdlp, rcode )

      call ncvpt( ncid, pgfxdnpid, s2, c2, pgfxdnp, rcode )
      call ncvpt( ncid, pgfydnpid, s2, c2, pgfydnp, rcode )
      call ncvpt( ncid, pgfzdnpid, s2, c2, pgfzdnp, rcode )
    endif

    !call ncvpt( ncid, diabid, s2, c2, diabp, rcode ) CJN commented out
    call ncvpt( ncid, tkepid, s2, c2, tkep,  rcode )

    call ncvpt( ncid, thrpid, s2, c2, thrp, rcode )
    call ncvpt( ncid, rhodid, s2, c2, rhodp, rcode )
    call ncvpt( ncid, thpid,  s2, c2, thp, rcode )
    call ncvpt( ncid, pipid,  s2, c2, pip, rcode )

    call ncvpt( ncid, circid, 1, ntimes, circ, rcode )
    call ncvpt( ncid, vdotl1id, s2, c2, vdotl1, rcode )
    call ncvpt( ncid, vdotl2id, s2, c2, vdotl2, rcode )
    call ncvpt( ncid, intBdzid, 1, ntimes, intBdz, rcode )

    if(output_vor.eq.1)then
      call ncvpt( ncid, xvortpid, s2, c2, xvortp, rcode )
      call ncvpt( ncid, yvortpid, s2, c2, yvortp, rcode )
      call ncvpt( ncid, zvortpid, s2, c2, zvortp, rcode )
      call ncvpt( ncid, svortpid, s2, c2, svortp, rcode )
      call ncvpt( ncid, cvortpid, s2, c2, cvortp, rcode )

      call ncvpt( ncid, psipid, s2, c2, psip, rcode )

      call ncvpt( ncid, zvort_tiltid, s2, c2, zvort_tiltp, rcode )
      call ncvpt( ncid, zvort_stretchid, s2, c2, zvort_stretchp, rcode )

      call ncvpt( ncid, svort_tiltid, s2, c2, svort_tiltp, rcode )
      call ncvpt( ncid, svort_stretchid, s2, c2, svort_stretchp, rcode )
      call ncvpt( ncid, svort_bclnid, s2, c2, svort_bclnp, rcode )
      call ncvpt( ncid, svort_xchgid, s2, c2, svort_xchgp, rcode )

      call ncvpt( ncid, dpsidtid, s2, c2, dpsidt, rcode )

      call ncvpt( ncid, cvort_tiltid, s2, c2, cvort_tiltp, rcode )
      call ncvpt( ncid, cvort_stretchid, s2, c2, cvort_stretchp, rcode )
      call ncvpt( ncid, cvort_bclnid, s2, c2, cvort_bclnp, rcode )
      call ncvpt( ncid, cvort_xchgid, s2, c2, cvort_xchgp, rcode )

    endif

    call ncclos( ncid, rcode )

  end subroutine    

    !===========================================================================!
    !###########################################################################!
    !#                         SUBROUTINE NC_DYN_WRITEOUT                      #!
    !#                                                                         #!
    !# Write out the variables from DYN to a netCDF file.  Besides buoyancy    #!
    !# (BUOY), there are the perturbation pressures and forces.                #!
    !#-------------------------------------------------------------------------#!
    !# Ryan Hastings, 9 Sep 2011                                               #!
    !#-------------------------------------------------------------------------#!
    !# v1.0, rewrote comments, eliminated pi_spin, and hydrostatic pressure    #!
    !# stuff, Ryan Hastings 14 Dec 2016                                        #!
    !###########################################################################!
    subroutine nc_dyn_writeout( outname, thr0, thrpert, &
        pi_dl, p_dl, pi_b, p_b, pi_dn, p_dn, &
        pgfb, pgfdl, pgfdn )

      !#########################################################################!
      !############################## VARIABLES ################################!

      !-------------------------- input variables ------------------------------!
      character(len=100)            outname
      real,dimension(ni,nj,nk)   :: thrpert
      real,dimension(nk)         :: thr0
      real,dimension(ni,nj,nk)   :: pi_dl, p_dl, pi_b, p_b, pi_dn, p_dn
      real,intent(in),dimension(ni,nj,nk+1,3) :: pgfb, pgfdl, pgfdn

      !-------------------------- local variables ------------------------------!

      integer                  ncid,  rcode, i
      integer  nidid, nip1did, njdid, njp1did, nkdid, nkp1did, timedid, itdid
      integer  xhid, xfid, yhid, yfid, zhid, zfid, timeid
      integer  thrpertid, thr0id, theid
      integer  pibid, pbid
      integer  pidlid, pdlid
      integer  pidnid, pdnid
      integer  pgfxbid, pgfybid, pgfzbid
      integer  pgfxdlid, pgfydlid, pgfzdlid
      integer  pgfxdnid, pgfydnid, pgfzdnid
      integer,dimension(3)  :: dimids, strt3, dimidu, dimidv, dimidw, dimidsp1

      !#########################################################################!
      !############################## MAIN BODY ################################!

      write(*,*)
      write(*,*) '-----------------------------------------------------------'
      write(*,*) 'opening ',outname
      ncid = nccre( outname, NCNOCLOB, rcode )

      !-------------------------------------------------------------------------!
      ! create dimensions

      write(*,*) 'creating dimensions'
      nidid = ncddef(ncid, 'ni', ni, rcode)
      nip1did = ncddef(ncid, 'nip1', ni+1, rcode )
      njdid = ncddef(ncid, 'nj', nj, rcode)
      njp1did = ncddef(ncid, 'njp1', nj+1, rcode )
      nkdid = ncddef(ncid, 'nk', nk, rcode)
      nkp1did = ncddef(ncid,'nkp1',nk+1,rcode)
      timedid = ncddef(ncid, 'time', 1, rcode)

      !-------------------------------------------------------------------------!
      ! define variables

      ! count array
      dimids = (/nidid,njdid,nkdid/)
      dimidsp1 = (/nidid,njdid,nkp1did/)

      write(*,*) 'defining spatial arrays'

      xhid = ncvdef(ncid, 'xh', nf_real, 1, nidid, rcode)
      xfid = ncvdef(ncid, 'xf', nf_real, 1, nip1did, rcode)
      yhid = ncvdef(ncid, 'yh', nf_real, 1, njdid, rcode)
      yfid = ncvdef(ncid, 'yf', nf_real, 1, njp1did, rcode )
      zhid = ncvdef(ncid, 'zh', nf_real, 1, nkdid, rcode)
      zfid = ncvdef(ncid, 'zf', nf_real, 1, nkp1did, rcode )
      timeid = ncvdef(ncid, 'time', nf_real, 1, timedid, rcode)

      thr0id = ncvdef( ncid, 'thr0', nf_real, 1, nkdid, rcode )
      thrpertid = ncvdef( ncid, 'thrpert', nf_real, 3, dimids, rcode )

      ! linear dynamic pressure perturbation
      if(output_dl.eq.1)then

        write(*,*) 'defining linear dynamic perturbation pressure variables'

        pdlid = ncvdef( ncid, 'p_dl', nf_real, 3, dimids, rcode )
        pidlid = ncvdef( ncid, 'pi_dl', nf_real, 3, dimids, rcode )
        pgfxdlid = ncvdef( ncid, 'pgfx_dl',nf_real,3,dimids,rcode )
        pgfydlid = ncvdef( ncid, 'pgfy_dl',nf_real,3,dimids,rcode )
        pgfzdlid = ncvdef( ncid, 'pgfz_dl',nf_real,3,dimidsp1,rcode )
  
      endif

      if(output_dn.eq.1)then

        write(*,*) 'defining nonlinear dynamic perturbation pressure variables'

        pdnid = ncvdef( ncid, 'p_dn', nf_real, 3, dimids, rcode )
        pidnid = ncvdef( ncid, 'pi_dn', nf_real, 3, dimids, rcode )
        pgfxdnid = ncvdef( ncid, 'pgfx_dn', nf_real, 3, dimids, rcode )
        pgfydnid = ncvdef( ncid, 'pgfy_dn', nf_real, 3, dimids, rcode )
        pgfzdnid = ncvdef( ncid, 'pgfz_dn', nf_real, 3, dimidsp1, rcode )

      endif


      ! buoyancy pressure perturbation

      if(output_b.eq.1)then

        write(*,*) 'defining buoyancy perturbation pressure variables'

        pbid = ncvdef( ncid, 'p_b', nf_real, 3, dimids, rcode )
        pibid = ncvdef( ncid, 'pi_b', nf_real, 3, dimids, rcode )
        pgfxbid = ncvdef( ncid, 'pgfx_b', nf_real, 3, dimids, rcode )
        pgfybid = ncvdef( ncid, 'pgfy_b', nf_real, 3, dimids, rcode )
        pgfzbid = ncvdef( ncid, 'pgfz_b', nf_real, 3, dimidsp1, rcode )

      endif

      call ncendf(ncid,rcode)
      write(*,*) 'done defining variables'

      !---------------------------------------------------------------------------!
      ! writing variables

      write(*,*) 'writing variables'
      write(*,*) '..space and time..'
      ! first convert back to km
      xh = xh/1000
      yh = yh/1000
      zh = zh/1000

      ! counting arrays
      strt3 = (/1,1,1/)
      dimids = (/ni,nj,nk/)
      dimidu = (/nip1,nj,nk/)
      dimidv = (/ni,njp1,nk/)
      dimidw = (/ni,nj,nkp1/)

      ! put spatial variables

      call ncvpt( ncid,   xhid, 1, ni, xh(1:ni), rcode )
      call ncvpt( ncid,   xfid, 1, ni+1, xf(1:ni+1), rcode )
      call ncvpt( ncid,   yhid, 1, nj, yh(1:nj), rcode )
      call ncvpt( ncid,   yfid, 1, nj+1, yf(1:nj+1), rcode )
      call ncvpt( ncid,   zhid, 1, nk, zh(1:nk), rcode )
      call ncvpt( ncid, zfid, 1, nk+1, zf(1:nk+1), rcode )
      call ncvpt( ncid, timeid, 1,  1,     time, rcode )

      call ncvpt( ncid, thr0id, 1, nk, thr0, rcode )
      call ncvpt( ncid, thrpertid, strt3, dimids, thrpert, rcode )

      if(output_dl.eq.1)then

        write(*,*) '..linear dynamic pressure..'

        call ncvpt( ncid, pidlid, strt3, dimids, pi_dl, rcode )
        call ncvpt( ncid, pdlid, strt3, dimids, p_dl, rcode )
        call ncvpt( ncid, pgfxdlid, strt3, dimids, pgfdl(1:ni,1:nj,1:nk,1), rcode )
        call ncvpt( ncid, pgfydlid, strt3, dimids, pgfdl(1:ni,1:nj,1:nk,2), rcode )
        call ncvpt( ncid, pgfzdlid, strt3, dimids, pgfdl(1:ni,1:nj,1:nk+1,3), rcode )

      endif

      if(output_dn.eq.1)then

        write(*,*) '..nonlinear dynamic pressure..'

        call ncvpt( ncid, pidnid, strt3, dimids, pi_dn, rcode )
        call ncvpt( ncid, pdnid, strt3, dimids, p_dn, rcode )
        call ncvpt( ncid, pgfxdnid, strt3, dimids, pgfdn(1:ni,1:nj,1:nk,1), rcode )
        call ncvpt( ncid, pgfydnid, strt3, dimids, pgfdn(1:ni,1:nj,1:nk,2), rcode )
        call ncvpt( ncid, pgfzdnid, strt3, dimids, pgfdn(1:ni,1:nj,1:nk+1,3), rcode )

       ! call ncvpt( ncid, pispinid, strt3, dimids, pi_spin, rcode ) commented
       ! out by CJN

      endif

      if(output_b.eq.1)then

        write(*,*) '..buoyancy pressure..'

        call ncvpt( ncid, pibid, strt3, dimids, pi_b, rcode )
        call ncvpt( ncid, pbid, strt3, dimids, p_b, rcode )
        call ncvpt( ncid, pgfxbid, strt3, dimids, pgfb(1:ni,1:nj,1:nk,1), rcode )
        call ncvpt( ncid, pgfybid, strt3, dimids, pgfb(1:ni,1:nj,1:nk,2), rcode )
        call ncvpt( ncid, pgfzbid, strt3, dimids, pgfb(1:ni,1:nj,1:nk+1,3), rcode )

      endif

      !-------------------------------------------------------------------------!
      ! closing file
      write(*,*) 'closing ',outname
      call ncclos( ncid, rcode )

    end subroutine


end module
