module traj_utils
  !#################################################################!
  !#################################################################!
  !####                                                         ####!
  !####                     MODULE TRAJ_UTILS                   ####!
  !####                                                         ####!
  !####  Contains procedures used by traj and circuit, and are  ####!
  !####  available for any other code that looks to do traj-    ####!
  !####  ectories.                                              ####!
  !####                                                         ####!
  !####  CONTAINS:                                              ####!
  !####    SUBROUTINE INTERP_ALL : interpolates several values  ####!
  !####       to the parcel positions from the grid             ####!
  !####    SUBROUTINE INTERP_WINDS : Interpolates winds from    ####!
  !####       grid to parcel position                           ####!
  !####    SUBROUTINE FIND_SURROUNDING_GRID_POINTS : Finds grid ####!
  !####       points surrounding parcel                         ####!
  !####    SUBROUTINE INTERP_U : Interpolates from u-momentum   ####!
  !####      grid to parcel position                            ####!
  !####    SUBROUTINE INTERP_V : As interp_u, for v-grid.       ####!
  !####    SUBROUTINE INTERP_W : As interp_u, for w-grid.       ####!
  !####    SUBROUTINE INTERP_SCALARS : Interpolates from scalar ####!
  !####      grid to parcel position                            ####!
  !####    SUBROUTINE INTERP_DYNZ : Interpolates PGFs to parcel ####!
  !####      position                                           ####!
  !####    SUBROUTINE RK4 : Advances parcel with 4th-order      ####!
  !####      Runke-Kutta method                                 ####!
  !####    SUBROUTINE COMPUTE_VORTS : Computes values related   ####!
  !####      to vorticity and vortex dynamics                   ####!
  !####    SUBROUTINE INTERP_VORT3 : Interpolates vorticity     ####!
  !####      values to parcel position                          ####!
  !####    SUBROUTINE INTERP_XVORT : Interpolates x-vorticity   ####!
  !####    SUBROUTINE INTERP_YVORT : Interpolates y-vorticity   ####!
  !####    SUBROUTINE INTERP_ZVORT : Interpolates z-vorticity   ####!
  !####      to parcel position                                 ####!
  !####    REAL FUNCTION COMPUTE_DPSI2 : Computes difference    ####!
  !####      between two angles                                 ####!
  !####    SUBROUTINE COMPUTE_PSI : Computes wind angle         ####!
  !####=========================================================####!
  !#### v1.0, Ryan Hastings 13 December 2016                    ####!
  !#################################################################!
  !#################################################################!
  use globals
  implicit none

contains

   subroutine interpintime_winds( u1, u2, u, v1, v2, v, w1, w2, w, t1, t2, ptime )
     !8 mar 2018
     implicit none
     real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: u1, u2
     real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: v1, v2
     real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: w1, w2

     real,intent(in) :: t1, t2, ptime

     real,intent(out),dimension(1:ni+1,1:nj,0:nk) :: u
     real,intent(out),dimension(1:ni,1:nj+1,0:nk) :: v
     real,intent(out),dimension(1:ni,1:nj,1:nk+1) :: w

     call interp1_u(u1,u2,t1,t2,ptime,u)
     call interp1_v(v1,v1,t1,t2,ptime,v)
     call interp1_w(w1,w2,t1,t2,ptime,w)

  end subroutine

  subroutine interpintime_scalars( buoy1, buoy2, buoy, tke1, tke2, tke, thr1, thr2, thr, &
        qv1, qv2, qv, qc1, qc2, qc, qr1, qr2, qr, qi1, qi2, qi, qs1, qs2, qs, &
        qg1, qg2, qg, qi1, qi2, qi, qs1, qs2, qs, qg1, qg2, qg, qhl1, qhl2, qhl, &
        ncc1, ncc2, ncc, ncr1, ncr2, ncr, nci1, nci2, nci, ncs1, ncs2, ncs, &
        ngi1, ngi1, ngi, rhod1, rhod2, rhod, th1, th2, th, ppi1, ppi1, ppi, t1, t2, ptime )

    implicit none
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: buoy1, buoy2, tke1, tke2, thr1, thr2
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: qv1, qv2, qc1, qc2, qi1, qi2, qs1, qs2
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: qg1, qg2, qr1, qr2, qhl1, qhl2
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: ncc1, ncc2, ncr1, ncr2, nci1, nci2
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: ncs1, ncs1, ncg1, ncg2, rhod1, rhod2
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: th1, th2, ppi1, ppi2
    real,intent(in) :: t1, t2, ptime

    real,intent(out),dimension(1:ni,1:nj,0:nk) :: buoy, tke, thr
    real,intent(out),dimension(1:ni,1:nj,0:nk) :: qv, qc, qi, qs
    real,intent(out),dimension(1:ni,1:nj,0:nk) :: qg, qr, qhl
    real,intent(out),dimension(1:ni,1:nj,0:nk) :: ncc, ncr, nci
    real,intent(out),dimension(1:ni,1:nj,0:nk) :: ncs, ncg, rhod
    real,intent(out),dimension(1:ni,1:nj,0:nk) :: th, ppi

    call interp1_scalar(buoy1,buoy2,t1,t2,ptime,buoy)
    call interp1_scalar(qv1,qv2,t1,t2,ptime,qv)
    call interp1_scalar(qc1,qc2,t1,t2,ptime,qc)
    call interp1_scalar(qr1,qr2,t1,t2,ptime,qr)
    call interp1_scalar(qi1,qi2,t1,t2,ptime,qi)
    call interp1_scalar(qs1,qs2,t1,t2,ptime,qs)
    call interp1_scalar(qg1,qg2,t1,t2,ptime,qg)
    call interp1_scalar(qhl1,qhl2,t1,t2,ptime,qhl)

    call interp1_scalar(ncv1,ncv2,t1,t2,ptime,ncv)
    call interp1_scalar(ncc1,ncc2,t1,t2,ptime,ncc)
    call interp1_scalar(ncr1,ncr2,t1,t2,ptime,ncr)
    call interp1_scalar(nci1,nci2,t1,t2,ptime,nci)
    call interp1_scalar(ncs1,ncs2,t1,t2,ptime,ncs)
    call interp1_scalar(ncg1,ncg2,t1,t2,ptime,ncg)

    call interp1_scalar(rhod1,rhod2,t1,t2,ptime,rhod)
    call interp1_scalar(th1,th2,t1,t2,ptime,th)
    call interp1_scalar(ppi1,ppi2,t1,t2,ptime,ppi)

  end subroutine
 

  !################################################################!
  !#            SUBROUTINE INTERP_ALL                             #!
  !#                                                              #!
  !#  Interpolates all values from grids to parcel position.      #!
  !#  Passed variables:                                           #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    real,dimension(:) zzh : extended scalar vertical points   #!
  !#    real,dimension(:,:,:) u, v, w : wind fields               #!
  !#    real,dimension(:,:,:) buoy : buoyancy field               #!
  !#    real,dimension(:,:,:) pgfx, pgfy, pgfz : total PGFs       #!
  !#    real,dimension(:,:,:) tke : turbulent kinetic energy      #!
  !#    real,dimension(:) thr0 : base state density potential     #!
  !#      temperature                                             #!
  !#    real,dimension(:,:,:) thr : total density potential temp  #!
  !#    real,dimension(:,:,:) qv, qc, qr, qi, qs, qg, qhl :       #!
  !#      mixing ratios                                           #!
  !#    real,dimension(:,:,:) ncc, ncr, nci, ncs, ncg :           #!
  !#      number concentrations                                   #!
  !#    real,dimension(:,:,:) rhod : dry air density              #!
  !#    real,dimension(:,:,:) th : potential temperature          #!
  !#    real,dimension(:,:,:) ppi : nondimensionalized pressure   #!
  !#  Returned variables:                                         #!
  !#    real up, vp, wp : winds at parcel position                #!
  !#    real buoyp : buoyancy of parcel                           #!
  !#    real pgfxp, pgfyp, pgfzp : PGF at parcel position         #!
  !#    real tkep : TKE at parcel position                        #!
  !#    real thrp : perturbation density potential temperature of #!
  !#      parcel                                                  #!
  !#    real qvp, qcp, qrp, qip, qsp, qgp, qhl : mixing ratios of #!
  !#      parcel                                                  #!
  !#    real ncp, nrp, nip, nsp, ngp : number concentration at    #!
  !#      parcel position                                         #!
  !#    real rhodp : dry air density of parcel                    #!
  !#    real thp : potential temperature of parcel                #!
  !#    real pip : nondimensionalized pressure of parcel          #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings 13 December 2016                         #!
  !################################################################!
  subroutine interp_all( xpc, ypc, zpc, zzh, u, up, v, vp,&
    w, wp, buoy, buoyp, &
    pgfx, pgfxp, pgfy, pgfyp, pgfz, pgfzp, &
    tke, tkep, thr0, thr, thrp, &
    qv, qvp, qc, qcp, qr, qrp, qi, qip, &
    qs, qsp, qg, qgp, qhl, qhlp, &
    ncc, ncp, ncr, nrp, nci, nip, ncs, nsp,&
    ncg, ngp, rhod, rhodp, &
    th, thp, ppi, pip )

    implicit none
    !######################## VARIABLES ###########################!
    ! Passed variables
    real,intent(inout) :: xpc, ypc, zpc

    real,dimension(1:ni+1,1:nj,0:nk) :: u, pgfx
    real,dimension(1:ni,1:nj+1,0:nk) :: v, pgfy
    real,dimension(1:ni,1:nj,1:nk+1) :: w, tke, pgfz

    real,dimension(1:ni,1:nj,0:nk) :: buoy, thr, th, ppi
    real,dimension(1:ni,1:nj,0:nk) :: qv, qc, qr, qi, qs, qg, qhl, rhod
    real,dimension(1:ni,1:nj,0:nk) :: ncc, ncr, nci, ncs, ncg

    real,intent(in),dimension(0:nk) :: zzh

    real,intent(inout) :: xpc, ypc, zpc

    ! Returned variables
    real,intent(out) :: up, vp, wp, buoyp, pgfxp, pgfyp, pgfzp, tkep, thrp
    real,intent(out) :: qvp, qcp, qrp, qip, qsp, qgp, qhlp
    real,intent(out) :: ncp, nrp, nip, nsp, ngp, rhodp, thp, pip

    ! Local variables
    integer :: ih, iif, jh, jf, kh, kf
    real :: dum_thr, dum_thr0
    real,dimension(0:nk):: thr0_2

    !########################### MAIN BODY ########################!

    call find_surrounding_grid_points( xpc, ypc, zpc, ih, iif, jh, jf, kh, kf )

         IF((Iif.NE.9999).and.(jf.ne.9999).and.(kf.ne.9999))THEN ! the parcel has not gone off the edge

    call interp_u( u, xpc, ypc, zpc, iif, jh, kh, zzh, up )
    call interp_u( pgfx, xpc, ypc, zpc, iif, jh, kh, zzh, pgfxp )

    call interp_v( v, xpc, ypc, zpc,  ih, jf, kh, zzh, vp )
    call interp_v( pgfy, xpc, ypc, zpc,  ih, jf, kh, zzh, pgfyp )

    call interp_w( w, xpc, ypc, zpc,  ih, jh, kf, wp )
    call interp_w( pgfz, xpc, ypc, zpc, ih, jh, kf, pgfzp )

    call interp_scalars( buoy, xpc, ypc, zpc, ih, jh, kh, zzh, buoyp )
    call interp_w( tke, xpc, ypc, zpc, ih, jh, kf, tkep )
    call interp_scalars( thr, xpc, ypc, zpc, ih, jh, kh, zzh, dum_thr )

    call interp_scalars( qv, xpc, ypc, zpc, ih, jh, kh, zzh, qvp )
    call interp_scalars( qc, xpc, ypc, zpc, ih, jh, kh, zzh, qcp )
    call interp_scalars( qr, xpc, ypc, zpc, ih, jh, kh, zzh, qrp )
    call interp_scalars( qi, xpc, ypc, zpc, ih, jh, kh, zzh, qip )
    call interp_scalars( qg, xpc, ypc, zpc, ih, jh, kh, zzh, qgp )
    if(ptype.eq.27)then
      call interp_scalars( qhl, xpc, ypc, zpc, ih, jh, kh, zzh, qhlp )
    endif
    call interp_scalars( qs, xpc, ypc, zpc, ih, jh, kh, zzh, qsp )

    call interp_scalars( ncc, xpc, ypc, zpc, ih, jh, kh, zzh, ncp )
    call interp_scalars( ncr, xpc, ypc, zpc, ih, jh, kh, zzh, nrp )
    call interp_scalars( nci, xpc, ypc, zpc, ih, jh, kh, zzh, nip )
    call interp_scalars( ncs, xpc, ypc, zpc, ih, jh, kh, zzh, nsp )
    call interp_scalars( ncg, xpc, ypc, zpc, ih, jh, kh, zzh, ngp )

    call interp_scalars( rhod, xpc, ypc, zpc, ih, jh, kh, zzh, rhodp )
    call interp_scalars( th, xpc, ypc, zpc, ih, jh, kh, zzh, thp )
    call interp_scalars( ppi, xpc, ypc, zpc, ih, jh, kh, zzh, pip )

    thr0_2(1:nk) = thr0
    thr0_2(0) = thr0(1)

    dum_thr0 = thr0_2(kh) + (zpc-zzh(kh))*(thr0_2(kh+1)-thr0_2(kh))/(zzh(kh+1)-zzh(kh))

    thrp=dum_thr-dum_thr0

        ELSE ! parcel has gone off the edge

    xpc=-9.e9
    ypc=-9.e9
    zpc=-9.e9

    up=-9.e9
    pgfxp=-9.e9
    vp=-9.e9
    pgfyp=-9.e9
    wp=-9.e9
    pgfzp=-9.e9

    buoyp=-9.e9
    tkep=-9.e9
    thrp=-9.e9
    qvp=-9.e9
    qcp=-9.e9
    qrp=-9.e9
    qip=-9.e9
    qsp=-9.e9
    qgp=-9.e9
    qhlp=-9.9e9

        ENDIF

  end subroutine

  !################################################################!
  !#            SUBROUTINE INTERP_WINDS                           #!
  !#                                                              #!
  !#  Interpolates winds to parcel position.  Mainly used by RK4. #!
  !#  Passed variables:                                           #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    real,dimension(:) zzh : extended vertical scalar points   #!
  !#    real,dimension(:,:,:) u, v, w : wind fields               #!
  !#  Returned variables:                                         #!
  !#    real up, vp, wp : winds at parcel position                #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_winds( xpc, ypc, zpc, zzh, u, up, v, vp, w, wp )

    implicit none

    real,intent(inout) :: xpc, ypc, zpc

    real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: u
    real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: v
    real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: w

    real,intent(in),dimension(0:nk) :: zzh

    real,intent(out) :: up, vp, wp

    integer :: ih, iif, jh, jf, kh, kf

    call find_surrounding_grid_points( xpc, ypc, zpc, ih, iif, jh, jf, kh, kf )

    call interp_u( u, xpc, ypc, zpc, iif, jh, kh, zzh, up )
    call interp_v( v, xpc, ypc, zpc,  ih, jf, kh, zzh, vp )
    call interp_w( w, xpc, ypc, zpc,  ih, jh, kf, wp )

  end subroutine

  !################################################################!
  !#            SUBROUTINE FINE_SURROUNDING_GRID_POINTS           #!
  !#                                                              #!
  !#  Given a parcel position, find the indices of the grid point #!
  !#  in the lower southwest corner on each grid.                 #!
  !#  Passed variables:                                           #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#  Returned variables:                                         #!
  !#    integer ih, jh, kh : indices of lower southwest corner on #!
  !#      scalar grid                                             #!
  !#    integer iif, jf, kf : indices of lower southwest corner   #!
  !#      on u-, v-, and w-grids respectively                     #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 December 2016                        #!
  !################################################################!
  subroutine find_surrounding_grid_points( xpc, ypc, zpc, ih, iif, jh, jf, kh, kf )

    real,intent(inout) :: xpc, ypc, zpc
    integer,intent(out) :: ih, iif, jh, jf, kh, kf

    integer,dimension(1) :: tmp

    tmp=minloc(abs(xh(:)-xpc))
    ih=tmp(1)
    if(xpc.lt.xh(1))then
      ih=1
    elseif(xpc.gt.xh(ni))then
      ih=ni-1
    else
      if(xh(ih).gt.xpc)then
        ih=tmp(1)-1
      else
        ih=tmp(1)
      endif
    endif

    if(nj.eq.1)then
      jh=1
      jf=1
    else
      tmp=minloc(abs(yh(:)-ypc))
      jh=tmp(1)
      if(ypc.lt.yh(1))then
        jh=1
      elseif(ypc.gt.yh(nj))then
        jh=nj-1
      else
        if(yh(jh).gt.ypc)then
          jh=tmp(1)-1
        else
          jh=tmp(1)
        endif
      endif
    endif

    tmp=minloc(abs(zh(:)-zpc))
    kh=tmp(1)
    if(zpc.lt.zh(1))then
      kh=0
    elseif(zpc.gt.zh(nk))then
      kh=nk-1
    else
      if(zh(kh).gt.zpc)then
        kh=tmp(1)-1
      else
        kh=tmp(1)
      endif
    endif

    tmp=minloc(abs(xf(:)-xpc))
    iif=tmp(1)
    if(xpc.lt.xf(1))then
      iif=1
      xpc=xf(1)
    elseif(xpc.gt.xf(ni+1))then
      iif=ni
      xpc=xf(ni)
    else
      if(xf(iif).gt.xpc)then
        iif=tmp(1)-1
      else
        iif=tmp(1)
      endif
    endif

    tmp=minloc(abs(yf(:)-ypc))
    jf=tmp(1)
    if(ypc.lt.yf(1))then
      jf=1
      ypc=yf(1)
    elseif(ypc.gt.yf(nj+1))then
      jf=nj
      ypc=yf(nj)
    else
      if(yf(jf).gt.ypc)then
        jf=tmp(1)-1
      else
        jf=tmp(1)
      endif
    endif

    tmp=minloc(abs(zf(:)-zpc))
    kf=tmp(1)
    if(zpc.lt.0)then
      kf=1
      zpc=0
    elseif(zpc.gt.zf(nk+1))then
      kf=nk
      zpc=zf(nk)
    else
      if(zf(kf).gt.zpc)then
        kf=tmp(1)-1
      else
        kf=tmp(1)
      endif
    endif

  end subroutine

  subroutine interp1_scalar(s1,s2,t1,t2,ptime,s)
    implicit none
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: s1, s2
    real,intent(in) :: t1, t2, ptime
    real,intent(out),dimension(1:ni,1:nj,0:nk) :: s

    s=s1+(t2-ptime)*(s2-s1)/(t2-t1)
  endif

  subroutine interp1_u(s1,s2,t1,t2,ptime,s)
    implicit none
    real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: s1, s2
    real,intent(in) :: t1, t2, ptime
    real,intent(out),dimension(1:ni+1,1:nj,0:nk) :: s
    
    s=s1+(t2-ptime)*(s2-s1)/(t2-t1)
  endif

  subroutine interp1_v(s1,s2,t1,t2,ptime,s)
    implicit none
    real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: s1, s2
    real,intent(in) :: t1, t2, ptime
    real,intent(out),dimension(1:ni,1:nj+1,0:nk) :: s
    
    s=s1+(t2-ptime)*(s2-s1)/(t2-t1)
  endif

  subroutine interp1_w(s1,s2,t1,t2,ptime,s)
    implicit none
    real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: s1, s2
    real,intent(in) :: t1, t2, ptime
    real,intent(out),dimension(1:ni,1:nj,1:nk+1) :: s
    
    s=s1+(t2-ptime)*(s2-s1)/(t2-t1)
  endif

    
  !################################################################!
  !#                  SUBROUTINE INTERP_U                         #!
  !#                                                              #!
  !#  Interpolates value from u-grid (e.g., u, pgfx) to parcel    #!
  !#  position.                                                   #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) u : value                           #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    integer iif, jh, kh : indices of lower southwest corner   #!
  !#      of grid box surrounding parcel                          #!
  !#  Returned variable:                                          #!
  !#    real up : value at parcel position                        #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_u( u, xpc, ypc, zpc, iif, jh, kh, zzh, up )

    implicit none

    real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: u
    real,intent(in) :: xpc, ypc, zpc
    integer,intent(in) :: iif, jh, kh
    real,intent(in),dimension(0:nk) :: zzh

    real,intent(out) :: up

    real :: dx, dx1, dx2, dy, dy1, dy2, dz_, dz1, dz2

    if(nj.eq.1)then
      dx =   xf(iif+1)-xf(iif)
      dx1 = (xf(iif+1)-xpc)    /dx
      dx2 = (xpc-      xf(iif))/dx

      dz_ = zzh(kh+1)-zzh(kh)
      dz1 = (zzh(kh+1)-zpc)/dz_
      dz2 = (zpc-zzh(kh))/dz_

      up = u(iif  ,jh  ,kh  )*dx1*dz1 + &
           u(iif  ,jh  ,kh+1)*dx1*dz2 + &
           u(iif+1,jh  ,kh  )*dx2*dz1 + &
           u(iif+1,jh  ,kh+1)*dx2*dz2 + &
    else
      dx =   xf(iif+1)-xf(iif)
      dx1 = (xf(iif+1)-xpc)    /dx
      dx2 = (xpc-      xf(iif))/dx

      dy =   yh(jh+1)-yh(jh)
      dy1 = (yh(jh+1)-ypc)   /dy
      dy2 = (ypc     -yh(jh))/dy

      dz_ = zzh(kh+1)-zzh(kh)
      dz1 = (zzh(kh+1)-zpc)/dz_
      dz2 = (zpc-zzh(kh))/dz_

      up = u(iif  ,jh  ,kh  )*dx1*dy1*dz1 + &
           u(iif  ,jh  ,kh+1)*dx1*dy1*dz2 + &
           u(iif  ,jh+1,kh  )*dx1*dy2*dz1 + &
           u(iif  ,jh+1,kh+1)*dx1*dy2*dz2 + &
           u(iif+1,jh  ,kh  )*dx2*dy1*dz1 + &
           u(iif+1,jh  ,kh+1)*dx2*dy1*dz2 + &
           u(iif+1,jh+1,kh  )*dx2*dy2*dz1 + &
           u(iif+1,jh+1,kh+1)*dx2*dy2*dz2
    endif

  end subroutine

  !################################################################!
  !#               SUBROUTINE INTERP_V                            #!
  !#                                                              #!
  !#  Interpolates from v-grid to parcel position.                #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) v : value on v-grid                 #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    integer ih, jf, kh : indices of lower southwest corner of #!
  !#      box surrounding parcel on grid                          #!
  !#  Returned variable:                                          #!
  !#    real vp : interpolated to parcel position                 #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_v( v, xpc, ypc, zpc, ih, jf, kh, zzh, vp )

    implicit none

    real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: v
    real,intent(in) :: xpc, ypc, zpc
    integer,intent(in) :: ih, jf, kh
    real,intent(in),dimension(0:nk) :: zzh

    real,intent(out) :: vp

    real :: dx, dx1, dx2, dy, dy1, dy2, dz_, dz1, dz2

    if(nj.eq.1)then
      dx =   xh(ih+1)-xh(ih)
      dx1 = (xh(ih+1)-xpc)    /dx
      dx2 = (xpc-      xh(ih))/dx

      dz_ = zzh(kh+1)-zzh(kh)
      dz1 = (zzh(kh+1)-zpc)/dz_
      dz2 = (zpc-zzh(kh))/dz_

      vp = v(ih  ,jf  ,kh  )*dx1*dz1 + &
           v(ih  ,jf  ,kh+1)*dx1*dz2 + &
           v(ih+1,jf  ,kh  )*dx2*dz1 + &
           v(ih+1,jf  ,kh+1)*dx2*dz2 + &

    else
      dx =   xh(ih+1)-xh(ih)
      dx1 = (xh(ih+1)-xpc)    /dx
      dx2 = (xpc-      xh(ih))/dx

      dy =   yf(jf+1)-yf(jf)
      dy1 = (yf(jf+1)-ypc)   /dy
      dy2 = (ypc     -yf(jf))/dy

      dz_ = zzh(kh+1)-zzh(kh)
      dz1 = (zzh(kh+1)-zpc)/dz_
      dz2 = (zpc-zzh(kh))/dz_

      vp = v(ih  ,jf  ,kh  )*dx1*dy1*dz1 + &
           v(ih  ,jf  ,kh+1)*dx1*dy1*dz2 + &
           v(ih  ,jf+1,kh  )*dx1*dy2*dz1 + &
           v(ih  ,jf+1,kh+1)*dx1*dy2*dz2 + &
           v(ih+1,jf  ,kh  )*dx2*dy1*dz1 + &
           v(ih+1,jf  ,kh+1)*dx2*dy1*dz2 + &
           v(ih+1,jf+1,kh  )*dx2*dy2*dz1 + &
           v(ih+1,jf+1,kh+1)*dx2*dy2*dz2
    endif

  end subroutine

  !################################################################!
  !#               SUBROUTINE INTERP_W                            #!
  !#                                                              #!
  !#  Interpolates from w-grid to parcel position.                #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) w : value on w-grid                 #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    integer ih, jh, kf : indices of lower southwest corner of #!
  !#      box on grid surrounding parcel                          #!
  !#  Returned variable:                                          #!
  !#    real wp : value interpolated to parcel position           #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_w( w, xpc, ypc, zpc, ih, jh, kf, wp )

    implicit none

    real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: w
    real,intent(in) :: xpc, ypc, zpc
    integer,intent(in) :: ih, jh, kf
    real,intent(out) :: wp

    real :: dx, dx1, dx2, dy, dy1, dy2, dz_, dz1, dz2


    if(nj.eq.1)then
      dx =   xh(ih+1)-xh(ih)
      dx1 = (xh(ih+1)-xpc)    /dx
      dx2 = (xpc-      xh(ih))/dx

      dz_ =   zf(kf+1)-zf(kf)
      dz1 = (zf(kf+1)-zpc)/dz_
      dz2 = (zpc-zf(kf))/dz_

      wp = w(ih  ,jh  ,kf  )*dx1*dz1 + &
           w(ih  ,jh  ,kf+1)*dx1*dz2 + &
           w(ih+1,jh  ,kf  )*dx2*dz1 + &
           w(ih+1,jh  ,kf+1)*dx2*dz2 + &

    else

      dx =   xh(ih+1)-xh(ih)
      dx1 = (xh(ih+1)-xpc)    /dx
      dx2 = (xpc-      xh(ih))/dx

      dy =   yh(jh+1)-yh(jh)
      dy1 = (yh(jh+1)-ypc)   /dy
      dy2 = (ypc     -yh(jh))/dy

      dz_ =   zf(kf+1)-zf(kf)
      dz1 = (zf(kf+1)-zpc)/dz_
      dz2 = (zpc-zf(kf))/dz_

      wp = w(ih  ,jh  ,kf  )*dx1*dy1*dz1 + &
           w(ih  ,jh  ,kf+1)*dx1*dy1*dz2 + &
           w(ih  ,jh+1,kf  )*dx1*dy2*dz1 + &
           w(ih  ,jh+1,kf+1)*dx1*dy2*dz2 + &
           w(ih+1,jh  ,kf  )*dx2*dy1*dz1 + &
           w(ih+1,jh  ,kf+1)*dx2*dy1*dz2 + &
           w(ih+1,jh+1,kf  )*dx2*dy2*dz1 + &
           w(ih+1,jh+1,kf+1)*dx2*dy2*dz2
    endif


  end subroutine
  !################################################################!
  !#                SUBROUTINE INTERP_SCALARS                     #!
  !#                                                              #!
  !#  Interpolates from scalar grid to parcel position.           #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) scalar : value on scalar grid       #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    integer ih, jh, kh : indices of lower southwest corner of #!
  !#      box surrounding parcel on scalar grid                   #!
  !#    real,dimension(:) zzh : extended vertical scalar grid     #!
  !#  Returned variable:                                          #!
  !#    real scalarp : value interpolated to parcel position      #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_scalars( scalar, xpc, ypc, zpc, ih, jh, kh, zzh, scalarp )


    implicit none

    real,intent(in),dimension(1:ni,1:nj,0:nk) :: scalar
    real,intent(in) :: xpc, ypc, zpc
    integer,intent(in) :: ih, jh, kh
    real,intent(in),dimension(0:nk) :: zzh
    real,intent(out) :: scalarp

    real :: dx, dx1, dx2, dy, dy1, dy2, dz_, dz1, dz2

    if(nj.eq.1)then

      dx =   xh(ih+1)-xh(ih)
      dx1 = (xh(ih+1)-xpc)    /dx
      dx2 = (xpc-      xh(ih))/dx

      dz_ =   zzh(kh+1)-zzh(kh)
      dz1 = (zzh(kh+1)-zpc)/dz_
      dz2 = (zpc-zzh(kh))/dz_

      scalarp = scalar(ih  ,jh  ,kh  )*dx1*dz1 + &
                scalar(ih  ,jh  ,kh+1)*dx1*dz2 + &
                scalar(ih+1,jh  ,kh  )*dx2*dz1 + &
                scalar(ih+1,jh  ,kh+1)*dx2*dz2 + &

    else

      dx =   xh(ih+1)-xh(ih)
      dx1 = (xh(ih+1)-xpc)    /dx
      dx2 = (xpc-      xh(ih))/dx

      dy =   yh(jh+1)-yh(jh)
      dy1 = (yh(jh+1)-ypc)   /dy
      dy2 = (ypc     -yh(jh))/dy

      dz_ =   zzh(kh+1)-zzh(kh)
      dz1 = (zzh(kh+1)-zpc)/dz_
      dz2 = (zpc-zzh(kh))/dz_

      scalarp = scalar(ih  ,jh  ,kh  )*dx1*dy1*dz1 + &
                scalar(ih  ,jh  ,kh+1)*dx1*dy1*dz2 + &
                scalar(ih  ,jh+1,kh  )*dx1*dy2*dz1 + &
                scalar(ih  ,jh+1,kh+1)*dx1*dy2*dz2 + &
                scalar(ih+1,jh  ,kh  )*dx2*dy1*dz1 + &
                scalar(ih+1,jh  ,kh+1)*dx2*dy1*dz2 + &
                scalar(ih+1,jh+1,kh  )*dx2*dy2*dz1 + &
                scalar(ih+1,jh+1,kh+1)*dx2*dy2*dz2
    endif
  end subroutine

  !################################################################!
  !#                 SUBROUTINE INTERP_DYNZ                       #!
  !#                                                              #!
  !#  Interpolates values of dynamic interest (output from dyn)   #!
  !#  to parcel position.                                         #!
  !#  Passed variables:                                           #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    real,dimension(:,:,:) pgfx_b, pgfy_b, pgfz_b : buoyancy   #!
  !#      PPGFs                                                   #!
  !#    real,dimension(:,:,:) pgfx_dl, pgfy_dl, pgfz_dl : linear  #!
  !#      dynamic PPGFs                                           #!
  !#    real,dimension(:,:,:) pgfx_dn, pgfy_dn, pgfz_dn : non-    #!
  !#      linear dynamic PPGFs                                    #!
  !#  Returned variables:                                         #!
  !#    real pgfxbp, pgfybp, pgfzpg : buoyancy PPGFs at parcel    #!
  !#    real pgfxdlp, pgfydlp, pgfzdlp : linear dynamic PPGFs at  #!
  !#      parcel position                                         #!
  !#    real pgfxdnp, pgfydnp, pgfzdnp : nonlinear dynamic PPGFs  #!
  !#      at parcel position                                      #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_dynz( xpc, ypc, zpc, pgfx_b, pgfxbp, pgfy_b, pgfybp, pgfz_b, pgfzbp, &
    pgfx_dl, pgfxdlp, pgfy_dl, pgfydlp, pgfz_dl, pgfzdlp, &
    pgfx_dn, pgfxdnp, pgfy_dn, pgfydnp, pgfz_dn, pgfzdnp )

    implicit none

    real,intent(inout) :: xpc, ypc, zpc
    real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: pgfx_b, pgfx_dl, pgfx_dn
    real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: pgfy_b, pgfy_dl, pgfy_dn
    real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: pgfz_b, pgfz_dl, pgfz_dn

    real,intent(out) :: pgfxbp, pgfxdlp, pgfxdnp
    real,intent(out) :: pgfybp, pgfydlp, pgfydnp
    real,intent(out) :: pgfzbp, pgfzdlp, pgfzdnp

    integer :: ih, iif, jh, jf, kh, kf
    real,dimension(0:nk) :: zzh
    zzh(0)=-zh(1)
    zzh(1:nk) = zh(1:nk)

    call find_surrounding_grid_points( xpc, ypc, zpc, ih, iif, jh, jf, kh, kf )

      if((iif.ne.9999).and.(jf.ne.9999).and.(kf.ne.9999))then

    call interp_u( pgfx_b,  xpc, ypc, zpc, iif, jh, kh, zzh, pgfxbp )
    call interp_u( pgfx_dl, xpc, ypc, zpc, iif, jh, kh, zzh, pgfxdlp )
    call interp_u( pgfx_dn, xpc, ypc, zpc, iif, jh, kh, zzh, pgfxdnp )

    call interp_v( pgfy_b,  xpc, ypc, zpc, ih, jf, kh, zzh, pgfybp )
    call interp_v( pgfy_dl, xpc, ypc, zpc, ih, jf, kh, zzh, pgfydlp )
    call interp_v( pgfy_dn, xpc, ypc, zpc, ih, jf, kh, zzh, pgfydnp )

    call interp_w( pgfz_b,  xpc, ypc, zpc, ih, jh, kf, pgfzbp )
    call interp_w( pgfz_dl, xpc, ypc, zpc, ih, jh, kf, pgfzdlp )
    call interp_w( pgfz_dn, xpc, ypc, zpc, ih, jh, kf, pgfzdnp )

      else

    pgfxbp=0.
    pgfxdlp=0.
    pgfxdnp=0.

    pgfybp=0.
    pgfydlp=0.
    pgfydnp=0.

    pgfzbp=0.
    pgfzdlp=0.
    pgfzdnp=0.
      endif

  end subroutine

  !################################################################!
  !#               SUBROUTINE RK4                                 #!
  !#                                                              #!
  !#  Uses 4th-order Runge-Kutta method to advance parcel, given  #!
  !#  wind fields.                                                #!
  !#  Passed variables:                                           #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    real up, vp, wp : winds at parcel position                #!
  !#    real,dimension(:,:,:) u, v, w, utem, vtem, wtem : wind    #!
  !#      fields at two different time steps                      #!
  !#    real dt : time step                                       #!
  !#  Returned variables:                                         #!
  !#    real xpc, ypc, zpc : updated parcel position              #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine rk4( xpc, ypc, zpc, up, vp, wp, u, v, w, utem, vtem, wtem, dt )

    implicit none

    real,intent(inout),dimension(2) :: xpc, ypc, zpc
    real,intent(in) :: up, vp, wp
    real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: u, utem
    real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: v, vtem
    real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: w, wtem
    real,intent(in) :: dt

    real,dimension(0:nk) :: zzh

    real x1, x2, x3, x4
    real y1, y2, y3, y4
    real z1, z2, z3, z4

    real k1u, k1ua, k1ub
    real k2u, k2ua, k2ub
    real k3u, k3ua, k3ub
    real k4u, k4ua, k4ub

    real k1v, k1va, k1vb
    real k2v, k2va, k2vb
    real k3v, k3va, k3vb
    real k4v, k4va, k4vb

    real k1w, k1wa, k1wb
    real k2w, k2wa, k2wb
    real k3w, k3wa, k3wb
    real k4w, k4wa, k4wb

        IF((XPC(1).NE.-9.E9).AND.(YPC(1).NE.-9.E9).AND.(ZPC(1).NE.-9.E9).and. &
          & (xpc(2).ne.-9.e9).and.(ypc(2).ne.-9.e9).and.(zpc(2).ne.-9.e9))THEN

    zzh(0)=-zh(1)
    zzh(1:nk) = zh(1:nk)

    if (dt.gt.0) then
      x1=xpc(1)
      y1=ypc(1)
      z1=zpc(1)
    else
      x1=xpc(2)
      y1=ypc(2)
      z1=zpc(2)
    endif

    k1u=up
    k1v=vp
    k1w=wp

    x2 = x1 + 0.5*dt*k1u
    y2 = y1 + 0.5*dt*k1v
    z2 = z1 + 0.5*dt*k1w

    call interp_winds( x2, y2, z2, zzh, u, k2ua, v, k2va, w, k2wa )
    call interp_winds( x2, y2, z2, zzh, utem, k2ub, vtem, k2vb, wtem, k2wb )

    k2u=0.5*(k2ua+k2ub)
    k2v=0.5*(k2va+k2vb)
    k2w=0.5*(k2wa+k2wb)

    x3 = x1 + 0.5*dt*k2u
    y3 = y1 + 0.5*dt*k2v
    z3 = z1 + 0.5*dt*k2w

    call interp_winds( x3, y3, z3, zzh, u, k3ua, v, k3va, w, k3wa )
    call interp_winds( x3, y3, z3, zzh, utem, k3ub, vtem, k3vb, wtem, k3wb )

    k3u=0.5*(k3ua+k3ub)
    k3v=0.5*(k3va+k3vb)
    k3w=0.5*(k3wa+k3wb)

    x4 = x1 + dt*k3u
    y4 = y1 + dt*k3v
    z4 = z1 + dt*k3w

    if(dt.gt.0)then

      call interp_winds( x4, y4, z4, zzh, utem, k4u, vtem, k4v, wtem, k4w )

      xpc(2) = xpc(1) + dt*( k1u + 2*k2u + 2*k3u + k4u )/6
      ypc(2) = ypc(1) + dt*( k1v + 2*k2v + 2*k3v + k4v )/6
      zpc(2) = max( (zpc(1)+dt*(k1w+2*k2w+2*k3w+k4w)/6), 0.0 )

    else

      call interp_winds( x4, y4, z4, zzh, u, k4u, v, k4v, w, k4w )

      xpc(1) = xpc(2) + dt*( k1u + 2*k2u + 2*k3u + k4u )/6
      ypc(1) = ypc(2) + dt*( k1v + 2*k2v + 2*k3v + k4v )/6
      zpc(1) = max( (zpc(2)+dt*(k1w+2*k2w+2*k3w+k4w)/6), 0.0 )

    endif

        ELSE

    xpc=-9.e9
    ypc=-9.e9
    zpc=-9.e9

         ENDIF

  end subroutine
  !################################################################!
  !#                 SUBROUTINE COMPUTE_VORTS                     #!
  !#                                                              #!
  !#  Given buoyancy and wind fields, this subroutine computes    #!
  !#  several quantities relevant to vortex dynamics.  Note that  #!
  !#  the vorticity fields are computed on their "natural" grid,  #!
  !#  based on the derivatives.  E.g., z-vorticity is on a grid   #!
  !#  with dimensions (ni+1,nj+1,nk).                             #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) buoy : buoyancy                     #!
  !#    real,dimension(:,:,:) u, v, w : wind fields               #!
  !#  Returned variables:                                         #!
  !#    real,dimension(:,:,:) xvort, yvort, zvort : vorticity     #!
  !#    real,dimension(:,:,:) dpsi_dx, dpsi_dy, dpsi_dz :         #!
  !#      gradient of psi [tan(v/u)]                              #!
  !#    real,dimension(:,:,:) dbdx, dbdy : horizontal gradient of #!
  !#      buoyancy                                                #!
  !#    real,dimension(:,:,:) dvhdx, dvhdy, dvhdz : gradient of   #!
  !#      wind speed                                              #!
  !#    real,dimension(:,:,:) dwdx, dwdy, dwdz : gradient of w    #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine compute_vorts( buoy, u, v, w, xvort, yvort, zvort, &
    dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
    dwdy, dwdx, dwdz )

    !####################### VARIABLES ############################!
    ! Passed variables
    real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: u
    real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: v
    real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: w
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: buoy

    ! Returned variables
    real,intent(out),dimension(1:ni+1,1:nj,0:nk) :: dpsi_dx, dvhdx, dbdx
    real,intent(out),dimension(1:ni,1:nj+1,0:nk) :: dpsi_dy, dvhdy, dbdy
    real,intent(out),dimension(1:ni,1:nj,1:nk+1) :: dpsi_dz, dvhdz

    real,intent(out),dimension(1:ni,1:nj+1,1:nk+1) :: dwdy, xvort
    real,intent(out),dimension(1:ni+1,1:nj,1:nk+1) :: dwdx, yvort
    real,intent(out),dimension(1:ni+1,1:nj+1,0:nk) :: zvort
    real,intent(out),dimension(1:ni,1:nj,0:nk) :: dwdz

    ! Local variables
    real,dimension(1:ni,1:nj,0:nk) :: uinterp, vinterp, vh

    real,dimension(1:ni,1:nj+1,1:nk+1) :: dvdz
    real,dimension(1:ni+1,1:nj,1:nk+1) :: dudz
    real,dimension(1:ni+1,1:nj+1,0:nk) :: dudy, dvdx

    real,dimension(1:ni,1:nj,0:nk) :: psi

    integer i,j,k
    real,dimension(0:nk) :: zzh

    !##################### MAIN BODY ##############################!

    ! Set up array for extended scalar vertical coordinates
    zzh(1:nk) = zh(1:nk)
    zzh(0) = -zh(1)

    ! Compute shear and vorticity over grid
    write(*,*) 'computing shear and vorticity over grid'

    dvdz(1:ni,1:nj+1,1:nk+1) = 0.0
    do k=1,nk
    do j=1,nj+1
    do i=1,ni
      dvdz(i,j,k) = (v(i,j,k)-v(i,j,k-1))/(zzh(k)-zzh(k-1))
    enddo
    enddo
    enddo

    dwdy(1:ni,1:nj+1,1:nk+1) = 0.0
    do k=1,nk+1
    do j=2,nj
    do i=1,ni
      dwdy(i,j,k) = (w(i,j,k)-w(i,j-1,k))*rdz
    enddo
    enddo
    enddo

    xvort = dwdy-dvdz

    dudz(1:ni+1,1:nj,1:nk+1) = 0.0
    do k=1,nk
    do j=1,nj
    do i=1,ni+1
      dudz(i,j,k)  = (u(i,j,k)-u(i,j,k-1))/(zzh(k)-zzh(k-1))
    enddo
    enddo
    enddo

    dwdx(1:ni+1,1:nj,1:nk+1) = 0.0
    do k=1,nk+1
    do j=1,nj
    do i=2,ni
      dwdx(i,j,k) = (w(i,j,k)-w(i-1,j,k))*rdz
    enddo
    enddo
    enddo

    yvort = dudz - dwdx

    dudy(1:ni+1,1:nj+1,0:nk) = 0.0
    do k=0,nk
    do j=2,nj
    do i=1,ni+1
      dudy(i,j,k) = (u(i,j,k)-u(i,j-1,k))*rdz
    enddo
    enddo
    enddo

    dvdx(1:ni+1,1:nj+1,0:nk) = 0.0
    do k=0,nk
    do j=1,nj
    do i=2,ni
      dvdx(i,j,k) = (v(i,j,k)-v(i-1,j,k))*rdz
    enddo
    enddo
    enddo

    zvort = dvdx - dudy

    ! Compute psi and vh and derivatives
    do k=0,nk ! start by interpolating u and v to scalar grid
    do j=1,nj
    do i=1,ni
      uinterp(i,j,k) = 0.5*( u(i+1,j,k) + u(i,j,k) )
      vinterp(i,j,k) = 0.5*( v(i,j+1,k) + v(i,j,k) )
    enddo
    enddo
    enddo

    do k=0,nk
    do j=1,nj
    do i=1,ni
      call compute_psi( uinterp(i,j,k), vinterp(i,j,k), psi(i,j,k) )
    enddo
    enddo
    enddo

    vh = sqrt( uinterp**2 + vinterp**2 )

    do k=1,nk
    do j=1,nj
    do i=1,ni
      dwdz(i,j,k) = (w(i,j,k+1)-w(i,j,k))/(zf(k+1)-zf(k))
    enddo
    enddo
    enddo
    dwdz(:,:,0) = dwdz(:,:,1)

    dpsi_dx(:,:,:) = 0.0
    do k=0,nk
    do j=1,nj
    do i=2,ni-1
       dpsi_dx(i,j,k) = compute_dpsi2( psi(i-1,j,k), psi(i,j,k), &
         psi(i+1,j,k), 2*dz )
    enddo
    enddo
    enddo

    dpsi_dy(:,:,:) = 0.0
    do k=0,nk
    do j=2,nj-1
    do i=1,ni
      dpsi_dy(i,j,k) = compute_dpsi2( psi(i,j-1,k), psi(i,j,k), &
        psi(i,j+1,k), 2*dz )
    enddo
    enddo
    enddo

    dpsi_dz(:,:,:) = 0.0
    do k=2,nk-1
    do j=1,nj
    do i=1,ni
      dpsi_dz(i,j,k) = compute_dpsi2( psi(i,j,k-1), psi(i,j,k), &
        psi(i,j,k+1), zh(k+1)-zh(k-1) )
    enddo
    enddo
    enddo

    dvhdx(:,:,:) = 0.0
    do k=0,nk
    do j=1,nj
    do i=2,ni
      dvhdx(i,j,k) = (vh(i,j,k)-vh(i-1,j,k))*rdz
    enddo
    enddo
    enddo

    dvhdy(:,:,:) = 0.0
    do k=0,nk
    do j=2,nj
    do i=1,ni
      dvhdy(i,j,k) = (vh(i,j,k)-vh(i,j-1,k))*rdz
    enddo
    enddo
    enddo

    dvhdz(:,:,:) = 0.0
    do k=1,nk
    do j=1,nj
    do i=1,ni
      dvhdz(i,j,k) = (vh(i,j,k)-vh(i,j,k-1))*rdz
    enddo
    enddo
    enddo

    ! Compute horizontal gradient of buoyancy
    dbdx(:,:,:) = 0.0
    do k=0,nk
    do j=1,nj
    do i=2,ni
      dbdx(i,j,k) = (buoy(i,j,k)-buoy(i-1,j,k))*rdz
    enddo
    enddo
    enddo

    dbdy(:,:,:) = 0.0
    do k=0,nk
    do j=2,nj
    do i=1,ni
      dbdy(i,j,k) = (buoy(i,j,k)-buoy(i,j-1,k))*rdz
    enddo
    enddo
    enddo

  end subroutine

  !################################################################!
  !#                SUBROUTINE INTERP_VORT3                       #!
  !#                                                              #!
  !#  Computes vorticity dynamics terms at position of parcel.    #!
  !#  Based on code written by Edwin Adlerman.                    #!
  !#  Passed variables:                                           #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    real up, vp, wp : parcel winds                            #!
  !#    real,dimension(:,:,:) xvort, yvort, zvort : vorticity     #!
  !#    real,dimension(:,:,:) dpsi_dx, dpsi_dy, dpsi_dz :         #!
  !#      gradient of psi [tan(v/u)]                              #!
  !#    real,dimension(:,:,:) dvhdx, dvhdy, dvhdz: gradient of    #!
  !#      wind speed                                              #!
  !#    real,dimension(:,:,:) dwdx, dwdy, dwdz : gradient of w    #!
  !#  Returned variables:                                         #!
  !#    real xvortp, yvortp, zvortp : x-, y-, and z-vorticity     #!
  !#    real svortp, cvortp : streamwise and crosswise vorticity  #!
  !#    real psip: atan(v/u)                                      #!
  !#    real zvort_tiltp : tilting of z-vorticity                 #!
  !#    real zvort_stretchp : stretching of z-vorticity           #!
  !#    real svort_tiltp: tilting of streamwise vorticity         #!
  !#    real svort_stretchp : stretching of streamwise vorticity  #!
  !#    real svort_bclnp : baroclinic generation of streamwise    #!
  !#      vorticity                                               #!
  !#    real cvort_stretchp : stretching of crosswise vorticity   #!
  !#    real cvort_tiltp : tilting of crosswise vorticity         #!
  !#    real cvort_bclnp : baroclinic generation of crosswise     #!
  !#      vorticity                                               #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings 13 Dec 2016                              #!
  !################################################################!
  subroutine interp_vort3( xpc, ypc, zpc, up, vp, wp, xvort, yvort, zvort, &
    dpsi_dx, dvhdx, dbdx, dpsi_dy, dvhdy, dbdy, dpsi_dz, dvhdz, &
    dwdy, dwdx, dwdz, xvortp, yvortp, zvortp, svortp, &
    cvortp, psip, zvort_tiltp, zvort_stretchp, svort_tiltp, svort_stretchp, &
    svort_bclnp, cvort_stretchp, cvort_tiltp, cvort_bclnp )

    real,intent(in),dimension(1:ni+1,1:nj,0:nk) :: dpsi_dx, dvhdx, dbdx
    real,intent(in),dimension(1:ni,1:nj+1,0:nk) :: dpsi_dy, dvhdy, dbdy
    real,intent(in),dimension(1:ni,1:nj,1:nk+1) :: dpsi_dz, dvhdz
    real,intent(in),dimension(1:ni,1:nj+1,1:nk+1) :: dwdy, xvort
    real,intent(in),dimension(1:ni+1,1:nj,1:nk+1) :: dwdx, yvort
    real,intent(in),dimension(1:ni+1,1:nj+1,0:nk) :: zvort
    real,intent(in),dimension(1:ni,1:nj,0:nk) :: dwdz
    real,intent(inout) :: xpc, ypc, zpc, up, vp, wp

    real,intent(out) :: xvortp, yvortp, zvortp, svortp, cvortp, psip
    real,intent(out) :: zvort_tiltp, zvort_stretchp
    real,intent(out) :: svort_stretchp, svort_tiltp, svort_bclnp
    real,intent(out) :: cvort_stretchp, cvort_tiltp, cvort_bclnp

    integer iif, ih, jf, jh, kf, kh, i, j, k
    real,dimension(0:nk) :: zzh

    real sxp, syp, vhp, dwdxp, dwdyp, dwdzp, dpsidxp, dpsidyp, dpsidzp
    real dpsidsp, dpsidnp, dbdxp, dbdyp, dbdsp, dbdnp, dwdsp, dwdnp
    real dvhdxp, dvhdyp, dvhdzp, dvhdnp, dvhdsp

    zzh(1:nk) = zh(1:nk)
    zzh(0) = -zh(1)

    call find_surrounding_grid_points( xpc, ypc, zpc, ih, iif, jh, jf, kh, kf )
    call interp_xvort( xvort, xpc, ypc, zpc, ih, jf, kf, xvortp )
    call interp_yvort( yvort, xpc, ypc, zpc, iif, jh, kf, yvortp )
    call interp_zvort( zvort, xpc, ypc, zpc, iif, jf, kh, zzh, zvortp )

    vhp = sqrt( up**2 + vp**2 )

    sxp = up/vhp
    syp = vp/vhp

    svortp = sxp*xvortp+syp*yvortp
    cvortp = (-syp)*xvortp+sxp*yvortp

    ! vertical vorticity equation

    call interp_yvort( dwdx, xpc, ypc, zpc, iif, jh, kf, dwdxp )
    call interp_xvort( dwdy, xpc, ypc, zpc, ih, jf, kf, dwdyp )

    zvort_tiltp = dwdxp*xvortp + dwdyp*yvortp

    call interp_scalars( dwdz, xpc, ypc, zpc, ih, jh, kh, zzh, dwdzp )
    zvort_stretchp = dwdzp*zvortp

    ! streamwise vorticity equation

    call interp_u( dpsi_dx, xpc, ypc, zpc, iif, jh, kh, zzh, dpsidxp )
    call interp_v( dpsi_dy, xpc, ypc, zpc, ih, jf, kh, zzh, dpsidyp )
    call interp_w( dpsi_dz, xpc, ypc, zpc, ih, jh, kf, dpsidzp )

    call interp_u( dvhdx, xpc, ypc, zpc, iif, jh, kh, zzh, dvhdxp )
    call interp_v( dvhdy, xpc, ypc, zpc, ih, jf, kh, zzh, dvhdyp )
    call interp_w( dvhdz, xpc, ypc, zpc, ih, jh, kf, dvhdzp )

    call interp_u( dbdx, xpc, ypc, zpc, iif, jh, kh, zzh, dbdxp )
    call interp_v( dbdy, xpc, ypc, zpc, ih, jf, kh, zzh, dbdyp )

    svort_tiltp = ( xvortp*dvhdxp + yvortp*dvhdyp + zvortp*dvhdzp )
    svort_stretchp = 0.0
    svort_bclnp = ( dbdyp*sxp - dbdxp*syp )

    ! crosswise vorticity equation
    cvort_tiltp = xvortp*vhp*dpsidxp + yvortp*vhp*dpsidyp + zvortp*vhp*dpsidzp
    cvort_stretchp = 0.0

    cvort_bclnp = -(dbdxp*sxp+dbdyp*syp)

   end subroutine
  !################################################################!
  !#                SUBROUTINE INTERP_XVORT                       #!
  !#                                                              #!
  !#  Interpolates from (ni,nj+1,nk+1) grid to parcel position.   #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) xvort : value on grid               #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    integer ih, jf, kf : indices for lower southwest corner   #!
  !#      of box surrounding parcel on grid                       #!
  !#  Returned variable:                                          #!
  !#    real, xvortp : value interpolated to parcel position      #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_xvort( xvort, xpc, ypc, zpc, ih, jf, kf, xvortp )

    implicit none

    real,intent(in),dimension(1:ni,1:nj+1,1:nk+1) :: xvort
    real,intent(in) :: xpc, ypc, zpc
    integer,intent(in) :: ih, jf, kf
    real,intent(out) :: xvortp

    real :: dx, dx1, dx2, dy, dy1, dy2, dz_, dz1, dz2

      IF(IH.NE.9999)THEN

    dx = xh(ih+1)-xh(ih)
    dx1 = (xh(ih+1)-xpc)/dx
    dx2 = (xpc-xh(ih))/dx

    dy = yf(jf+1)-yf(jf)
    dy1 = (yf(jf+1)-ypc)/dy
    dy2 = (ypc-yf(jf))/dx

    dz_ = zf(kf+1)-zf(kf)
    dz1 = (zf(kf+1)-zpc)/dz_
    dz2 = (zpc-zf(kf))/dz_

    xvortp = xvort(ih  ,jf  ,kf  )*dx1*dy1*dz1 + &
             xvort(ih  ,jf  ,kf+1)*dx1*dy1*dz2 + &
             xvort(ih  ,jf+1,kf  )*dx1*dy2*dz1 + &
             xvort(ih  ,jf+1,kf+1)*dx1*dy2*dz2 + &
             xvort(ih+1,jf  ,kf  )*dx2*dy1*dz1 + &
             xvort(ih+1,jf  ,kf+1)*dx2*dy1*dz2 + &
             xvort(ih+1,jf+1,kf  )*dx2*dy2*dz1 + &
             xvort(ih+1,jf+1,kf+1)*dx2*dy2*dz2

      ELSE

   xvortp = 9.e9

      ENDIF

  END SUBROUTINE

  !################################################################!
  !#                     SUBROUTINE INTERP_YVORT                  #!
  !#                                                              #!
  !#  Interpolates from a (ni+1,nj,nk+1) grid to parcel position. #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) yvort : value on grid               #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    integer iif, jh, kf : indices for lower southwest corner  #!
  !#      of box on grid surrounding parcel position              #!
  !#  Returned variable:                                          #!
  !#    real yvortp : value interpolated to parcel position       #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine interp_yvort( yvort, xpc, ypc, zpc, iif, jh, kf, yvortp )

    implicit none

    real,intent(in),dimension(1:ni+1,1:nj,1:nk+1) :: yvort
    real,intent(in) :: xpc, ypc, zpc
    integer,intent(in) :: iif, jh, kf
    real,intent(out) :: yvortp

    real :: dx, dx1, dx2, dy, dy1, dy2, dz_, dz1, dz2


    dx = xf(iif+1)-xf(iif)
    dx1 = (xf(iif+1)-xpc)/dx
    dx2 = (xpc-xf(iif))/dx

    dy = yh(jh+1)-yh(jh)
    dy1 = (yh(jh+1)-ypc)/dy
    dy2 = (ypc-yh(jh))/dx

    dz_ = zf(kf+1)-zf(kf)
    dz1 = (zf(kf+1)-zpc)/dz_
    dz2 = (zpc-zf(kf))/dz_

    yvortp = yvort(iif  ,jh  ,kf  )*dx1*dy1*dz1 + &
             yvort(iif  ,jh  ,kf+1)*dx1*dy1*dz2 + &
             yvort(iif  ,jh+1,kf  )*dx1*dy2*dz1 + &
             yvort(iif  ,jh+1,kf+1)*dx1*dy2*dz2 + &
             yvort(iif+1,jh  ,kf  )*dx2*dy1*dz1 + &
             yvort(iif+1,jh  ,kf+1)*dx2*dy1*dz2 + &
             yvort(iif+1,jh+1,kf  )*dx2*dy2*dz1 + &
             yvort(iif+1,jh+1,kf+1)*dx2*dy2*dz2


  END SUBROUTINE

  !################################################################!
  !#                  SUBROUTINE INTERP_ZVORT                     #!
  !#                                                              #!
  !#  Interpolates from (ni+1,nj+1,nk) grid to parcel position.   #!
  !#  Passed variables:                                           #!
  !#    real,dimension(:,:,:) zvort : value on grid               #!
  !#    real xpc, ypc, zpc : parcel position                      #!
  !#    integer iif, jf, kh : indices for lower southwest corner  #!
  !#      of box on grid surrounding parcel                       #!
  !#    real,dimension(:) zzh : extended scalar vertical points   #!
  !#  Returned variable:                                          #!
  !#    real zvortp : value interpolated to parcel position       #!
  !#==============================================================#!
  !# v1.0 Ryan Hastings, 13 Dec 2016                              #!
  !################################################################!
  subroutine interp_zvort( zvort, xpc, ypc, zpc, iif, jf, kh, zzh, zvortp )

    implicit none

    real,intent(in),dimension(1:ni+1,1:nj+1,0:nk) :: zvort
    real,intent(inout) :: xpc, ypc, zpc
    real,intent(in),dimension(0:nk) :: zzh
    integer,intent(in) :: iif, jf, kh
    real,intent(out) :: zvortp
    integer ih

    real :: dx, dx1, dx2, dy, dy1, dy2, dz_, dz1, dz2

      IF(IH.NE.9999)THEN

    dx = xf(iif+1)-xf(iif)
    dx1 = (xf(iif+1)-xpc)/dx
    dx2 = (xpc-xf(iif))/dx

    dy = yf(jf+1)-yf(jf)
    dy1 = (yf(jf+1)-ypc)/dy
    dy2 = (ypc-yf(jf))/dx

    dz_ = zzh(kh+1)-zzh(kh)
    dz1 = (zzh(kh+1)-zpc)/dz_
    dz2 = (zpc-zzh(kh))/dz_

    zvortp = zvort(iif  ,jf  ,kh  )*dx1*dy1*dz1 + &
             zvort(iif  ,jf  ,kh+1)*dx1*dy1*dz2 + &
             zvort(iif  ,jf+1,kh  )*dx1*dy2*dz1 + &
             zvort(iif  ,jf+1,kh+1)*dx1*dy2*dz2 + &
             zvort(iif+1,jf  ,kh  )*dx2*dy1*dz1 + &
             zvort(iif+1,jf  ,kh+1)*dx2*dy1*dz2 + &
             zvort(iif+1,jf+1,kh  )*dx2*dy2*dz1 + &
             zvort(iif+1,jf+1,kh+1)*dx2*dy2*dz2

      ELSE

   zvortp = 9.e9

      ENDIF

  END SUBROUTINE

  !################################################################!
  !#                 REAL FUNCTION COMPUTE_DPSI2                  #!
  !#                                                              #!
  !#  Compute difference between angles over distance.  Based on  #!
  !#  code from Edwin Adlerman.                                   #!
  !#  Passed variables:                                           #!
  !#    real psi1, psi, psi2: angles                              #!
  !#    real dx : distance                                        #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  real function compute_dpsi2( psi1, psi, psi2, dx )
    real :: psi1, psi, psi2, dx
    real :: dpsi

    dpsi = psi2-psi1
    if(psi2.gt.psi1)then
      if((psi.le.psi2).and.(psi.ge.psi1))then
      else
        dpsi=psi2-2*pi+psi1
      endif
    elseif(psi2.lt.psi1)then
      if((psi.le.psi1).and.(psi.ge.psi2))then
      else
        dpsi = psi2+2*pi-psi1
      endif
    endif
    if(abs(dpsi).gt.pi*1.25)then
      if(dpsi.gt.0.0)then
        dpsi=dpsi-2*pi
      else
        dpsi=dpsi+2*pi
      endif
    endif

    compute_dpsi2 = dpsi/dx

  end function

  !################################################################!
  !#                 SUBROUTINE COMPUTE_PSI                       #!
  !#                                                              #!
  !#  Computes atan(v/u)                                          #!
  !#  Passed variables:                                           #!
  !#    real up, vp : parcel winds                                #!
  !#  Returned variable:                                          #!
  !#    real psip : angle at parcel position                      #!
  !#==============================================================#!
  !# v1.0, Ryan Hastings, 13 Dec 2016                             #!
  !################################################################!
  subroutine compute_psi(up,vp,psip)

    real,intent(in) :: up
    real,intent(in) :: vp
    real,intent(out) :: psip

    if((up.ge.0).and.(vp.eq.0))then
      psip=0.0
    elseif((up.gt.0).and.(vp.gt.0))then
      psip=atan(vp/up)
    elseif((up.eq.0).and.(vp.gt.0))then
      psip=0.5*pi
    elseif((up.lt.0).and.(vp.gt.0))then
      psip=pi+atan(vp/up)
    elseif((up.lt.0).and.(vp.eq.0))then
      psip=pi
    elseif((up.lt.0).and.(vp.lt.0))then
      psip=pi+atan(vp/up)
    elseif((up.eq.0).and.(vp.lt.0))then
      psip=1.5*pi
    elseif((up.gt.0).and.(vp.lt.0))then
      psip=2*pi+atan(vp/up)
    endif

  end subroutine



end module
