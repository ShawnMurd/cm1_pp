module globals
  !#############################################################################!
  !#############################################################################!
  !####                                                                     ####!
  !####                            MODULE GLOBALS                           ####!
  !####                                                                     ####!
  !####  Module GLOBALS contains variables to be used by all programs and   ####!
  !####  subroutines in the CM1_PP postprocessing package for CM1.  This is ####!
  !####  used in lieu of a include file, mostly because I read on some web- ####!
  !####  site that it was better because modules are more intrinsic to F90  ####!
  !####  than the 'include' construction.  It's all magic to me, like fire  ####!
  !####  or the light of the daystar.                                       ####!
  !####=====================================================================####!
  !####  Ryan Hastings, 12 Sep 2011                                         ####!
  !####---------------------------------------------------------------------####!
  !#### v1.0, eliminated Poisson solver variables (they were for an itera-  ####!
  !####   tive solver which is no longer used) and added comments           ####!
  !####   14 Dec 2016                                                       ####!
  !#############################################################################!
  !#############################################################################!
  !-----------------------------------------------------------------------------!
  ! spatial parameters

  integer,public :: ni, nip1, nj, njp1, nk, nkp1, testq, ntimes

  real,public :: dz, rdz, dt
  real,public :: time

  real,dimension(:),allocatable,public :: xh, xf, xp, yh, yf, yp, zh, zf, zp
  real,dimension(:),allocatable,public :: mzh, mzf

  !--------------------------------------------------------------------!
  ! physical constants..all units are SI units

  real,public,parameter :: g=9.80665 ! gravitational constant
  real,public,parameter :: pi=3.14159265 ! pi..duh
  real,public,parameter :: rd = 287.04 ! gas constant for dry air
  real,public,parameter :: rv = 461.5  ! gas constant for water vapor
  real,public,parameter :: eps = rd/rv
  real,public,parameter :: reps = rv/rd
  real,public,parameter :: cp = 1005.7 ! heat capacity for dry air
  real,public,parameter :: rcp = 1.0/cp
  real,public,parameter :: kappa = rd/cp
  real,public,parameter :: rkappa = cp/rd
  real,public,parameter :: cpl = 4190.0 ! heat capacity for liquid water
  real,public,parameter :: cpv = 1870.0 ! heat capacity for water vapor
  real,public,parameter :: t00 = 273.15 ! freezing point of water
  real,public,parameter :: xlv = 2501000.0 ! latent heat of water vapor
  real,public,parameter :: lv1 = xlv+(cpl-cpv)*t00 ! for computing theta-e
  real,public,parameter :: lv2 = cpl-cpv ! for computing theta-e
  real,public,parameter :: p00 = 1.e5 ! standard pressure

  !--------------------------------------------------------------------!
  ! outputs
  integer,public :: output_dl, output_dn, output_b ! for dyn.input
  integer,public :: output_dyn, output_vor, output_z_only, noncq ! for traj.input

  integer,public :: missing_flag
  real,public,parameter :: missing_value = -9e9



  integer,public :: ptype ! microphysics type

end module
