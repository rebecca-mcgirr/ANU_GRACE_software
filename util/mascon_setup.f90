    module mascon_setup

! mascon_setup
!
! main arrays and variables used to build and test mascon files:
! --------------------------------------------------------------
!
! pm/sm/tm - primary/secondary/ternary mascon prefix
! (for pm below, read also sm...,tm...)
! pmlat - primary mascon centre latitude
! pmlon - primary mascon centre longitude
! pmr - primary mascon geocentric radius
! pma - primary mascon area
! pmh - primary mascon altitude
! pmg - primary mascon geoid height
! pmden - primary mascon water density (salt/fresh)
! pcland - percent land in primary mascon
! np,ns,nt - number of p/s/t mascons
! npex,nsex,ntex - expected number of p/s/t mascons
! sap,tap,tas - secondary associated primary etc.
! nsip,ntip,ntis - number of secondaries in a primary, etc.
! maxsip,maxtip,maxtis - maximum # of secondaries in a primary mascon, etc.
! sscol,tpcol,tscol - secondary secondary colour, etc.
! (p/s/t)region - primary region, etc.
! (p/s/t)type - primary type, etc.

    real (kind=8), parameter    :: pi = 3.141592653589793d0
    integer, parameter      :: MHR = 1000

    real (kind=8), dimension(:), allocatable    :: pmlat,pmlon,pmr,pma,pmh,pmg,pmden,pcland
    real (kind=8), dimension(:), allocatable    :: smlat,smlon,smr,sma,smh,smg,smden
    real (kind=8), dimension(:), allocatable    :: tmlat,tmlon,tmr,tma,tmh,tmg,tmden
    real (kind=8), dimension(:), allocatable    :: dspmax,dtpmax,dtsmax
    integer, dimension(:), allocatable          :: nsip,ntip,tidemask
    integer, dimension(:), allocatable          :: ntis,sap,sscol
    integer, dimension(:), allocatable          :: tap,tas,tpcol,tscol
    integer, dimension(:), allocatable          :: isip,itip,kips,kpsc
    integer, dimension(:), allocatable          :: itis,kist,kstc
    integer, dimension(:), allocatable          :: hsind,hpind
    integer, dimension(:), allocatable          :: scheck,tcheck

    character*15, dimension(:), allocatable	   :: pregion
    character*15, dimension(:), allocatable	   :: sregion
    character*15, dimension(:), allocatable	   :: tregion
    character(6), dimension(:), allocatable    :: ptype
    character(6), dimension(:), allocatable    :: stype
    character(6), dimension(:), allocatable    :: ttype

    character(150), dimension(1:MHR)  :: headrec

    end module mascon_setup
