! $Id: jv_cmn.h,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
!
!----jv_cmn.h---COMMON BLOCKS for new FAST-J code (wild/prather 7/99)
!
!  Parameters
!  ----------
!
!     NB    Number of levels in CTM plus one for above model top
!     NC    Number of levels in the fundamental Fast-J grid
!     NS    Maximum number of species which require J-values calculating
!     NW    Maximum number of wavelength bins that can be used
!     NP    Maximum number of aerosol/cloud types that can be used
!     MX    Number of aerosol/cloud types supplied from CTM
!     NDUST Number of mineral dust categories
!
!                                        Note: THETA(NL) no longer used
!
!  NOTES for CTM Interface (bmy, 10/27/99, 3/23/03)
!  =====================================================================
!  (1) Change JPNL and JPPJ from parameters to variables, which are 
!      set in "inphot.f".  This allows the user to switch the number 
!      of levels at run-time via the CTM inputs. 
!
!  (2) Now make RAD, ZZHT, DTAUMAX, DTAUSUB, DSUBDIV, SZAMAX into
!      parameters instead of holding them in common blocks.  
!
!  (3) Create new common blocks /WLLOC/ and /JVLOC/ to hold certain
!      quantities -Xlocal for parallel code (ppm, 4/98, bmy, 9/21/99)
!
!  (4) The common blocks that must be held -Xlocal are:
!         /ATMOS/, /JVSUB/, /WLLOC/, /JVLOC/ 
! 
!  (4a) Declare the above commons THREADPRIVATE for the Compaq
!       Alpha platform (bmy, 7/10/01)
!
!  (5) Break MIEDX off from the WLLOC common block, since it must
!      not be declared LOCAL for the parallelization. (bmy, 5/2/00)
!  
!  (6) For including aerosol optical depths: (rvm, bmy, 9/30/00)
!      (a) Increase MX from 3 to 10 .  
!      (c) Add ODMDUST(IPAR,JPAR,LPAR,NDUST) to common block /CLIM/
! 
!  (7) Move NDUST to CMN_SIZE to avoid conflicts (bmy, 11/15/01)
!  
!  (8) For updating aerosol optical depths again (rvm, bmy, 2/27/02):
!      (a) Change NP from 21 to 56
!      (b) Change MX from 10 to 35
!      (c) Add ODAER(IPAR,JPAR,LPAR,NAER*NRH) to common block /CLIM/
!  
!  (9) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!       compilation.  Also added & continuation characters in column 73
!       to allow header files to be included in F90 freeform files.
!       Also changed comment character from "C" to "!" to allow this
!       file to be inlined into freeform source code. (bmy, 6/25/02)
!
!  (10) Renamed cpp switch from DEC_COMPAQ to COMPAQ.  Also declare common
!        blocks ATMOS, JVLOC, WLLOC, JVSUB as !$OMP THREADPRIVATE for
!        all platforms. (bmy, 3/23/03)
!  (11) Added new pressure denpendencies algorithm parameters 
!         for MGLY. (tmf, 1/7/09)
!  (12) Added 'pdepf' as pressure dependancy function selector. (tmf, 1/31/06)
!-----------------------------------------------------------------------------
      INTEGER      NB, NC, NS, NW, NP, MX
      PARAMETER   (NB=LPAR+1, NC=2*NB, NS=51, NW=15, NP=56, MX=35)
      CHARACTER*20 TITLEA(NP)
      CHARACTER*78 TITLE0
      CHARACTER*7  TITLEJ(3,NS), jlabel(JPMAX)
      INTEGER jind(JPMAX),jadsub(nc)
      INTEGER NJVAL,NW1,NW2,MIEDX,NAA,NLBATM,npdep,jpdep(NS)
      TYPE (XPLEX) TJ,PJ,DM,DO3,Z,AER,AMF,RAD,RFLECT,SZA,U0,TANHT,ZZHT
      TYPE (XPLEX) WBIN,WL,FL,QO2,QO3,Q1D,QQQ,QRAYL,TQQ,FFF,VALJ,WAA,
     &           QAA,
     &           PAA
      REAL*8 D_WBIN,D_WL,D_FL,D_QRAYL,D_QBC,D_Q1D,D_QQQ,D_QO2,D_QO3,
     & D_TQQ,D_zpdep(NW,7),D_QAA,D_RAA,D_SSA,D_PAA,D_WAA,
     & D_mglypdep(NW,3),D_TREF,D_OREF
      TYPE (XPLEX) RAA,SSA,TREF,OREF,BREF,QBC,DBC,zpdep(NW,7)
      TYPE (XPLEX) dtaumax,szamax,zj(LPAR,JPMAX),jfacta(JPMAX)
      TYPE (XPLEX) dtausub,dsubdiv
      TYPE (XPLEX) ODMDUST,ODAER
      INTEGER PDEPF(7)
      TYPE (XPLEX) MGLYPDEP(NW, 3)

!-----------------------------------------------------------------------
! These common blocks MUST NOT be held local (bmy, 5/2/00)
      COMMON /TITLS/TITLE0,TITLEJ,TITLEA
      COMMON /CCWVL/WBIN(NW+1),WL(NW),FL(NW),QO2(NW,3),QO3(NW,3),       &
     &              Q1D(NW,3),QQQ(NW,2,NS-3),QRAYL(NW),TQQ(3,NS),       &
     &              WAA(4,NP),QAA(4,NP),                                & 
     &              PAA(8,4,NP),RAA(4,NP),SSA(4,NP),QBC(NW),            &
     &              D_WBIN(NW+1),D_WL(NW),D_FL(NW),D_QO2(NW,3),         &
     &              D_QO3(NW,3),D_Q1D(NW,3),D_QQQ(NW,2,NS-3),           &
     &              D_QRAYL(NW),D_TQQ(3,NS),D_WAA(4,NP),D_QAA(4,NP),    &
     &              D_PAA(8,4,NP),D_RAA(4,NP),D_SSA(4,NP),D_QBC(NW),    &
     &              NJVAL,NW1,NW2,NAA,NLBATM
      COMMON /CLIM/ TREF(51,18,12),OREF(51,18,12),BREF(51),             &
     &              D_TREF(51,18,12),D_OREF(51,18,12),                  &
     &              ODMDUST(IPAR,JPAR,LPAR,NDUST),                      &
     &              ODAER(IPAR,JPAR,LPAR,NAER*NRH)

      COMMON /JVALS/jfacta,zpdep,npdep,jpdep,jind,jlabel,               &
     &              pdepf,mglypdep,D_zpdep,D_mglypdep

      COMMON /JVIDX/MIEDX(MX)
!-----------------------------------------------------------------------
! These common blocks MUST be held local for the parallelization (bmy, 5/2/00)
      COMMON /ATMOS/TJ(NB),PJ(NB+1),DM(NB),DO3(NB),DBC(NB),Z(NB),       &
     &              AER(MX,NB),AMF(NB,NB),RFLECT,SZA,U0,TANHT
      COMMON /JVLOC/zj
      COMMON /WLLOC/FFF(NW,lpar),VALJ(NS)
      COMMON /JVSUB/jadsub

      !=================================================================
      ! Declare the following common blocks as THREADPRIVATE for the
      ! OpenMP parallelization on all platforms (bmy, 3/23/03)
      !=================================================================
!$OMP THREADPRIVATE( /ATMOS/ )
!$OMP THREADPRIVATE( /JVLOC/ )
!$OMP THREADPRIVATE( /WLLOC/ )
!$OMP THREADPRIVATE( /JVSUB/ )
!-----------------------------------------------------------------------
! Parameters for FAST-J 
      PARAMETER ( RAD     = xplex(6375.d5,0d0) )
      PARAMETER ( ZZHT    = xplex(5.d5,0d0)    )
      PARAMETER ( dtaumax = xplex(1.d0,0d0)    )
      PARAMETER ( dtausub = xplex(1.d0,0d0)    )
      PARAMETER ( dsubdiv = xplex(10.d0,0d0)   )
      PARAMETER ( szamax  = xplex(98.0d0,0d0)  )
!-----------------------------------------------------------------------
