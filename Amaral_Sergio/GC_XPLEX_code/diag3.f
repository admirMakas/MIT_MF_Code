! $Id: diag3.f,v 1.3 2012/03/01 22:00:26 daven Exp $
      SUBROUTINE DIAG3                                                      
! 
!******************************************************************************
!  Subroutine DIAG3 prints out diagnostics to the BINARY format punch file 
!  (bmy, bey, mgs, rvm, 5/27/99, 12/15/08)
!
!  NOTES: 
!  (40) Bug fix: Save levels 1:LD13 for ND13 diagnostic for diagnostic
!        categories "SO2-AC-$" and "SO2-EV-$".  Now reference F90 module
!        "tracerid_mod.f".  Now reference NUMDEP from "drydep_mod.f".
!        Now save anthro, biofuel, biomass NH3 in ND13; also fixed ND13
!        tracer numbers.  For ND13, change scale factor from SCALESRCE to 1.
!        Now references "wetscav_mod.f".  Now also save true tracer numbers 
!        for ND38 and ND39 diagnostic.  Now also write out biomass SO2.
!        Now convert ND01, ND02, ND44 diagnostics for Rn/Pb/Be from kg to 
!        kg/s here. (bmy, 1/24/03)
!  (41) Now save out natural NH3 in ND13 as "NH3-NATU" (rjp, bmy, 3/23/03)
!  (42) Now replace DXYP(JREF) by routine GET_AREA_M2, GET_XOFFSET, and
!        GET_YOFFSET of "grid_mod.f".  Now references "time_mod.f".
!        DIAGb, DIAGe are now local variables.  Now remove obsolete statements
!        IF ( LBPNCH > 0 ).  Removed SCALE1, replaced with SCALEDYN. 
!        (bmy, 2/24/03)
!  (43) Added TSKIN, PARDF, PARDR, GWET to ND67 diagnostic.  For GEOS-4/fvDAS,
!        UWND, VWND, TMPU, SPHU are A-6 fields.  Adjust the ND66 scale factors 
!        accordingly.  Delete KZZ from ND66.  Updated comments. (bmy, 6/23/03)
!  (44) Bug fix: use LD68 instead of ND68 in DO-loop to avoid out-of-bounds 
!        error. (bec, bmy, 7/15/03)
!  (45) Now print out NTRACE drydep fluxes for tagged Ox.  Also tagged Ox 
!        now saves drydep in molec/cm2/s.  Now print out Kr85 prod/loss in 
!        ND03. (bmy, 8/20/03)
!  (46) Now use actual tracer number for ND37 diagnostic. (bmy, 1/21/04)
!  (47) Now loop over the actual # of soluble tracers for ND17, ND18.  
!        (bmy, 3/19/04)
!  (48) Now use the actual tracer # for ND17 and ND18 diagnostics. 
!        Rearrange ND44 code for clarity. (bmy, 3/23/04)
!  (49) Added ND06 (dust aerosol) and ND07 (carbon aerosol) diagnostics.
!        Now scale online dust optical depths by SCALECHEM in ND21 diagnostic.
!        (rjp, tdf, bmy, 4/5/04)
!  (50) Added ND08 (seasalt aerosol) diagnostic (rjp, bec, bmy, 4/20/04)
!  (51) Now save out SO2 from ships (if LSHIPSO2=T) (bec, bmy, 5/20/04)
!  (52) Added NVOC source diagnostics for ND07 (rjp, bmy, 7/13/04)
!  (53) Now reference "logical_mod.f", "tracer_mod.f", and "diag_pl_mod.f".
!        Bug fix in write to DMS_BIOG. (bmy, 7/20/04)
!  (54) Comment out ND27 for GEOS-4.  It isn't working 100% right.  If you
!        examine the flux at 200 hPa, you get the same info. (bmy, 10/15/04)
!  (55) Added biofuel SO4 to the bpch file under ND13.  Bug fix: replace ND68 
!        with LD68 in call to BPCH2 (auvray, bmy, 11/17/04)
!  (56) Now save ND03 mercury diagnostic arrays to bpch file.  Also updated
!        ND44 for tagged Hg tracers (eck, bmy, 12/14/04)
!  (57) Now print out extra ND21 diagnostics for crystalline sulfur tracers.  
!        Also now save total oceanic mass of Hg0 and Hg2.  Now call 
!        WRITE_DIAG03 from "diag03_mod.f" (bmy, 1/21/05)
!  (58) Now call WRITE_DIAG41 from "diag41_mod.f" (bmy, 2/17/05)
!  (59) Add P(SO4s) to row 8 of ND05 diagnostic.  Also remove special tracer
!        numbers for the ND67 diagnostic.  Now do not save CLDMAS for ND67
!        for GEOS-4, since GEOS-4 convection uses different met fields.
!        (bec, bmy, 5/3/05)
!  (60) Bug fix in ND68 diagnostic: use LD68 instead of ND68 in call to BPCH2.
!        Now modified for GEOS-5 and GCAP met fields.  Remove references to
!        CO-OH param simulation.  Also remove references to TRCOFFSET since
!        that is always zero now.  Now call GET_HALFPOLAR from "bpch2_mod.f" 
!        to get the HALFPOLAR value for GEOS or GCAP grids. (swu, bmy, 6/24/05)
!  (61) References ND04, WRITE_DIAG04 from "diag04_mod.f".  Also now updated
!        ND30 diagnostic for land/water/ice flags.  Also remove reference
!        to LWI array. (bmy, 8/18/05)
!  (62) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (63) Added MBO as tracer #5 in ND46 diagnostic (tmf, bmy, 10/20/05)
!  (64) Removed duplicate variable declarations.  Now remove restriction on 
!        printing out cloud mass flux in GEOS-4 for the ND66 diagnostic. 
!        (bmy, 3/14/06)
!  (65) References ND56, WRITE_DIAG56 from "diag56_mod.f" (ltm, bmy, 5/5/06)
!  (66) Now remove TRCOFFSET; it's obsolete.  References ND42, WRITE_DIAG42 
!        from "diag42_mod.f" (dkh, bmy, 5/22/06)
!  (67) Updated ND36 diagnostic for CH3I (bmy, 7/25/06)
!  (68) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (69) Replace TINY(1d0) with 1d-32 to avoid problems on SUN 4100 platform
!        (bmy, 9/5/06)
!  (70) Now write diag 54 (time in the troposphere) if asked for (phs, 9/22/06)
!  (71) Now use new time counters for ND43 & ND45,  Also now average between
!        0 and 24 UT for ND47.  Bug fix in ND36. (phs, bmy, 3/5/07)
!  (72) Bug fix in ND65: use 3-D counter array (phs, bmy, 3/6/07)
!  (73) Bug fix in ND07: now save out IDTSOA4 tracer.  Modifications for H2/HD
!        diagnostics (ND10, ND27, ND44) (tmf, phs, bmy, 9/18/07)
!  (74) Now save out true pressure at 3-D level edges for ND31.  Change ND31
!        diagnostic category name to "PEDGE-$". Bug fix in ND28 diagnostic to 
!        allow you to print out individual biomass tracers w/o having to print 
!        all of them. (bmy, dkh, 1/24/08)
!  (75) Bug fix: Now divide ALBEDO in ND67 by SCALE_I6 for GEOS-3 met, but
!        by SCALE_A3 for all other met types (phs, bmy, 10/7/08)
!  (76) Fix ND65, ND47, and ozone case in ND45. Now only ND45 depends
!        on LD45 (phs, 11/17/08)
!  (77) Bug fix: Select the right index of AD34 to write.  Pick the right 
!         tracer field from AD22 if only a subset of tracers are requested 
!         to be printed out. (ccc, 12/15/08)
!  (78) Added ND52 for gamma(HO2) (jaegle, 02/26/09)
!  (79) Updated test on ship emissions flag for AD13 (phs, 3/3/09)     
!  (80) Add AD07_SOAGM for dicarbonyl SOA formation (tmf, 3/6/09)
!  (81) Add output in AD22 for dicarbonyl photolysis J values (tmf, 3/6/09)
!  (82) Add output in AD46 for biogenic C2H4 emissions (tmf, 3/6/09)
!  (87) Add diagnostics 19, 58 and 60 for methane. (kjw, 8/18/09, adj32_023)
!******************************************************************************
! 
      ! References to F90 modules
      USE BPCH2_MOD
      USE BIOMASS_MOD,  ONLY : BIOTRCE,     NBIOMAX
      USE BIOFUEL_MOD,  ONLY : NBFTRACE,    BFTRACE
      USE DIAG_MOD,     ONLY : AD01,        AD02,        AD05    
      USE DIAG_MOD,     ONLY : AD06,        AD07,        AD07_BC
      USE DIAG_MOD,     ONLY : AD07_SOAGM
      USE DIAG_MOD,     ONLY : AD07_OC,     AD07_HC,     AD08
      USE DIAG_MOD,     ONLY : AD09,        AD09_em,     AD11
      USE DIAG_MOD,     ONLY : AD12,        AD13_DMS,    AD13_SO2_ac 
      USE DIAG_MOD,     ONLY : AD13_SO2_an, AD13_SO2_bb, AD13_SO2_bf
      USE DIAG_MOD,     ONLY : AD13_SO2_ev, AD13_SO2_nv, AD13_SO4_an
      USE DIAG_MOD,     ONLY : AD13_SO4_bf, AD13_SO2_sh, AD13_NH3_an
      USE DIAG_MOD,     ONLY : AD13_NH3_na, AD13_NH3_bb, AD13_NH3_bf
      USE DIAG_MOD,     ONLY : CONVFLUP,    TURBFLUP,    AD16
      USE DIAG_MOD,     ONLY : CT16,        AD17,        CT17
      USE DIAG_MOD,     ONLY : AD18,        CT18,        AD21
      USE DIAG_MOD,     ONLY : AD21_cr,     AD22,        LTJV
      USE DIAG_MOD,     ONLY : CTJV,        MASSFLEW,    MASSFLNS
      USE DIAG_MOD,     ONLY : MASSFLUP,    AD28,        AD29
      USE DIAG_MOD,     ONLY : AD30,        AD31
      USE DIAG_MOD,     ONLY : AD32_ac,     AD32_an,     AD32_bb
      USE DIAG_MOD,     ONLY : AD32_bf,     AD32_fe,     AD32_li
      USE DIAG_MOD,     ONLY : AD32_so,     AD32_ub,     AD33
      USE DIAG_MOD,     ONLY : AD34,        AD35,        AD36
      USE DIAG_MOD,     ONLY : AD37,        AD38,        AD39
      USE DIAG_MOD,     ONLY : AD43,        LTNO
      USE DIAG_MOD,     ONLY : CTNO,        LTOH,        CTOH
      USE DIAG_MOD,     ONLY : LTHO2,       CTHO2,       LTNO2
      USE DIAG_MOD,     ONLY : CTNO2,       LTNO3,       CTNO3
      USE DIAG_MOD,     ONLY : AD44,        AD45,        LTOTH
      USE DIAG_MOD,     ONLY : CTOTH,       AD46,        AD47
      USE DIAG_MOD,     ONLY : AD52
      USE DIAG_MOD,     ONLY : AD54,        CTO3,        CTO3_24h
      USE DIAG_MOD,     ONLY : AD19,        AD58,        AD60
      USE DIAG_MOD,     ONLY : AD55,        AD66,        AD67
      USE DIAG_MOD,     ONLY : AD68,        AD69
      USE DIAG_MOD,     ONLY : AD10,        AD10em
      USE DIAG03_MOD,   ONLY : ND03,        WRITE_DIAG03
      USE DIAG04_MOD,   ONLY : ND04,        WRITE_DIAG04
      USE DIAG41_MOD,   ONLY : ND41,        WRITE_DIAG41
      USE DIAG42_MOD,   ONLY : ND42,        WRITE_DIAG42
      USE DIAG56_MOD,   ONLY : ND56,        WRITE_DIAG56
!     diag59 added, (lz,10/11/10) 
      USE DIAG59_MOD,   ONLY : ND59,        WRITE_DIAG59
      USE DIAG_PL_MOD,  ONLY : AD65
      USE DRYDEP_MOD,   ONLY : NUMDEP,      NTRAIND
      USE FILE_MOD,     ONLY : IU_BPCH
      USE GRID_MOD,     ONLY : GET_AREA_M2, GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,  ONLY : LCARB,       LCRYST,      LDUST    
      USE LOGICAL_MOD,  ONLY : LSHIPSO2,    LSOA,        LSSALT
      USE LOGICAL_MOD,  ONLY : LEDGARSHIP,  LARCSHIP,    LEMEPSHIP
      USE TIME_MOD,     ONLY : GET_DIAGb,   GET_DIAGe,   GET_CT_A3   
      USE TIME_MOD,     ONLY : GET_CT_A6,   GET_CT_CHEM, GET_CT_CONV 
      USE TIME_MOD,     ONLY : GET_CT_DYN,  GET_CT_EMIS, GET_CT_I6   
      USE TRACER_MOD,   ONLY : N_TRACERS,   STT,         TRACER_MW_G
      USE TRACER_MOD,   ONLY : TRACER_NAME
      USE TRACER_MOD,   ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,   ONLY : ITS_A_CH3I_SIM
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,   ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,   ONLY : ITS_A_MERCURY_SIM
      USE TRACER_MOD,   ONLY : ITS_A_RnPbBe_SIM
      USE TRACER_MOD,   ONLY : ITS_A_TAGOX_SIM
      USE TRACERID_MOD, ONLY : IDTPB,       IDTDST1,     IDTDST2 
      USE TRACERID_MOD, ONLY : IDTDST3,     IDTDST4,     IDTBCPI 
      USE TRACERID_MOD, ONLY : IDTOCPI,     IDTALPH,     IDTLIMO 
      USE TRACERID_MOD, ONLY : IDTSOA1,     IDTSOA2,     IDTSOA3 
      USE TRACERID_MOD, ONLY : IDTSALA,     IDTSALC,     IDTDMS 
      USE TRACERID_MOD, ONLY : IDTSO2,      IDTSO4,      IDTNH3 
      USE TRACERID_MOD, ONLY : IDTOX,       IDTNOX,      IDTHNO3 
      USE TRACERID_MOD, ONLY : IDTISOP,     IDTACET,     IDTPRPE 
      USE TRACERID_MOD, ONLY : IDTH2,       IDTHD
      USE TRACERID_MOD, ONLY : NEMANTHRO ,  IDTSOA4
      USE TRACERID_MOD, ONLY : IDTSOAG,     IDTSOAM
      USE TRACERID_MOD, ONLY : IDTMONX,     IDTMBO, IDTC2H4
      USE WETSCAV_MOD,  ONLY : GET_WETDEP_NSOL
      USE WETSCAV_MOD,  ONLY : GET_WETDEP_IDWETD  

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! IFLX, LPAUSE
#     include "CMN_DIAG"   ! Diagnostic switches & arrays
#     include "CMN_O3"     ! FMOL, XNUMOL
#     include "comode.h"   ! IDEMS

      ! Local variables
      INTEGER            :: I, IREF, J, JREF, L, M, MM, MMB, LMAX
      INTEGER            :: N, NN, NMAX, NTEST
      INTEGER            :: IE, IN, IS, IW, ITEMP(3)
      TYPE (XPLEX)         :: SCALE_TMP(IIPAR,JJPAR)
      TYPE (XPLEX)         :: SCALE_I6,  SCALE_A6,  SCALE_A3,  SCALED    
      TYPE (XPLEX)         :: SCALEDYN,  SCALECONV, SCALESRCE, SCALECHEM 
      TYPE (XPLEX)         :: SCALEX,    SECONDS,   PMASS,     PRESSX
      TYPE (XPLEX)         :: FDTT,      AREA_M2,   DIAGb,     DIAGe
      
      ! For binary punch file, version 2.0
      CHARACTER (LEN=40) :: CATEGORY 
      TYPE (XPLEX)             :: ARRAY(IIPAR,JJPAR,LLPAR+1)
      TYPE (XPLEX)             :: LONRES, LATRES
      INTEGER            :: IFIRST, JFIRST, LFIRST
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180 = 1
      CHARACTER (LEN=20) :: MODELNAME 
      CHARACTER (LEN=40) :: UNIT
      CHARACTER (LEN=40) :: RESERVED = ''
!
!******************************************************************************
!  DIAG3 begins here!
!
!  Define scale factors for division.  
!  Add a small number (e.g. 1d-32) to prevent division by zero errors.
!******************************************************************************
!
      ! Now use counter variables from "time_mod.f" (bmy, 3/27/03)
      DIAGb     = GET_DIAGb()
      DIAGe     = GET_DIAGe()
      SECONDS   = ( DIAGe - DIAGb ) * 3600d0
      SCALED    = 1d0
      SCALEDYN  = XPLX( GET_CT_DYN()  ) + 1d-32
      SCALECONV = XPLX( GET_CT_CONV() ) + 1d-32
      SCALESRCE = XPLX( GET_CT_EMIS() ) + 1d-32
      SCALECHEM = XPLX( GET_CT_CHEM() ) + 1d-32
      SCALE_A3  = XPLX( GET_CT_A3()   ) + 1d-32
      SCALE_A6  = XPLX( GET_CT_A6()   ) + 1d-32
      SCALE_I6  = XPLX( GET_CT_I6()   ) + 1d-32
!
!******************************************************************************
!  Setup for binary punch file:
!
!  IFIRST, JFIRST, LFIRST = I, J, L indices of the starting grid box 
!  LONRES                 = DISIZE, cast to TYPE (XPLEX)
!  LATRES                 = DJSIZE, cast to TYPE (XPLEX)
!******************************************************************************
!
      IFIRST = GET_XOFFSET( GLOBAL=.TRUE. ) + 1
      JFIRST = GET_YOFFSET( GLOBAL=.TRUE. ) + 1
      LFIRST = 1
      LONRES = DISIZE
      LATRES = DJSIZE

      ! Get the proper model name and HALFPOLAR setting for the bpch file
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
!
!******************************************************************************
!  ND01: Rn, Pb, Be emissions (Category: "RN--SRCE")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Emissions of 222Rn             : kg/s       : SCALESRCE
!  (2)  Pb210  : Emissions of 210Pb             : kg/s       : SCALECHEM
!  (3)  Be7    : Emissions of 7Be               : kg/s       : SCALESRCE
!
!  and  Rn, Pb, Be lost to radioactive decay (Category: "RN-DECAY")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Loss of 222Rn                  : kg/s       : SCALECHEM
!  (2)  Pb210  : Loss of 210Pb                  : kg/s       : SCALECHEM
!  (3)  Be7    : Loss of 7Be                    : kg/s       : SCALECHEM
!******************************************************************************
!
      IF ( ND01 > 0 ) THEN
         CATEGORY = 'RN--SRCE'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(1)
            N  = TINDEX(1,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N
               
            ! Pb "emission" comes from chemical decay of Rn, which happens 
            ! in the chemistry routine, so use SCALECHEM (bmy, 1/27/03)
            IF ( N == IDTPB ) THEN
               SCALEX = SCALECHEM
            ELSE
               SCALEX = SCALESRCE
            ENDIF
                 
            ! Divide by # of emission timesteps
            DO L = 1, LD01
               ARRAY(:,:,L) = AD01(:,:,L,N) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD01,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD01) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND02: Rn, Pb, Be lost to radioactive decay (Category: "RN-DECAY")
!
!   # : Field  : Description                    : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  Rn222  : Loss of 222Rn                  : kg/s       : SCALECHEM
!  (2)  Pb210  : Loss of 210Pb                  : kg/s       : SCALECHEM
!  (3)  Be7    : Loss of 7Be                    : kg/s       : SCALECHEM
!******************************************************************************
!
      IF ( ND02 > 0 ) THEN
         CATEGORY = 'RN-DECAY'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(2)
            N  = TINDEX(2,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N

            ! Divide by # of chemistry timesteps
            DO L = 1, LD02
               ARRAY(:,:,L) = AD02(:,:,L,N) / SCALECHEM
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD02,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD02) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND03: Diagnostics from Hg0/Hg2/HgP offline simulation (eck, bmy, 1/20/05)
!******************************************************************************
!
      IF ( ND03 > 0 ) CALL WRITE_DIAG03
!
!******************************************************************************
!  ND04: Diagnostics from CO2 simulation (pns, bmy, 7/26/05)
!******************************************************************************
!
      IF ( ND04 > 0 ) CALL WRITE_DIAG04
!
!******************************************************************************
!  ND05: Production/Loss for coupled fullchem/aerosol runs (NSRCX==3) or
!        offline sulfate chemistry runs (NSRCX==10).      
!
!   # : Field  : Description                   : Units        : Scale factor
!  ----------------------------------------------------------------------------
!  (1 ) SO2dms : P(SO2) from DMS + OH          : kg S         : SCALEX
!  (2 ) SO2no3 : P(SO2) from DMS + NO3         : kg S         : SCALEX
!  (3 ) SO2    : Total P(SO2)                  : kg S         : SCALEX
!  (4 ) MSAdms : P(MSA) from DMS               : kg S         : SCALEX
!  (5 ) SO4gas : P(SO4) gas phase              : kg S         : SCALEX
!  (6 ) SO4aq  : P(SO4) aqueous phase          : kg S         : SCALEX
!  (7 ) PSO4   : Total P(SO4)                  : kg S         : SCALEX
!  (8 ) PSO4s  : Total P(SO4 from seasalt)     : kg S         : SCALEX
!  (9 ) LOH    : L(OH) by DMS                  : kg OH        : SCALEX
!  (10) LNO3   : L(NO3) by DMS                 : kg NO3       : SCALEX
!******************************************************************************
!
      IF ( ND05 > 0 ) THEN
         CATEGORY = 'PL-SUL=$'

         DO M = 1, TMAX(5)
            N = TINDEX(5,M)

            ! Tracers 9, 10 are OH, NO3
            ! and are in [kg] instead of [kg S]
            IF ( N < 9 ) THEN 
               UNIT = 'kg S'
            ELSE
               UNIT = 'kg'
            ENDIF

            NN     = N
            SCALEX = 1.d0

            DO L = 1, LD05
               ARRAY(:,:,L) = AD05(:,:,L,N) / SCALEX
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD05,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD05) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND06: Dust aerosol emissions
!
!   # : Field    : Description                     : Units      : Scale factor
!  --------------------------------------------------------------------------
!  (1)  DUST     : Soil dust (4 different classes) : kg         : 1
!******************************************************************************
!
      IF ( ND06 > 0 .and. LDUST ) THEN

         ! Category & unit string
         UNIT     = 'kg'
         CATEGORY = 'DUSTSRCE'

         ! Loop over # of dust bins
         DO N = 1, NDSTBIN 

            ! At present we have 4 dust bins
            IF ( N == 1 ) NN = IDTDST1 
            IF ( N == 2 ) NN = IDTDST2 
            IF ( N == 3 ) NN = IDTDST3 
            IF ( N == 4 ) NN = IDTDST4 

            ! Save dust into ARRAY
            ARRAY(:,:,1) = AD06(:,:,N) 

            ! Write to BPCH file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF     
!
!******************************************************************************
!  ND07: Emissions of BC and OC aerosols
!
!   # : Field    : Description                  : Units        : Scale factor
!  --------------------------------------------------------------------------
!  (1)  Carbon   : Carbonaceous aerosols        : kg           : 1
!******************************************************************************
!
      IF ( ND07 > 0 .and. LCARB ) THEN

         ! Unit
         UNIT = 'kg'
         
         !-------------------
         ! BC ANTHRO source
         !-------------------
         CATEGORY     = 'BC-ANTH'
         N            = IDTBCPI
         ARRAY(:,:,1) = AD07(:,:,1) 

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !-------------------
         ! BC BIOMASS source
         !-------------------
         CATEGORY     = 'BC-BIOB'
         N            = IDTBCPI
         ARRAY(:,:,1) = AD07(:,:,2) 
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !------------------- 
         ! BC BIOFUEL source
         !------------------- 
         CATEGORY     = 'BC-BIOF'
         N            = IDTBCPI
         ARRAY(:,:,1) = AD07(:,:,3) 
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !------------------------------ 
         ! H-philic BC from H-phobic BC
         !------------------------------ 
         CATEGORY     = 'PL-BC=$'
         N            = IDTBCPI

         DO L = 1, LD07
            ARRAY(:,:,L) = AD07_BC(:,:,L) 
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     LD07,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD07) )

         !------------------------------ 
         ! OC ANTHRO source
         !------------------------------ 
         CATEGORY     = 'OC-ANTH'
         N            = IDTOCPI
         ARRAY(:,:,1) = AD07(:,:,4) 
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !------------------------------ 
         ! OC BIOMASS source
         !------------------------------
         CATEGORY     = 'OC-BIOB'
         N            = IDTOCPI
         ARRAY(:,:,1) = AD07(:,:,5) 
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !------------------------------ 
         ! OC BIOFUEL source
         !------------------------------
         CATEGORY     = 'OC-BIOF'
         N            = IDTOCPI
         ARRAY(:,:,1) = AD07(:,:,6) 
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !------------------------------ 
         ! OC BIOGENIC source
         !------------------------------
         CATEGORY     = 'OC-BIOG'
         N            = IDTOCPI
         ARRAY(:,:,1) = AD07(:,:,7) 
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !------------------------------ 
         ! H-philic OC from H-phobic OC
         !------------------------------
         CATEGORY     = 'PL-OC=$'
         N            = IDTOCPI

         DO L = 1, LD07
            ARRAY(:,:,L) = AD07_OC(:,:,L) 
         ENDDO
         
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     LD07,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD07) )

         ! Only save extra SOA diagnostics if LSOA=T
         IF ( LSOA ) THEN

            !------------------------------
            ! NVOC SOURCE diagnostics
            !------------------------------
            DO N = 8, 12

               SELECT CASE ( N )

                  ! ALPH
                  CASE ( 8 )
                     CATEGORY = 'OC-ALPH'
                     NN       = IDTALPH

                  ! LIMO
                  CASE ( 9 )
                     CATEGORY = 'OC-LIMO'
                     NN       = IDTLIMO

                  ! TERP
                  CASE ( 10 )
                     CATEGORY = 'OC-TERP'
                     NN       = IDTLIMO + 1

                  ! ALCO
                  CASE ( 11 )
                     CATEGORY = 'OC-ALCO'
                     NN       = IDTLIMO + 2

                  ! SESQ
                  CASE ( 12 )
                     CATEGORY = 'OC-SESQ'
                     NN       = IDTLIMO + 3

               END SELECT

               ARRAY(:,:,1) = AD07(:,:,N)
                  
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     1,        IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )

            ENDDO

            !-----------------------------------------------
            ! SOA Production from NVOC oxidation [kg]
            ! 1:ALPH+LIMO+TERP, 2:ALCO, 3:SESQ, 4:ISOP
            !-----------------------------------------------
            CATEGORY = 'PL-OC=$'

            DO N = 1, 4

               IF ( N == 1 ) NN = IDTSOA1
               IF ( N == 2 ) NN = IDTSOA2
               IF ( N == 3 ) NN = IDTSOA3
               IF ( N == 4 ) NN = IDTSOA4  ! (tmf, bmy, 3/20/07)

               DO L = 1, LD07
                  ARRAY(:,:,L) = AD07_HC(:,:,L,N)
               ENDDO

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     LD07,     IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD07) )
            
            ENDDO

            !-----------------------------------------------
            ! SOA Production from GLYX and MGLY [kg]
            ! 1: SOAG <- GLYX;  2: SOAM <- MGLY   IN AEROSOL
            ! 3: SOAG <- GLYX;  4: SOAM <- MGLY   INCLOUD
            ! (tmf, 1/7/09)
	    ! Test if SOAG and SOAM tracers are valid before 
            ! saving them. (ccc, 1/7/09)
            !-----------------------------------------------
	    IF ( IDTSOAG /= 0 .AND. IDTSOAM /= 0 ) THEN
               CATEGORY     = 'SOAGM=$'

               DO N = 1, 4

                  IF ( N == 1 ) NN = 91
                  IF ( N == 2 ) NN = 92
                  IF ( N == 3 ) NN = 93
                  IF ( N == 4 ) NN = 94

                  DO L = 1, LD07
                     ARRAY(:,:,L) = AD07_SOAGM(:,:,L,N) 
                  ENDDO

                  CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                        HALFPOLAR, CENTER180, CATEGORY, NN,
     &                        UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                        IIPAR,     JJPAR,     LD07,     IFIRST,
     &                        JFIRST,    LFIRST,    ARRAY(:,:,1:LD07) )
            
               ENDDO
            ENDIF
         ENDIF
      ENDIF   
!
!******************************************************************************
!  ND08: Sea salt aerosol emissions
!
!   # : Field    : Description                     : Units      : Scale factor
!  --------------------------------------------------------------------------
!  (1)  SALA     : Accumulation mode seasalt       : kg         : 1
!  (2)  SALC     : Coarse mode seasalt             : kg         : 1
!******************************************************************************
!
      IF ( ND08 > 0 .and. LSSALT ) THEN

         ! Category & unit string
         UNIT     = 'kg'
         CATEGORY = 'SALTSRCE'

         ! Loop over seasalt tracers
         DO N = 1, 2

            ! At present we have 2 seasalts
            IF ( N == 1 ) NN = IDTSALA
            IF ( N == 2 ) NN = IDTSALC

            ! Save seasalts into ARRAY
            ARRAY(:,:,1) = AD08(:,:,N) 

            ! Write to BPCH file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF   
!
!******************************************************************************
!  ND09: HCN/CH3CN sources/sinks (Categories: "HCN-PL-$", "HCN-SRCE")
!
!  # : Field    : Description                     : Units       : Scale factor
! ----------------------------------------------------------------------------
! (1:N) sink    : Loss of tagged tracer to OH     : kg
! (N+1) HCNbb   : HCN   from biomass burning      : molec/cm2/s : SCALESRCE
! (N+2) CH3CNbb : CH3CN from biomass burning      : molec/cm2/s : SCALESRCE
! (N+3) HCNdf   : HCN   from domestic fossil fuel : molec/cm2/s : SCALESRCE
! (N+4) CH3CNdf : CH3CN from domestic fossil fuel : molec/cm2/s : SCALESRCE
! (N+5) HCNoc   : HCN   loss to ocean uptake      : molec/cm2/s : SCALECHEM
! (N+6) CH3CNoc : CH3CN loss to ocean uptake      : molec/cm2/s : SCALECHEM
!******************************************************************************
!
      IF ( ND09 > 0 ) THEN

         ! Binary punch file
         DO M = 1, TMAX(9)
            N  = TINDEX(9,M)
            IF ( N > N_TRACERS+6 ) CYCLE

            ! Test tracer number
            IF ( N <= N_TRACERS ) THEN

               !---------------------------
               ! HCN/CH3CN sinks
               !---------------------------              
               CATEGORY  = 'HCN-PL-$'
               UNIT      = 'kg'
               NN        = N
              
               DO L = 1, LD09
                  ARRAY(:,:,L) = AD09(:,:,L,N)
               ENDDO

               ! Save to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     LD09,     IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD09) )

            ELSE 

               !---------------------------
               ! HCN/CH3CN sources
               !---------------------------
               CATEGORY     = 'HCN-SRCE'
               UNIT         = 'molec/cm2/s'
               NN           = N - N_TRACERS

               ! Pick proper scale
               IF ( NN <= 4 ) THEN
                  SCALEX = SCALESRCE
               ELSE
                  SCALEX = SCALECHEM
               ENDIF

               ! Scale data
               ARRAY(:,:,1) = AD09_em(:,:,NN) / SCALEX
 
               ! Write to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     1,        IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDIF
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND10: H2/HD source diagnostics, prod and loss (phs, 9/18/07)
!
!   #  Field   : Description             : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) H2oh   : H2 Loss by OH           : mol/cm3/s : SCALECHEM
!  (2 ) H2iso  : H2 Prod from isoprene   : mol/cm3/s : SCALECHEM
!  (3 ) H2ch4  : H2 Prod from CH4        : mol/cm3/s : SCALECHEM
!  (4 ) H2ch3oh: H2 Prod from CH3OH      : mol/cm3/s : SCALECHEM
!  (5 ) H2mono : H2 Prod from monoprene  : mol/cm3/s : SCALECHEM
!  (6 ) H2acet : H2 Prod from acetone    : mol/cm3/s : SCALECHEM
!  (7 ) H2o1d  : H2 Loss by strat O1D    : mol/cm3/s : SCALECHEM
!
!  (8 ) HDoh   : H2 Loss by OH           : mol/cm3/s : SCALECHEM
!  (9 ) HDiso  : H2 Prod from isoprene   : mol/cm3/s : SCALECHEM
!  (10) HDch4  : H2 Prod from CH4        : mol/cm3/s : SCALECHEM
!  (11) HDch3oh: H2 Prod from CH3OH      : mol/cm3/s : SCALECHEM
!  (12) HDmono : H2 Prod from monoprene  : mol/cm3/s : SCALECHEM
!  (13) HDacet : H2 Prod from acetone    : mol/cm3/s : SCALECHEM
!  (14) HDo1d  : H2 Loss by strat O1D    : mol/cm3/s : SCALECHEM
!
!  (15) ALPHA  : OH k rates kHD/kH2 ratio: unitless  : SCALECHEM
!
!  (16) H2anth : H2 from Anthro Sources  : mol/cm2/s : SCALESRCE
!  (17) H2bb   : H2 from Biomass Burning : mol/cm2/s : SCALESRCE
!  (18) H2bf   : H2 from Biofuel Burning : mol/cm2/s : SCALESRCE
!  (19) H2ocean: H2 from Ocean           : mol/cm2/s : SCALESRCE
!  (19) HDocean: HD from Ocean           : mol/cm2/s : SCALESRCE
!
!  NOTES:
!  (1 ) Non zero only if ND10>0 and it is a H2/HD offline simulation
!  (2 ) 
!******************************************************************************
!
      IF ( ND10 > 0 ) THEN
         DO M = 1, TMAX(10)
            N = TINDEX(10,M)
            IF ( N > PD10 ) CYCLE

            ! Test tracer number (NEMISS=5, see "ndxx_setup.f" )
            IF ( N <= ( PD10 - 5 ) ) THEN

               !---------------------------
               ! H2/HD Prod-Loss
               !---------------------------              
               CATEGORY  = 'PL-H2HD-'
               UNIT      = 'molec/cm3/s'
               NN        = N
              
               DO L = 1, LD10
                  ARRAY(:,:,L) = AD10(:,:,L,N) / SCALECHEM
               ENDDO

               ! Save to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     LD10,     IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD10) )

            ELSE 

               !---------------------------
               ! H2/HD sources
               !---------------------------
               CATEGORY     = 'H2HD-SRC'
               UNIT         = 'molec/cm2/s'
               NN           = N - 15

               ! Scale data
               ARRAY(:,:,1) = AD10em(:,:,NN) / SCALESRCE
 
               ! Write to disk
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     1,        IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDIF
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND11: Acetone source & sink diagnostic (Category: "ACETSRCE")
!
!   # : Field  : Description                        : Units      : Scale factor
!  ----------------------------------------------------------------------------
!  (1)  ACETmo : Acetone source from MONOTERPENES   : at C/cm2/s : SCALESRCE
!  (2)  ACETmb : Acetone source from METHYL BUTENOL : at C/cm2/s : SCALESRCE 
!  (3)  ACETbg : Acetone source from DIRECT EMISSION: at C/cm2/s : SCALESRCE 
!  (4)  ACETdl : Acetone source from DRY LEAF MATTER: at C/cm2/s : SCALESRCE 
!  (5)  ACETgr : Acetone source from GRASSLANDS     : at C/cm2/s : SCALESRCE 
!  (6)  ACETop : Acetone source from OCEANS         : at C/cm2/s : SCALESRCE 
!  (7)  ACETol : Acetone sink   from OCEANS         : at C/cm2/s : SCALECHEM
!******************************************************************************
!
      IF ( ND11 > 0 ) THEN
         CATEGORY = 'ACETSRCE'
         UNIT     = 'atoms C/cm2/s'

         DO M = 1, TMAX(11)
            N  = TINDEX(11,M)
            IF ( N > PD11 ) CYCLE
            NN = N 
               
            ! Acetone ocean sink is on the chemistry timestep
            ! but acetone sources are all on the emission timestep
            IF ( N == 7 ) THEN
               SCALEX = SCALECHEM
            ELSE
               SCALEX = SCALESRCE
            ENDIF

            ARRAY(:,:,1) = AD11(:,:,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND12: distribution of suface emissions in the boundry layer: [fraction]
!
!   # : Field   : Description                         : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1) EMDIS-BL : Fraction of BL occupied by level L  : unitless : SCALECHEM
!******************************************************************************
!
      IF ( ND12 > 0 ) THEN
         UNIT     = 'unitless'
         CATEGORY = 'EMDIS-BL'

         DO L = 1, LD12
            ARRAY(:,:,L) = AD12(:,:,L) / SCALECHEM
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLTROP,   IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD12) )
      ENDIF
!
!******************************************************************************
!  ND13: Sulfur emissions (for DMS/SO2/SO4/MSA/NH3/NH4/NIT chemistry)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) DMS-BIOG : Biogenic DMS emission           : kg S     : 1
!  (2 ) SO2-AC-$ : Aircraft SO2 emission           : kg S     : 1
!  (3 ) SO2-AN-$ : Anthropogenic SO2 emission      : kg S     : 1
!  (4 ) SO2-BIOB : Biomass SO2 emission            : kg S     : 1
!  (5 ) SO2-BIOF : Biofuel SO2 emission            : kg S     : 1
!  (6 ) SO2-NV-$ : Non-eruptive volcano SO2 em.    : kg S     : 1
!  (7 ) SO2-EV-$ : Eruptive volcano SO2 emissions  : kg S     : 1
!  (8 ) SO4-AN-$ : Anthropogenic SO4 emission      : kg S     : 1
!  (9 ) NH3-ANTH : Anthropogenic NH3 emission      : kg NH3   : 1
!  (10) NH3-NATU : Natural source NH3 emission     : kg NH3   : 1
!  (11) NH3-BIOB : Biomass burning NH3 emission    : kg NH3   : 1
!  (12) NH3-BIOF : Biofuel burning NH3 emission    : kg NH3   : 1
!******************************************************************************
!
      IF ( ND13 > 0 .and. 
     &   ( ITS_A_FULLCHEM_SIM() .or. ITS_AN_AEROSOL_SIM() ) ) THEN
         UNIT = 'kg S'

         !==============================================================
         ! Biogenic DMS 
         !==============================================================
         CATEGORY     = 'DMS-BIOG'
         ARRAY(:,:,1) = AD13_DMS(:,:)
         N            = IDTDMS
         
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Aircraft SO2 
         !==============================================================
         CATEGORY = 'SO2-AC-$'
         N        = IDTSO2

         DO L = 1, LD13
            ARRAY(:,:,L) = AD13_SO2_ac(:,:,L)
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N, 
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD13,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD13) )

         !==============================================================
         ! Anthropogenic SO2 
         !==============================================================
         CATEGORY = 'SO2-AN-$'
         N        = IDTSO2
         
         DO L = 1, NOXEXTENT 
            ARRAY(:,:,L) = AD13_SO2_an(:,:,L)
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &               IIPAR,     JJPAR,     NOXEXTENT, IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:NOXEXTENT) )

         !==============================================================
         ! Biomass SO2 
         !==============================================================
         CATEGORY     = 'SO2-BIOB'
         ARRAY(:,:,1) = AD13_SO2_bb(:,:)
         N            = IDTSO2

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )


         !==============================================================
         ! Biofuel SO2
         !==============================================================
         CATEGORY     = 'SO2-BIOF'
         ARRAY(:,:,1) = AD13_SO2_bf(:,:)
         N            = IDTSO2

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Eruptive volcano SO2 
         !==============================================================
         CATEGORY = 'SO2-EV-$'
         N        = IDTSO2

         DO L = 1, LD13
            ARRAY(:,:,L) = AD13_SO2_ev(:,:,L)
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD13,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD13) )

         !==============================================================
         ! Non-eruptive volcano SO2 
         !==============================================================
         CATEGORY = 'SO2-NV-$'
         N        = IDTSO2

         DO L = 1, LD13 
            ARRAY(:,:,L) = AD13_SO2_nv(:,:,L)
         ENDDO
               
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD13,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD13) )


         !==============================================================
         ! Ship SO2     bec (5/17/04)
         ! New test on logical flag (phs, 3/2/09)
         !==============================================================
         IF ( LSHIPSO2 .OR. LEDGARSHIP .OR. LARCSHIP .OR.
     $        LEMEPSHIP ) THEN
            
            CATEGORY     = 'SO2-SHIP'
            ARRAY(:,:,1) = AD13_SO2_sh(:,:)
            N            = IDTSO2

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, N,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDIF

         !==============================================================
         ! Anthropogenic SO4 
         !==============================================================
         CATEGORY = 'SO4-AN-$'
         N        = IDTSO4

         DO L = 1, NOXEXTENT 
            ARRAY(:,:,L) = AD13_SO4_an(:,:,L)
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  N,
     &               UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &               IIPAR,     JJPAR,     NOXEXTENT, IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:NOXEXTENT) )

         !==============================================================
         ! Biofuel SO4
         !==============================================================
         CATEGORY     = 'SO4-BIOF'
         ARRAY(:,:,1) = AD13_SO4_bf(:,:)
         N            = IDTSO4

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )


         !==============================================================
         ! Anthropogenic NH3 
         !==============================================================
         UNIT         = 'kg'
         CATEGORY     = 'NH3-ANTH'
         ARRAY(:,:,1) = AD13_NH3_an(:,:) 
         N            = IDTNH3

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Natural source NH3
         !==============================================================
         CATEGORY     = 'NH3-NATU'
         ARRAY(:,:,1) = AD13_NH3_na(:,:) 
         N            = IDTNH3

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Biomass NH3
         !==============================================================
         CATEGORY     = 'NH3-BIOB'
         ARRAY(:,:,1) = AD13_NH3_bb(:,:)
         N            = IDTNH3

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Biofuel NH3 
         !==============================================================
         CATEGORY     = 'NH3-BIOF'
         ARRAY(:,:,1) = AD13_NH3_bf(:,:)
         N            = IDTNH3

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF    
!
!******************************************************************************
!  ND14: Upward mass flux from wet convection (NFCLDMX)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  CONVFLUP : Upward mass flux from wet conv  : kg/s     : SCALECONV
!
!  NOTES:
!  (1) Bug fix -- only write LD14 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND14 > 0 ) THEN
         CATEGORY = 'CV-FLX-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(14)
            N  = TINDEX(14,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N
               
            ARRAY(:,:,1:LD14) = CONVFLUP(:,:,1:LD14,N) / SCALECONV

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD14,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD14) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND15: Upward mass flux from boundary layer mixing (TURBDAY)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  TURBFLUX : Upward mass flux from BL mixing : kg/s     : SCALECONV
!
!  NOTES:
!  (1) Bug fix -- only write LD15 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND15 > 0 ) THEN
         CATEGORY = 'TURBMC-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(15)
            N  = TINDEX(15,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N

            ARRAY(:,:,1:LD15) = TURBFLUP(:,:,1:LD15,N) / SCALECONV
               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD15,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD15) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND16: Fraction of grid box experiencing LS or convective precipitation
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-FLS-$ : LS precip fraction          : unitless  : CT16(:,:,:,1)
!  (2) WD-FCV-$ : Convective precip fraction  : unitless  : CT16(:,:,:,2)
!******************************************************************************
!
      IF ( ND16 > 0 ) THEN

         ! Large-scale area of precipitation
         CATEGORY = 'WD-FRC-$'
         UNIT     = 'unitless'

         DO M = 1, TMAX(16)
            N  = TINDEX(16,M)
            IF ( N > PD16 ) CYCLE
            NN = N 
              
            DO L = 1, LD16
               SCALE_TMP(:,:) = XPLX( CT16(:,:,L,N) ) + 1d-20
               ARRAY(:,:,L)   = AD16(:,:,L,N) / SCALE_TMP(:,:)
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD16,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD16) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND17: Fraction of tracer lost rainout in LS and convective precip
!
!   #  Field    : Description                  : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-LSR-$ : Rainout fraction/LS Precip   : unitless  : CT17(:,:,:,1)
!  (2) WD-CVR-$ : Rainout fraction/conv precip : unitless  : CT17(:,:,:,2)
!
!  NOTES:
!  (1) Now loop over all soluble tracers (bmy, 3/19/04)
!  (2) Now use actual tracer number (bmy, 3/23/04)
!******************************************************************************
!
      IF ( ND17 > 0 ) THEN
         UNIT = 'unitless'

         ! Get max # of soluble tracers for this simulation
         NMAX = GET_WETDEP_NSOL()

         ! Loop over soluble tracers
         DO N = 1, NMAX

            ! Tracer number 
            NN = GET_WETDEP_IDWETD( N )

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(17) )
               MMB = TINDEX(17,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= NN ) CYCLE

            ! Large-scale rainout/washout fractions
            CATEGORY = 'WD-LSR-$'
               
            DO L = 1, LD17
               SCALE_TMP(:,:) = XPLX( CT17(:,:,L,1) ) + 1d-20
               ARRAY(:,:,L)   = AD17(:,:,L,N,1) / SCALE_TMP(:,:) 
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD17,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD17) )


            ! Convective rainout/washout fractions
            CATEGORY = 'WD-CVR-$'

            DO L = 1, LD17
               SCALE_TMP(:,:) = XPLX( CT17(:,:,L,2) ) + 1d-20
               ARRAY(:,:,L)   = AD17(:,:,L,N,2) / SCALE_TMP(:,:) 
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD17,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD17) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND18: Fraction of tracer lost to washout in LS or convective precip
!
!   #  Field    : Description                  : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WD-LSW-$ : Washout fraction/LS precip   : unitless  : CT18(:,:,:,1)
!  (2) WD-CVW-$ : Washout fraction/conv precip : unitless  : CT18(:,:,:,2)
!
!  NOTES:
!  (1) Now loop over all soluble tracers (bmy, 3/19/04)
!  (2) Now use actual tracer number (bmy, 3/23/04)
!******************************************************************************
!
      IF ( ND18 > 0 ) THEN
         UNIT = 'unitless'

         ! Get max # of soluble tracers for this simulation
         NMAX = GET_WETDEP_NSOL()

         DO N = 1, NMAX

            ! Tracer number
            NN = GET_WETDEP_IDWETD( N )

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(18) )
               MMB = TINDEX(18,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= NN ) CYCLE

            ! Large-scale rainout/washout fractions
            CATEGORY = 'WD-LSW-$'

            DO L = 1, LD18
               SCALE_TMP(:,:) = XPLX( CT18(:,:,L,1) ) + 1d-20
               ARRAY(:,:,L)   = AD18(:,:,L,N,1) / SCALE_TMP(:,:) 
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD18,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD18) )


            ! Convective washout fractions
            CATEGORY = 'WD-CVW-$'

            DO L = 1, LD18
               SCALE_TMP(:,:) = XPLX( CT18(:,:,L,2) ) + 1d-20
               ARRAY(:,:,L)   = AD18(:,:,L,N,2) / SCALE_TMP(:,:) 
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD18,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD18) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND19: CH4 loss
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) CH4-LOSS : CH4 removing by OH              : kg CH4   : 1
!******************************************************************************
!
      IF ( ND19 > 0 ) THEN

         UNIT = 'kg'

         !==============================================================
         ! CH4 Loss
         !==============================================================
         CATEGORY     = 'CH4-LOSS'
         N            = 1
        
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     LLPAR,    IFIRST,
     &               JFIRST,    LFIRST,    AD19(:,:,:) )

      ENDIF
!
!******************************************************************************
!  ND21: Optical depth diagnostics
!
!   # : Field : Description                        : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) OPTD   Cloud Optical Depth                 : unitless : SCALECHEM
!  (2 ) CLMO   Max Overlap Cloud Fraction (GEOS1,S): unitless : SCALECHEM
!   or  CLDF   3-D Total Cloud fraction   (GEOS3,4): unitless : SCALECHEM
!  (3 ) CLRO   Random  Overlap Cloud Fraction      : unitless : SCALECHEM
!  (4 ) OPD    Mineral Dust Optical Depth (400 nm) : unitless : none
!  (5 ) SD     Mineral Dust Surface Area           : cm2/cm3  : none
!  (6 ) OPSO4  Sulfate Optical Depth (400 nm)      : unitless : SCALECHEM
!  (7 ) HGSO4  Hygroscopic growth of SO4           : unitless : SCALECHEM
!  (8 ) SSO4   Sulfate Surface Area                : cm2/cm3  : SCALECHEM
!  (9 ) OPBC   Black Carbon Optical Depth (400 nm) : unitless : SCALECHEM
!  (10) HGBC   Hygroscopic growth of BC            : unitless : SCALECHEM
!  (11) SBC    Black Carbon Surface Area           : cm2/cm3  : SCALECHEM
!  (12) OPOC   Organic C Optical Depth (400 nm)    : unitless : SCALECHEM
!  (13) HGOC   Hygroscopic growth of OC            : unitless : SCALECHEM
!  (14) SOC    Organic Carbon Surface Area         : cm2/cm3  : SCALECHEM
!  (15) OPSSa  Sea Salt (accum) Opt Depth (400 nm) : unitless : SCALECHEM
!  (16) HGSSa  Hygroscopic growth of SSa           : unitless : SCALECHEM
!  (17) SSSa   Sea Salt (accum) Surface Area       : cm2/cm3  : SCALECHEM
!  (18) OPSSc  Sea Salt (coarse) Opt Depth(400 nm) : unitless : SCALECHEM
!  (19) HGSSc  Hygroscopic growth of SSc           : unitless : SCALECHEM
!  (20) SSSc   Sea Salt (coarse) Surface Area      : cm2/cm3  : SCALECHEM  
!
!  NOTES:
!  (1 ) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2 ) Don't divide monthly mean AOD by SCALECHEM (rvm, bmy, 12/8/00)
!  (3 ) Use SCALE_A6 for GEOS-2, GEOS-3 fields, since optical depths are read
!        in from disk every 6 hours.  Use SCALECHEM for GEOS-1, GEOS-STRAT
!        fields, since optical depths are computed every chemistry timestep.
!        Use SCALEDYN for CO-OH parameterization simulation. (bmy, 4/23/01)
!  (4 ) Now GEOS-2, GEOS-3 use SCALECHEM for ND21 (bmy, 8/13/01)
!  (5 ) Updated tracers for new aerosols from Mian Chin (rvm, bmy, 3/1/02)
!  (6 ) Now scale online dust fields by SCALECHEM (bmy, 4/9/04)
!  (7 ) Also save out extra diagnostics for cryst sulfur tracers (bmy, 1/5/05)
!******************************************************************************
!
      IF ( ND21 > 0 ) THEN
         CATEGORY = 'OD-MAP-$'

         ! ND21 is updated every chem timestep 
         SCALEX = SCALECHEM

         DO M = 1, TMAX(21)
            N  = TINDEX(21,M)
            IF ( N > PD21 ) CYCLE
            NN = N 
               
            ! Select proper unit string (cf list above)
            SELECT CASE( N ) 
               CASE ( 5, 8, 11, 14, 17, 20 )
                  UNIT = 'cm2/cm3'
               CASE DEFAULT
                  UNIT = 'unitless'
            END SELECT
               
            IF ( N > 3 .AND. N < 6 ) THEN

               ! Online or offline dust fields?
               IF ( LDUST ) THEN

                  ! If LDUST=T, then we are using online dust fields,
                  ! so we must scale by the chemistry timestep. (4/9/04)
                  ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N) / SCALEX
      
               ELSE
                  
                  ! If LDUST=F, then we are using offline monthly-mean
                  ! dust fields.  These don't have to be scaled by
                  ! the chemistry timestep. (bmy, 4/9/04)
                  ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N)

               ENDIF

            ELSE

               ! For all other types of optical depths, we need 
               ! to scale by the chemistry timestep (bmy, 4/9/04)
               ARRAY(:,:,1:LD21) = AD21(:,:,1:LD21,N) / SCALEX

            ENDIF
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD21,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD21) )
         ENDDO    

         !==============================================================
         ! If we are using the crystalline sulfate tracers (LCRYST=T),
         ! then also save out the extra ND21 diagnostics:
         !
         ! #21: Opt depth for HYSTERESIS CASE            [unitless]
         ! #22: Opt depth for SOLID CASE                 [unitless]
         ! #23: Opt depth for LIQUID CASE                [unitless]
         ! #24: Opt depth HYSTERESIS - Opt depth SOLID   [unitless] 
         ! #25: Opt depth HYSTERESIS - Opt depth LIQUID  [unitless]
         ! #26: Radiative forcing                        [W/m2    ]    
         !==============================================================
         IF ( LCRYST ) THEN
            
            ! Category
            CATEGORY = 'OD-MAP-$'

            ! Loop over extra 
            DO N = 1, 6               

               ! Define unit string
               IF ( N == 6 ) THEN
                  UNIT = 'W/m2'
               ELSE
                  UNIT = 'unitless'
               ENDIF

               ! Scale by chemistry timestep
               ARRAY(:,:,1) = AD21_cr(:,:,N) / SCALECHEM

               ! Save to BPCH file
               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, N+PD21,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     1,        IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDDO
         ENDIF
      ENDIF
!
!******************************************************************************
!  ND22: J-value diagnostics
!
!   #  : Field : Description                   : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1  ) JNO2  : NO2   J-Value                 : s-1   : SCALE_JV
!  (2  ) JHNO3 : HNO3  J-Value                 : s-1   : SCALE_JV
!  (3  ) JH2O2 : H2O2  J-Value                 : s-1   : SCALE_JV
!  (4  ) JCH2O : CH2O  J-Value                 : s-1   : SCALE_JV
!  (5  ) JO3   : O3    J-Value                 : s-1   : SCALE_JV
!  (6  ) POH   : OH-source from O3 photolysis  : s-1   : SCALE_JV
!  (7  ) JGLYX : GLYX  J-Value                 : s-1   : SCALE
!  (8  ) JMGLY : MGLY  J-Value                 : s-1   : SCALE
!  (71 ) JCH3I : CH3I  J-value (s^-1)          : s-1   : SCALE_JV
!  (81 ) JHCN  : HCN   J-value (s^-1)          : s-1   : SCALE_JV
!
!  NOTES:
!  (1) We must add TRCOFFSET for CH3I and HCN runs, so that GAMAP can
!       recognize those photo rates as distinct from the NO2, HNO3,
!       H2O2, CH2O, O3, and POH photo rates.
!  (2) Pick the right tracer field from AD22 if only a subset of tracers
!       are requested to be printed out. (ccc, 12/15/08)
!  (3) Add GLYX and MGLY tracers (tmf, 3/6/09)
!******************************************************************************
!
      IF ( ND22 > 0 ) THEN
         CATEGORY  = 'JV-MAP-$'
         SCALE_TMP = XPLX( CTJV ) + 1d-20
         UNIT      = 's-1'

         DO M = 1, TMAX(22)
            N  = TINDEX(22,M)
            !-----------------------------------------------------------------
            ! Prior to 12/15/08:
            !IF ( N > PD22 ) CYCLE
            !-----------------------------------------------------------------
            NN = N

            !-----------------------------------------------------------------
            ! NOTE: We can no longer select "all" in "input.geos", but we
            ! must specify the tracer #'s for ND22 explicitly:
            !
            ! Fullchem:                CH3I         HCN
            ! 1           = NOx        1 = CH3I     1 = HCN
            ! 7           = HNO3
            ! 8           = H2O2
            ! 20          = CH2O
            ! 55          = GLYX
            ! 56          = MGLY
            ! N_TRACERS+1 = O3 & OH
            !
            ! (ccc, bmy, 12/15/08)
            !-----------------------------------------------------------------
            IF ( NN >= N_TRACERS+1 ) THEN
               MM = 5    ! Write 'O3' and 'OH'
            ELSE
               SELECT CASE ( TRIM( TRACER_NAME(NN) ) )
                  CASE ( 'NOx', 'HCN', 'CH3I' )
                     MM = 1
                  CASE ( 'HNO3' )
                     MM = 2
                  CASE ( 'H2O2' )
                     MM = 3
                  CASE ( 'CH2O' )
                     MM = 4
                  CASE ( 'GLYX' )
                     MM = 7
                  CASE ( 'MGLY' )
                     MM = 8
                  CASE DEFAULT
                     MM = 0
               END SELECT
            ENDIF

            ! Skip if not a valid index
            IF ( MM == 0 ) CYCLE

            DO L = 1, LD22
               !---------------------------------------------------------------
               ! Prior to 12/15/08:
               !ARRAY(:,:,L) = AD22(:,:,L,N) / SCALE_TMP(:,:)
               !---------------------------------------------------------------
               ARRAY(:,:,L) = AD22(:,:,L,MM) / SCALE_TMP(:,:)
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, MM,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD22,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD22) )

            ! If we have just written out O3, then write out OH
            ! (ccc, bmy, 12/15/08)
            IF ( MM == 5 ) THEN
               MMB = 6
               DO L = 1, LD22
                  ARRAY(:,:,L) = AD22(:,:,L,MMB) / SCALE_TMP(:,:)
               ENDDO

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, MMB,    
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     LD22,     IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD22) )
            ENDIF
         ENDDO    
      ENDIF     
!
!******************************************************************************
!  ND24: Eastward mass flux from transport (TPCORE, XTP)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLEW : Eastward mass flux - transport  : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLEW is TYPE (XPLEX)...store to ARRAY, which is COMPLEX*16
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD24 levels out to the bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND24 > 0 ) THEN
         CATEGORY = 'EW-FLX-$'
         UNIT = 'kg/s'

         DO M = 1, TMAX(24)
            N  = TINDEX(24,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N
            
            ! (dkh, 02/09/12, adj32_022) 
            !ARRAY(:,:,1:LD24) = MASSFLEW(:,:,1:LD24,N) / SCALEDYN
            ARRAY(:,:,1:LD24) = MASSFLEW(:,:,LD24:1:-1,N) / SCALEDYN


            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD24,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD24) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND25: Northward mass flux from transport (TPCORE, YTP)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLNS : Northward mass flux - transport : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLNS is TYPE (XPLEX)...store to ARRAY, which is COMPLEX*16
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD25 levels out to the bpch file (bmy, 12/7/00)
!******************************************************************************
!  
      IF ( ND25 > 0 ) THEN
         CATEGORY = 'NS-FLX-$'
         UNIT = 'kg/s'

         DO M = 1, TMAX(25)
            N  = TINDEX(25,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N

            ! (dkh, 02/09/12, adj32_022) 
            !ARRAY(:,:,1:LD25) = MASSFLNS(:,:,1:LD25,N) / SCALEDYN
            ARRAY(:,:,1:LD25) = MASSFLNS(:,:,LD25:1:-1,N) / SCALEDYN


            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD25,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD25) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND26: Upward mass flux from transport (TPCORE, FZPPM)
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1)  MASSFLUP : Upward mass flux - transport    : kg/s     : SCALEDYN
!
!  NOTES:
!  (1) MASSFLNS is TYPE (XPLEX)...store to ARRAY, which is COMPLEX*16
!       before sending to BPCH or IJSCAL (bey, bmy, 4/23/99)
!  (2) Now only write LD26 levels to the bpch file (bmy, 12/7/00)
!******************************************************************************
!  
      IF ( ND26 > 0 ) THEN
         CATEGORY = 'UP-FLX-$'
         UNIT     = 'kg/s'

         DO M = 1, TMAX(26)
            N  = TINDEX(26,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N
            
            ! (dkh, 02/09/12, adj32_022) 
            !ARRAY(:,:,1:LD26) = MASSFLUP(:,:,1:LD26,N) / SCALEDYN
            ARRAY(:,:,1:LD26) = MASSFLUP(:,:,LD26:1:-1,N) / SCALEDYN

               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES, 
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD26,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD26) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND27: Cross-tropopause Stratospheric Influx of Ox 
!
!   #  : Field : Description                   : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1) : Ox    : Ox from the stratosphere      : kg/s  : SCALEDYN
!
!  NOTES:
!  (1) Only print out if we are doing a NOx-Ox-HC run (NSRCX == 3)
!       or a single tracer Ox run (NSRCX == 6). (bey, bmy, 11/10/99)
!  (2) Now consider the cross-tropopause stratospheric influx of ozone, 
!       which, in some grid boxes, includes horizontal influxes as well as 
!       up(down)ward flux. (qli, 1/5/2000) 
!  (3) Now error check for N > NTRACE (bmy, 10/23/01)
!  (4) NOTE: There is a problem with for ND27 with GEOS-4.  Djj says that 
!       the downward flux at the 200 hPa level should be more or less the
!       same as the ND27 diagnostic. (bmy, 10/15/04)
!  (5) Now provides stratrospheric flux of H2/HD if it is a H2/HD simulation
!       (lyj, phs, 9/18/07)
!******************************************************************************
!
#if   !defined( GEOS_4 ) 
      IF ( ND27 > 0 .and. IDTOX > 0 ) then
         IF ( ( IDTOX > 0 .and. 
     &        ( ITS_A_FULLCHEM_SIM() .or. ITS_A_TAGOX_SIM() ) ) .OR. 
     &        ( ITS_A_H2HD_SIM()                              ) ) THEN

            CATEGORY = 'STRT-FLX'
            UNIT     = 'kg/s'

            ! Full chemistry   -- compute NOx, Ox, HNO3 fluxes
            ! H2/HD            -- compute H2,  HD fluxes
            ! Single tracer Ox -- compute Ox flux only, hardwire
            !                     to tracer = 1 (bmy, 2/7/00)
            IF ( ITS_A_FULLCHEM_SIM() ) THEN
               ITEMP = (/ IDTNOX, IDTOX, IDTHNO3 /)
            ELSE IF ( ITS_A_H2HD_SIM() ) THEN
               ITEMP = (/ IDTH2, IDTHD, 0 /)
            ELSE
               ITEMP = (/ 1, 0, 0 /)
            ENDIF
            
            ! Loop over tracers
            DO M = 1, 3
               N = ITEMP(M)
               IF ( N == 0        ) CYCLE
               IF ( N > N_TRACERS ) CYCLE

               ! Loop over grid boxes
               DO J = 1, JJPAR
               DO I = 1, IIPAR

                  ! Get the level of the tropopause
                  L = LPAUSE(I,J)

                  ! Initialize integer flags
                  IS = 0
                  IN = 0
                  IW = 0
                  IE = 0

                  ! Set integer flags based on the value of each bit of IFLX
                  IF ( BTEST( IFLX(I,J), 0 ) ) IS = 1
                  IF ( BTEST( IFLX(I,J), 1 ) ) IN = 1
                  IF ( BTEST( IFLX(I,J), 2 ) ) IW = 1
                  IF ( BTEST( IFLX(I,J), 3 ) ) IE = 1

                  ! Add fluxes from the top, south, and west
                  ARRAY(I,J,1) = MASSFLUP(I,J,L,N)          +
     &                           ( MASSFLNS(I,J,L,N) * IS ) +  
     &                           ( MASSFLEW(I,J,L,N) * IW ) 
               
                  ! Add fluxes from the north 
                  ! (take poles into account !)
                  IF ( J < JJPAR ) THEN
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                              ( MASSFLNS(I,J+1,L,N) * IN )
                  ELSE 
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                              ( MASSFLNS(I,  1,L,N) * IN )
                  ENDIF

                  ! Add fluxes from the east 
                  !(wrap around dateline if necessary)
                  IF ( I < IIPAR ) THEN
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                              ( MASSFLEW(I+1,J,L,N) * IE )
                  ELSE 
                     ARRAY(I,J,1) = ARRAY(I,J,1) -
     &                             ( MASSFLEW(  1,J,L,N) * IE )
                  ENDIF
               ENDDO
               ENDDO

               UNIT = 'kg/s'
               NN   = N
                  
               ARRAY(:,:,1) = ARRAY(:,:,1) / SCALEDYN

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                     IIPAR,     JJPAR,     PD27,     IFIRST,
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1) )
            ENDDO
         ENDIF
      ENDIF
#endif
!
!******************************************************************************
!  ND28: Biomass burning diagnostic 
!
!   # : Field : Description   : Units            : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) NOx   : NOx            : molec NOx /cm2/s : SCALESRCE
!  (4 ) CO    : CO             : molec CO  /cm2/s : SCALESRCE
!  (9 ) ACET  : Acetone        : atoms C   /cm2/s : SCALESRCE
!  (10) MEK   : Ketones(>C3)   : atoms C   /cm2/s : SCALESRCE
!  (11) ALD2  : Acetaldehyde   : atoms C   /cm2/s : SCALESRCE
!  (18) PRPE  : Propene        : atoms C   /cm2/s : SCALESRCE
!  (19) C3H8  : Propane        : atoms C   /cm2/s : SCALESRCE
!  (20) C2HO  : Formaldehyde   : molec CH2O/cm2/s : SCALESRCE
!  (21) C2H6  : Ethane         : atoms C   /cm2/s : SCALESRCE
!  (26) SO2   : Sulfur dioxide : molec SO2 /cm2/s : SCALESRCE
!  (30) NH3   : Ammonia        : molec NH3 /cm2/s : SCALESRCE
!  (34) BCPO  : Black carbon   : atoms C   /cm2/s : SCALESRCE
!  (35) OCPO  : Organic carbon : atoms C   /cm2/s : SCALESRCE
!
!  NOTES:
!  (1) Use the F90 intrinsic "ANY" function to make sure that N 
!       corresponds to actual biomass burning tracers (bmy, 4/8/99)
!  (2) ND28 now uses allocatable array AD28 instead of AIJ. (bmy, 3/16/00)
!  (3) Now write biofuel burning tracers to the punch file in the same order 
!       as they are listed in "diag.dat". (bmy, 4/17/01)
!******************************************************************************
!
      IF ( ND28 > 0 ) THEN
         CATEGORY = 'BIOBSRCE'
         UNIT     = ''

         DO M = 1, TMAX(28)
            N  = TINDEX(28,M)
            IF ( .not. ANY( BIOTRCE == N ) ) CYCLE
            NN = N
            
            DO MM = 1, NBIOMAX
               IF ( BIOTRCE(MM) == NN ) THEN
                  MMB = MM
                  EXIT
               ENDIF
            ENDDO

            ARRAY(:,:,1) = AD28(:,:,MMB) / SCALESRCE 

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND29: CO source diagnostics
!
!   #  Field  : Description             : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) COanth : CO from Anthro Sources  : mol/cm2/s : SCALESRCE
!  (2) CObb   : CO from Biomass Burning : mol/cm2/s : SCALESRCE
!  (3) CObf   : CO from Biofuel Burning : mol/cm2/s : SCALESRCE
!  (4) COmeth : CO from Methanol        : mol/cm2/s : SCALESRCE
!  (5) COmono : CO from Monoterpenes    : mol/cm2/s : SCALESRCE
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) ND29 now uses allocatable array AD29 instead of AIJ. (bmy, 3/16/00)
!  (3) Added CO-sources from isoprene and monoterpenes (bnd, bmy, 1/2/01)
!******************************************************************************
!
      IF ( ND29 > 0 ) THEN
         CATEGORY ='CO--SRCE'
         UNIT    = 'mol/cm2/s'

         DO M = 1, TMAX(29)
            N  = TINDEX(29,M)
            IF ( N > PD29 ) CYCLE
            NN = N 

            ARRAY(:,:,1) = AD29(:,:,N) / SCALESRCE

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND30: Land map diagnostic
!
!   #  Field : Description             : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) LWI   : GMAO Land-Water indices  : unitless  : SCALED 
!
!  NOTES: 
!  (1) Values are: 0=water; 1=land; 2=ice (bmy, 8/18/05)
!******************************************************************************
!
      IF ( ND30 > 0 ) THEN
         CATEGORY = 'LANDMAP'
         UNIT     = 'unitless'
            
         ARRAY(:,:,1) = AD30(:,:) / SCALEDYN
         NN           = 1 

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF
!
!******************************************************************************
!  ND31: Surface pressure diagnostic
!
!   #  Field : Description                      : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) Pedge : Pressure at bot edge of level L  : mb        : SCALEDYN
!
!  NOTES: 
!  (1) The ASCII punch file was using SCALE2 instead of SCALE1.
!       This has now been fixed. (hyl, bmy, 12/21/99).
!  (2) Now use AD31 dynamically allocatable array (bmy, 2/17/00)
!  (3) Bug fix: write out 1 level to the bpch file (bmy, 12/7/00)
!  (4) Now remove SCALE1, replace with SCALEDYN (bmy, 2/24/03)
!  (5) Now save out true pressure at level edges.  Now   (bmy, 5/8/07)
!******************************************************************************
!   
      IF ( ND31 > 0 ) THEN
         CATEGORY          = 'PEDGE-$'
         UNIT              = 'mb'
         ARRAY(:,:,1:LD31) = AD31(:,:,1:LD31) / SCALEDYN
         NN                = 1

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, NN,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD31,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD31) )
      ENDIF
!
!******************************************************************************
!  ND32: NOx source diagnostic
!
!  Levels        : Field                  : Units       : Scale Factor
!  -------------------------------------------------------------------------
!  1 - LLTROP    : Aircraft NOx           : molec/cm2/s : SCALESRCE
!  1 - NOXEXTENT : Anthropogenic NOx      : molec/cm2/s : SCALESRCE
!  Surface       : Biomass Burning NOx    : molec/cm2/s : SCALESRCE
!  Surface       : Biofuel Burning NOx    : molec/cm2/s : SCALESRCE
!  Surface       : Fertilizer NOx         : molec/cm2/s : SCALESRCE
!  1 - LLCONVM   : Lightning NOx          : molec/cm2/s : SCALESRCE
!  Surface       : Soil NOx               : molec/cm2/s : SCALESRCE
!  Above TP      : NOx from upper boundary: molec/cm2/s : SCALEDYN
!
!  Print out all of the types of NOx, for all levels.
!
!  NOTES:
!  (1) Only print out ND32 if for an O3 chemistry run ( NSRCX == 3 ),
!       and if NOx is a defined tracer ( IDTNOX > 0 ). (bmy, 5/26/99)
!  (2) ND32 now uses allocatable arrays instead of AIJ. (bmy 3/16/00)
!  (3) Added biofuel burning to ND32 diagnostic (bmy, 9/12/00)
!******************************************************************************
!
      IF ( ND32 > 0 .and. IDTNOX > 0 .and. ITS_A_FULLCHEM_SIM() ) THEN

         ! All categories of NOx are in molec/cm2/s
         UNIT = 'molec/cm2/s'

         !==============================================================
         ! Aircraft NOx
         !==============================================================
         CATEGORY = 'NOX-AC-$'
            
         DO L = 1, LLTROP
            ARRAY(:,:,L) = AD32_ac(:,:,L) / SCALESRCE               
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLTROP,   IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LLTROP) )

         !==============================================================
         ! Anthropogenic NOx
         !==============================================================
         CATEGORY = 'NOX-AN-$'

         DO L = 1, NOXEXTENT 
            ARRAY(:,:,L) = AD32_an(:,:,L) / SCALESRCE
         ENDDO
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,    LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY,  IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,     RESERVED,   
     &               IIPAR,     JJPAR,     NOXEXTENT, IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:NOXEXTENT) )

         !==============================================================
         ! Biomass Burning NOx
         !==============================================================
         CATEGORY     = 'NOX-BIOB'
         ARRAY(:,:,1) = AD32_bb(:,:) / SCALESRCE

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Binary punch file: NOx from Biofuel
         !==============================================================
         CATEGORY     = 'NOX-BIOF'
         ARRAY(:,:,1) = AD32_bf(:,:) / SCALESRCE

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &                  HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Fertilizer NOx
         !==============================================================
         CATEGORY     = 'NOX-FERT'
         ARRAY(:,:,1) = AD32_fe(:,:) / SCALESRCE
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Lightning NOx
         !==============================================================
         CATEGORY = 'NOX-LI-$'

         DO L = 1, LLCONVM 
            ARRAY(:,:,L) = AD32_li(:,:,L) / SCALESRCE
         ENDDO
               
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX, 
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LLCONVM,  IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LLCONVM) )

         !==============================================================
         ! Soil NOx
         !==============================================================
         CATEGORY     = 'NOX-SOIL'
         ARRAY(:,:,1) = AD32_so(:,:) / SCALESRCE
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

         !==============================================================
         ! Stratospheric NOx (boundary condition)
         !==============================================================
         CATEGORY     = 'NOX-STRT'
         ARRAY(:,:,1) = AD32_ub(:,:) / SCALEDYN
            
         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &               HALFPOLAR, CENTER180, CATEGORY, IDTNOX,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )
      ENDIF
!
!******************************************************************************
!  ND33: Atmospheric column sum of Tracer
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) COLUMN-T : Trop. Column Sum of Tracer  : kg        : SCALEDYN
!
!  NOTES: 
!  (1) Now use dynamically allocatable array AD33 (bmy, 2/17/00)
!  (2) Rename category to COLUMN-T, since this is a column sum of tracer over
!       the entire atmosphere, not just the troposphere. (bmy, 4/3/02)
!  (3) Now replace SCALE1 with SCALEDYN (bmy, 3/27/03)
!******************************************************************************
!
      IF ( ND33 > 0 ) THEN
         CATEGORY = 'COLUMN-T'
         UNIT     = 'kg'

         DO M = 1, TMAX(33)
            N  = TINDEX(33,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N
            
            ARRAY(:,:,1) = AD33(:,:,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,     
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF  
!
!******************************************************************************
!  ND34: Biofuel burning diagnostic 
!
!   # : Field : Description         : Units            : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) NOx   : NOx                 : molec NOx /cm2/s : SCALESRCE
!  (4 ) CO    : CO                  : molec CO  /cm2/s : SCALESRCE
!  (5 ) ALK4  : Alkanes(>C4)        : atoms C   /cm2/s : SCALESRCE
!  (9 ) ACET  : Acetone             : atoms C   /cm2/s : SCALESRCE
!  (10) MEK   : Metyl Ethyl Ketone  : atoms C   /cm2/s : SCALESRCE
!  (11) ALD2  : Acetaldehyde        : atoms C   /cm2/s : SCALESRCE
!  (18) PRPE  : Alkenes(>=C3)       : atoms C   /cm2/s : SCALESRCE
!  (19) C3H8  : Propane             : atoms C   /cm2/s : SCALESRCE
!  (20) CH2O  : Formaldehyde        : molec CH2O/cm2/s : SCALESRCE
!  (21) C2H6  : Ethane              : atoms C   /cm2/s : SCALESRCE
!
!  NOTES:
!  (1) Use the F90 intrinsic "ANY" function to make sure that N 
!       corresponds to actual biofuel burning tracers (bmy, 3/15/01)
!  (3) Now write biofuel burning tracers to the punch file in the same order 
!       as they are listed in "diag.dat". (bmy, 4/17/01)
!  (4) Use BFTRACE and NBFTRACE to get the right index for AD34. 
!      (ccc, 12/8/2008)
!******************************************************************************
!
      IF ( ND34 > 0 ) THEN
         CATEGORY = 'BIOFSRCE'
         UNIT     = ''
         
         DO M = 1, TMAX(34)
            N  = TINDEX(34,M)
            IF ( .not. ANY( BFTRACE == N ) ) CYCLE
            NN = N

            DO MM = 1, NBFTRACE
               IF ( BFTRACE(MM) == NN ) THEN
                  MMB = MM
                  EXIT
               ENDIF
            ENDDO

            ARRAY(:,:,1) = AD34(:,:,MMB) / SCALESRCE

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND35: Tracer concentration at 500 mb 
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) 500-AVRG : Tracer at 500 mb            : v/v       : SCALEDYN
!
!  NOTES:
!  (1) Now use dynamically allocatable array AD35 (bmy, 2/17/00)
!  (2) Now replace SCALE1 with SCALEDYN (bmy, 2/24/03)
!******************************************************************************
!
      IF ( ND35 > 0 ) THEN
         CATEGORY = '500-AVRG'        
         UNIT     = ''

         DO M = 1, TMAX(35)
            N  = TINDEX(35,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N
               
            ARRAY(:,:,1) = AD35(:,:,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF  
!
!******************************************************************************
!  ND36: Anthropogenic source diagnostic
!
!   #   Field  : Description                     : Units         : S. Factor
!  ---------------------------------------------------------------------------
!  (1 ) NOx    : NOx                             : mol/cm2/s     : SCALE3
!  (4 ) CO     : CO                              : mol/cm2/s     : SCALE3
!  (5 ) ALK4   : Alkanes(>C4)                    : atoms C/cm2/s : SCALE3
!  (9 ) ACET   : Acetone                         : atoms C/cm2/s : SCALE3
!  (10) MEK    : Ketones(>C3)                    : atoms C/cm2/s : SCALE3
!  (18) PRPE   : Propene                         : atoms C/cm2/s : SCALE3 
!  (19) C3H8   : Propane                         : atoms C/cm2/s : SCALE3 
!  (21) C2H6   : Ethane                          : atoms C/cm2/s : SCALE3
!  (71) CH3Ioc : Methyl Iodide (oceanic source)  : ng/m2/s       : SCALE3
!  (72) CH3Ibb : Methyl Iodide (biomass burning) : ng/m2/s       : SCALE3
!  (73) CH3Iwb : Methyl Iodide (wood burning)    : ng/m2/s       : SCALE3
!  (74) CH3Irc : Methyl Iodide (rice paddies)    : ng/m2/s       : SCALE3
!  (75) CH3Iwl : Methyl Iodide (wetlands)        : ng/m2/s       : SCALE3
!
!  NOTES:
!  (1) ND36 is also used for CH3I emissions diagnostics when NSRCX=2.
!  (2) For an O3 run (NSRCX = 3, the "default" run) make sure that the 
!       tracer number N matches an entry in the IDEMS emission index 
!       array (bmy, 4/9/99)  
!  (3) Write the tracers out to the punch file in the same order as
!       they are listed in the IDEMS array.  Thus, we have to re-assign
!       N = IDEMS(M) after we test to make sure it is a valid tracer
!       number (bmy, 4/16/99)
!  (4) For a CH3I run, make sure that the tracer number N is not larger
!       than NTRACE (bmy, 4/9/99) 
!  (5) ND36 now uses the AD36 array instead of AIJ. (bmy, 3/16/00)
!  (6) Rewritten for clarity; also fixed for CH3I (bmy, 7/25/06)
!  (7) Bug fix: given the tracer number, now search for entry in IDEMS
!       to jive with historical baggage (bmy, 3/6/07)
!******************************************************************************
!                     
      IF ( ND36 > 0 ) THEN

         ! Loop over # of tracers
         DO M = 1, TMAX(36)
            
            ! Get the tracer # from input.geos
            N  = TINDEX(36,M)

            IF ( ITS_A_CH3I_SIM() ) THEN

               !--------------------------------------------------------
               ! For CH3I simulation only 
               !--------------------------------------------------------
               CATEGORY = 'CH3ISRCE'
               UNIT     = 'ng/m2/s'
               IF ( N > NEMANTHRO ) CYCLE 
               
               ! Tracer number
               NN = N
               
               ! Index for AD36 array
               MM = M

            ELSE 

               !--------------------------------------------------------
               ! For full-chemistry.  Note, due to historical baggage,
               ! the order of the tracers in AD36 array corresponds to
               ! the order as given in IDEMS.  Therefore, for the given
               ! tracer number N, we must find the corresponding entry
               ! in IDEMS. (bmy, 3/5/07)
               !--------------------------------------------------------
               CATEGORY = 'ANTHSRCE'
               UNIT     = ''

               ! reset these
               MM = 0
               NN = 0

               ! Given the tracer number N, find the proper entry in the 
               ! IDEMS array and select that for output (bmy, 3/5/07)
               DO NMAX = 1, NEMANTHRO
                  IF ( N == IDEMS(NMAX) ) THEN
                     MM = NMAX
                     NN = N
                     EXIT
                  ENDIF
               ENDDO

               ! We haven't found a match, skip to next tracer
               IF ( MM == 0 ) CYCLE
                  
            ENDIF

            ! Divide by seconds
            ARRAY(:,:,1) = AD36(:,:,MM) / SECONDS

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND37: Fraction of tracer scavenged in convective cloud updrafts
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETCVF-$ : Scavenging fraction         : unitless  : SCALECONV
!******************************************************************************
!
      IF ( ND37 > 0 ) THEN
         CATEGORY = 'MC-FRC-$'
         UNIT     = 'unitless'

         ! Get actual # of soluble tracers
         NMAX = GET_WETDEP_NSOL()

         ! Loop over soluble tracers
         DO N = 1, NMAX

            ! Tracer number 
            NN = GET_WETDEP_IDWETD( N )

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(37) )
               MMB = TINDEX(37,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= NN ) CYCLE

            DO L = 1, LD37
               ARRAY(:,:,L) = AD37(:,:,L,N) / SCALECONV
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD37,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD37) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND38: Rainout loss of tracer in convective updrafts
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETDCV-$ : Rainout loss of tracer      : kg/s      : SCALECONV
!
!  NOTES:
!  (1) Now write only LD38 levels to bpch file (bmy, 12/7/00)
!******************************************************************************
!
      IF ( ND38 > 0 ) THEN
         CATEGORY = 'WETDCV-$'
         UNIT     = 'kg/s'

         ! Get actual # of soluble tracers
         M = GET_WETDEP_NSOL()

         ! Loop over soluble tracers
         DO N = 1, M

            ! Tracer number
            NN = GET_WETDEP_IDWETD( N )

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(38) )
               MMB = TINDEX(38,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= NN ) CYCLE

            ! Divide by # of convective timesteps
            DO L = 1, LD38
               ARRAY(:,:,L) = AD38(:,:,L,N) / SCALECONV
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD38,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD38) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND39: Rainout loss of tracer in large scale rains 
!
!   #  Field    : Description                 : Units     : Scale factor
!  --------------------------------------------------------------------------
!  (1) WETDLS-$ : Large-scale loss of tracer  : kg/s      : SCALEDYN
!******************************************************************************
!
      IF ( ND39 > 0 ) THEN
         CATEGORY = 'WETDLS-$'
         UNIT     = 'kg/s'

         ! Get actual # of soluble tracers
         M = GET_WETDEP_NSOL()
            
         ! Loop over soluble tracers
         DO N = 1, M
               
            ! Tracer number
            NN = GET_WETDEP_IDWETD( N )

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(39) )
               MMB = TINDEX(39,MM)
               MM  = MM + 1
            ENDDO
               
            IF ( MMB /= NN ) CYCLE

            ! Divide by # of wetdep (= dynamic) timesteps
            DO L = 1, LD39
               ARRAY(:,:,L) = AD39(:,:,L,N) / SCALEDYN
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     LD39,     IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD39) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND41: Afternoon boundary layer heights
!******************************************************************************
!
      IF ( ND41 > 0 ) CALL WRITE_DIAG41
!
!******************************************************************************
!  ND42: SOA concentrations [ug/m3]
!******************************************************************************
!
      IF ( ND42 > 0 ) CALL WRITE_DIAG42
!
!******************************************************************************
!  ND42: Free diagnostic as of 11/24/99
!
!  ND43: Chemical production of OH and NO
!
!   # : Field : Description             : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1)  OH    : OH  Chemical Diagnostic : mol/cm3 : CTOH
!  (2)  NO    : NO  Chemical Diagnostic : v/v     : CTNO
!  (3)  HO2   : HO2 Chemical Diagnostic : v/v     : CTHO2
!  (4)  NO2   : NO2 Chemical Diagnostic : v/v     : CTNO2
!  (5)  NO3   : NO3 Chemical Diagnostic : v/v     : CTNO3
!
!  NOTES:
!  (1) Print output for either a NOx-Ox-HC run (NSRCX == 3), or a CO run
!       with parameterized OH (NSRCX== 5).  (bmy, 4/17/00)
!  (2) Add parentheses in IF test since .AND. has higher precedence
!       than .OR. (jsw, bmy, 12/5/00)
!  (3) Added HO2, NO2 to ND43 (rvm, bmy, 2/27/02)
!  (4) Added NO3 to ND43 (bmy, 1/16/03)
!  (5) Now uses 3D counters (phs, 1/24/07)
!  (6) Now assume that LD43 can't be higher than LD45 (phs, 1/24/07)
!  (7) Check that CTxx are not zero, instead of adding 1e-20 (phs, 11/13/07)
!******************************************************************************
!
      IF ( ND43 > 0 .and. ITS_A_FULLCHEM_SIM() ) THEN

         CATEGORY = 'CHEM-L=$' 

         DO M = 1, TMAX(43)
            N  = TINDEX(43,M)
            NN = N 

            ! default units
            UNIT = 'v/v'


            SELECT CASE ( N )

               ! OH
               CASE ( 1 )
                  WHERE( CTOH /= 0 )
                     ARRAY(:,:,1:LD43)%r =AD43(:,:,1:LD43,N)%r / CTOH
                     ARRAY(:,:,1:LD43)%i = AD43(:,:,1:LD43,N)%i /CTOH 
                  ELSEWHERE
                     ARRAY(:,:,1:LD43)%r = 0.d0
                     ARRAY(:,:,1:LD43)%i = 0.d0
                  ENDWHERE

                  UNIT = 'molec/cm3'
                  
               ! NO
               CASE ( 2 ) 
                  WHERE( CTNO /= 0 )
                     ARRAY(:,:,1:LD43)%r = AD43(:,:,1:LD43,N)%r /
     $                   ( CTNO )
                     ARRAY(:,:,1:LD43)%i = AD43(:,:,1:LD43,N)%i /
     $                   ( CTNO )
                  ELSEWHERE
                     ARRAY(:,:,1:LD43)%r = 0.d0
                     ARRAY(:,:,1:LD43)%i = 0.d0
                  ENDWHERE
                  

               ! HO2 (rvm, bmy, 2/27/02)
               CASE ( 3 )
                  WHERE( CTHO2 /= 0 )
                     ARRAY(:,:,1:LD43)%r = AD43(:,:,1:LD43,N)%r /
     $                    ( CTHO2 )
                     ARRAY(:,:,1:LD43)%i = AD43(:,:,1:LD43,N)%i /
     $                    ( CTHO2 )
                  ELSEWHERE
                     ARRAY(:,:,1:LD43)%r = 0.d0
                     ARRAY(:,:,1:LD43)%i = 0.d0
                  ENDWHERE


               ! NO2 (rvm, bmy, 2/27/02)
               CASE ( 4 ) 
                  WHERE( CTNO2 /= 0 )
                     ARRAY(:,:,1:LD43)%r = AD43(:,:,1:LD43,N)%r /
     $                    ( CTNO2 )
                     ARRAY(:,:,1:LD43)%i = AD43(:,:,1:LD43,N)%i /
     $                    ( CTNO2 )
                  ELSEWHERE
                     ARRAY(:,:,1:LD43)%r = 0.d0
                     ARRAY(:,:,1:LD43)%i = 0.d0
                  ENDWHERE
                  

               ! NO3 (rjp, bmy, 1/16/03)
               CASE ( 5 )
                  WHERE( CTNO3 /= 0 )
                     ARRAY(:,:,1:LD43)%r = AD43(:,:,1:LD43,N)%r /
     $                    ( CTNO3 )
                     ARRAY(:,:,1:LD43)%i = AD43(:,:,1:LD43,N)%i /
     $                    ( CTNO3 )
                  ELSEWHERE
                     ARRAY(:,:,1:LD43)%r = 0.d0
                     ARRAY(:,:,1:LD43)%i = 0.d0
                  ENDWHERE


               CASE DEFAULT
                  CYCLE

            END SELECT


            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD43,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD43) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND44: Drydep flux (molec/cm2/s) and velocity (cm/s) diagnostics
!
!   #   : Field    : Quantity           : Units               : Scale factor
!  -------------------------------------------------------------------------
!  (1 ) : DRYD-FLX : drydep fluxes      : molec/cm2/s or kg/s : SCALECHEM
!  (2 ) : DRYD-VEL : drydep velocities  : cm/s                : SCALECHEM
!
!  NOTES: 
!  (1 ) Remove diagnostics for wetdep HNO3, H2O2 from ND44.
!  (2 ) For NSRCX == 1 (Rn-Pb-Be), save the actual tracer number 
!        instead of the dry deposition index.  Add TRCOFFSET to N.
!  (3 ) For NSRCX == 6 (single tracer Ox), drydep fluxes are in kg/s.
!  (4 ) ND44 now uses allocatable array AD44 instead of AIJ. (bmy, 3/16/00)
!  (5 ) Add code from amf for multi-tracer Ox (bmy, 7/3/01)
!  (6 ) Now divide by SCALECHEM since DRYFLX is only called after the
!        chemistry routines for all relevant simulations (bmy, 1/27/03)
!  (7 ) Now print out NTRACE drydep fluxes for tagged Ox.  Also tagged Ox 
!        now saves drydep in molec/cm2/s. (bmy, 8/19/03)
!  (8 ) Rearrange ND44 code for clarity (bmy, 3/24/04)
!  (9 ) Add code for H2/HD simulation (phs, 5/8/07)
!******************************************************************************
!
      IF ( ND44 > 0 ) THEN
         
         !==============================================================
         ! Drydep fluxes
         !==============================================================

         ! Category name
         CATEGORY = 'DRYD-FLX'

         ! # of drydep flux tracers
         IF ( ITS_A_TAGOX_SIM() .or. ITS_A_MERCURY_SIM() ) THEN
            M = N_TRACERS
         ELSE
            M = NUMDEP
         ENDIF

         ! Loop over drydep tracers
         DO N = 1, M

            IF ( ITS_A_RnPbBe_SIM() .or. ITS_A_H2HD_SIM() ) THEN

               ! Radon or H2/HD
               UNIT = 'kg/s'       
               NN   = NTRAIND(N)
 
            ELSE IF ( ITS_A_TAGOX_SIM() .or. ITS_A_MERCURY_SIM() ) THEN

               ! Tagged Ox or Tagged Hg
               UNIT = 'molec/cm2/s'
               NN   = N
  
            ELSE 
     
               ! Other simulations
               UNIT = 'molec/cm2/s'
               NN   = NTRAIND(N)
 
            ENDIF

            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(44) )
               MMB = TINDEX(44,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= NN ) CYCLE

            ! Save into ARRAY
            ARRAY(:,:,1) = ( AD44(:,:,N,1) / SCALECHEM )

            ! Write to file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO

         !==============================================================
         ! Drydep velocities
         !==============================================================

         ! Category and Unit 
         CATEGORY  = 'DRYD-VEL'
         UNIT      = 'cm/s'

         ! # of drydep velocity tracers
         IF ( ITS_A_TAGOX_SIM() ) THEN
            M = 1
         ELSE IF ( ITS_A_MERCURY_SIM() ) THEN
            M = 2
         ELSE
            M = NUMDEP
         ENDIF

         ! Loop over drydep tracers
         DO N = 1, M

            NN           = NTRAIND(N)
            ! To output only the species asked in input.geos 
            ! (ccc, 5/15/09)
            MM  = 1
            MMB = 0
            DO WHILE ( MMB /= NN .AND. MM <= TMAX(44) )
               MMB = TINDEX(44,MM)
               MM  = MM + 1
            ENDDO
            
            IF ( MMB /= NN ) CYCLE

            ! Tracer number plus GAMAP offset
            ARRAY(:,:,1) = AD44(:,:,N,2) / SCALESRCE

            ! Write to file
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND45: Tracer Mixing Ratio (v/v) for Levels L = 1, LD45
!        averaged between hours OTH_HR1 and OTH_HR2
!
!   # : Field   : Description            : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) IJ-AVG-$ : Tracer mixing ratio    : v/v     : CTOTH
!
!  NOTES:
!  (1) For NSRCX == 3 (NOx-Ox-HC run), store pure O3 with index NTRACE+1.
!  (2) Now store pure O3 as NNPAR+1 (now tracer #32). (bmy, 1/10/03)
!  (3) Now uses CTO3 instead of CTOH for pure O3 (phs, 1/24/07)
!  (4) Better handling of O3 case (phs, 11/17/08)
!******************************************************************************
!
      IF ( ND45 > 0 ) THEN
         CATEGORY    = 'IJ-AVG-$'
         SCALE_TMP   = XPLX( CTOTH ) + 1d-20
         UNIT        = ''   

         DO M = 1, TMAX(45)
            N  = TINDEX(45,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N

            DO L = 1, LD45
               ARRAY(:,:,L) = AD45(:,:,L,N) / SCALE_TMP(:,:)
            ENDDO
            
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD45,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD45) )

            ! Store pure O3 as NNPAR+1 (bmy, 1/10/03)
            IF ( ITS_A_FULLCHEM_SIM() .and. NN == IDTOX ) THEN 

                  WHERE( CTO3 /= 0 )
                     ARRAY(:,:,1:LD45)%r=AD45(:,:,1:LD45,N_TRACERS+1)%r/
     $                    ( CTO3 )
                     ARRAY(:,:,1:LD45)%i=AD45(:,:,1:LD45,N_TRACERS+1)%i/
     $                    ( CTO3 )
                  ELSEWHERE
                     ARRAY(:,:,1:LD45)%r = 0.d0
                     ARRAY(:,:,1:LD45)%i = 0.d0
                  ENDWHERE
                  
               NN = N_TRACERS + 1

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     LD45,     IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD45) )
            ENDIF
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND46: Biogenic source diagnostic
!
!   # : Field : Description    : Units         : Scale Factor
!  ---------------------------------------------------------------------------
!  (1)  ISOP  : Isoprene       : atoms C/cm2/s : SCALE3
!  (2)  ACET  : Acetone        : atoms C/cm2/s : SCALE3
!  (3)  PRPE  : Propene        : atoms C/cm2/s : SCALE3
!  (4)  MONOT : Monoterpenes   : atoms C/cm2/s : SCALE3
!  (5)  MBO   : Methyl Butenol : atoms C/cm2/s : SCALE3
!  (6)  C2H4  : Ethene         : atoms C/cm2/s : SCALE3
!  
!  NOTES:
!  (1) ND46 now uses allocatable array AD46 instead of AIJ (bmy, 3/16/00)
!  (2) Also write out PRPE for CO-OH run (NSRCX == 5), regardless of
!       the setting of IDTPRPE.  This is to print out monterpene 
!       diagnostics. (bnd, bmy, 4/18/00)
!  (3) Added monoterpenes as tracer #4.  This requires updated versions
!       of "tracerinfo.dat" and "diaginfo.dat" for GAMAP. (bmy, 1/2/01)
!  (4) Added MBO as tracer #5. (tmf, bmy, 10/20/05)
!  (5) Added C2H4 as tracer #6. (tmf, 1/20/09)
!******************************************************************************
!
      IF ( ND46 > 0 ) THEN
         CATEGORY = 'BIOGSRCE'
         UNIT     = ''

         DO M = 1, TMAX(46)
            N  = TINDEX(46,M)
            IF ( N > PD46 ) CYCLE
            NN = N

            ! Skip if ISOP, ACET, PRPE are not tracers
            IF ( N == 1 .and. IDTISOP == 0 ) CYCLE
            IF ( N == 2 .and. IDTACET == 0 ) CYCLE
            IF ( N == 3 .and. IDTPRPE == 0 ) CYCLE
            IF ( N == 6 .and. IDTC2H4 == 0 ) CYCLE
            
            ARRAY(:,:,1) = AD46(:,:,N) / SCALESRCE
               
            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND47: Tracer Mixing Ratio (v/v) for Levels L = 1, LD47
!        *always* averaged between 0000 and 2400 Local time.
!
!   # : Field   : Description                 : Units   : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) IJ-24H-$ : 24h avg Tracer mixing ratio : v/v     : SCALEDYN
!
!  NOTES:
!  (1) For NSRCX == 3 (NOx-Ox-HC run), store pure O3 with index NTRACE+1.
!  (2) Now store pure O3 as NNPAR+1 (now tracer #32). (bmy, 1/10/03) 
!  (3) Now replace SCALE1 with SCALEDYN
!  (4) Now averaged between 0 and 24 UT.  Replace SCALEDYN with CTOH and
!       CTO3 (phs, 1/24/07)
!  (5) Revert to SCALEDYN for all species, except O3, which uses new 
!       CTO3_24h counter (phs, 11/17/08)
!******************************************************************************
!
      IF ( ND47 > 0 ) THEN
         CATEGORY = 'IJ-24H-$'
         UNIT     = ''   

         DO M = 1, TMAX(47)
            N  = TINDEX(47,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N
            
            DO L = 1, LD47
               ARRAY(:,:,L) = AD47(:,:,L,N) / SCALEDYN
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD47,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD47) )

            ! Store pure O3 as NNPAR+1 (bmy, 1/10/03)
            IF ( ITS_A_FULLCHEM_SIM() .and. NN == IDTOX ) THEN 

               WHERE( CTO3_24h(:,:,1:LD47) /= 0 )
                  ARRAY(:,:,1:LD47)%r=AD47(:,:,1:LD47,N_TRACERS+1)%r /
     $                                ( CTO3_24h(:,:,1:LD47) )
                  ARRAY(:,:,1:LD47)%i=AD47(:,:,1:LD47,N_TRACERS+1)%i /
     $                                ( CTO3_24h(:,:,1:LD47) )
               ELSEWHERE
                  ARRAY(:,:,1:LD47)%r = 0.d0
                  ARRAY(:,:,1:LD47)%i = 0.d0
               ENDWHERE
               
               NN = N_TRACERS + 1

               CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                     HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                     UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                     IIPAR,     JJPAR,     LD47,     IFIRST,     
     &                     JFIRST,    LFIRST,    ARRAY(:,:,1:LD47) )
            ENDIF
         ENDDO
      ENDIF

!******************************************************************************
!  ND52: gamma HO2 and aerosol radius (jaegle 02/26/09)
!   # Category 
!   # : Field     : Description                      : Units     : Scale factor
!  ----------------------------------------------------------------------------
!  (1): GAMMAHO2  : Uptake coef for HO2              : unitless  : SCALECHEM
!
!******************************************************************************
      IF ( ND52 > 0 ) THEN
        CATEGORY  = 'GAMMA'
        UNIT      = 'unitless'

        DO L = 1, LD52
            ARRAY(:,:,L) = AD52(:,:,L) / SCALECHEM
        ENDDO

        ! Save to disk
        CALL BPCH2( IU_BPCH,    MODELNAME, LONRES,   LATRES,
     &                    HALFPOLAR,  CENTER180, CATEGORY, 1,
     &                    UNIT,       DIAGb,     DIAGe,    RESERVED,
     &                    IIPAR,      JJPAR,     LD52,     IFIRST,
     &                    JFIRST,     LFIRST,    ARRAY(:,:,1:LD52) )

      ENDIF


!
!******************************************************************************
!  ND54: Time-in-the-Troposphere diagnostic
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) TIME-TPS : Time spend in troposphere   : fraction  : SCALEDYN
!******************************************************************************
!
      IF ( ND54 > 0 ) THEN
         CATEGORY = 'TIME-TPS'
         UNIT = 'unitless'

         DO L = 1, LD54
            ARRAY(:,:,L) = AD54(:,:,L) / SCALEDYN
         ENDDO

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     LD54,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,1:LD54) )
      ENDIF
!
!******************************************************************************
!  ND55: Tropopause diagnostics
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) TP-LEVEL : Tropopause level            : unitless  : SCALEDYN
!  (2) TP-HGHT  : Tropopause height           : km        : SCALEDYN
!  (3) TP-PRESS : Tropopause pressure         : mb        : SCALEDYN
!******************************************************************************
!
      IF ( ND55 > 0 ) THEN
         CATEGORY = 'TR-PAUSE'

         DO M = 1, TMAX(55)
            N  = TINDEX(55,M)
            IF ( N > PD55 ) CYCLE
            NN = N

            ! Pick the appropriate unit string
            SELECT CASE ( N )
               CASE ( 1 )
                  UNIT = 'unitless'
               CASE ( 2 )
                  UNIT = 'km'
               CASE ( 3 )
                  UNIT = 'mb'
            END SELECT

            ARRAY(:,:,1) = AD55(:,:,N) / SCALEDYN

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND56: Lightning flash rate diagnostics (ltm, bmy, 5/5/06))
!******************************************************************************
!
      IF ( ND56 > 0 ) CALL WRITE_DIAG56
!
!******************************************************************************
!  ND58: CH4 Emission Diagnostics
!
!   #  : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1 ) CH4-TOT  : CH4 Emissions total(w/o sab): kg        : 1
!  (2 ) CH4-GAO  : CH4 Emissions gas & oil     : kg        : 1
!  (3 ) CH4-COL  : CH4 Emissions coal          : kg        : 1
!  (4 ) CH4-LIV  : CH4 Emissions livestock     : kg        : 1
!  (5 ) CH4-WST  : CH4 Emissions waste         : kg        : 1
!  (6 ) CH4-BFL  : CH4 Emissions biofuel       : kg        : 1
!  (7 ) CH4-RIC  : CH4 Emissions rice          : kg        : 1
!  (8 ) CH4-OTA  : CH4 Emissions other anthro  : kg        : 1
!  (9 ) CH4-BBN  : CH4 Emissions bioburn       : kg        : 1
!  (10) CH4-WTL  : CH4 Emissions wetlands      : kg        : 1
!  (11) CH4-SAB  : CH4 Emissions soil abs      : kg        : 1
!  (12) CH4-OTN  : CH4 Emissions other natural : kg        : 1
!******************************************************************************
!
      IF ( ND58 > 0 ) THEN
         CATEGORY = 'CH4-EMIS'

         DO M = 1, TMAX(58)
            N  = TINDEX(58,M)
            IF ( N > PD58 ) CYCLE
            NN = N

            UNIT = 'kg'

            ARRAY(:,:,1) = AD58(:,:,N)

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,
     &                  IIPAR,     JJPAR,     1,        IFIRST,
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND59: NH3 concentrations [ug/m3] (diag59 added, lz,10/11/10)
!******************************************************************************
!
      IF ( ND59 > 0 ) CALL WRITE_DIAG59
!
!******************************************************************************
!  ND60: Wetland fraction
!
!   # : Field    : Description                     : Units    : Scale factor
!  --------------------------------------------------------------------------
!  (1 ) WET-FRAC : WETLAND FRACTION                : unitless : 1
!******************************************************************************
!
      IF ( ND60 > 0 ) THEN

         UNIT = 'unitless'
         CATEGORY     = 'WET-FRAC'
         ARRAY(:,:,1) = AD60(:,:)
         N            = 1

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, N,
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,
     &               IIPAR,     JJPAR,     1,        IFIRST,
     &               JFIRST,    LFIRST,    ARRAY(:,:,1) )

      ENDIF
!
!*****************************************************************************
!  ND62: I-J Instantaneous Column Maps for Tracers (molec/cm^2)  
!
!  The unit conversion is as follows:
!
!    STT (kg) | 6.022e23 molec |   mole   | 1000 g |    1        |   m^2
!    ---------+----------------+----------+--------+-------------+----------
!             |     mole       |  MOLWT g |  kg    | AREA_M2 m^2 | 10^4 cm^2
!
!
!  which is equivalent to
!
!   ( STT * 6.022e22 ) / ( MOLWT * AREA_M2 )
!*****************************************************************************
!
      IF ( ND62 > 0 ) THEN
         CATEGORY = 'INST-MAP'

         DO M = 1, TMAX(62)
            N  = TINDEX(62,M)
            IF ( N > N_TRACERS ) CYCLE
            NN = N

            DO J = 1, JJPAR

               ! Grid box surface area [cm2]
               AREA_M2 = GET_AREA_M2( J )

               DO I = 1, IIPAR
                  ARRAY(I,J,1) = ( SUM( STT(I,J,:,N) ) * 6.022d22 )
     &                         / ( TRACER_MW_G(N)       * AREA_M2  )
               ENDDO
            ENDDO

            ! Write the proper unit string
            IF ( TRACER_MW_G(N) > 12d0 ) THEN
               UNIT = 'molec/cm2'
            ELSE
               UNIT = 'atoms C/cm2'
            ENDIF

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )

         ENDDO
      ENDIF
!
!******************************************************************************
!  ND65: Production/Loss of specified chemical families
!
!   # : Field   : Description                 : Units     : Scale Factor
!  ---------------------------------------------------------------------------
!  (1) PORL-L=$ : Chemical family P-L rates   : mol/cm3/s : SCALECHEM
!
!  NOTES:
!  (1 ) Make sure the units for NSRCX == 6 (single tracer O3) P-L 
!        coincide with those in "chemo3.f".  
!  (2 ) ND65 now uses allocatable array AD65 instead of AIJ. (bmy, 3/16/00)
!  (3 ) Add L(CH3I) to the ND65 diagnostic -- do not take the average 
!        but instead compute the total sum of L(CH3I) (nad, bmy, 3/20/01)
!  (4 ) Add updates for multi-tracer Ox run from amf (bmy, 7/3/01)
!  (5 ) Now account for time in troposphere for full chemistry. It is
!        assumed that LD45 >= LD65 in using CTO3 (phs, 3/6/07)
!  (6 ) Do not use CTO3 anymore, but the new CTO3_24h, which is the 3D
!        tropospheric chemistry counter (phs, 11/17/08)
!******************************************************************************
!
      IF ( ND65 > 0 ) THEN     
         CATEGORY = 'PORL-L=$'

         ! Loop over ND65 families
         DO N = 1, NFAMILIES

            ! Don't add TRCOFFSET for single tracer Ox
            ! Also select proper unit string
            IF ( ITS_A_CH3I_SIM() ) THEN
               NN          = N
               UNIT        = 'kg/s'

               DO L = 1, LD65
                  ARRAY(:,:,L) = AD65(:,:,L,N)
               ENDDO

            ELSE IF ( ITS_A_TAGOX_SIM() ) THEN
               NN          = N
               UNIT        = 'kg/s'
               
               WHERE( CTO3_24h(:,:,1:LD65) /= 0 )
                  ARRAY(:,:,1:LD65)%r = AD65(:,:,1:LD65,N)%r /
     $                                ( CTO3_24h(:,:,1:LD65) )
                  ARRAY(:,:,1:LD65)%i = AD65(:,:,1:LD65,N)%i /
     $                                ( CTO3_24h(:,:,1:LD65) )
               ELSEWHERE
                  ARRAY(:,:,1:LD65)%r = 0.d0
                  ARRAY(:,:,1:LD65)%i = 0.d0
               ENDWHERE

            ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
               NN          = N
               UNIT        = 'mol/cm3/s'

               DO L = 1, LD65
                  ARRAY(:,:,L) = AD65(:,:,L,N) / SCALECHEM
               ENDDO

            ELSE
               NN     = N
               UNIT   = 'mol/cm3/s'

               WHERE( CTO3_24h(:,:,1:LD65) /= 0 )
                  ARRAY(:,:,1:LD65)%r = AD65(:,:,1:LD65,N)%r /
     $                                ( CTO3_24h(:,:,1:LD65) )
                  ARRAY(:,:,1:LD65)%i = AD65(:,:,1:LD65,N)%i /
     $                                ( CTO3_24h(:,:,1:LD65) )
               ELSEWHERE
                  ARRAY(:,:,1:LD65)%r = 0.d0
                  ARRAY(:,:,1:LD65)%i = 0.d0
               ENDWHERE
               
            ENDIF

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD65,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD65) )
         ENDDO       
      ENDIF
!
!******************************************************************************
!  ND66: GMAO 3-D fields 
!
!   # : Field  : Description                       : Units   : Scale factor
!  --------------------------------------------------------------------------
!  (1)  UWND   : GMAO Zonal Winds                  : m/s     : SCALE_I6 or _A6
!  (2)  VWND   : GMAO Meridional Winds             : m/s     : SCALE_I6 or _A6
!  (3)  TMPU   : GMAO Temperatures                 : K       : SCALE_I6 or _A6
!  (4)  SPHU   : GMAO Specific Humidity            : g/kg    : SCALE_I6 or _A6
!  (5)  CLDMAS : GMAO Cloud Mass Flux              : kg/m2/s : SCALE_A6 or _A6
!  (6)  DTRAIN : GMAO Detrainment flux             : kg/m2/s : SCALE_A6
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) Add CLDMAS to ND66 diagnostic as field #6, but with tracer index
!       #7 (for compatibility with the existing GAMAP).  (rvm, bmy, 9/8/00)
!  (3) For GEOS-4/fvDAS, UWND, VWND, TMPU, SPHU are A-6 fields.  Adjust
!       the scale factors accordingly.  Also delete KZZ. (bmy, 6/23/03)
!  (4) Modified for GEOS-5 and GCAP (bmy, 6/9/05)
!******************************************************************************
!
      IF ( ND66 > 0 ) THEN
         CATEGORY = 'DAO-3D-$'

         DO M = 1, TMAX(66)
            N  = TINDEX(66,M)
            NN = N 
            
            SELECT CASE ( N )

               ! UWND, VWND
               CASE ( 1,2 )
#if   defined( GEOS_3 )
                  SCALEX = SCALE_I6
#else
                  SCALEX = SCALE_A6
#endif
                  UNIT   = 'm/s'

               ! TMPU
               CASE ( 3 )
#if   defined( GEOS_3 )
                  SCALEX = SCALE_I6
#else
                  SCALEX = SCALE_A6
#endif
                  UNIT   = 'K'

               ! SPHU
               CASE ( 4 )
#if   defined( GEOS_3 )
                  SCALEX = SCALE_I6
#else
                  SCALEX = SCALE_A6
#endif
                  UNIT   = 'g/kg'

               ! CLDMAS, DTRAIN
               CASE( 5, 6 )
                  SCALEX = SCALE_A6
                  UNIT   = 'kg/m2/s'

               CASE DEFAULT
                  CYCLE
            END SELECT

            ARRAY(:,:,1:LD66) = AD66(:,:,1:LD66,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD66,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD66) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND67: GMAO surface fields
!
!   # : Field  : Description                      : Units    : Scale factor
!  -----------------------------------------------------------------------
!  (1 ) HFLUX  : GMAO Sensible Heat Flux          : W/m2     : SCALE_A3
!  (2 ) RADSWG : GMAO Insolation @ Surface        : W/m2     : SCALE_A3
!  (3 ) PREACC : GMAO Accum Precip @ Surface      : mm/day   : SCALE_A3
!  (4 ) PRECON : GMAO Conv Precip @ Surface       : mm/day   : SCALE_A3
!  (5 ) TS     : GMAO Surface Air Temperature     : K        : SCALE_A3
!  (6 ) RADSWT : GMAO Insolation @ Top of Atm     : W/m2     : SCALE_A3
!  (7 ) USTAR  : GMAO Friction Velocity           : m/s      : SCALE_A3
!  (8 ) Z0     : GMAO Roughness Height            : m        : SCALE_A3
!  (9 ) PBL    : GMAO PBL depth                   : mb       : SCALE_A3
!  (10) CLDFRC : GMAO Cloud Fraction              : unitless : SCALE_A3
!  (11) U10M   : GMAO U-wind @ 10 m               : m/s      : SCALE_A3
!  (12) V10M   : GMAO V-wind @ 10 m               : m/s      : SCALE_A3
!  (13) PS-PBL : GMAO Boundary Layer Top Pressure : mb       : SCALEDYN
!  (14) ALBD   : GMAO Surface Albedo              : unitless : SCALE_I6 
!  (15) PHIS   : GMAO Geopotential Heights        : m        : SCALED 
!  (16) CLTOP  : GMAO Cloud Top Height            : levels   : SCALE_A6
!  (17) TROPP  : GMAO Tropopause pressure         : mb       : SCALE_I6
!  (18) SLP    : GMAO Sea Level pressure          : mb       : SCALE_I6
!  (19) TSKIN  : Ground/sea surface temp.         : hPa      : SCALE_A3
!  (20) PARDF  : Photosyn active diffuse rad.     : W/m2     : SCALE_A3
!  (21) PARDR  : Photosyn active direct  rad.     : W/m2     : SCALE_A3
!  (22) GWET   : Top soil wetness                 : unitless : SCALE_A3
!
!  NOTES:
!  (1 ) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2 ) Now use AD67 allocatable array (bmy, 2/17/00)
!  (3 ) Add TROPP as tracer #17 and SLP as tracer #18 (bmy, 10/11/00)
!  (4 ) Now replace SCALE1 with SCALEDYN (bmy, 3/27/03)
!  (5 ) Added TSKIN, PARDF, PARDR, GWET for GEOS-4 (bmy, 6/23/03)
!  (6 ) Fix SCALEX for ALBEDO: use I6 for GEOS-3 only, and A3 for other
!     models (phs, 9/3/08)
!******************************************************************************
!
      IF ( ND67 > 0 ) THEN
         CATEGORY = 'DAO-FLDS'

         ! Binary punch file
         DO M = 1, TMAX(67)
            N  = TINDEX(67,M)
            NN = N 

            SELECT CASE ( N ) 
               CASE ( 1, 2, 6 )
                  SCALEX = SCALE_A3
                  UNIT   = 'W/m2'
               CASE ( 3, 4 )
                  SCALEX = SCALE_A3
                  UNIT   = 'mm/day'
               CASE ( 5 )
                  SCALEX = SCALE_A3
                  UNIT   = 'K'
               CASE ( 7, 11, 12 )
                  SCALEX = SCALE_A3
                  UNIT   = 'm/s'
               CASE ( 8 )
                  SCALEX = SCALE_A3
                  UNIT   = 'm'
               CASE ( 9 )
                  SCALEX = SCALE_A3
                  UNIT   = 'hPa'
               CASE ( 10 )
                  SCALEX = SCALE_A3
                  UNIT   = 'unitless'

#if   defined( GCAP )
                  ! CLDFRC is a 6-hr field in GCAP, GEOS-STRAT 
                  ! (swu, bmy, 6/9/05)
                  SCALEX = SCALE_A6
#endif

               CASE ( 13 )
                  SCALEX = SCALEDYN
                  UNIT   = 'hPa'
               CASE ( 14 ) 
                  ! Bug fix: For GEOS-3, ALBEDO is an I-6 field, but
                  ! for GEOS-4, GEOS-5, GCAP, it is an A-3 field.
                  ! (lyj, phs, bmy, 10/7/08)
#if   defined( GEOS_3 )
                  SCALEX = SCALE_I6 
#else
                  SCALEX = SCALE_A3
#endif
                  UNIT   = 'unitless'
               CASE ( 15 )
                  SCALEX = SCALED
                  UNIT   = 'm'
               CASE ( 16 ) 
                  SCALEX = SCALE_A6
                  UNIT   = 'levels'
               CASE ( 17 ) 
                  SCALEX = SCALE_I6
                  UNIT   = 'hPa'
               CASE ( 18 ) 
                  SCALEX = SCALE_I6
                  UNIT   = 'hPa'
               CASE ( 19 )
                  SCALEX = SCALE_A3
                  UNIT   = 'K'
               CASE ( 20 ) 
                  SCALEX = SCALE_A3
                  UNIT   = 'W/m2'
               CASE ( 21 ) 
                  SCALEX = SCALE_A3
                  UNIT   = 'W/m2'
               CASE ( 22 ) 
                  SCALEX = SCALE_A3
                  UNIT   = 'unitless'
               CASE DEFAULT
                  CYCLE
            END SELECT
                  
            ARRAY(:,:,1) = AD67(:,:,N) / SCALEX

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     1,        IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND68: Grid box diagnostics
!
!   # : Field   : Description                       : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1) BXHEIGHT : Grid box height                   : m     : SCALEDYN
!  (2) AD       : Air mass in grid box              : kg    : SCALEDYN
!  (3) AVGW     : Mixing ratio of water vapor       : v/v   : SCALEDYN
!  (4) N(AIR)   : Number density of air             : m^-3  : SCALEDYN
!
!  NOTES:
!  (1) We don't need to add TRCOFFSET to N.  These are not CTM tracers.
!  (2) Now replaced SCALE1 with SCALEDYN (bmy, 2/24/03)
!  (3) Bug fix: replace ND68 with LD68 in call to BPCH2 (swu, bmy, 6/9/05)
!******************************************************************************
!
      IF ( ND68 > 0 ) THEN
         CATEGORY = 'BXHGHT-$'
         UNIT     = ''

         DO M = 1, TMAX(68)
            N  = TINDEX(68,M)
            IF ( N > PD68 ) CYCLE
            NN = N 

            DO L = 1, LD68
               ARRAY(:,:,L) = AD68(:,:,L,N) / SCALEDYN
            ENDDO

            CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY, NN,    
     &                  UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &                  IIPAR,     JJPAR,     LD68,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,1:LD68) )
         ENDDO
      ENDIF
!
!******************************************************************************
!  ND69: Grid Box Surface Areas
!
!   # : Field : Description                       : Units : Scale factor
!  --------------------------------------------------------------------------
!  (1) DXYP   : Surface area of grid box          : m^2   : SCALED = 1.0
!
!  NOTES:
!  (1) Only print DXYP for the first timestep, as it is an invariant field.
!  (2) We don't need to add TRCOFFSET to N.  This is not a CTM tracer.
!  (3) Now use the AD69 dynamically allocatable array (bmy, 2/17/00)
!******************************************************************************
!
      IF ( ND69 > 0 ) THEN 
         CATEGORY = 'DXYP'
         UNIT     = 'm2'

         CALL BPCH2( IU_BPCH,   MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, 1,    
     &               UNIT,      DIAGb,     DIAGe,    RESERVED,   
     &               IIPAR,     JJPAR,     1,        IFIRST,     
     &               JFIRST,    LFIRST,    AD69(:,:,1) )

         ! Set ND69 = 0 so we won't print it out again
         ND69 = 0
      ENDIF

      ! Echo output
      WRITE( 6, '(a)' ) '     - DIAG3: Diagnostics written to bpch!'

      ! Return to calling program
      END SUBROUTINE DIAG3    
