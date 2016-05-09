! $Id: diag_mod.f,v 1.2 2012/03/01 22:00:26 daven Exp $
      MODULE DIAG_MOD 
!
!******************************************************************************
!  Module DIAG_MOD contains declarations for allocatable arrays for use with 
!  GEOS-CHEM diagnostics. (amf, bdf, bmy, 11/30/99, 11/18/08)
!
!  Module Routines:
!  ============================================================================
!  (1 ) CLEANUP_DIAG : Deallocates all module arrays  
!
!  GEOS-CHEM modules referenced by diag_mod.f
!  ============================================================================
!  none
!
!  NOTES:
!  (1 ) DIAG_MOD is written in Fixed-Format F90.
!  (2 ) Call subroutine CLEANUP at the end of the MAIN program to deallocate
!        the memory before the run stops.  It is always good style to free
!        any memory we have dynamically allocated when we don't need it
!        anymoren
!  (3 ) Added ND13 arrays for sulfur emissions (bmy, 6/6/00)
!  (4 ) Moved ND51 arrays to "diag51_mod.f" (bmy, 11/29/00)
!  (5 ) Added AD34 array for biofuel burning emissions (bmy, 3/15/01)
!  (6 ) Eliminated old commented-out code (bmy, 4/20/01)
!  (7 ) Added AD12 array for boundary layer emissions in routine "setemis.f".
!        (bdf, bmy, 6/15/01)
!  (8 ) Added CHEML24, DRYDL24, CTCHDD for archiving daily mean chemical
!        and drydep loss in chemo3 and chemo3.f (amf, bmy, 7/2/01)
!  (9 ) Add ND43 arrays LTNO2, CTNO2, LTHO2, CTHO2 (rvm, bmy, 2/27/02)
!  (10) Add AD01, AD02 arrays for Rn-Pb-Be simulation (hyl, bmy, 8/7/02)
!  (11) Add AD05 array for sulfate P-L diagnostic (rjp, bdf, bmy, 9/20/02)
!  (12) Added subroutine CLEANUP_DIAG...moved code here from "cleanup.f", 
!        so that it is internal to "diag_mod.f".  Added arrays AD13_NH3_bb,
!        AD13_NH3_bf, AD13_NH3_an for NH3 emissons in ND13.  Deleted obsolete
!        allocatable arrays CHEML24, DRYDL24, CTCHDD.  Now also added LTNO3
!        and CTNO3 arrays for ND43 diagnostic.  Added AD13_SO2_bf array for
!        SO2 biofuel. (bmy, 1/16/03)
!  (13) Added array AD13_NH3_na for ND13 diagnostic (rjp, bmy, 3/23/03)
!  (14) Removed P24H and L24H -- these are now defined w/in "tagged_ox_mod.f"
!        Also added AD03 array for Kr85 prod/loss diag. (jsw, bmy, 8/20/03)
!  (15) Added ND06 (dust emission) and ND07 (carbon aerosol emission) 
!        diagnostic arrays (rjp, tdf, bmy, 4/5/04)
!  (16) Added AD13_SO2_sh diagnostic array for ND13 (bec, bmy, 5/20/04)
!  (17) Added AD07_HC diagnostic array for ND07 (rjp, bmy, 7/13/04)
!  (18) Moved AD65 & FAMPL to "diag65_mod.f" (bmy, 7/20/04)
!  (19) Added array AD13_SO4_bf (bmy, 11/17/04)!
!  (20) Added extra arrays for ND03 mercury diagnostics (eck, bmy, 12/7/04)
!  (21) Added extra ND21 array for crystalline sulfur tracers.  Also remove
!        ND03 and ND48 arrays; they are obsolete (bmy, 1/21/05)
!  (22) Removed AD41 and AFTTOT arrays; they're obsolete (bmy, 2/17/05)
!  (23) Added AD09, AD09_em arrays for HCN/CH3CN simulation (xyp, bmy, 6/27/05)
!  (24) Added AD30 array for land/water/ice output (bmy, 8/18/05)
!  (25) Added AD54 array for time spend in the troposphere (phs, 9/22/06)
!  (26) Added CTO3 counter. Convert ND43 counter arrays from 2D to 3D, for
!        the variable tropopause. (phs, 1/19/07)
!  (27) Added AD10 and AD10em arrays for ND10 H2-HD-sim diag (phs, 9/18/07)
!  (28) Added CTO3_24h to account for time in the troposphere for O3 in
!        ND47 (phs, 11/17/08)
!  (29) Added AD52 for Gamma HO2 diagnostic. (jaegle, ccc, 2/26/09)
!  (30) Updated to save out GLYX production of SOAG in ND07.
!       (tmf, 3/6/09)
!******************************************************************************
!     
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      USE MYTYPE
      ! For ND01 -- Rn, Pb, Be emissions
      TYPE (XPLEX),  ALLOCATABLE :: AD01(:,:,:,:)

      ! For ND02 -- Rn, Pb, Be decay
      TYPE (XPLEX),  ALLOCATABLE :: AD02(:,:,:,:)

      !--------------------------------------------
      !! For ND03 -- Kr85 prod/loss
      !TYPE (XPLEX),  ALLOCATABLE :: AD03(:,:,:,:)
      !--------------------------------------------

      ! For ND05 -- Sulfate prod/loss diagnostics
      TYPE (XPLEX),  ALLOCATABLE :: AD05(:,:,:,:)

      ! For ND06 -- Dust aerosol emission
      TYPE (XPLEX),  ALLOCATABLE :: AD06(:,:,:)

      ! For ND07 -- Carbon aerosol emission
      TYPE (XPLEX),  ALLOCATABLE :: AD07(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD07_BC(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD07_OC(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD07_HC(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD07_SOAGM(:,:,:,:)

      ! For ND08 -- seasalt emission
      TYPE (XPLEX),  ALLOCATABLE :: AD08(:,:,:)

      ! For ND09 -- HCN / CH3CN simulation
      TYPE (XPLEX),  ALLOCATABLE :: AD09(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD09_em(:,:,:)

      ! For ND10 -- H2/HD prod, loss, & emiss diagnostics
      TYPE (XPLEX),  ALLOCATABLE :: AD10(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD10em(:,:,:)

      ! For ND12 -- boundary layer multiplication factor
      TYPE (XPLEX),  ALLOCATABLE :: AD11(:,:,:)

      ! For ND12 -- boundary layer multiplication factor
      TYPE (XPLEX),  ALLOCATABLE :: AD12(:,:,:)

      ! For ND13 -- Sulfur emissions
      TYPE (XPLEX),  ALLOCATABLE :: AD13_DMS(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO2_ac(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO2_an(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO2_bb(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO2_bf(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO2_nv(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO2_ev(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO2_sh(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO4_an(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_SO4_bf(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_NH3_an(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_NH3_na(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_NH3_bb(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD13_NH3_bf(:,:)

      ! For ND14 -- wet convection mass flux diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: CONVFLUP(:,:,:,:)

      ! For ND15 -- BL mixing mass flux diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: TURBFLUP(:,:,:,:)

      ! For ND16 -- Fraction of grid box that is precipitating
      TYPE (XPLEX),  ALLOCATABLE :: AD16(:,:,:,:)  
      INTEGER, ALLOCATABLE :: CT16(:,:,:,:)
      
      ! For ND17 -- Fraction of tracer lost to rainout 
      TYPE (XPLEX),  ALLOCATABLE :: AD17(:,:,:,:,:)   
      INTEGER, ALLOCATABLE :: CT17(:,:,:,:)

      ! For ND18 -- Fraction of tracer lost to washout
      TYPE (XPLEX),  ALLOCATABLE :: AD18(:,:,:,:,:)   
      INTEGER, ALLOCATABLE :: CT18(:,:,:,:)

      ! For ND21 -- Optical Depth diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD21(:,:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD21_cr(:,:,:)

      ! For ND22 -- J-value diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD22(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTJV(:,:)
      INTEGER, ALLOCATABLE :: CTJV(:,:) 

      ! For ND23 -- CH3CCl3 lifetime diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: DIAGCHLORO(:,:,:,:)

      ! For ND24 -- E/W transport mass flux diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: MASSFLEW(:,:,:,:)

      ! For ND25 -- N/S transport mass flux diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: MASSFLNS(:,:,:,:)

      ! For ND26 -- UP/DOWN transport mass flux diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: MASSFLUP(:,:,:,:)

      ! For ND28 -- Biomass burning diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD28(:,:,:)

      ! For ND29 -- CO source diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD29(:,:,:)

      ! For ND30 -- land / water / ice flags
      TYPE (XPLEX),  ALLOCATABLE :: AD30(:,:)

      ! For ND31 -- surface pressures
      TYPE (XPLEX),  ALLOCATABLE :: AD31(:,:,:)

      ! For ND32 -- NOx sources 
      TYPE (XPLEX),  ALLOCATABLE :: AD32_ac(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD32_an(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD32_bb(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD32_bf(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD32_fe(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD32_li(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD32_so(:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD32_ub(:,:)

      ! For ND33 -- tropopsheric sum of tracer
      TYPE (XPLEX),  ALLOCATABLE :: AD33(:,:,:)

      ! For ND34 -- biofuel emissions
      TYPE (XPLEX),  ALLOCATABLE :: AD34(:,:,:)

      ! For ND35 -- 500 mb tracer
      TYPE (XPLEX),  ALLOCATABLE :: AD35(:,:,:)

      ! For ND36 -- Anthropogenic source diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD36(:,:,:)

      ! For ND37 -- Fraction of tracer scavenged in cloud updrafts
      TYPE (XPLEX),  ALLOCATABLE :: AD37(:,:,:,:)      

      ! For ND38 -- Rainout in moist convection diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD38(:,:,:,:)      

      ! For ND39 -- Washout in aerosol wet deposition diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD39(:,:,:,:)      

      ! For ND43 -- OH, NO, NO2, HO2 chemical diagnostics
      TYPE (XPLEX),  ALLOCATABLE :: AD43(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTNO(:,:)
      INTEGER, ALLOCATABLE :: CTNO(:,:,:)
      INTEGER, ALLOCATABLE :: LTOH(:,:)
      INTEGER, ALLOCATABLE :: CTOH(:,:,:)
      INTEGER, ALLOCATABLE :: LTNO2(:,:)
      INTEGER, ALLOCATABLE :: CTNO2(:,:,:)
      INTEGER, ALLOCATABLE :: LTHO2(:,:)
      INTEGER, ALLOCATABLE :: CTHO2(:,:,:)
      INTEGER, ALLOCATABLE :: LTNO3(:,:)
      INTEGER, ALLOCATABLE :: CTNO3(:,:,:)

      ! For ND44 -- Dry deposition fluxes & velocities
      TYPE (XPLEX),  ALLOCATABLE :: AD44(:,:,:,:)

      ! For ND45 -- Tracer concentration diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD45(:,:,:,:)      
      INTEGER, ALLOCATABLE :: LTOTH(:,:)
      INTEGER, ALLOCATABLE :: CTOTH(:,:)
      INTEGER, ALLOCATABLE :: CTO3(:,:,:)

      ! For ND46 -- Tracer concentration diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD46(:,:,:)      

      ! For ND47 -- 24-h tracer concentration diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD47(:,:,:,:)      

      ! For ND47(O3) / ND65 -- 24-h tracer diagnostic
      INTEGER, ALLOCATABLE :: CTO3_24h(:,:,:)

      ! Dynamically allocatable array -- local only to DIAG50.F
      TYPE (XPLEX),  ALLOCATABLE :: STT_TEMPO2(:,:,:,:)

      ! For ND52 -- gamma HO2 diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD52(:,:,:)

      ! For ND54 -- tropopause diagnostics
      TYPE (XPLEX),  ALLOCATABLE :: AD54(:,:,:)

      ! For ND55 -- tropopause diagnostics
      TYPE (XPLEX),  ALLOCATABLE :: AD55(:,:,:)

      ! -- for methane simulation diagnostics
      ! (kjw, dkh, 02/12/12, adj32_023) 
      TYPE (XPLEX),  ALLOCATABLE :: AD19(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD58(:,:,:)
      TYPE (XPLEX),  ALLOCATABLE :: AD60(:,:)

      ! For ND66 -- I-6 fields diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD66(:,:,:,:)      

      ! For ND67 -- DAO surface fields diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD67(:,:,:)      

      ! For ND68 -- BXHEIGHT, AD, AVGW diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD68(:,:,:,:)      

      ! For ND69 -- DXYP diagnostic
      TYPE (XPLEX),  ALLOCATABLE :: AD69(:,:,:)      

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG
!
!******************************************************************************
!  Subroutine CLEANUP_DIAG deallocates all module arrays.
!  (bmy, 12/13/02, 9/18/07)
!
!  NOTES:
!  (1 ) Now also deallocate AD13_NH3_an, AD13_NH3_bb, AD13_NH3_bf arrays
!        for the ND13 diagnostic.  (bmy, 12/13/02)
!  (2 ) Now also deallocate AD13_NH3_na array for ND13 (rjp, bmy, 3/23/03)
!  (3 ) Removed P24H and L24H, these are now defined within "tagged_ox_mod.f".
!       Now also deallocate AD03 array for Kr85 prod/loss (jsw, bmy, 8/20/03)
!  (4 ) Now also deallocate AD06 and AD07* arrays (rjp, bdf, bmy, 4/5/04)
!  (5 ) Now also deallocate AD08 array (rjp, bec, bmy, 4/20/04)
!  (6 ) Now also deallocaes AD13_SO2_sh array (bec, bmy, 5/20/04)
!  (7 ) Now also deallocates AD07_HC array (rjp, bmy, 7/13/04)
!  (8 ) Now also deallocate AD13_SO4_bf array (bmy, 11/17/04)
!  (9 ) Now deallocate extra arrays for ND03 diagnostics (eck, bmy, 12/7/04)
!  (10) Now deallocates AD21_cr array.  Remove reference to arrays for ND03
!        and ND48 diagnostics, they're obsolete. (cas, sas, bmy, 1/21/05)
!  (11) Removed AD41 and AFTTOT arrays; they're obsolete (bmy, 2/17/05)
!  (12) Now also deallocate AD09 and AD09_em (bmy, 6/27/05)
!  (13) Now deallocate AD30 (bmy, 8/18/05)
!  (14) Now deallocate CTO3, AD10, AD10em arrays (phs, 9/18/07)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_DIAG begins here!
      !=================================================================
      IF ( ALLOCATED( AD01        ) ) DEALLOCATE( AD01        )
      IF ( ALLOCATED( AD02        ) ) DEALLOCATE( AD02        )
      IF ( ALLOCATED( AD06        ) ) DEALLOCATE( AD06        )
      IF ( ALLOCATED( AD07        ) ) DEALLOCATE( AD07        )
      IF ( ALLOCATED( AD07_BC     ) ) DEALLOCATE( AD07_BC     )
      IF ( ALLOCATED( AD07_OC     ) ) DEALLOCATE( AD07_OC     )
      IF ( ALLOCATED( AD07_HC     ) ) DEALLOCATE( AD07_HC     )
      IF ( ALLOCATED( AD07_SOAGM  ) ) DEALLOCATE( AD07_SOAGM  )
      IF ( ALLOCATED( AD08        ) ) DEALLOCATE( AD08        )
      IF ( ALLOCATED( AD09        ) ) DEALLOCATE( AD09        )
      IF ( ALLOCATED( AD09_em     ) ) DEALLOCATE( AD09_em     )
      IF ( ALLOCATED( AD10        ) ) DEALLOCATE( AD10        )
      IF ( ALLOCATED( AD10em      ) ) DEALLOCATE( AD10em      )
      IF ( ALLOCATED( AD11        ) ) DEALLOCATE( AD11        )
      IF ( ALLOCATED( AD12        ) ) DEALLOCATE( AD12        )
      IF ( ALLOCATED( AD13_DMS    ) ) DEALLOCATE( AD13_DMS    )
      IF ( ALLOCATED( AD13_SO2_ac ) ) DEALLOCATE( AD13_SO2_ac )
      IF ( ALLOCATED( AD13_SO2_an ) ) DEALLOCATE( AD13_SO2_an )
      IF ( ALLOCATED( AD13_SO2_bb ) ) DEALLOCATE( AD13_SO2_bb )
      IF ( ALLOCATED( AD13_SO2_bf ) ) DEALLOCATE( AD13_SO2_bf )
      IF ( ALLOCATED( AD13_SO2_nv ) ) DEALLOCATE( AD13_SO2_nv )
      IF ( ALLOCATED( AD13_SO2_ev ) ) DEALLOCATE( AD13_SO2_ev )
      IF ( ALLOCATED( AD13_SO2_sh ) ) DEALLOCATE( AD13_SO2_sh )
      IF ( ALLOCATED( AD13_SO4_an ) ) DEALLOCATE( AD13_SO4_an )
      IF ( ALLOCATED( AD13_SO4_bf ) ) DEALLOCATE( AD13_SO4_bf )
      IF ( ALLOCATED( AD13_NH3_an ) ) DEALLOCATE( AD13_NH3_an )
      IF ( ALLOCATED( AD13_NH3_na ) ) DEALLOCATE( AD13_NH3_na )
      IF ( ALLOCATED( AD13_NH3_bb ) ) DEALLOCATE( AD13_NH3_bb )
      IF ( ALLOCATED( AD13_NH3_bf ) ) DEALLOCATE( AD13_NH3_bf )
      IF ( ALLOCATED( AD16        ) ) DEALLOCATE( AD16        )
      IF ( ALLOCATED( AD17        ) ) DEALLOCATE( AD17        )
      IF ( ALLOCATED( AD18        ) ) DEALLOCATE( AD18        )
      IF ( ALLOCATED( AD21        ) ) DEALLOCATE( AD21        )
      IF ( ALLOCATED( AD21_cr     ) ) DEALLOCATE( AD21_cr     )
      IF ( ALLOCATED( AD22        ) ) DEALLOCATE( AD22        ) 
      IF ( ALLOCATED( AD28        ) ) DEALLOCATE( AD28        ) 
      IF ( ALLOCATED( AD29        ) ) DEALLOCATE( AD29        ) 
      IF ( ALLOCATED( AD30        ) ) DEALLOCATE( AD30        ) 
      IF ( ALLOCATED( AD31        ) ) DEALLOCATE( AD31        ) 
      IF ( ALLOCATED( AD32_ac     ) ) DEALLOCATE( AD32_ac     ) 
      IF ( ALLOCATED( AD32_an     ) ) DEALLOCATE( AD32_an     )
      IF ( ALLOCATED( AD32_bb     ) ) DEALLOCATE( AD32_bb     )
      IF ( ALLOCATED( AD32_bf     ) ) DEALLOCATE( AD32_bf     )
      IF ( ALLOCATED( AD32_fe     ) ) DEALLOCATE( AD32_fe     )
      IF ( ALLOCATED( AD32_li     ) ) DEALLOCATE( AD32_li     )
      IF ( ALLOCATED( AD32_so     ) ) DEALLOCATE( AD32_so     )
      IF ( ALLOCATED( AD32_ub     ) ) DEALLOCATE( AD32_ub     )
      IF ( ALLOCATED( AD33        ) ) DEALLOCATE( AD33        )
      IF ( ALLOCATED( AD34        ) ) DEALLOCATE( AD34        )
      IF ( ALLOCATED( AD35        ) ) DEALLOCATE( AD35        )
      IF ( ALLOCATED( AD36        ) ) DEALLOCATE( AD36        )
      IF ( ALLOCATED( AD37        ) ) DEALLOCATE( AD37        )
      IF ( ALLOCATED( AD38        ) ) DEALLOCATE( AD38        )  
      IF ( ALLOCATED( AD39        ) ) DEALLOCATE( AD39        )
      IF ( ALLOCATED( AD43        ) ) DEALLOCATE( AD43        )
      IF ( ALLOCATED( AD44        ) ) DEALLOCATE( AD44        )
      IF ( ALLOCATED( AD45        ) ) DEALLOCATE( AD45        )
      IF ( ALLOCATED( AD46        ) ) DEALLOCATE( AD46        )
      IF ( ALLOCATED( AD47        ) ) DEALLOCATE( AD47        )
      IF ( ALLOCATED( AD52        ) ) DEALLOCATE( AD52        )
      IF ( ALLOCATED( AD54        ) ) DEALLOCATE( AD54        )
      IF ( ALLOCATED( AD55        ) ) DEALLOCATE( AD55        )
      IF ( ALLOCATED( AD19        ) ) DEALLOCATE( AD19        )
      IF ( ALLOCATED( AD58        ) ) DEALLOCATE( AD58        )
      IF ( ALLOCATED( AD60        ) ) DEALLOCATE( AD60        )
      IF ( ALLOCATED( AD66        ) ) DEALLOCATE( AD66        )
      IF ( ALLOCATED( AD68        ) ) DEALLOCATE( AD68        )
      IF ( ALLOCATED( AD69        ) ) DEALLOCATE( AD69        )
      IF ( ALLOCATED( CONVFLUP    ) ) DEALLOCATE( CONVFLUP    )
      IF ( ALLOCATED( CT16        ) ) DEALLOCATE( CT16        )
      IF ( ALLOCATED( CT17        ) ) DEALLOCATE( CT17        )
      IF ( ALLOCATED( CT18        ) ) DEALLOCATE( CT18        )
      IF ( ALLOCATED( CTJV        ) ) DEALLOCATE( CTJV        )
      IF ( ALLOCATED( CTNO        ) ) DEALLOCATE( CTNO        )
      IF ( ALLOCATED( CTO3        ) ) DEALLOCATE( CTO3        )
      IF ( ALLOCATED( CTO3_24h    ) ) DEALLOCATE( CTO3_24h    )
      IF ( ALLOCATED( CTOH        ) ) DEALLOCATE( CTOH        )
      IF ( ALLOCATED( CTNO2       ) ) DEALLOCATE( CTNO2       )
      IF ( ALLOCATED( CTNO3       ) ) DEALLOCATE( CTNO3       )
      IF ( ALLOCATED( CTHO2       ) ) DEALLOCATE( CTHO2       )
      IF ( ALLOCATED( CTOTH       ) ) DEALLOCATE( CTOTH       )
      IF ( ALLOCATED( DIAGCHLORO  ) ) DEALLOCATE( DIAGCHLORO  )
      IF ( ALLOCATED( LTJV        ) ) DEALLOCATE( LTJV        )
      IF ( ALLOCATED( LTNO        ) ) DEALLOCATE( LTNO        )
      IF ( ALLOCATED( LTOH        ) ) DEALLOCATE( LTOH        )
      IF ( ALLOCATED( LTNO2       ) ) DEALLOCATE( LTNO2       )
      IF ( ALLOCATED( LTNO3       ) ) DEALLOCATE( LTNO3       )
      IF ( ALLOCATED( LTHO2       ) ) DEALLOCATE( LTHO2       )
      IF ( ALLOCATED( LTOTH       ) ) DEALLOCATE( LTOTH       )
      IF ( ALLOCATED( MASSFLEW    ) ) DEALLOCATE( MASSFLEW    )
      IF ( ALLOCATED( MASSFLNS    ) ) DEALLOCATE( MASSFLNS    )
      IF ( ALLOCATED( MASSFLUP    ) ) DEALLOCATE( MASSFLUP    )
      IF ( ALLOCATED( TURBFLUP    ) ) DEALLOCATE( TURBFLUP    )
      IF ( ALLOCATED( STT_TEMPO2  ) ) DEALLOCATE( STT_TEMPO2  ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG
      
!------------------------------------------------------------------------------

      END MODULE DIAG_MOD 

