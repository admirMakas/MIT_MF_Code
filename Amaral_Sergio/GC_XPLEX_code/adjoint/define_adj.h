!$Id: define_adj.h,v 1.11 2012/04/13 22:29:08 nicolas Exp $
!******************************************************************************
!  Include file "define_adj.h" specifies C-preprocessor "switches" that are 
!  used to include or exclude certain sections of ADJOINT code, mostly for
!  controlling observation datasets used. The reason they are pre-processor
!  switches instead of logical flags is so that we can omit the code which
!  requires installation of hdf libaries and such. All are independent of
!  each other, but not of simulation and tracer type  
!  (adj_group, 6/08/09)
!
!  List of "Switches"
!  ===========================================================================
!  (1 ) TES_NH3_OBS         : Use NH3 data from TES
!  (2 ) PM_ATTAINMENT       : Compute PM attainment
!  (3 ) SOMO35_ATTAINMENT   : Compute ozone attainment
!  (4 ) SCIA_KNMI_NO2_OBS   : Use NO2 obs from SCIA KNMI retrieval
!  (5 ) IMPROVE_SO4_NIT_OBS : Use sulfate-nitrate from IMPROVE network
!  (6 ) CASTNET_NH4_OBS     : Use amonia from CASTNET network
!  (7 ) TES_O3_OBS          : Use O3 obs from TES
!  (8 ) SCIA_DAL_NO2_OBS    : Use NO2 obs from SCIA Dalhousie retrieval
!  (9 ) SCIA_DAL_SO2_OBS    : Use SO2 obs from SCIA Dalhousie retrieval
!  (10) MOPITT_V3_CO_OBS    : Use v3 CO obs from MOPITT 
!  (11) MOPITT_V4_CO_OBS    : Use v4 CO obs from MOPITT 
!  (12) SCIA_BRE_CO_OBS     : Use CO obs from SCIA Bremen retrieval
!  (13) AIRS_CO_OBS         : Use CO obs from AIRS (UMBC) retrieval
!  (14) PSEUDO_OBS          : Generate pseudo obs if no data selected
!  (15) LOG_OPT             : Optimized log of scaling factors
!  (16) SOMO35_ATTAINMENT   : Ozone attainment 
!  (17) PM_ATTAINMENT       : PM    attainment 
!  (18) LIDORT              : Online radiative forcing calculations 
!  (19) GOSAT_CO2_OBS       : Use CO2 obs from GOSAT retrieval 
!  (20) MODIS_AOD_OBS       : Use AOD obs from MODIS 
!  (21) IMPROVE_BC_OC_OBS   : Use BC and OC aerosol obs from IMPROVE  
!  (22) MOPITT_V5_CO_OBS    : Use v5 CO obs form MOPITT 
!  (23) TES_O3_IRK          : Use radiative kernels for TES O3
!
! NOTES:
! (1 )  Replace MOPITT_IR_CO_OBS with MOPITT_V3_CO_OBS and MOPITT_V4_CO_OBS
!       (zhe, dkh, 02/04/11) 
! (2 ) Add MODIS_AOD_OBS (xxu, dkh, 01/09/12, adj32_011) 
! (3 ) Addd IMPROVE_BC_OC_OBS (yhmao, dkh/ 01/16/12, adj32_013) 
! (4 ) Add MOPITT_V5_CO_OBS (zhej, dkh, 01/16/12, adj32_016) 
! (5 ) Add CH4 obs operators (kjw, dkh, 02/12/12, adj32_023) 
!******************************************************************************
!
!==============================================================================
! Undefine all "switches" so that they cannot be accidentally reset  
!==============================================================================

#undef TES_NH3_OBS
#undef PM_ATTAINMENT
#undef SOMO35_ATTAINMENT
#undef SCIA_KNMI_NO2_OBS
#undef IMPROVE_SO4_NIT_OBS
#undef CASTNET_NH4_OBS
#undef TES_O3_OBS
#undef TES_O3_IRK
#undef SCIA_DAL_NO2_OBS
#undef SCIA_DAL_SO2_OBS
#undef MOPITT_V3_CO_OBS
#undef MOPITT_V4_CO_OBS
#undef SCIA_BRE_CO_OBS
#undef AIRS_CO_OBS
#undef GOSAT_CO2_OBS
#undef PM_ATTAINMENT
#undef SOMO35_ATTAINMNET
#undef PSEUDO_OBS
#undef LOG_OPT
#undef LIDORT 
! (xxu, dkh, 01/09/12, adj32_011) 
#undef MODIS_AOD_OBS
! (yhmao, dkh, 01/13/12, adj32_013) 
#undef IMPROVE_BC_OC_OBS
! (zhej, dkh, 01/16/12, adj32_016)
#undef MOPITT_V5_CO_OBS
! (kjw, dkh, 02/12/12, adj32_023) 
#undef TES_CH4_OBS
#undef SCIA_CH4_OBS
#undef MEM_CH4_OBS
#undef LEO_CH4_OBS
#undef GEOCAPE_CH4_OBS


!----------CO observations------
! pick any combination
! => MOPITT CO                  
!    => MOPITT v3
!    => MOPITT v4
!    => MOPITT v5
! => AIRS CO                    
! => SCIA Bremen CO    
!#define  MOPITT_V3_CO_OBS 'MOPITT_V3_CO_OBS'
!#define  MOPITT_V4_CO_OBS 'MOPITT_V4_CO_OBS'
!#define  MOPITT_V5_CO_OBS 'MOPITT_V5_CO_OBS'
!#define  AIRS_CO_OBS     'AIRS_CO_OBS'
!#define  SCIA_BRE_CO_OBS   'SCIA_BRE_CO_OBS'

!---------aerosol-related----------
!NH3 observations               
! = > TES_NH3_OBS               
!SO2 observations
! => SCIA_DAL_SO2_OBS           
!Aerosol observations           
! => PM_ATTAINMENT              
! => IMPROVE_SO4_NIT_OBS        
! => IMPROVE_BC_OC_OBS        
! => CASTNET_NH4_OBS             
!#define TES_NH3_OBS         'TES_NH3_OBS'
!#define SCIA_DAL_SO2_OBS    'SCIA_DAL_SO2_OBS'
!#define PM_ATTAINMENT       'PM_ATTAINMENT' 
!#define IMPROVE_SO4_NIT_OBS 'IMPROVE_SO4_NIT_OBS'
!#define IMPROVE_BC_OC_OBS   'IMPROVE_BC_OC_OBS'   
!#define CASTNET_NH4_OBS     'CASTNET_NH4_OBS'
!#define MODIS_AOD_OBS       'MODIS_AOD_OBS'

!--------ozone-related--------------
! => SOMO35_ATTAINMENT          
! => TES O3
! => TES O3 IRKs 
!#define SOMO35_ATTAINMENT 'SOMO35_ATTAINMENT' 
!#define TES_O3_OBS    'TES_O3_OBS'
!#define TES_O3_IRK    'TES_O3_IRK'

!-------CH4 Observations------------
! => TES CH4
! => SCIA CH4
! => MEM CH4
! => Generic LEO instrument CH4
! => GEOCAPE CH4
!#define TES_CH4_OBS         'TES_CH4_OBS'
!#define SCIA_CH4_OBS        'SCIA_CH4_OBS'
!#define MEM_CH4_OBS         'MEM_CH4_OBS'
!#define LEO_CH4_OBS         'LEO_CH4_OBS'
!#define GEOCAPE_CH4_OBS     'GEOCAPE_CH4_OBS'

!-------NO2 observations------------
! => SCIA_KNMI_NO2_OBS          
! => SCIA_DAL_NO2_OBS           
!#define SCIA_KNMI_NO2_OBS   'SCIA_KNMI_NO2_OBS'
!#define SCIA_DAL_NO2_OBS    'SCIA_DAL_NO2_OBS'

!-------CO2 observations------------
! => GOSAT_CO2_OBS
!#define GOSAT_CO2_OBS


!------other options-----------------
!#define PSEUDO_OBS    'PSEUDO_OBS'
!#define LOG_OPT       'LOG_OPT' 
!#define LIDORT        'LIDORT'
 
