!$Id: logical_adj_mod.f,v 1.5 2012/03/01 22:00:26 daven Exp $

      MODULE LOGICAL_ADJ_MOD
!
!******************************************************************************
!  Module LOGICAL_ADJ_MOD contains all of the logical switches used by 
!  adjoint GEOS-CHEM.
!  (adj_group, 6/07/09)
!
!  Module Variables:
!  ============================================================================
!  (1 ) LADJ       (LOGICAL) : ON/OFF switch for adjoint run
!  (2 ) LADJ_TRAN  (LOGICAL) : ON/OFF switch for adjoint transport
!  (3 ) LADJ_CHEM  (LOGICAL) : ON/OFF switch for adj chemistry
!  (4 ) LAERO_THERM(LOGICAL) : ON/OFF switch for aerosol thermo
!  (5 ) LFD_SPOT   (LOGICAL) : ON/OFF switch for FD in 1 box
!  (6 ) LFD_GLOB   (LOGICAL) : ON/OFF switch for FD in 3d
!  (7 ) LSENS      (LOGICAL) : ON/OFF switch for sensitivity run
!  (8 ) L4DVAR     (LOGICAL) : ON/OFF switch for 4dvar run
!  (9 ) L3DVAR     (LOGICAL) : ON/OFF switch for 3dvar run
!  (10) LAPSRC     (LOGICAL) : ON/OFF switch for adding 2nd part of J
!  (11) LBKCOV     (LOGICAL) : ON/OFF switch for computing cov. matrix
!  (12) LINVH      (LOGICAL) : ON/OFF switch for computing inv. Hessian
!  (13) LLINOZ     (LOGICAL) : ON/OFF switch for LINOZ (fwd AND adj)
!  (14) LFDTEST    (LOGICAL) : ON/OFF switch for FD test SPOT or GLOBAL
!  (15) LADJ_EMS   (LOGICAL) : ON/OFF switch for emission optimization
!  (16) LICS       (LOGICAL) : ON/OFF switch for initial condi. optimization
!  (17) LRXNR      (LOGICAL) : ON/OFF switch for RXN rates as control vars
!  (18) LADJDIAG   (LOGICAL) : ON/OFF switch for saving adj diagnostics
!  (19) LJSAVE     (LOGICAL) : ON/OFF switch for saving .save and .save2 files
!  (20) LADJ_TRAJ  (LOGICAL) : ON/OFF switch for saving trajectory files
!  (21) LDCOSAT    (LOGICAL) : ON/OFF switch for saving CO satellite diagnostics
!  (22) LHMOD      (LOGICAL) : ON/OFF switch for saving H(model)
!  (23) LHOBS      (LOGICAL) : ON/OFF switch for saving h(obs)
!  (24) LHMODIFF   (LOGICAL) : ON/OFF switch for saving H(model)-h(obs)
!  (25) LADJ_FORCE (LOGICAL) : ON/OFF switch for saving adjoint forcing
!  (26) LMODBIAS   (LOGICAL) : ON/OFF switch for saving model bias
!  (27) LOBS_COUNT (LOGICAL) : ON/OFF switch for saving obs count/gridbox
!  (28) LDOFS      (LOGICAL) : ON/OFF switch for saving gridded DOFs
!  (29) LPRINTFD   (LOGICAL) : ON/OFF switch for printing adj debug
!  (30) LDEL_CHKPT (LOGICAL) : ON/OFF switch for deleting checkpoint files
!  (31) LITR       (LOGICAL) : ON/OFF switch for saving iteration diagnostics
!  (32) LDEVOC     (LOGICAL) : ON/OFF switch for including CO VOCs 
!  (33) LATF       (LOGICAL) : ON/OFF switch for saving iteration diagnostics
!  (34) LMAX_OBS   (LOGICAL) : ON/OFF switch for capping number of obs
!  (35) LTRAJ_SCALE(LOGICAL) : ON/OFF switch for scaling STT_ADJ in *.adj.* files
!  (36) LKGBOX     (LOGICAL) : ON/OFF switch for cost function in kg/box
!  (37) LUGM3      (LOGICAL) : ON/OFF switch for cost function in kg/box
!  (38) LSTT_PPB   (LOGICAL) : ON/OFF switch for cost function in ppb
!  (39) LSTT_PPB_TROP_PPM (L): ON/OFF switch for cost function in ppm
!  (40) LCSPEC_PPB (LOGICAL) : ON/OFF switch for cost function in cspec ppb
!  (41) LCSPEC_OBS (LOGICAL) : ON/OFF switch for observing any cspec species 
!  (42) LTES_BLVMR (LOGICAL) : ON/OFF switch for looking at TES BLVMR
!  (43) LFILL_ADJ  (LOGICAL) : ON/OFF switch for filling during adj advenction
!  (44) LEMS_ABS   (LOGICAL) : ON/OFF switch for sense w.r.t abosulte emissions
!  (45) LTES_PSO   (LOGICAL) : ON/OFF switch for CH4
!  (46) LPOP_UGM3  (LOGICAL) : ON/OFF switch for sense w.r.t pop weighted conc
!  (47) LAD_STRAT  (LOGICAL) : ON/OFF switch for strat chem adjoint
! 
!  NOTES:
!  (1 ) Added LITR, LDEVOC and LATF (zhe, dkh, 02/04/11) 
!  (2 ) Added lots of new flags:  (dkh, 02/09/11) 
!         LMAX_OBS,   LTRAJ_SCALE, LKGBOX,     LUGM3,     LSTT_PPB, LSTT_TROP_PPM,
!         LCSPEC_PPB, LCSPEC_OBS,  LTES_BLVMR, LFILL_ADJ, LEMS_ABS
!  (3 ) Added LTES_PSO (kjw, dkh, 02/12/12, adj32_023)
!  (4 ) Added LPOP_UGM3 (sev, dkh, 02/13/12, adj32_024) 
!  (5 ) Added LADJ_STRAT (hml, dkh, 02/14/12, adj32_025) 
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      LOGICAL :: LADJ
      LOGICAL :: LADJ_ONLY !(jkoo, 09/26/2011)
      LOGICAL :: LADJ_TRAN
      LOGICAL :: LADJ_CHEM
      LOGICAL :: LAERO_THERM 
      LOGICAL :: LFD_SPOT
      LOGICAL :: LFD_GLOB
      LOGICAL :: LSENS
      LOGICAL :: L4DVAR
      LOGICAL :: L3DVAR
      LOGICAL :: LAPSRC
      LOGICAL :: LBKCOV
      LOGICAL :: LINVH
      LOGICAL :: LLINOZ
      LOGICAL :: LFDTEST
      LOGICAL :: LADJ_EMS
      LOGICAL :: LICS
      LOGICAL :: LRXNR
      LOGICAL :: LADJDIAG
      LOGICAL :: LJSAVE
      LOGICAL :: LADJ_TRAJ
      LOGICAL :: LDCOSAT
      LOGICAL :: LHMOD
      LOGICAL :: LHOBS
      LOGICAL :: LHMODIFF
      LOGICAL :: LADJ_FORCE
      LOGICAL :: LMODBIAS
      LOGICAL :: LOBS_COUNT
      LOGICAL :: LDOFS
      LOGICAL :: LPRINTFD
      LOGICAL :: LDEL_CHKPT
      LOGICAL :: LITR
      LOGICAL :: LDEVOC
      LOGICAL :: LATF
      LOGICAL :: LMAX_OBS
      LOGICAL :: LKGBOX
      LOGICAL :: LUGM3
      LOGICAL :: LSTT_PPB   
      LOGICAL :: LSTT_TROP_PPM
      LOGICAL :: LCSPEC_PPB
      LOGICAL :: LTRAJ_SCALE
      LOGICAL :: LCSPEC_OBS
      LOGICAL :: LTES_BLVMR
      LOGICAL :: LFILL_ADJ
      LOGICAL :: LEMS_ABS
      LOGICAL :: LTES_PSO
      LOGICAL :: LPOP_UGM3      
      LOGICAL :: LADJ_STRAT


      ! End of module
      END MODULE LOGICAL_ADJ_MOD
