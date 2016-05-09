!$Id: inverse_driver.f,v 1.32 2012/03/01 22:00:26 daven Exp $
!$Log: inverse_driver.f,v $
!Revision 1.32  2012/03/01 22:00:26  daven
!Submit beta_v32_026 for testing (dkh, 03/01/12)
!
!Revision 1.31  2011/02/23 00:08:47  daven
!UPDATES in forward model:
! - add diag59 (lz, 11/18/10)
! - GCv8-02-04: add EPA/NEI05.
!   - completely update scale_anthro_mod.f to GCv9-01-01
!   - completely update epa_nei_mod.f to GCv9-01-01
!   - add nei2005_anthro_mod.f to GCv9-01-01
! - completely update error_mod.f to GCv9-01-01
!
!BUG FIXES in forward model:
! - GCv8-02-03: Corrected_Bond_et_al_BC.2FOC_emissions
! - GCv8-02-04: Bug_fix_in_emfossil.f_for_0.5_x_0.666_nested_grid_tagged-CO_option
! - GCv8-03-02: Fix_for_EPA.2FNEI_2005_emissions
! - GCv8-03-02: Minor_fixes_in_gamap_mod.f
! - GCv9-01-01: Bug_fix_for_biofuels_in_EPA.2FNEI05
! - GCv?-??-??: Add scaling of aromatic emissions over the US. (hotp, 11/23/09)
! - GCv9-01-01: Important_bug_fixes_for_ship_emissions
! - GCv9-01-01: Fix_to_prevent_div-by-zero_in_sulfate_mod.f
! - GCv9-01-02: Double_counting_of_biofuel_emissions_over_Asia
! - GCv9-01-02: fix SET_TINDEX for ND17, 18, 38, 39 so that all wet diagnostics get written out (dkh, 02/16/11)
!
!UPDATES in adjoint model:
! - add LITR for iteration diagnostics (zhe, dkh, 02/04/11)
! - Make sure ITS_A_NEW_MONTH is true only once per month during adjoint,
!     which minimized i/o.  (dbm, 02/10/11)
! - now make run script copy the executable rather than move it, thus
!     avoiding exessive recompilation (dbm, 02/10/11)
! - update MOPITT obs operators to support v3 and v4 (zhe, 02/04/11)
! - now add support for nested grid with offline CO (zhe, 02/04/11)
! - now emit biomass burning emissions for offline CO throughout
!     the boundary layer (dbm, 02/10/11)
! - updated input.gcadj (dkh, 02/10/11)
!    - better distinction between tracers and species
!    - better distinction between observations and control parameters
!    - additional input flags and parameters to replace hard wired options
!       - ICS_SF_DEFAULT, EMS_SF_DEFAULT, EMS_ERROR, ICS_ERROR
!       - LTRAJ_SCALE, LITR, NSPAN, LMAX_OBS, LEMS_ABS
! - replace OPT_THIS_SPECIES with OPT_THIS_TRACER (dkh, 02/10/11)
! - allow for flux filling during adjoint advection LFILL_ADJ (jkoo, dkh, 02/11/11)
! - add TES_BLVRM flag for tes NH3 observation  (dkh, 02/14/11)
! - in lidort_mod (dkh, 01/27/11)
!    - use dry diameter of BC to estimate number concentration
!    - add BC mass absorption enhancement factor ABS_FAC
!    - use growth curve for sulfate wet size rather than H2O from rpmares
! - implent LEMS_ABS option to output sensitivities w.r.t emissions
!    rather than emissions scaling factors (dkh, 02/17/11)
! - enforce LMAX_OBS = T and NSPAN = 1 for FD_GLOB (dkh, 02/21/11)
!
!BUG FIXES in  adjoint model:
! - Missing a factor of 1d6 for the cspec_ppb case in CALC_ADJ_FORCE_FOR_SENSE (fgap, dkh, 02/03/11)
! - LVARTROP treated correctly (dkh, 01/26/11)
! - For LUNZIP = T, don't delete met files during forward run (zj, dkh, 07/30/10)
! - Convert units before and after transport to account for discrete <--> continuous
!     adjoint (jkoo, dkh,  02/14/11)
! - Set the min value of CSPEC checkpt arrays to be SMAL2 (dkh, 02/19/11)
! - update sulfate_adj_mod to account for Fix_to_prevent_div-by-zero_in_sulfate_mod.f (dkh, 02/19/11)
! - fix the FD_SPOT test (dkh, 02/21/11)
!   - now only make fdglob files if FD_GLOB, not if FD_SPOT
!   - now evaluate the adjoint on 1st and 2nd iterations, and halt model if users
!       asks for a third iteration in MAYBE_DO_GEOS_CHEM_ADJ
! - Force DAYS to be at least 1 to allow for simulations less than 1 day. (dkh, 02/22/11)
!
!Revision 1.30  2010/11/19 07:05:24  daven
!BUG FIXES:
! - fix use of ICSFD and EMSFD in SET_SF and SET_LOG_SF (lz, dkh)
! - add 0.001 to diag of S_OER in tes_o3_mod before inverting (ks,mm,dkh)
! - fix bug with concnox in partition_adj.f (jkoo)
! - implement SDFLAG flag in inverse_mod.f (zhe)
!UPDATES:
! - update run script to run on prospero by default (dkh)
!   - no backup
!   - append iteration to ctm.bpch
!   - use echo instead of ex. *
! - make sure that if LADJ = F, FD_GLOB = F (dkh)
! - update tes_nh3_mod.f (dkh)
! - update CALC_APRIORI to include option for TES_NH3_OBS (dkh)
! - update comments in adj_arrays about units of EMS_SF_ADJ and ICS_SF_ADJ  (ajt)
! - add GOSAT co2 obs operator. gosat_co2_mod.f, Makefile.ifort.netcdf, adj_arrays_mod.f,
!    geos_chem_adj_mod.f, input_adj_mod.f, define_adj.h (dkh)
!
!Revision 1.29  2010/07/30 23:47:04  daven
!Patch several bugs:
! - BUG FIX: update co2 fwd model, see Ray's email 5/18
! - BUG FIX: enforce defualt scaling factors before using SF_tmp
! - BUF FIX: declare QTMP and FTMP thread private in fvdas_convect_adj_mod.f
! - BUG FIX: if an obs operator is defined, don't crash with No observations!
! - BUG FIX: use CHK_STT in MAKE_ADJ_FILE for scale factor instead of STT
! - BUF FIX: now declare BL_FRAC thread private in subroutine CHEM_OCPI_ADJ
! - BUF FIXES from Zhe Jiang, see email 5/17
!    - tagged_co_adj-mod.f (STT(I,J,L,1))
!    - define_adj.h (MOPITT_IR_CO_OBS)
!    - mopitt_obs_mod.f (MOP_COL_GRID)
!    - don't deleted unzipped met fields (geos_chem_mod)
!    - reset EMS_SF_ADJ each iteration to prevent buildup
!Cleanup and enhancements
! - don't need to call MAKE_PRESSURE_CHKFILE -
! - Replace CSPEC_O3_FORCE with CSPEC_ADJ_FORCE
! - remove 'ddd fwd' debug printout
! - change format line 112 in tes_nh3_mod.f to match that in tes_o3_mod.f
! - now allocate CHK_STT_BEFCHEM for LCHEM or LWETD
! - add define.h to checkpoint_mod.f in dep list makefiles
!New features
! - add online LIDORT and MIE code
!
!Revision 1.28  2010/05/07 20:39:47  daven
!General cleanup and streamining
! - update checkpoint_mod.f to be cleaner, remove files after used
! - remove unused directories (code_adj_emis, changsub, monika)
! - update comments at top of geos_chem_adj_mod
!Add stratospheric chemistry adjoint
! - add schem_adj.f and CO_strat_pl_adj.f
! - reinstate call to SCHEM in fwd model
!Add CO2 adjoint
! - implement fwd model updates from Ray:
!   - co2_mod.f
!   - dag04_mod.f
!   - gamap_mod.f
!   - input_mod.f
!   - logical_mod.f
! - move co2_mod to code/modified
! - update makefiles
! - add CO2 emissions IDs to adj_arrays_mod
! - add 'ppm_free_trop' as sensitivity option
! - add normalized gradients, IJ-GDEN$
!Add TES O3 obs operator
! - update Makefile.ifort.netcdf
!   - add tes_o3_mod.f
!   - link to LAPACK libraries
! - save strat O3 profile from SET_PROF in O3_PROF_SAV
! - always call SET_PROF in photoj.f
!Update TES NH3 obs operator
!
!Revision 1.27  2010/04/28 21:00:00  daven
!Now support adjoint runs spanning multiple months / years (dkh, 04/28/10)
! - update ITS_A_NEW_MONTH and ITS_A_NEW_YEAR
! - move DIRECTION to time_mod.f
!
!Revision 1.26  2010/04/25 17:18:58  daven
!BUG FIX: correctly reset adjoints in GEOS-5 convection (dkh, 04/21/10)
!BUG FIX: fix directory for cleaning *.adj.* files (jk, dkh, 04/24/10)
!Now updated support for LADJ = F (dkh, 04/25/10)
! - works with X=0 and XSTOP=0
! - updated input_adj_mod and soilnox_mod to check for LADJ
! - now use HSAVE from commode_mod instead of checkpt_mod
!Now make running with LINOZE and UPBD on as the default
!
!Revision 1.25  2010/04/01 07:09:43  daven
!Add adjoint of deposition and emissions in gas solver (dkh, 04/01/10)
! - add calcrate_adj.f, setemis_adj.f
! - apply emission scaling factors in setemis.f, move to code/modified
! - update Makefiles *
! - add to adj_arrays_mod.f: DEPSAV_ADJ, REMIS_ADJ
! - for KPP, create DMAP to speed up calculation of V_R. Saves > 10% time.
!
!Revision 1.24  2010/03/09 15:03:46  daven
!General updates and fixes
! - add define.h to dep list for inverse_mod.f in Makefiles
! - GFED2 2008 monthly data is now available (gfed2_biomass_mod)
! - upgrade to the newer bpch2_mod.f from v8-02-04
! - now only checkpt XYLAI if LCHEM in checkpt_mod (for read and write)
! - BUG FIX: correct typo in thread private pramas in SRCNH3_ADJ
! - remove obsolete 4dvar_driver.f, and references to
!   - MAKE_IMIX_CHKFILE, READ_IMIX_CHKFILE
!   - MAKE_FPBL_CHKFILE, READ_FPBL_CHKFILE
!Now include adjoint of acetone oceean sink
! - now call OCEAN_ACET_SINK in chemistry_mod
! - now call OCEAN_ACET_SINK (self-adjoint) in chemistry_adj_mod
! - update the forward model OCEAN_ACET_SINK to be more stable and
!   more precisely self-adjoint.
!Now include adjoint of UPBDFLX_NOY
! - reinstate fluxes in forward model
! - add routine UPBDFLX_NOY_ADJ
!Correct the following fwd model BUG FIXES from v8-02-04
! - update reactions in sulfate_mod.f
! - Bug fix for EMEP ship emissions
! - Minor bug fix in gamap_mod.f
! - Fixes and updates in seasalt_mod.f
! - Add EFLUX to ND67 (this actually from an earlier code update)
! - Bug fix in DIAG20 (diag_pl_mod.f)
! - Div-by-zero error encountered in arsl1k.f (just update the whole file)
! - Fix for diagnostic arrays in TPCORE
!Correct the following fwd model BUG FIXES from v8-02-05
! - make STREETS thread private in READ_ANTHRO_NH3
! - Fix for initialization of EMEP ship emissions
!Now support LADJ_TRAJ diagnostic option
! - update MAKE_ADJ_FILE
! - update gamap_mod to include IJ-ADJ-$
!Now support Tagged Ox simulation (Lin Zhang et al., GRL 2009)
! - update chemistry_adj_mod.f, geos_chem_adj_mod.f, adj_arrays_mod.f,
!   input_adj_mod.f, tagged_ox_mod.f, add tagged_ox_adj_mod.f
! - update Makefiles
! - treat it as an LADJ_EMS options as the sensitivities are w.r.t. sources
! - works OK but not exact yet.  Still needs some debugging.
!
!Revision 1.23  2010/02/10 06:25:03  daven
!Updates for additional features (dkh, 02/09/10)
! - update lightning NOx with patches from 7/10/09 from v8-02-03
! - update SO2 emissions adjoints
! - comment out IDADJ_ENOxso in adj_arrays_mod.f for now
! - now include adjoint output in tracerinfor.dat, diaginfo.dat
!   - move gama_mod.f to code/modified and update makefiles
!
!Revision 1.22  2010/01/28 17:37:21  daven
!Update for additional emissions and a few bug fixes (dkh, 01/28/10)
!- Add checkpointing to support use of MEGAN emissions
!  - checkpoint T_15_AVG and T_DAY
!  - move megan_mod.f to modified/
!- Add checkpointing to support use of lightning NOx emissions
!  - move lightning_mod.f to modified
!  - now checkpoint SLBASE
!  - move lightning_mod.o to after checkpt_mod.f in Makefiles
!  - take out the temp hack by Lee Murray to use specieal reprocessed OTD fields
!- Turn on all the standard emissions in geos5 input.geos
!- move ITS_TIME_FOR_(some met field)_ADJ functions to time_mod.f
!- BUG FIX: now use NSECb from geos_chem_mod
!  - always readin in met files, even if in the 'turn around' zone.
!- To be safe, add some constraints on the KPP <--> SMVGEAR mapping
!  of active species following recomendations of Claire Carouge
!- make geos5 benchmark use Makefile.ifort instead of Makefile.ifort.netcdf
!- update use of DIRECTION in chemdr_adj and chemdr.f
!- decrease bufsize to 4000 in gckpp_adj_Integrator.f90
!
!Revision 1.21  2010/01/06 23:05:04  daven
!Several small bug fixes and updates (dkh, 01/06/10)
!- fix hardwiring of QC_SO2 allocation, should be NSTEP (mak, 11/19/09)
!- read/wring XYLAI in checkpt_mod.f -- only do this for fullchem (mak, 11/19/09)
!- reinstate OMP pragmas in fvdas_convect_adj_mod (mak, 12/09)
!- added a prior constraint for full chem LOG_OPT (dkh, 12/14/09)
!- decrease bufsize to 4000 in gckpp_adj_Integrator (dkh, 01/06/10)
!- add define.h to dep list for adj_arrays_mod.f in all the Makefiles* (dkh, 01/06/10)
!
!Revision 1.20  2009/11/18 07:09:33  daven
!Fix several bugs in the forward model that have been found since
! the release of v8-02-01 (dkh, 11/17/09)
!- apply patch for forward model bug in
!  biomass_mod.f (mak, 11/17/09)
!From the list of bugs fixed in v8-02-02,
!http://wiki.seas.harvard.edu/geos-chem/index.php/Bugs_and_fixes
!- Bug with ND52 diagnostic
!- EPA/NEI inventory: reset other species to zero
!- Scale factor for oceanic acetone for GEOS5 2x2.5
!- Bug with PRIVATE declaration in sulfate_mod.f
!- Bug with online 2ndary aerosol (this was already fixed)
!- Bug for dust in ND48
!From the list of bugs fixed in v8-02-03
!http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_v8-02-03
!- Several bug fixes in sulfate_mod.f
!- Missing NOx data in S.E.-Asia
!- ( Mis-calculation of Courant numbers in tpcore_fvdas_mod.f90 was
!    fixed in a previous update to adjoint code )
!- Format problem in planeflight_mod.f
!- Minor fixes in wet deposition
!- Minor fixes for IBM XLF compiler
!- Don't apply 'Minor fixes in gamap_mod.f' becuase I think they
!   got it wrong (switched ALD2 with PRPE)
!- Avoiding the "Too many levels in photolysis code" error
!Implement newer BC/OC emissions from v8-02-02.
!- add USE_BOND_BIOBURN to carbon_mod.f
!- add LCOOKE to input_mod.f, logical_mod.f and switch in input.geos
!
!Revision 1.19  2009/11/12 00:45:48  daven
!Updates to emissions adjoints and general performance (dkh, 11/11/09)
! - Update TES NH3
!   - change to 4D Var mode instead of sensitivity forcing
!   - switch to TES_v4
!     - skip a few more NT
!     - now read in QFLAG, DFLAG
!   - switch to TES(NT)%VAR from TES%VAR(NT) format
!   - more diagnostic output, including doubled NH3
! - Adjoint of NH3 emissions
!   - Now include emission scaling factors in SRCNH3
!   - Add EMSSULFATE_ADJ
!   - Add SRCNH3_ADJ
!   - Update emissions adjoint IDs to ADJ_ARRAYS_MOD for fullchem
!      and get rid of old hard-coded IDs.
! - BUG FIX:  OBJSc should be OBJSe in Makefile.ifort and Makefile.ifort.netcdf
! - Update SET_LOG_SF to use ICS_SF_tmp and EMS_SF_tmp for PSEUDO_OBS
! - Update input_adj_mod.f to stop if using unsupported option LBKCOV (mak)
! - Cleanup and update the RESCALE and LOG_RESCALE routines. Now move
!   all regularization / apriori / penalty stuff elsewhere.
! - Print NSECb to prevent corruption in LOG_OPT
! - Adjoint of BC and OC emissions
!   - define ID #'s in adj_arrays_mod
!   - update carbon_mod.f to include scaling factors
!   - update carbon_adj_mod.f to include EMISSCARBON_ADJ and
!     EMITHIGH_ADJ
! - Now include Soil NOx
!   - modify soilnoxems.f to checkpt emissions and include
!     scaling factors for adjoint.
!   - add IDADJ_ENOxso to adj_arrays_mod.f
!   - move soilnoxems.f to modified/ and update Makefiles *
!   - Make DIRECTION a module variable in ADJ_ARRAYS_MOD and
!     add new routines GET_DIRECTION and SET_DIRECTION
! - Add counting of active emissions for groups of species
!   - Now include N_CARB_EMS_ADJ and N_SULF_EMS_ADJ
!   - Now include IS_CARB_EMS_ADJ and IS_SULF_EMS_ADJ
!
!Revision 1.18  2009/10/26 18:54:15  daven
!BUG FIX: recalculation of isoprene emissions (dkh, 10/26/09)
! - change conditions in geos_chem_mod for determining NEW_DAY
! - now checkpoint XYLAI -- it is tricky to recalculate
! TURN CHEM BACK ON
!Update TES NH3 operater (dkh, 10/26/09)
! - Add Makefile.ifort.netcdf
! - add new i/o diagnostics
! - update to v3 retrievals
! - BUG FIX: need to include define_adj.h at the top of
!   geos_chem_mod.f
! - move EXPAND_NAME to adj_arrays_mod.f so that it is
!   more widely accessible
! - add CMN_DEP to checkpt_mod.f in Makefiles
! - when use LIBS, take out the *.o before the -o
! - now move tes_nh3_mod to after checkpt_mod in compile list
!
!Revision 1.17  2009/10/12 18:08:52  daven
!Add TES NH3 operator (dkh, 10/12/09)
! - add to project and Makefiles
! - add GET_IJ to grid_mod.f
!Debug and test WETDEP adjoint (dkh, 10/12/09)
! - update AD_WASHOUT to match GCv8
! - BUG FIX: recalculate ALPHA correctly
! - BUG FIX: checkpoint MCHK values correclty at L = 1
! - update loops in adj wetdep routines to be parallel
!    - now make RAINFRAC_0 and WASHFRAC_0 local variable in
!      soubroutine ADJ_SO2_WETDEP
!    - declare SO2_MCHK as THREADPRIVATE in ADJ_SO2_WETDEP
! - Reset STT(SO2) and STT(SO4) at the end of DO_WETDEP_ADJ
!
!Revision 1.16  2009/10/05 01:25:15  daven
!Several imporant updates and fixes (dkh, mak 10/04/09)
! - Update Makefile (mak)
!   - move HDF pieces to new file, Makefile.ifort.hdf
!   - add ErrorModule and sciabr_co_obs_mod to Makefile.ifort
! - Now include sulfate chemistry (dkh)
!   - add sulfate_adj_mod.f
!   - move sulfate_mod.f to modified/
!     - remove sea salt interaction with SO4, NIT
!     - add checkpointing
!   - call INIT_WETSCAV_ADJ in geos_chem_adj if
!     LSULF and LCHEM to allocate SO2s_ADJ and H2O2s_ADJ
! - Now include full chem wet deposition (dkh, but not tested!)
!    - move wetscav_mod to /modified
!    - add wetscav_adj_mod to Makefile.ifort
!    - update adjoint routines
!     - change ADJ_STT --> STT_ADJ
!     - ADJ_SO2s, ADJ_H2O2s --> SO2s_ADJ, H2O2s_ADJ
!     - change IDADJxxx --> IDTxxx
!     - now just cyle past N = IDTSO2 for adjoint
!       of wetdep for non-SO2 species (old method
!       was to set RAINFRAC, WASHFRAC for SO2 = 0 )
!     - change LINUX_EFC --> SGI_MIPS for preproc
!       directtives around the parallel do
! - BUG FIX: apply tpcore_fvdas patch, no longer need
!   to set va = 0d0 (dkh)
!   - keep LFILL as an argument so that we can set it to
!     .FALSE. for adjoint transport.
! - Add CALC_APRIORI (mak)
! - Now get first guess of scaling factors from input.gcadj (mak)
!   - Update input_adj_mod and input.gcadj to input guesses
!   - Add ICS_SF_tmp and EMS_SF_tmp to ADJ_ARRAYS_MOD
!   - Update SET_SF in inverse_mod
!      - only apply to EMSFD or ICSFD
!   - Update and fix SET_LOG_SF similarly (dkh)
! - Update the default GEOS-4 tagged CO simulation
!   - 50% error for MOPITT
!   - LAERO_THEM = F
!   - turn on anthro emissions
!     - EMEP, BRAVO, STREETS, NEI99, CAC
!
!Revision 1.15  2009/09/21 01:54:19  daven
!Debug GEOS-4 convection adjoint (dkh, 09/20/09)
! - remove obsolete MAKE_CONVECTION_CHKFILE from geos_chem_mod
! - BUG FIX: now use CHK_STT_CON in DO_GEOS4_CONVECT_ADJ instead of STT
!   - and make the Q array TYPE (XPLEX) in fvdas_convect_adj_mod
! - kludge: set N_SPEC = 1 for TAGCO sim at the top of gfed2_biomass_mod
!Debug GEOS-5 advection adjoint
! - BUG FIX. During the adjoint call to GEOS-5 transport, the array "va" sometimes
!    ends up with random values, say in locations like va(71,2), which are never
!    inititialized or explicitly defined. Shouldn't they be defined somewhere?
!    That could be a bug in fwd model... but initializing va to 0d0 at the
!    start of TPCORE fixes the problem.  Note that the symptom is:
!     forrtl: severe (408): fort: (3): Subscript #2 of the array QQUWK has
!     value -2 which is less than the lower bound of -1
!General
! - now only call DO_EMISSIONS_ADJ if LADJ_EMS
! - cleanup inverse_mod.f a bit
! - make the repository geos4 simulation a tagged CO inverse ICS test,
!    while the geos5 simulation is full chem global FD test.
! - simplify INIT_WEIGHT to prevent it from crashing
!
!Revision 1.14  2009/09/15 16:10:28  daven
!Update input.gcadj (mak, dkh, 09/15/09)
! - now can specify IFD, JFD directly
!
!Revision 1.13  2009/09/15 05:33:02  daven
!Implement het chem adjoints (dkh, 09/14/09)
! - turnon CHEMCARBON in chemistry_mod.f
! - turnon CHEMCARBON_ADJ in chemistry_adj_mod.f
! - add carbon_adj_mod.f
! - move carbon_mod.f to code/modified/ and update Makefile
! - make DRYxxx public in carbon_mod
!Add adjoint of aerosol thermodynamics (dkh, 09/09/09)
! - implement LAERO_THERM flag in do_chemistry and do_chemistry_adj
! - make RECOMP_RPMARES for recalculating intermediate values
! - make rpamres_adj_mod
! - make the following routine in rpmares_mod public so that
!    they can be used in rpmares_adj_mod:
!      - POLY4, POLY6, CUBIC, AWATER, ACTCOF
!Unrelated
! - Don't stop the simulation if VAR in fwd is 1.0003d-99 and
!    the recalculated value is 1.0000d-99 in CINSPECT
!
!Revision 1.12  2009/09/08 04:18:25  daven
!Update CO emissions adjoint (mak, dkh, 09/07/09)
! - rename tagged_adj_co_mod --> tagged_co_adj_mod
!   * ( did this in Makefile, need to actually do it ) *
! - add emissions_adj_mod
! - don't recalculate forward emissions during adjoint if
!    its a tagged co simulation
! - BUG_FIX: now set initial guess scaling factors to perturbed
!    value every time passing through N_CALC = 1
!Add aerosol thermodynamics to forward code (dkh, 09/07/09)
! - move rmpares_mod.f to code/modified/
! - reinstate CALL RPMARES in chemistry_mod
! - make RPAMRES_FORADJ
!
!Revision 1.11  2009/09/07 20:12:47  daven
!Updated convection adjoint (dkh, 08/25/09)
! - modify convection_mod.f to checkpoint arrays
! - move convection_mod.f to modified/convection_mod.f
! - delete extra copy of wetscav_mod.f in code/
! - updated NFCLDMX_ADJ to support GEOS-5
!
!Revision 1.10  2009/08/17 03:59:52  daven
!Turn on chemistry for tagged CO (dkh, 07/27/09)
! - Remove the tagged_co_mod.f file from /code, as we
!    use the one in /code/modified NEED TO DO THIS
! - Remove the bpch2_mod.f file from /code, as we
!    use the one in /code/modified NEED TO DO THIS
! - GEOS-5 tagged CO with chemistry turned on will
!    crash if not using 72 vertical levels because the
!    geos5 OH file in GEOS_MEAN hasn't been reduced
!    to 30 vertical levels.  So as a temporary hack,
!    force the simulation to use the GEOS_4 OH fields.
! - same for GEOS-5 P/L fields
! - Add export OMP_NUM_THREADS=8 to run script
! - Add more informative printout to run script
! - Now run script checks for gctm.sf.* at the
!    end of each iteration
! - Fix line overflow on 1225 of input_adj_mod.f
! - Only print out the optimizaiton header if
!    ITERATE = T
! - Move call to CLEAN_FILE_DIRS to input_adj_mod
!    to avoid deleting *.obs.* files for pseudo tests
!    and to remove old gctm.sf.* files before calling
!    ARE_FLAGS_VALID (helps inform run script of crash)
! - Now check to ensure that 1 < N_CALC < 3 for FDTEST
!    in subroutine ARE_FLAGS_VALID
! - Verified FDTEST using LOG_OPT
!    - implemented LOG_RESCALE for LOG_OPT, LICS
!    - implemented LOG_OPT in APPLY_IC_SCALING,
!       adding include define_adj.h
!    - make sure call SET_LOG_SF when
!       N_CALC == N_CALC_STOP == 1
! - Use STT_ORIG in RESCALE_ADJOINT and LOG_RESCALE_ADJOINT
!    rather than reading the restart file again.
!
!Implement the full chem simulation (dkh, 08/16/09)
! - Update input.geos and input.gcadj
! - Update INIT_CHECKPT (use N_TRACERS instead of NOBS)
! - Remove chemistry_mod.f and chemdr.f from /code, as they
!    are in /code/modified
! - Remove restart_mod.f from /code; it is in /code/modified
! - Remove physproc.f    from /code; it is in /code/modified
! - Use GCKPP_ADJ_DRIVER from dkh GCv6 adjoint for both
!    forward and backward integration.
! - Loop N up to N_TRACERS instead of NOBS in INIT_WEIGHT
! - Pass back the value of IERR into ISTATUS in gckpp_adj_Integartor.f90
! - Minimal rescaling for LFDTEST or LSENS
! - IMPLEMENT OMP -- change the makefile to use parallel F90
!    compile command.
!   - Take out USE GCKPP_ADJ_Model in GCKPP_ADJ_DRIVER, so reference
!      everything explicitly from GCKPP_ADJ_GLOBAL
!   - Add gckpp_* files to dependency list for chemistry_mod.f
!      and chemistry_adj_mod.f in the Makefile.ifort
!   - Update dependancy list for all gckpp_adj_* files in Makefile.ifort
!   - Declare THREADPRIVATE in gckpp_adj_Global:
!     - JLOOP, C, VAR, VAR_ADJ, FIX, V_CSPEC, V_CSPEC_ADJ, TIME,
!        VAR_R_ADJ, RCONST
!     - stack_ptr (moved here from gckpp_adj_Integrator.f90)
!   - Declare THREADPRIVATE in gckpp_adj_Function:
!     - A
!   - Remove EQIVALENCE statment for C, VAR, FIX
!      in gckpp_adj_Global.  To compensate, define C
!      from VAR and FIX in gckpp_adj_Initialize
!   - SET VAR(15) (LISOPOH) to 1d-99.
!    - Manually set VAR(13) and VAR(14) (CO2 and DRYDEP) to be zero as well.
! - Add to Makefile.ifort and CVS
!   - chemdr_adj.f
!   - lump_adj.f
!   - partition_adj.f
! - Add INIT_KPP to gckpp_adj_Util.f90 to initialize JCOEFF.
! - Move partition.f to modified/partition.f, add PART_CASE
! - Remove CSPEC_ADJ_FOR_KPP and CSPEC_FOR_KPP_ADJ.  Just use CSPEC_ADJ.
!     The reason we have CSPEC_FOR_KPP is so that
!     you can run KPP and SMVGEAR side-by-side in the forward model.  Since
!     KPP is the only way to calculate CSPEC_ADJ, don't need to make an
!     extra copy.
! - Reset NEMIS and NNADDV before calling READCHEM when FIRSTCHEM in
!    chemrd_adj.f.  Otherwise these will get double counted, causing
!    segfault crashes in calcrate. Same for NNADDA, NNADDN, NNADDC,
!    NNADDD, NNADDF, NNADDH, NNADDG
! - Now use NTLOOP instead of ITLOOP when checkpointing CSPEC arrays
! - Call SAVE_FULL_TROP before GASCONC in CHEMDR_ADJ.  Otherwise, CSPEC
!    could get overwritten with old values in CSPEC_FULL
! - Use HSAVE (dkh) instead of HSAVE_KPP (ks) as HSAVE is set up to
!    rotate and checkpt properly.
! - Now call DO_DRYDEP and DO_EMISSIONS in GEOS_CHEM_ADJ right
!    before the call to adjoint of chemistry
! - Now call CINSPECT to check for consistancy between the forward
!    and backward values of RCONST and VAR
! - Now save CHECK_STT_BEFCHEM before DO_EMISSIONS so that SO2, SO4
!    and DMS are correct in the adjoint gas-phase chemistry. May
!    want to change this once the adjoint of the emissions
!    and sulfate chemistry are in place.
! - Now call DO_PBL_MIX(.FALSE.) up top in geos_chem_adj in
!    order for FPBL to be calculated for subsequent processes
! - Now save CSPEC_PRIOR the first time through gasconc so
!    that partition_adj works at NHMSb (otherwise CSPEC_PRIOR
!    will be zero). Make IX, IY, IZ threadprivate.
! - Can' run with soil NOx on yet until we make it recalculate
!    emissions in adj mode.  Disable it for now.
! - Reinstate CALL OPTDEPTH in DO_CHEMISTRY_ADJ so that
!    photolysis rates get recalculated correctly
!
!Revision 1.9  2009/07/14 23:51:27  daven
!Updated to run with GEOS-5 and PBL mixing (dkh, 07/14/09)
! - add support for GEOS-5
!   - if using GEOS-5, make sure that IN_CLOUD_OD is defined in define.h
!   - update GET_A3_TIME_ADJ to treat GEOS_5 the same as GEOS_4
!   - implement printout for GEOS_5 in DISPLAY_MET
! - turn on LNEI99 (or else the code bombs when trying to access USA_MASK
!    which is not allocated). This is a bug in the standard forward code.
! - move call to CALC_ADJ_FORCE to after interpolation of the I-6 fields
! - add constraint for FDTEST in ITS_TIME_FOR_OBS that it must be a TS_CHEM
! - implement TURBDAY_ADJ, now recalculate IMIX and FPBL rather than checkpoint
!   - this means atting GET_IMIX and GET_FPBL to pbs_mix_mod
!   - modify adjoint code to reference these routines
!   - remove unused stuff for checkpointing
!   - remove old CLEANUP_PBL_MIX_ADJ from cleanup and cleanup_adj
! - apply FD diff to ICSFD in SET_SF_FORFD, FD_SPOT, LICS
! - ensure that UNITS of COST_FUNC are the default for FD_GLOB
!
!Revision 1.8  2009/06/26 03:57:58  daven
!Updated, CO FD_GLOB LICS now works w/o any processes (dkh, 06/25/09)
! - remove ctm.bpch and geos from repository to speedup
!    checkout
! - have gctm.model.* gctm.costfn* and gctm.obs* be written
!    to DIAGADJ_DIR instead of OPT_DIR
! - Take out the LADJ_TRAN flag. Would this ever
!    differ from LTRAN?
! - Now for FDTEST force FLAG to TRUE on the first attempt during
!    the adjoint integration and false otherwise.
! - clean out old *.fd.* files on N_CALC == 1
! - replace restart file with Monikas from 20040501
! - add ICSFD.  This selects the denominator of the sensitivities
!    for an initial conditions finite difference test independently
!    of the species being included in the numerator.
! - add more comments to FD menu in input.gcadj, swap order of
!    NFD and MFD. move FD_SPOT and FD_GLOB to FD menu, and
!    move definition of MMSCL to control variables menu
! - uypdate ARE_FLAGS_VALID to be more rigorous
!   - make sure LADJ_CHEM and LCHEM match
!   - make sure 1 and only 1 type of simulation selected
!   - make sure some obs are selected for 3D or 4DVar
!   - make sure no obs are selected for FDTEST
!   - make sure that 1 and only 1 of LADJ_EMS and LICS
!      is included for FDTEST
! - Re-reading the restart file gave negative STT values. Instead,
!    let's go back to using STT_ORIG.
! - fixed inconsistancies in lots of preprocessor tags:
!    MOPITT_OBS    --> MOPITT_IR_CO_OBS
!    O3_ATTAINMENT --> SOMO35_ATTAINMENT
!    CASTNET_OBS   --> CASTNET_NH4_OBS
!    IMPROVE_OBS   --> IMPROVE_SO4_NIT_OBS
! - don't leave PSEUDO_OBS on by default, as it can mess up an FD test
! - add define_adj.h to dep list in Makefile.ifort for
!    - adj_arrays_mod.o
!    - input_adj_mod.o
!    - inverse_driver.o
!    - geos_chem_mod.o
! - in input_mod.f, make sure that TS_DYN always stays the same
!    value regardless of which processes are turned on or off.
!    Move this file to code/modified and update Makefile accordingly
!
!Revision 1.7  2009/06/23 06:47:07  daven
!Updates (mak, dkh, 06/23/09)
! - move tagged_co_mod from code/modified/monika to
!    code/modified (mak)
! - add background error (mak, untested)
! - reinstate LDEL_CHKPT flag (dkh)
! - update CO obs operator defs in define_adj.h (mak)
! - distinguish between LEMS and LADJ_EMS (mak)
! - added CO diagnostics (mak)
! - updated Makefile, but still mostly commented out (mak)
!
!Revision 1.6  2009/06/19 07:05:23  daven
!Updates (mak, dkh, 06/19/09)
! - switch the test simulation to
!    - 4DVar
!    - use a real restart file
!    - not LOG_OPT in define_adj.h
!    - 30 levels (in define.h
!    - change flags in input.gcadj
!    - turn off biogenic emissions in input.geos
! - add LADJ_EMS flag to many places
! - add RESCALE routine in geos_chem_adj_mod
! - now call INIT_WEIGHT at the beginning of
!    geos_chem_adj_mod
!
!Revision 1.5  2009/06/17 07:39:04  daven
!Update met field i/o for GEOS_4 adj integration (dkh, 06/17/09)
! - decrement adjoint time before reading checkpt files
! - update routines in i6_read_mod, dao_mod
! - add SLP_TMP, LWI_TMP, TO3_TMP and TTO3_TMP
! - test w/ and /wo transport on.  Adjust call to INTERP
!    accordingly.
! - add dependencies to geos_chem_mod in Makefile
! - move get_read_mod.f to modified/
!
!Revision 1.4  2009/06/15 06:44:12  daven
!Updates  (dkh, 06/15/09)
! - New Makefile layout with folders and dependency (ks)
! - move calcrate.f to modified (ks)
! - update transport_mod (ks)
! - add run script to run directory (dkh)
!
!
      PROGRAM INVERSE
!
!*****************************************************************************************
! Program inverse is the master driver for the inverse and adjoint modeling capabilities
!  of the GEOS-Chem chemical transport model.  (dkh, ks, mak, cs  06/07/09) 
! 
! NOTES
! (1 ) Add support for inverse Hessian LINVH (dkh, 01/13/12, adj32_012) 
! (2 ) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/15/12, adj32_055) 
!*****************************************************************************************
!
      ! Reference to f90 modules
      USE A3_READ_MOD,          ONLY : UNZIP_A3_FIELDS
      USE A6_READ_MOD,          ONLY : UNZIP_A6_FIELDS
      USE ADJ_ARRAYS_MOD,       ONLY : COST_FUNC
      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC_STOP
      USE ADJ_ARRAYS_MOD,       ONLY : NOPT, IFD, JFD, LFD, MFD
      USE ADJ_ARRAYS_MOD,       ONLY : NFD,EMSFD
      USE ADJ_ARRAYS_MOD,       ONLY : STT_ADJ_FD, EMS_SF_ADJ,ICS_SF_ADJ
      USE ADJ_ARRAYS_MOD,       ONLY : INIT_ADJ_ARRAYS
      USE ADJ_ARRAYS_MOD,       ONLY : COST_FUNC_SAV
      USE ADJ_ARRAYS_MOD,       ONLY : LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,       ONLY : STRFD
      USE CHECKPT_MOD,          ONLY : MAKE_EMS_ADJ_FILE 
      USE ERROR_MOD,            ONLY : DEBUG_MSG
      USE FILE_MOD,             ONLY : CLOSE_FILES
      USE GEOS_CHEM_MOD,        ONLY : DO_GEOS_CHEM
      USE GEOS_CHEM_ADJ_MOD,    ONLY : DO_GEOS_CHEM_ADJ
      USE GWET_READ_MOD,        ONLY : UNZIP_GWET_FIELDS
      USE I6_READ_MOD,          ONLY : UNZIP_I6_FIELDS
      USE INPUT_MOD,            ONLY : READ_INPUT_FILE
      USE INPUT_ADJ_MOD,        ONLY : READ_INPUT_ADJ_FILE 
      USE INV_HESSIAN_MOD,      ONLY : UPDATE_HESSIAN
      USE INVERSE_MOD,          ONLY : SET_SF,        SET_LOG_SF
      USE INVERSE_MOD,          ONLY : SET_SF_FORFD
      USE INVERSE_MOD,          ONLY : MAKE_SF_FILE
      USE INVERSE_MOD,          ONLY : MAKE_GDT_FILE
      USE INVERSE_MOD,          ONLY : MAKE_CFN_FILE
      USE INVERSE_MOD,          ONLY : READ_GDT_FILE
      USE INVERSE_MOD,          ONLY : READ_CFN_FILE
      USE INVERSE_MOD,          ONLY : SET_OPT_RANGE
      USE INVERSE_MOD,          ONLY : INIT_INVERSE
      USE INVERSE_MOD,          ONLY : GET_X_FROM_SF
      USE INVERSE_MOD,          ONLY : GET_SF_FROM_X
      USE INVERSE_MOD,          ONLY : GET_GRADNT_FROM_ADJ
      USE INVERSE_MOD,          ONLY : X
      USE INVERSE_MOD,          ONLY : GRADNT
      USE INVERSE_MOD,          ONLY : CALC_NOPT
      USE INVERSE_MOD,          ONLY : DISPLAY_STUFF
      USE INVERSE_MOD,          ONLY : MAKE_SAT_DIAG_FILE
      USE INVERSE_MOD,          ONLY : ITER_CONDITION  
      USE INVERSE_MOD,          ONLY : MAYBE_DO_GEOS_CHEM_ADJ
      USE LOGICAL_MOD,          ONLY : LPRT
      USE LOGICAL_MOD,          ONLY : LUNZIP
      USE LOGICAL_ADJ_MOD,      ONLY : LFDTEST 
      USE LOGICAL_ADJ_MOD,      ONLY : LINVH
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ_EMS
      USE LOGICAL_ADJ_MOD,      ONLY : LICS
      USE LOGICAL_ADJ_MOD,      ONLY : LDCOSAT
      USE LOGICAL_ADJ_MOD,      ONLY : LATF
      USE LOGICAL_ADJ_MOD,      ONLY : LITR
      USE LOGICAL_ADJ_MOD,      ONLY : LEMS_ABS
      USE LOGICAL_ADJ_MOD,      ONLY : LADJ_STRAT

      !USE PHIS_READ_MOD,        ONLY : UNZIP_PHIS_FIELD

#     include "define_adj.h"   

      ! Force all variables to be declared explicitly
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      ! Program variables 
      LOGICAL                  :: ITERATE = .TRUE.

      ! Variables and parameters for optimization -- see setulb.f for
      ! definitions of these.
      INTEGER                  ::   iprint, isave(44)
      CHARACTER*60             ::   task, csave
      TYPE (XPLEX)         ::   factr, pgtol, dsave(29)
      TYPE (XPLEX)         ::   F
      LOGICAL                  ::   lsave(4)
      INTEGER, PARAMETER       ::   MMAX = 17
      INTEGER, PARAMETER       ::   MM   = 5 
      INTEGER                  ::   LENWA
      INTEGER, ALLOCATABLE     ::   nbd(:)
      INTEGER, ALLOCATABLE     ::   iwa(:)
      TYPE (XPLEX),  ALLOCATABLE     ::   llim(:)
      TYPE (XPLEX),  ALLOCATABLE     ::   u(:) 
      TYPE (XPLEX),  ALLOCATABLE     ::   wa(:) 
      INTEGER                  ::   IOPT 

      !=================================================================
      ! INVERSE starts here!
      !=================================================================
 
      ! Read forward model input file and call init routines from 
      ! other modules
      CALL READ_INPUT_FILE
      IF ( LPRT ) CALL DEBUG_MSG( '### INVERSE: a READ_INPUT_FILE' )
      
      ! Read final iteration number from file 
      OPEN(  65, file = 'ITER' )
      READ(  65,*) N_CALC_STOP 
      CLOSE( 65 )

      ! Read input file for adjoint model 
      CALL READ_INPUT_ADJ_FILE

      ! Initialize arrays for optimization
      CALL INIT_SETULB 

      ! Initialize inverse modeling module
      CALL INIT_INVERSE

      ! Curent iteration 
      N_CALC = 0 

      ! Initialize adjoint arrays
      CALL INIT_ADJ_ARRAYS 
      IF ( LPRT ) CALL DEBUG_MSG( '### INVERSE: a INIT_ADJ_ARRAYS' )

      ! Now do this in input_adj_mod.f (dkh, 07/28/09) 
      !! Clean out file directories (rm *.chk.* , *.adj.* , *.ics.* and 
      !!  *.gdt.* files )
      !CALL CLEAN_FILE_DIRS

      !=================================================================
      !           ***** R E F E R E N C E   C A L C U L A T I O N ***** 
      !                    for generating pseudo observations 
      !=================================================================
      IF ( N_CALC_STOP == 0 ) THEN 
      
         ! Now only call this once above (dkh, 07/27/09) 
         !! Remove files from previous runs 
         !CALL CLEAN_FILE_DIRS
      
         ! Set IC's to their reference values
#if   defined ( LOG_OPT ) 
            CALL SET_LOG_SF
#else
            CALL SET_SF
#endif
      
         ! Call GEOS-CHEM 
         CALL DO_GEOS_CHEM
      
         ! Make SF file
         CALL MAKE_SF_FILE 

         ! EXIT
         ITERATE = .FALSE. 
      
      ENDIF


      ! Allow for use of this driver to run only the forward model as
      ! a reference calculation. 
      IF ( .not. LADJ ) ITERATE  = .FALSE. 

 
      !=================================================================
      !           ***** S E T   S C A L I N G   F A C T O R S ***** 
      !=================================================================

      ! Now only call this once above (dkh, 07/27/09) 
      !! this call was deleting obs files! need to either delete it or
      !! replace some options inside (mak, 6/18/09)
      !CALL CLEAN_FILE_DIRS

      ! Perturb the initial conditions
      IF ( ITERATE ) THEN 
#if   defined ( LOG_OPT ) 
         CALL SET_LOG_SF
#else
         CALL SET_SF
#endif
      ENDIF 

      !=================================================================
      !           ***** O P T I M I Z A T I O N *****
      !=================================================================

      ! Set parameters for optimization. See setulb.f for definitions.
      ! Let PGTOL be very small for FDTEST, as we're not actually doing 
      ! an optimization in this case. 
      IPRINT    = 1
      FACTR     = 1.0D01
      IF ( LFDTEST ) THEN 
         PGTOL  = 1.0D-12
      ELSE
         PGTOL  = 1.0D-05
      ENDIF
#if   defined ( LOG_OPT ) 
      DO IOPT = 1, NOPT
         NBD(IOPT)    = 0  ! 0 = no bounds
      ENDDO
#else
      DO IOPT = 1, NOPT
         NBD(IOPT)    = 1
         LLIM(IOPT)   = 0.0D0
      ENDDO
#endif

      task   = 'START'
 
      ! Mare array of scaling factors into a vector for optimization
      CALL GET_X_FROM_SF

      !=================================================================
      ! OPTIMIZATION loop starts here!
      !================================================================= 
        
      IF ( ITERATE ) THEN 
         ! Echo some input to the screen
         WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
         WRITE( 6, '(a,/)'   ) 'S T A R T   O P T I M I Z A T I O N'
         WRITE (6,16)
   16    format(/ 5x, 'Solving GEOS-Chem Adjoint.'
     +          / 5x, ' (f = 0.0 at the optimal solution.)' /)
      ENDIF 

      ! Beginning of the loop
      DO WHILE( ITERATE )

         print*, ' do setulb ' 
         ! Call the L-BFGS-code
         CALL SETULB( NOPT,   MM,     X,      LLIM,   U,     NBD,
     &                F,      GRADNT, FACTR,  PGTOL,  WA,    IWA,
     &                TASK,   IPRINT, CSAVE,  LSAVE,  ISAVE, DSAVE )

         print*, ' done setulb ' , TASK(1:2)

         ! Force it to continue for FD tests, as cost func or gradients 
         ! may be very small or zero (dkh, 02/11/11) 
         IF ( LFDTEST ) TASK(1:2) = 'FG'

         IF ( TASK(1:2) == 'FG' ) THEN

            ! Iteration diagnostics (zhe 11/28/10)
            IF ( LITR ) THEN 
               IF ( N_CALC .GT. 0 )  CALL ITER_CONDITION( N_CALC ) 
               LATF = .FALSE.
            ENDIF 

            ! The minimization routine has returned to request the
            ! function f and gradient g values at the current x

            ! Update iteration count
            N_CALC = N_CALC + 1

            ! Resent cost function for this iteration
            COST_FUNC = 0.D0

            IF ( N_CALC < N_CALC_STOP ) THEN

               WRITE(6,*)  'READING SAVED DATA for N_CALC = ', N_CALC

               ! Read scaling factor values from disk
               CALL GET_SF_FROM_X

               CALL DISPLAY_STUFF( 1 )

               ! Read gradients from disk
               CALL READ_GDT_FILE

               ! Read cost function from disk
               CALL READ_CFN_FILE

               ! Put adjoints into GRADNT vector
               CALL GET_GRADNT_FROM_ADJ

               !Save the current adjoint in the finite difference test cell
               ! Initial conditions test 
               IF ( LFDTEST .AND. LICS .AND. LADJ_EMS) THEN

                  PRINT*, 'WE HAVE A PROBLEM WITH STT_ADJ_FD when LICS &
     &                    LADJ_EMS are both TRUE'

               ELSEIF ( LFDTEST .AND. LICS ) THEN

                  STT_ADJ_FD(N_CALC) = ICS_SF_ADJ(IFD,JFD,LFD,NFD) 

               ELSEIF ( LFDTEST .AND. LADJ_EMS ) THEN

                  
                  ! Emissions test 
                  IF ( .NOT. LADJ_STRAT ) THEN

                     STT_ADJ_FD(N_CALC) = EMS_SF_ADJ(IFD,JFD,MFD,EMSFD)

                  ! Strat prod and loss sense (hml, adj32_025)
                  ELSEIF ( LADJ_STRAT ) THEN

                     STT_ADJ_FD(N_CALC) = LOSS_SF_ADJ
     &                                    (IFD,JFD,MFD,STRFD)

                  ENDIF

               ENDIF

               ! Copy value of COST_FUNC to the optimization variable F 
               F = COST_FUNC

               ! Save current cost function
               COST_FUNC_SAV(N_CALC) = COST_FUNC

               CALL DISPLAY_STUFF( 2 )

               ! to estimate inverse Hessian (offline) (dkh, 01/12/12, adj32_012) 
               IF ( N_CALC == 1 .and. LINVH ) CALL UPDATE_HESSIAN

               ! Return to beginning of loop

            ELSEIF ( N_CALC == N_CALC_STOP ) THEN

               ! Done if we are just estimating inverse Hessian (dkh, 01/12/12, adj32_012) 
               IF ( LINVH ) THEN 
                  WRITE(6,*) ' Force quit'
                  STOP      
               ENDIF 

               ! UPDATE THE INITIAL CONDITIONS 

               ! If we're doing a finite difference test, reset to the orginal
               ! SF and augment by amount FD_DIFF. Don't use X in this case.
               ! old:
               !IF ( ACTIVE_VARS == 'FDTEST' .AND. N_CALC == 2 ) THEN
               ! new: now support 2nd order FDTEST (MAKE_SAVE_FILE_2)
               IF ( LFDTEST .AND. N_CALC > 1 ) THEN

                  CALL SET_SF_FORFD

               ELSEIF ( N_CALC == 1 ) THEN
                  
                  ! don't need to call this again ??
                  !CALL SET_SF
#if   defined ( LOG_OPT ) 
                  CALL SET_LOG_SF
#else
                  CALL SET_SF
#endif

               ELSE

                  ! Update the scaling factors to the current X
                  CALL GET_SF_FROM_X

               ENDIF

               CALL DISPLAY_STUFF( 3 )


               !==============================================================
               ! OPTIONAL: uncomment to use scaling factors from another run
               !==============================================================
               !CALL READ_SF_FILE

               !==============================================================
               ! FORWARD RUN 
               !==============================================================
               CALL DO_GEOS_CHEM

               !==============================================================
               ! ADJOINT CALCULATION 
               !==============================================================
               IF ( .not. LFDTEST ) THEN 
                   
                   CALL DO_GEOS_CHEM_ADJ
                 
               ! For finite difference test, we may or may not do adjoint 
               ELSE 

                   CALL MAYBE_DO_GEOS_CHEM_ADJ

               ENDIF 
                 
               !==============================================================
               ! SAVE RESULTS TO DISK and EXIT OPTIMIZATION LOOP
               !==============================================================

               ! Zero the gradients of the species that we do not wish to optimize
               ! or in places that you don't want optimized 
               CALL SET_OPT_RANGE
                   ! Add to this Kumaresh's spatial filter 

               ! Write gradients 
               CALL MAKE_GDT_FILE

               ! Write scaling factors
               CALL MAKE_SF_FILE

               ! Write cost function
               CALL MAKE_CFN_FILE

               !==============================================================
               ! Diagnostics 
               !=================================================

               ! store satellite diagnostics
               ! for now CO, but subroutines all general, just need linking
               ! (mak 6/19/09)
               IF ( LDCOSAT ) THEN
                  !Store FORCING, MOP_MOD_DIFF and MODEL_BIAS
                  !CALL MAKE_FORCING_FILE
                  !CALL MAKE_MOPMOD_FILE
                  ! store model, mopitt and model bias to files
                  ! model
                  CALL MAKE_SAT_DIAG_FILE(  1 )
              
                  ! obs and DOFs
                  IF( N_CALC_STOP == 1) THEN
                     CALL MAKE_SAT_DIAG_FILE(  2 )
                  ENDIF
                  CALL MAKE_SAT_DIAG_FILE(  6 )
                 
                  CALL MAKE_SAT_DIAG_FILE(  7 )
 
                  ! model bias (wrt satellite data)
                  CALL MAKE_SAT_DIAG_FILE(  3 )
              
                  ! store COST_ARRAY, OBS_COUNT, OBS_HOUR*
                  CALL MAKE_SAT_DIAG_FILE(  5 )
              
               ENDIF

               IF ( LEMS_ABS ) CALL MAKE_EMS_ADJ_FILE 
           
               ! Write results to screen 
               CALL DISPLAY_STUFF( 4 )

               ! Exit loop
               ITERATE = .FALSE.

            ENDIF 


         ELSEIF ( TASK(1:5) == 'NEW_X' ) THEN

            ! The minimization routine has returned with a new iterate, 
            ! and we have opted to continue the interation

            ! Update the inverse hessian approximation (dkh, 01/12/12, adj32_012) 
            IF ( LINVH ) THEN 
               CALL UPDATE_HESSIAN  
            ENDIF 

         ELSE

            ! We terminate execution when TASK is neither FG nor NEW_X.
            ! We print the information contained in the string TASK
            ! if the default output is not used and the execution is 
            ! not stopped intentionally by the user.
            IF ( IPRINT == -1 .AND. TASK(1:4) /= 'STOP' )
     &         WRITE(6,*) TASK

            WRITE(6,*) TASK

            ! Exit loop
            ITERATE = .FALSE.

         ENDIF

      !=================================================================
      ! OPTIMIZATION loop ends here!
      !=================================================================
      ENDDO
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)'   ) 'F I N I S H   I T E R A T I O N '       

      ! Clean up and quit
      CALL CLOSE_FILES
      CALL CLEANUP
      CALL CLEANUP_ADJ

      ! Remove all met files from temporary directory
      IF ( LUNZIP ) THEN
         CALL UNZIP_A3_FIELDS(  'remove all' )
         CALL UNZIP_A6_FIELDS(  'remove all' )
         CALL UNZIP_I6_FIELDS(  'remove all' )
         !CALL UNZIP_PHIS_FIELD( 'remove all' )

#if   defined( GEOS_3 )
         ! We only need to remove the GWET fields if we are
         ! using the online dust simulation (bmy, 4/1/04)
         IF ( LDUST ) THEN
            CALL UNZIP_GWET_FIELDS( 'remove all' )
         ENDIF
#endif

      ENDIF
 
      ! Write the final iteration number for the next iteration to file 
      OPEN(  65, file = 'ITER' )
      WRITE(  65,*) N_CALC_STOP + 1
      CLOSE( 65 )


      WRITE( 6, '(a,/)'   ) 'D O N E'
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
 
      CONTAINS 
!------------------------------------------------------------------------------

      SUBROUTINE INIT_SETULB( ) 
!
!******************************************************************************
!  Subroutine INIT_SETULB initializes arrays used by the optimization routine,
!  setulb, whose size depends upon the model simulation type and resolution. 
!  (dkh, 06/07/09) 
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : NOPT 
      USE ERROR_MOD,          ONLY : ALLOC_ERR 
      USE INVERSE_MOD,        ONLY : CALC_NOPT

      ! Local variables
      INTEGER                :: AS

      !=================================================================
      ! INIT_SETULB begins here!
      !=================================================================

      ! Calculate the maximum number of control parameters that could
      ! be optimized, NOPT
      CALL CALC_NOPT 
     
      LENWA = 2 * MM * NOPT + 4 * NOPT + 11 * MM * MM + 8 * MM

      ALLOCATE( NBD( NOPT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NBD' )
      !NBD = 0

      ALLOCATE( IWA( 3*NOPT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IWA' )
      IWA = 0

      ALLOCATE( LLIM( NOPT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LLIM' )
      ! Don't set default bounds (dkh, 11/07/09) 
      !LLIM = 0

      ALLOCATE( U( NOPT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LLIM' )
      ! Don't set default bounds (dkh, 11/07/09) 
      !U = 0

      ALLOCATE( WA( LENWA ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WA' )
      WA = 0

      ! Return to calling program
      END SUBROUTINE INIT_SETULB 

!------------------------------------------------------------------------------
     
      END PROGRAM INVERSE

