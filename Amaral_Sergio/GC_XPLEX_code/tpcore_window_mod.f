! $Id: tpcore_window_mod.f,v 1.1.1.1 2009/06/09 21:51:53 daven Exp $
      MODULE TPCORE_WINDOW_MOD
!
!******************************************************************************
!  Module TPCORE_MOD contains the TPCORE transport subroutine package by
!  S-J Lin, version 7.1. (yxw, bmy, 12/2/03, 11/5/08)
!  
!  Module routines:
!  ============================================================================
!  (1 ) TPCORE_WINDOW   : TPCORE driver routine for nested-grid simulation
!  (2 ) COSA            : TPCORE internal subroutine
!  (3 ) COSC            : TPCORE internal subroutine           
!  (4 ) FCT3D           : TPCORE internal subroutine 
!  (5 ) FILEW           : TPCORE internal subroutine   
!  (6 ) FILNS           : TPCORE internal subroutine
!  (7 ) FXPPM           : TPCORE internal subroutine
!  (8 ) FYPPM           : TPCORE internal subroutine
!  (9 ) FZPPM           : TPCORE internal subroutine    
!  (10) HILO            : TPCORE internal subroutine
!  (11) HILO3D          : TPCORE internal subroutine
!  (12) QCKXYZ          : TPCORE internal subroutine
!  (13) LMTPPM_x        : TPCORE internal subroutine
!  (14) LMTPPM_y        : TPCORE internal subroutine
!  (15) LMTPPM_z        : TPCORE internal subroutine
!  (16) XADV            : TPCORE internal subroutine
!  (17) XMIST           : TPCORE internal subroutine
!  (18) XTP             : TPCORE internal subroutine
!  (19) YMIST           : TPCORE internal subroutine
!  (20) YTP             : TPCORE internal subroutine
!  (21) PRESS_FIX       : TPCORE pressure-fixer driver routine
!  (22) DYN0            : TPCORE pressure-fixer internal subroutine
!  (23) PFILTR          : TPCORE pressure-fixer internal subroutine
!  (24) LOCFLT          : TPCORE pressure-fixer internal subroutine       
!  (25) POLFLT          : TPCORE pressure-fixer internal subroutine
!  (26) DIAG_FLUX       : Computes ND24, ND25, ND26 mass flux diagnostics 
!  (27) POSITION_WINDOW : TPCORE internal subroutine
!
!  Reference Diagram: 
!  ============================================================================
!
!  <-------------------------------------- IGLOB ---------------------->
!
!  +-------------------------------------------------------------------+   ^
!  | GLOBAL REGION                                                     |   |
!  |                                                                   |   |
!  |                       <-------------- IIPAR ------------->        |   |
!  |                                                                   |   |
!  |                       +=================================[Y]  ^    |   |
!  |                       |  WINDOW REGION (met field size)  |   |    |   |
!  |                       |                                  |   |    |   |
!  |                       |      <------- IM_W ------->      |   |    |   |
!  |                       |      +--------------------+  ^   |   |    |   |
!  |                       |      |  TPCORE REGION     |  |   |   |    |   |
!  |                       |      |  (transport is     |  |   |   |    |   |
!  |<------- I0 ---------->|<---->|   done in this     | JM_W | JJPAR  | JGLOB
!  |                       | I0_W |   window!!!)       |  |   |   |    |   |
!  |                       |      |                    |  |   |   |    |   |
!  |                       |      +--------------------+  V   |   |    |   |
!  |                       |        ^                         |   |    |   |
!  |                       |        | J0_W                    |   |    |   |
!  |                       |        V                         |   |    |   |
!  |                      [X]=================================+   V    |   |
!  |                                ^                                  |   |
!  |                                | J0                               |   |
!  |                                V                                  |   |
! [1]------------------------------------------------------------------+   V
!
!  DIAGRAM NOTES:
!  (a) The outermost box ("Global Region") is the global grid size.  This 
!      region has IGLOB boxes in longitude and JGLOB boxes in latitude.  
!      The origin of the "Global Region" is at the south pole, at the 
!      lower left-hand corner (point [1]). 
!
!  (b) The next innermost box ("Window Region") is the nested-grid window.
!      This region has IIPAR boxes in longitude and JJPAR boxes in latitude.
!      This is the size of the trimmed met fields that will be used for
!      a 1 x 1 "nested-grid" simulation.  
!          
!  (c) The innermost region ("TPCORE Region") is the actual area in which
!      TPCORE transport will be performed.  Note that this region is smaller
!      than the "Window Region".  It is set up this way since a cushion of 
!      grid boxes is needed TPCORE Region for boundary conditions.
!
!  (d) I0 is the longitude offset (# of boxes) and J0 is the latitude offset
!      (# of boxes) which translate between the "Global Region" and the
!      "Window Region". 
!
!  (e) I0_W is the longitude offset (# of boxes), and J0_W is the latitude
!      offset (# of boxes) which translate between the "Window Region"
!      and the "TPCORE Region".  
!
!  (f) The lower left-hand corner of the "Window Region" (point [X]) has
!      longitude and latitude indices (I1_W, J1_W).  Similarly, the upper
!      right-hand corner (point [Y]) has longitude and latitude indices 
!      (I2_W, J2_W).
!
!  (g) Note that if I0=0, J0=0, I0_W=0, J0_W=0, IIPAR=IGLOB, JJPAR=JGLOB
!      specifies a global simulation.  In this case the "Window Region"
!      totally coincides with the "Global Region".  
!
!  (h) In order for the nested-grid to work we must save out concentrations
!      over the WINDOW REGION from a coarse model (e.g. 4x5) corresponding to
!      the same WINDOW REGION at 1x1.  These concentrations are copied along
!      the edges of the 1x1 WINDOW REGION and are thus used as boundary
!      conditions for TPCORE.  We assume that we will save out concentrations
!      from the 4x5 model since presently it takes too long to run at 2x25.
!
!  NOTES:
!  (1 ) Denote differences from "tpcore_mod.f" by !%%%.  Also assume that the
!        window region does not include the polar caps.   Now assume all 
!        platforms other than CRAY use OPENMP parallelization commands	
!        (yxw, bmy, 3/10/03)
!  (2 ) Updated information output depending on what type of machine it is.
!        (bmy, 12/2/03)
!  (3 ) Commented out call to FLUSH(6) (bmy, 1/26/04)
!  (4 ) Simplify PRIVATE definitions.  Also fixed bug in FZPPM which was
!        preventing the nested grid run from working on Altix (bmy, 11/9/04)
!  (5 ) Remove obsolete CO-OH code (bmy, 6/24/05)
!  (6 ) Now print output for IFORT compiler in "tpcore_window" (bmy, 10/18/05)
!  (7 ) Now do not parallelize DO loop 2500 in TPCORE_WINDOW.  For some reason 
!        this results in NaN's.  All other parallel loops may be left 
!        activated.  Also, now place all parallel loops in all routines w/in 
!        an #if defined block. (bmy, 11/5/08)
!******************************************************************************
!
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "tpcore_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
      PRIVATE

      ! ... except this routine
      PUBLIC :: TPCORE_WINDOW
      
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE TPCORE_WINDOW( IGD,  Q,    PS1,  PS2,  U,    V,   
     &                          W,    NDT,  IORD, JORD, KORD, NC,   
     &                          IM,   JM,   J1,   I0,   J0,   I0_W, 
     &                          J0_W, I1_W, J1_W, I2_W, J2_W, IM_W,
     &                          JM_W, IGZD, NL,   AP,   BP,   PT,
     &                          AE,   FILL, MFCT, Umax )
 
C TransPort module for Goddard Chemistry Transport Model (G-CTM), Goddard
C Earth Observing System General Circulation Model (GEOS-GCM), and Data
C Assimilation System (GEOS-DAS).
 
C Purpose: perform the transport of  3-D mixing ratio fields using 
C          externally specified winds on the hybrid Eta-coordinate.
C          One call to tpcore updates the 3-D mixing ratio
C          fields for one time step (NDT). [vertical mass flux is computed
C          internally using a center differenced hydrostatic mass
C          continuity equation].
 
C Schemes: Multi-dimensional Flux Form Semi-Lagrangian (FFSL) schemes
C          (Lin and Rood 1996, MWR) with a modified MFCT option (Zalesak 1979).
 
C Multitasking version: 7.1
C Last modified: Sept 2, 1999
C Changes from version 7.m: large-time-step bug in xtp fixed.
C Suggested compiler options:
C CRAY f77 compiler:  cf77 -Zp -c -Wd'-dec' -Wf' -a stack -exm'
C CRAY f90 compiler:  f90 -c -eZ -DCRAY -Dmultitask
C SGI Origin: f77 -c -DSGI -Dmultitask -r8 -64 -O3 -mips4 -mp
C             loader: f77 -64 -mp
C
C Send comments/suggestions to
C
C                 S.-J. Lin
C Address:
C                 Code 910.3, NASA/GSFC, Greenbelt, MD 20771
C                 Phone: 301-614-6161
C                 E-mail: slin@dao.gsfc.nasa.gov
C
C The algorithm is based on the following papers:
 
C 1. Lin, S.-J., and R. B. Rood, 1996: Multidimensional flux form semi-
C    Lagrangian transport schemes. Mon. Wea. Rev., 124, 2046-2070.
C
C 2. Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: A class of
C    the van Leer-type transport schemes and its applications to the moist-
C    ure transport in a General Circulation Model. Mon. Wea. Rev., 122,
C    1575-1593.
C
C 3. Lin, S.-J., and R. B. Rood, 1997: Multidimensional flux form semi-
C    Lagrangian transport schemes- MFCT option. To be submitted.
 
C ======
C INPUT:
C ======
 
C IGD: (horizontal) grid type on which winds are defined.
C IGD = 0  A-Grid  [all variables defined at the same point from south
C                   pole (j=1) to north pole (j=JM) ]

C IGD = 1  GEOS-GCM C-Grid (Max Suarez's center difference dynamical core)

C                                  [North]

C                                   V(i,j)
C                                      |
C                                      |
C                                      |
C                 [WEST]  U(i-1,j)---Q(i,j)---U(i,j) [EAST]
C                                      |
C                                      |
C                                      |
C                                   V(i,j-1)
 
C                                   [South]

C         U(i,   1) is defined at South Pole.
C         V(i,   1) is half grid north of the South Pole.
C         V(i,JM-1) is half grid south of the North Pole.
C
C         V must be defined at j=1 and j=JM-1 if IGD=1
C         V at JM need not be defined.
 
C Q(IM,JM,NL,NC): mixing ratios at current time (t)
C NC: total # of constituents
C IM: first (E-W) dimension; # of Grid intervals in E-W is IM
C JM: 2nd (N-S) dimension;   # of Grid intervals in N-S is JM-1
C NL: 3rd dimension (# of layers); vertical index increases from 1 at
C       the model top to NL near the surface (see fig. below).
C       It is assumed that NL > 5.
C
C PS1(IM,JM): surface pressure at current time (t)
C PS2(IM,JM): surface pressure at mid-time-level (t+NDT/2)
C PS2 is replaced by the predicted PS (at t+NDT) on output.
C Note: surface pressure can have any unit or can be multiplied by any
C       const.
C
C The hybrid ETA-coordinate:
C
C pressure at layer edges are defined as follows:
C
C        p(i,j,k) = AP(k)*PT  +  BP(k)*PS(i,j)          (1)
C
C Where PT is a constant having the same unit as PS.
C AP and BP are unitless constants given at layer edges.
C In all cases  BP(1) = 0., BP(NL+1) = 1.
C The pressure at the model top is PTOP = AP(1)*PT
C
C *********************
C For pure sigma system
C *********************
C AP(k) = 1 for all k, PT = PTOP,
C BP(k) = sige(k) (sigma at edges), PS = Psfc - PTOP, where Psfc
C is the true surface pressure.
C
C                  /////////////////////////////////
C              / \ ------ Model top P=PTOP ---------  AP(1), BP(1)
C               |
C    delp(1)    |  ........... Q(i,j,1) ............
C               |
C     W(k=1)   \ / ---------------------------------  AP(2), BP(2)
C
C
C
C     W(k-1)   / \ ---------------------------------  AP(k), BP(k)
C               |
C    delp(K)    |  ........... Q(i,j,k) ............
C               |
C      W(k)    \ / ---------------------------------  AP(k+1), BP(k+1)
C
C
C
C              / \ ---------------------------------  AP(NL), BP(NL)
C               |
C    delp(NL)   |  ........... Q(i,j,NL) .........
C               |
C     W(NL)=0  \ / -----Earth's surface P=Psfc ------ AP(NL+1), BP(NL+1)
C                 //////////////////////////////////
 
C U(IM,JM,NL) & V(IM,JM,NL):winds (m/s) at mid-time-level (t+NDT/2)
C Note that on return U and V are destroyed.
 
C NDT (integer): time step in seconds (need not be constant during the course of
C      the integration). Suggested value: 30 min. for 4x5, 15 min. for 2x2.5
C      (Lat-Lon) resolution. Smaller values maybe needed if the model
C      has a well-resolved stratosphere and Max(V) > 225 m/s
C
C J1 determines the size of the polar cap:
C    South polar cap edge is located at -90 + (j1-1.5)*180/(JM-1) deg.
C    North polar cap edge is located at  90 - (j1-1.5)*180/(JM-1) deg.
C There are currently only two choices (j1=2 or 3).
C IM must be an even integer if j1 = 2. Recommended value: J1=3.
C
C IORD, JORD, and KORD are integers controlling various options in E-W, N-S,
C and vertical transport, respectively. 
C
C
C  _ORD=
C        1: 1st order upstream scheme (too diffusive, not a TYPE (XPLEX) option; it
C           can be used for debugging purposes; this is THE only known "linear"
C           monotonic advection scheme.).
C        2: 2nd order van Leer (full monotonicity constraint;
C           see Lin et al 1994, MWR)
C        3: monotonic PPM* (Collela & Woodward 1984)
C        4: semi-monotonic PPM (same as 3, but overshoots are allowed)
C        5: positive-definite PPM (constraint on the subgrid distribution is
C           only strong enough to prevent generation of negative values;
C           both overshoots & undershootes are possible).
C        6: un-constrained PPM (nearly diffusion free; faster but
C           positivity of the subgrid distribution is not quaranteed. Use
C           this option only when the fields and winds are very smooth or
C           when MFCT=.true.)
C        7: Huynh/Van Leer/Lin full monotonicity constraint
C Only KORD can be set to 7 to enable the use of Huynh's 2nd monotonicity
C constraint for piece-wise parabolic distribution.
C
C *PPM: Piece-wise Parabolic Method
C
C Recommended values:
C IORD=JORD=3 for high horizontal resolution.
C KORD=6 or 7  if MFCT=.true.
C KORD=3 or 7  if MFCT=.false.
C
C The implicit numerical diffusion decreases as _ORD increases.
C DO not use option 4 or 5 for non-positive definite scalars
C (such as Ertel Potential Vorticity).
C
C If numerical diffusion is a problem (particularly at low horizontal
C resolution) then the following setup is recommended:
C IORD=JORD=KORD=6 and MFCT=.true.
C
C AE: Radius of the sphere (meters).
C     Recommended value for the planet earth: 6.371E6
C
C FILL (logical):   flag to do filling for negatives (see note below).
C MFCT (logical):   flag to do a Zalesak-type Multidimensional Flux
C                   correction. It shouldn't be necessary to call the
C                   filling routine when MFCT is true.
C
C Umax: Estimate (upper limit) of the maximum U-wind speed (m/s).
C (225 m/s is a good value for troposphere model; 300 m/s otherwise)
C
C  *****************************************************************
C  **Input added for window calculation (2x2.5) (yxw, 8/21/01)******
C  *****************************************************************
C
C  I0_W, J0_W:      window index offset
C  (I1_W, J1_W):    left-low corner index of the window
C  (I2_W, J2_W);    right-high corner index of the window
C  IM_W, JM_W:      maximum index of the window (in coarse grid)
C  
C
C ======
C Output
C ======
C
C Q: the updated mixing ratios at t+NDT (original values are over-written)
C W(;;NL): large-scale vertical mass flux as diagnosed from the hydrostatic
C          relationship. W will have the same unit as PS1 and PS2 (eg, mb).
C          W must be divided by NDT to get the correct mass-flux unit.
C          The vertical Courant number C = W/delp_UPWIND, where delp_UPWIND
C          is the pressure thickness in the "upwind" direction. For example,
C          C(k) = W(k)/delp(k)   if W(k) > 0;
C          C(k) = W(k)/delp(k+1) if W(k) < 0.
C              ( W > 0 is downward, ie, toward surface)
C PS2: predicted PS at t+NDT (original values are over-written)
C
C Memory usage:
C This code is optimized for speed. it requres 18 dynamically allocated
C 3D work arrays (IM,JM,NL) regardless of the value of NC.
C Older versions (version 4 or 4.5) use less memory if NC is small.

C =====
C NOTES:
C =====
C
C This forward-in-time upstream-biased transport scheme degenerates to
C the 2nd order center-in-time center-in-space mass continuity eqn.
C if Q = 1 (constant fields will remain constant). This degeneracy ensures
C that the computed vertical velocity to be identical to GEOS-1 GCM
C for on-line transport.
C
C A larger polar cap is used if j1=3 (recommended for C-Grid winds or when
C winds are noisy near poles).
C
C The user needs to change the parameter Jmax or Kmax if the resolution
C is greater than 0.25 deg in N-S or 500 layers in the vertical direction.
C (this TransPort Core is otherwise resolution independent and can be used
C as a library routine).
 
C PPM is 4th order accurate when grid spacing is uniform (x & y); 3rd
C order accurate for non-uniform grid (vertical sigma coord.).
 
C Time step is limitted only by transport in the meridional direction.
C (the FFSL scheme is not implemented in the meridional direction).
 
C Since only 1-D limiters are applied, negative values could
C potentially be generated when large time step is used and when the
C initial fields contain discontinuities.
C This does not necessarily imply the integration is unstable.
C These negatives are typically very small. A filling algorithm is
C activated if the user set "fill" to be true.
C Alternatively, one can use the MFCT option to enforce monotonicity.

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h"

C ****6***0*********0*********0*********0*********0*********0**********72
      INTEGER, PARAMETER :: Jmax = 721, kmax = 200 
C  add ghost zone depth here (yxw, 08/23/01)       
C ****6***0*********0*********0*********0*********0*********0**********72
 
C Input-Output arrays
 
      TYPE (XPLEX) Q(IM,JM,NL,NC),PS1(IM,JM),PS2(IM,JM),W(IM,JM,NL),
     &     U(IM,JM,NL),V(IM,JM,NL),AP(NL+1),BP(NL+1)
      LOGICAL  ZCROSS, FILL, MFCT, deform
 
C Local dynamic arrays

!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite array dimension limits and rename with _W suffix
!%%% 
      TYPE (XPLEX) CRX_W(1-IGZD:IM_W+IGZD+1,1-IGZD:JM_W+IGZD,NL),
     &   CRY_W(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD+1,NL),
     &   DELP_W(IM_W,JM_W,NL),PT,AE,UMAX,
     &   XMASS_W(IM_W+1,JM_W,NL),YMASS_W(IM_W,JM_W+1,NL),
     &   DELP2_W(0:IM_W+1,0:JM_W+1,NL),
     &   DG1_W(IM_W),DG2_W(IM_W,0:JM_W+1),DPI_W(IM_W,JM_W,NL),
     &   QLOW_W(0:IM_W+1,0:JM_W+1,NL), DG3_W(JM_W, IM_W), 
     &   WK_W(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD,NL),
     &   PU_W(IM_W+1,JM_W,NL), 
     &   DQ_W(IM_W,JM_W,NL),
     &   DELP1_W(IM_W,JM_W,NL),
     &   FX_W(IM_W+1,JM_W,NL),FY_W(IM_W,JM_W+1,NL),
     &   FZ_W(IM_W,JM_W,NL+1),
     &   QZ_W(IM_W,JM_W,NL),QMAX_W(IM,JM,NL),QMIN_W(IM,JM,NL),
     &   U_W(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD,NL), 
     &   V_W(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD,NL),
     &   W_W(IM_W,JM_W,NL) 
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Add new logical variable OUT
!%%%
      LOGICAL OUT       
 
! bey, 6/20/00. for mass-flux diagnostic
      TYPE (XPLEX) fx1_tp_w(IM_w,JM_w,NL), fy1_tp_w(IM_w,JM_w,NL), 
     &     fz1_tp_w(IM_w,JM_w,NL) 
 
      INTEGER JS_w(NL),JN_w(NL),j1_in,j2_in
 
C Local static  arrays
 
      TYPE (XPLEX) DTDX_w(-10:Jmax), DTDX5_w(-10:Jmax), 
     &     acosp_w(-10:Jmax),
     &     cosp_w(-10:Jmax), cose_w(-10:Jmax), 
     &     DAP(kmax), DBK(kmax)
      INTEGER NDT0,NSTEP
      DATA NDT0, NSTEP /0, 0/
      DATA ZCROSS /.true./

C Saved internal variables:
      SAVE DTDY_w, DTDY5_w, JS0, JN0,  DTDX_w, out, j1_in,j2_in,
     &     DTDX5_w, acosp_w, COSP_w, COSE_w, DAP,DBK

! New variables for TPCORE pressure fixer (bdf, bmy, 10/11/01)
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Modify array dimension limits
!%%% 
      TYPE (XPLEX) YMASS_PF_W(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD+1,NL), 
     &     XMASS_PF_W(1-IGZD:IM_W+IGZD+1,1-IGZD:JM_W+IGZD,NL), 
     &     TEMP_W(IM,JM,NL),PI,DL,DP,DT,CR1,MaxDT,ztc,DTDY_w,DTDY5_w,
     &     D5,sum1,sum2

      ! Other new variables for window TPCORE pressure fixer (yxw)
      TYPE (XPLEX) DELP2_P(-IGZD:IM_W+IGZD+1, -IGZD:JM_W+IGZD+1, NL),
     &     PU_P(1-IGZD: IM_W+IGZD+1, 1-IGZD:JM_W+IGZD+1, NL)
      INTEGER jm1,j1,i0,j0,i1_w,i2_w,j1_w,j2_w,i0_w,j0_w,iord,jord,kord,
     & NDT,k,igd,j,iml,js0,jn0,i,ic,js,jn,jt,nc,im,jm,igzd,im_w,jm_w,nl
      LOGICAL :: PRESSURE_FIX
      PRESSURE_FIX = .TRUE.

C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================ 
      deform = .false.
      JM1 = 181 -1
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Comment out variables below
!%%%      IMH = IM/2
!%%%      j2= JM - j1 + 1
!%%% 
      NSTEP = NSTEP + 1  
 
C****6***0*********0*********0*********0*********0*********0**********72
C Initialization
C****6***0*********0*********0*********0*********0*********0**********72
 
      ! For mass flux diagnostics (bey, 6/20/00)
      fx1_tp_w(:,:,:) = 0d0
      fy1_tp_w(:,:,:) = 0d0
      fz1_tp_w(:,:,:) = 0d0
      
      ! Also need to initialize these arrays, so that the flux diagnostics 
      ! will be identical for single or multi processor (bmy, 9/29/00)
      fx_w(:,:,:) = 0d0
      fy_w(:,:,:) = 0d0
      fz_w(:,:,:) = 0d0

      ! Need to initialize these arrays in order to avoid
      ! floating-point exceptions on Alpha (lyj, bmy, 4/19/02)
      YMASS_PF_w(:,:,:) = 0d0
      XMASS_PF_w(:,:,:) = 0d0

      if(NSTEP.eq.1) then
 
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) 'TPCORE_WINDOW -- FFSL TransPort Core v. 7.1'
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' ) 'Originally written by S-J Lin'
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' )
     & 'Window version created for GEOS-CHEM by Yuxuan Wang and'
      WRITE( 6, '(a)' ) 
     & 'Bob Yantosca, with the addition of flux diagnostics and the'
      WRITE( 6, '(a)' ) 'DYN0 pressure fixer from M. Prather'
      WRITE( 6, '(a)' )
      WRITE( 6, '(a)' ) 'Last Modification Date: 3/13/03'
      WRITE( 6, '(a)' )
 
#if   ( multitask )
      WRITE( 6, '(a)' ) 'TPCORE_WINDOW was compiled for multitasking'
#if   defined( CRAY )
      WRITE( 6, '(a)' ) 'for CRAY'
#elif defined( SGI_MIPS  )
      WRITE( 6, '(a)' ) 'for SGI Origin/Power Challenge machines'
#elif defined( COMPAQ  )
      WRITE( 6, '(a)' ) 'for COMPAQ/HP RISC Alpha machines'
#elif defined( LINUX_PGI )
      WRITE( 6, '(a)' ) 'for Linux environment w/ PGI compiler'
#elif defined( LINUX_IFORT )
      WRITE( 6, '(a)' ) 'for Linux environment w/ Intel IFORT compiler'
#elif defined( SPARC )
      WRITE( 6, '(a)' ) 'for SUN/Sparc machines'
#elif defined( IBM_AIX )
      WRITE( 6, '(a)' ) 'for IBM/AIX machines'
#endif
#endif

      ! Added output on the first time TPCORE is called (bmy, 10/11/01)
      IF ( PRESSURE_FIX ) THEN
         WRITE( 6, '(a)' ) 'TPCORE PRESSURE FIXER is turned ON!'
      ENDIF

      if( MFCT ) then
          write(6,*) '   MFCT option is on'
      endif
 
      WRITE(6,*) 'IM=',IM,' JM=',JM,' NL=',NL,' j1=',j1
      WRITE(6,*) 'I0=',I0,' J0=',J0
C
C     write window size information (yxw, 8/21/2001)

      WRITE(6,*) 'IM_W=', IM_W, ' JM_W=', JM_W
      WRITE(6,*) 'I1_W=', I1_W, ' I2_W=', I2_W
      WRITE(6,*) 'J1_W=', J1_W, ' J2_W=', J2_W
      WRITE(6,*) 'I0_W=', I0_W, ' J0_W=', J0_W
      WRITE(6,*) NC, IORD,JORD,KORD,NDT
      
      if(NL.LT.6) then
        write(6,*) 'stop in module tpcore'
        write(6,*) 'NL must be >=6'
        stop
      endif
 
      if(Jmax.lt.JM .or. Kmax.lt.NL) then
        write(6,*) 'stop in module tpcore'
        write(6,*) 'Jmax or Kmax is too small; see documentation'
        stop
      endif
 
      DO k=1,NL
         DAP(k) = (AP(k+1) - AP(k))*PT
         DBK(k) =  BP(k+1) - BP(k)
      ENDDO
      
      PI = 4. * ATAN(1.)
      DL = 2.*PI / xplx(360)
      DP =    PI / xplx(JM1)
C
C     for window calculation, we have to redefine DL and DP. But since it's for
C     2x2.5, we skip this step for simplicity. (yxw, 8/21/2001)
C
      IF(IGD.EQ.0) THEN
C Compute analytic cosine at cell edges
C
C also need to change cosa for window calculation (yxw, 8/21/2001)
C 
         CALL COSA(COSP_W,COSE_W,JM_W,J0_W,PI,DP,IGZD,J0,JMAX)
      ELSE
C Define cosine consistent with GEOS-GCM (using dycore2.0 or later)
         CALL COSC(COSP_W,COSE_W,JM_W,J1_W,J2_W,PI,DP,IGZD,JMAX)
      ENDIF

      DO J=-IGZD,JM_W+IGZD+1
         ACOSP_W(J) = 1./COSP_W(J)
      ENDDO
c15    write (6,*) 'cosp(',j+j0_w+j0,')=',cosp_w(j)
 
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% We don't need this code, it's for the polar cap 
!%%% ! Inverse of the Scaled polar cap area.
!%%%      agle = (float(j1)-1.5)*DP
!%%%      RCAP  = DP / ( float(IM)*(1.-COS(agle)) )
!%%%      acosp(1)  = RCAP
!%%%      acosp(JM) = RCAP
!%%%
      ENDIF
 
      if(NDT0 .ne. NDT) then
      DT   = NDT
      NDT0 = NDT

      CR1  = abs(Umax*DT)/(DL*AE)
      MaxDT = DP*AE / abs(Umax) + 0.5
      write(6,*)'Largest time step for max(V)=',Umax,' is ',MaxDT

      if(MaxDT .lt. abs(NDT)) then
         write(6,*) 'Warning!!! NDT maybe too large!'
      endif

      if(CR1.ge.0.95) then
         JS0 = J1_w-igzd
     	 JN0 = J2_w+igzd       !(yxw,eulerian)
         IML = IM-2
         ZTC = 0.
      else
         ZTC = acos(CR1) * (180./PI)
         JS0 = xplx(JM1)*(90.-ZTC)/180. + 2
         JS0 = max(JS0, J1+1)
         IML = min(6*JS0/(J1-1)+2, 4*360/5)
         JN0 = 181-JS0+1
      endif
 
      WRITE(6,*) 'IGZD = ', IGZD
      write(6,*) 'ZTC= ',ZTC,' JS= ',JS0,' JN= ',JN0,' IML= ',IML
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% determine the relationship of (JS0,JN0) and (J1_W, J2_W)
!%%% (ji_in,j2_in) are the part of the window inside the region (JS0,JN0)
!%%%
      CALL POSITION_WINDOW( JS0,       JN0, J1_W-IGZD,
     &                      J2_W+IGZD, OUT, J1_IN,    J2_IN ) 

      WRITE(6,*) 'J1_IN=', J1_IN,' ', 'J2_IN=', J2_IN

      DO J = 1-IGZD, JM_W+IGZD
         DTDX_W(J)  = DT / ( DL*AE*COSP_W(J) )
         DTDX5_W(J) = 0.5 * DTDX_W(J)
      ENDDO
      
      DTDY_w  = DT /(AE*DP)
      DTDY5_w = 0.5*DTDY_w

      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      ENDIF              ! END INITIALIZATION.
 
C****6***0*********0*********0*********0*********0*********0**********72
C Compute Courant number
C****6***0*********0*********0*********0*********0*********0**********72

      if(IGD.eq.0) then
 
C Convert winds on A-Grid to Courant # on C-Grid.

#if   defined( multitask  )
#if   defined( CRAY       ) 
CMIC$ do all shared(NL,im,jm1,jm,U,V,dtdx5_w,dtdy5_w,CRX_w,CRY_w,im_w,jm_w)
CMIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      do 900 k=1,NL
      DO J=1-IGZD,JM_W+IGZD
      DO I=1-IGZD,IM_W+1+IGZD
C     calculate Courant # at grid edge, using offsetted wind (yxw 08/23/01)
C
         CRX_w(i,j,k) = dtdx5_w(j)*(U(i+i0_w,j+j0_w,k)+
     &        U(i-1+i0_w,j+j0_w,k))
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Comment out
!%%%      if (CRX_w(i,j,k) .gt. 1) then
!%%%      write(6,555) i,j,k, CRX_w(i,j,k)
!%%%555   format('CRX is larger than 1 at grid ', 3(I3,1x), F12.5)
!%%%      endif
!%%%
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Replace w/ code below
!%%% ! for i=1
!%%%      do 48 j=2,JM1
!%%% 48    CRX(1,j,k) = dtdx5(j)*(U(1,j,k)+U(IM,j,k))
!%%%
      DO J = 1-IGZD, JM_W+1+IGZD
      DO I = 1-IGZD, IM_W+IGZD
       CRY_w(i,j,k) = DTDY5_w*(V(i+i0_w,j+j0_w,k)+V(i+i0_w,j-1+j0_w,k))
      ENDDO
      ENDDO

900   continue

      else

C Convert winds on C-grid to Courant #
C Beware of the index shifting!! (GEOS-GCM)

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all shared(NL,im,jm1,jm,U,V,dtdx_w,dtdy_w,CRX_w,CRY_w,jm_w,im_w)
CMIC$* private(i,j,k)
#else 
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

C  calculte Courant # using offsetted wind (yxw,08/23/01)
C
      DO 65 k=1,NL

      DO J =1-IGZD, JM_W+IGZD
      DO I =1-IGZD, IM_W+1+IGZD
         CRX_W(I,J,K) = DTDX_W(J) * U( I-1+I0_W, J+J0_W, K )
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Replace w/ code below
!%%%      do 55 j=2,JM1
!%%% 55    CRX(1,j,k) = dtdx(j)*U(IM,j,k)
!%%%
      DO J = 1-IGZD, JM_W+1+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         CRY_W(I,J,K) = DTDY_W * V( I+I0_W, J+1-J0_W, K )
      ENDDO
      ENDDO

 65   continue
      endif

      !=================================================================
      ! *****  T P C O R E   P R E S S U R E   F I X E R  *****
      !
      ! Run pressure fixer to fix mass conservation problem.  Pressure
      ! fixer routines PRESS_FIX, DYN0, PFILTR, LOCFLT, and POLFLT
      ! change the mass fluxes so they become consistant with met field
      ! pressures.  (bdf, bmy, 10/11/01)
      !
      ! NOTE: The pressure fixer is not 100% perfect; tracer mass will
      !       increase on the order of 0.5%/yr.  However, this is much
      !       better than w/o the pressure fixer, where the mass may
      !       increase by as much as 40%/yr.  (bdf, bmy, 10/22/01)
      !=================================================================
      IF ( PRESSURE_FIX ) THEN

         ! Loop over vertical levels --
         ! added parallel loop #if statements (bmy, 10/11/01)
#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all shared(NL,IM,JM,JM1,COSE,XMASS_PF,YMASS_PF,DELP2,CRX,CRY)
CMIC$* private(I,J,K,D5)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, D5 )
#endif
#endif
         DO K = 1, NL

            ! DELP = pressure thickness:
            ! the pseudo-density in a hydrostatic system.
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits 
!%%%            DO J = 1, JM
!%%%            DO I = 1, IM 
!%%%
            DO J = -IGZD, JM_W+IGZD+1
            DO I = -IGZD, IM_W+IGZD+1
               DELP2_P(I,J,K) = DAP(K) + DBK(K)*PS2(I+I0_W,J+J0_W)
            ENDDO
            ENDDO

            ! calculate mass fluxes for pressure fixer.

            ! N-S component
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%            DO J = J1, J2+1
!%%%
            DO J=1-IGZD, JM_W+IGZD+1
               D5 = 0.5 * COSE_W(J)
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%               DO I = 1, IM
!%%%
               DO I= 1-igzd, IM_W+igzd
                  YMASS_PF_w(I,J,K) =
     &                 CRY_w(I,J,K) * D5 * (DELP2_p(I,J,K)+
     &                  DELP2_p(I,J-1,K))
               ENDDO
            ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff is useless for nested simulation
!%%%            ! Enlarged polar cap.
!%%%            IF(J1.NE.2) THEN
!%%%               DO I=1,IM
!%%%                  YMASS_PF(I,1,K) = 0
!%%%                  YMASS_PF(I,JM1+1,K) = 0
!%%%               ENDDO
!%%%            ENDIF
!%%%
            ! E-W component
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%            DO J = J1, J2
!%%%            DO I =  2, IM
!%%%
            DO J = 1-IGZD, JM_W+IGZD+1
            DO I = 1-IGZD, IM_W+IGZD+1
               PU_P(I,J,K) = 0.5 * (DELP2_P(I,J,K) + DELP2_P(I-1,J,K))
            ENDDO
            ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap is useless for nested simulation 
!%%%            DO J = J1, J2
!%%%               PU(1,J,K) = 0.5 * (DELP2(1,J,K) + DELP2(IM,J,K))
!%%%            ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits 
!%%%            DO J = J1, J2
!%%%
            DO J = 1-IGZD, JM_W+IGZD
            DO I = 1-IGZD, IM_W+IGZD+1
               XMASS_PF_W(I,J,K) = PU_P(I,J,K) * CRX_W(I,J,K)
            ENDDO
            ENDDO

         ENDDO

         !==============================================================
         ! Call PRESS_FIX to apply the pressure fix to the mass fluxes
         ! XMASS_PF, YMASS_PF.  PRESS_FIX will call routine DYN0, etc.
         !==============================================================
         CALL PRESS_FIX( XMASS_PF_w, YMASS_PF_w, NDT,  ACOSP_w, Jmax,
     &                   I0_W,       J0_W,       IM_W, JM_W,    IGZD )

         ! Loop over vertical levels --
         ! added parallel loop #if statements (bmy, 10/11/01)
#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all shared(NL,IM,JM,XMASS_PF,PU,YMASS_PF,DELP2,COSE,CRX,CRY)
CMIC$* private(I,J,K,D5)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, D5 )
#endif
#endif
         DO K = 1, NL

            ! Recreate the CRX variable with the new values
            ! of XMASS_PF, which has been adjusted by DYN0
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%            DO J = J1, J2
!%%%            DO I =  1, IM
!%%%
            DO J = 1-IGZD, JM_W+IGZD
            DO I = 1-IGZD, IM_W+IGZD+1
               CRX_w(I,J,K) = XMASS_PF_w(I,J,K) / PU_p(I,J,K)
            ENDDO
            ENDDO

            ! Recreate the CRY variable with the new values
            ! of YMASS_PF, which has been adjusted by DYN0
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%            DO J = J1, J2+1
!%%%
            DO J = 1-IGZD, JM_W+IGZD+1
               D5 = 0.5 * COSE_W(J)

               DO I = 1-IGZD, IM_W+IGZD
                  CRY_W(I,J,K) = YMASS_PF_W(I,J,K) /
     &                 ( D5 * ( DELP2_P(I,J,K) + DELP2_P(I,J-1,K) ) )
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! End of TPCORE PRESSURE FIXER -- continue as usual
      !=================================================================


C**********************************************************************
C Check whether CRY_w is larger than one  (yxw, eulerian)
C********************************************************************** 
      IF ( maxval(CRY_w) .gt. 1.0) THEN
         write (6,*) 'CRY is larger than one!!'
         write (6,*) 'Decrease timestep NTDT!!'
         STOP
      ENDIF
 
C****6***0*********0*********0*********0*********0*********0**********72
C Find JN and JS
C****6***0*********0*********0*********0*********0*********0**********72
 
#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope shared(JS_w,JN_w,CRX_w,CRY_w,PS2,U,V,DPI_w,ymass_w,delp2_w,PU_W)
CMIC$* shared(xmass_w)
CMIC$* private(i,j,k,sum1,sum2,D5)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, SUM1, SUM2, D5 )
#endif
#endif

!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Find JS_w and JN_w, using CRX_w as criterion  
!%%% here, j1_in and j2_in have been calculated in the initialization part 
!%%%
      do 1000 k=1,NL
      if (.not. out) then 
      JS_w(k) = j1_w-igzd-1
      JN_w(k) = j2_w+igzd+1
 
      DO J = J1_IN, J1_W-IGZD,-1
      DO I = 1,     IM_W
         IF(ABS(CRX_W(I,J-J0_W-J0,K)) .GT. 1.) THEN
            JS_W(K) = J
            GO TO 112
         ENDIF
      ENDDO
      ENDDO

112   continue
 
      DO J = J2_IN, J2_W+IGZD
      DO I = 1,     IM_W
         IF(ABS(CRX_W(I,J-J0_W-J0,K)) .GT. 1.) THEN
            JN_W(K) = J
            GO TO 133
         ENDIF
      ENDDO
      ENDDO

133   continue
      
      else
      js_w(k)=-1 
      jn_w(k)=-1 

      DO J = J1_W, J2_W
      DO I = 1,    IM_W
         IF (ABS(CRX_W(I,J-J0_W-J0,K)) .LT. 1.) THEN
            JS_W(K)=J
            GO TO 134
         ENDIF
      ENDDO
      ENDDO

134   continue
    
      IF (JS_W(K) .NE. -1) THEN 

      DO J = JS_W(K), J2_W
      DO I = 1,       IM_W
         IF (ABS(CRX_W(I,J-J0_W-J0,K)) .GT. 1.) THEN
            JN_W(K)=J
            GO TO 146
         ENDIF
      ENDDO
      ENDDO

146   continue
      endif
      endif

C****6***0*********0*********0*********0*********0*********0**********72
C ***** Compute horizontal mass fluxes *****
C****6***0*********0*********0*********0*********0*********0**********72
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% now it is the window version,using offsetted pressure fields 
!%%% Rewrite DO-loop limits as necessary 
!%%% 
C delp = pressure thickness: the psudo-density in a hydrostatic system.
      DO J = 0, JM_W+1
      DO I = 0, IM_W+1
         DELP2_W(I,J,K) = DAP(K) + DBK(K)*PS2(I+I0_W,J+J0_W)
      ENDDO
      ENDDO

C N-S componenet
 
      do j=1,jm_w+1
      D5 = 0.5 * COSE_w(j)
      do i=1,IM_w
      ymass_w(i,j,k) = CRY_w(i,j,k)*D5*(delp2_w(i,j,k) + 
     &                   delp2_w(i,j-1,k))
      enddo
      enddo
 
      DO J = 1, JM_W
      DO I = 1, IM_W
         DPI_W(I,J,K) = (YMASS_W(I,J,K)-YMASS_W(I,J+1,K)) * ACOSP_W(J)
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%      if(j1.ne.2) then           ! Enlarged polar cap.
!%%%      do 95 i=1,IM
!%%%      DPI(i,  2,k) = 0.
!%%% 95    DPI(i,JM1,k) = 0.
!%%%      endif
!%%%
!%%% Poles
!%%%      sum1 = ymass(IM,j1  ,k)
!%%%     sum2 = ymass(IM,j2+1,k)
!%%%      do 98 i=1,IM-1
!%%%      sum1 = sum1 + ymass(i,j1  ,k)
!%%% 98    sum2 = sum2 + ymass(i,j2+1,k)
!%%% 
!%%%      sum1 = - sum1 * RCAP
!%%%      sum2 =   sum2 * RCAP
!%%%      do 100 i=1,IM
!%%%      DPI(i, 1,k) = sum1
!%%% 100   DPI(i,JM,k) = sum2
 
C E-W component
      do J = 1, JM_W
      DO I = 1, IM_W+1
         PU_W(I,J,K) = 0.5 * (DELP2_W(I,J,K) + DELP2_W(I-1,J,K))
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%      do j=j1,j2
!%%%      PU(1,j,k) = 0.5 * (delp2(1,j,k) + delp2(IM,j,k))
!%%%      enddo
!%%% 
      DO j = 1, jm_w
      DO i = 1, IM_w+1
         xmass_w(i,j,k) = PU_w(i,j,k)*CRX_w(i,j,k)
      ENDDO
      ENDDO

      DO j = 1, jm_w
      DO i = 1, IM_w
         DPI_w(i,j,k) = DPI_w(i,j,k) + xmass_w(i,j,k) - xmass_w(i+1,j,k)
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%
!%%%    DO 130 j=j1,j2
!%%% 130   DPI(IM,j,k) = DPI(IM,j,k) + xmass(IM,j,k) - xmass(1,j,k)
!%%% 
C****6***0*********0*********0*********0*********0*********0**********72
C Compute Courant number at cell center
C****6***0*********0*********0*********0*********0*********0**********72
 
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         IF(CRX_W(I,J,K)*CRX_W(I+1,J,K) .GT. 0.) THEN
            IF(CRX_W(I,J,K) .GT. 0.) THEN
               U_W(I,J,K) = CRX_W(I,J,K)
            ELSE
               U_W(I,J,K) = CRX_W(I+1,J,K)
            ENDIF
         ELSE
            U_W(I,J,K) = 0.
         ENDIF
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%      i=IM
!%%%      DO 136 j=2,JM1
!%%%      if(CRX(i,j,k)*CRX(1,j,k) .gt. 0.) then
!%%%         if(CRX(i,j,k) .gt. 0.) then
!%%%         U(i,j,k) = CRX(i,j,k)
!%%%         else
!%%%         U(i,j,k) = CRX(1,j,k)
!%%%         endif
!%%%      else
!%%%         U(i,j,k) = 0.
!%%%      endif
!%%% 136   continue
!%%% 
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         IF(CRY_W(I,J,K)*CRY_W(I,J+1,K) .GT. 0.) THEN
            IF(CRY_W(I,J,K) .GT. 0.) THEN
               V_W(I,J,K) = CRY_W(I,J,K)
            ELSE
               V_W(I,J,K) = CRY_W(I,J+1,K)
            ENDIF
         ELSE
            V_W(I,J,K) = 0.
         ENDIF
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%      do 139 i=1,IMH
!%%%      V(i,     1,k) = 0.5*(CRY(i,2,k)-CRY(i+IMH,2,k))
!%%%      V(i+IMH, 1,k) = -V(i,1,k)
!%%%      V(i,    JM,k) = 0.5*(CRY(i,JM,k)-CRY(i+IMH,JM1,k))
!%%% 139   V(i+IMH,JM,k) = -V(i,JM,k)
!%%%
1000  continue
!C****6***0*********0*********0*********0*********0*********0**********72
!C Compute vertical mass flux (same dimensional unit as PS)
!C****6***0*********0*********0*********0*********0*********0**********72
 
C compute total column mass CONVERGENCE.

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope shared(im,jm,DPI_w,PS1,PS2,W_w,DBK,jm_w,im_w)
CMIC$* shared(DPI_w,PS1,PS2,W_w,DBK)
CMIC$* private(i,j,k,DG1_w)
#else 
!$OMP PARALLEL DO PRIVATE( I, J, K, DG1_w )
#endif
#endif

!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits as necessary
!%%%
      do 395 j=1,jm_w

      do i=1,IM_w
         DG1_w(i) = DPI_w(i,j,1)
      ENDDO
 
      do k=2, NL
      do i=1, IM_w
         DG1_w(i)  = DG1_w(i) + DPI_w(i,j,k)
      ENDDO
      ENDDO

      do 360 i=1,IM_w
 
C Compute PS2 (PS at n+1) using the hydrostatic assumption.
C Changes (increases) to surface pressure = total column mass convergence
         PS2(i+i0_w,j+j0_w)  = PS1(i+i0_w,j+j0_w) + DG1_w(i)
 
C compute vertical mass flux from mass conservation principle.
         W_w(i,j,1) = DPI_w(i,j,1) - DBK(1)*DG1_w(i)
         W_w(i,j,NL) = 0.
360   continue
 
      DO K = 2, NL-1
      DO I = 1, IM_W
         W_W(I,J,K) = W_W(I,J,K-1) + DPI_W(I,J,K) - DBK(K)*DG1_W(I)
      ENDDO
      ENDDO

395   continue

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all
CMIC$* shared(deform,NL,im,jm,delp_w,delp1_w,delp2_w,DPI_w,DAP,DBK,PS1,PS2)
CMIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO 390 k=1,NL

      DO J =1, JM_W
      DO I =1, IM_W
         DELP1_W(I,J,K) = DAP(K) + DBK(K)*PS1(I+I0_W,J+J0_W)
         DELP2_W(I,J,K) = DAP(K) + DBK(K)*PS2(I+I0_W,J+J0_W)
         DELP_W (I,J,K) = DELP1_W(I,J,K) + DPI_W(I,J,K)
      ENDDO
      ENDDO
 
C Check deformation of the flow fields
      if(deform) then

!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Don't need this part 
!%%%      DO j=1,JM_w
!%%%      DO i=1,IM_w
!%%%      if(delp_w(i,j,k) .le. 0.) then
!%%%        write(6,*) k,'Noisy wind fields -> delp* is negative!'
!%%%        write(6,*) ' *** Smooth the wind fields or reduce NDT'
!%%%         stop
!%%%      endif
!%%%      ENDDO
!%%%      ENDDO

      endif
390   continue

C****6***0*********0*********0*********0*********0*********0**********72
C Do transport one tracer at a time.
C****6***0*********0*********0*********0*********0*********0**********72
 
      DO 5000 IC=1,NC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% MODIFICATIONS FOR NESTED GRID (bmy, 11/5/08)
!%%% Comment out the parallel loop statements for DO loop 2500, because
!%%% for some reason this causes NaN's.  We can leave all the other parallel
!%%% loops activated. (bmy, 11/5/08)
!%%%
!%%%#if   defined( multitask  )
!%%%#if   defined( CRAY       )
!%%%!CMIC$ do all autoscope
!%%%!CMIC$* shared(q,DQ_w,delp1_w,U_w,V_w,j1,JS_w,JN_w,im,jm,IML,IC,IORD,JORD,jm_w,im_w)
!%%%!CMIC$* shared(CRX_w,CRY_w,PU_w,xmass_w,ymass_w,fx_w,fy_w,acosp_w,qz_w)
!%%%!CMIC$* shared(fx1_tp_w, fy1_tp_w)
!%%%!CMIC$* private(i,j,k,jt,wk_w,DG2_w)
!%%%#else 
!%%%!$OMP PARALLEL DO PRIVATE( I, J, K, JT, WK_w, DG2_w )
!%%%#endif
!%%%#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      do 2500 k=1,NL
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%
!%%%      if(j1.ne.2) then
!%%%      DO 405 I=1,IM
!%%%      q(I,  2,k,IC) = q(I, 1,k,IC)
!%%% 405   q(I,JM1,k,IC) = q(I,JM,k,IC)
!%%%      endif
!%%%

C Initialize DQ
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% use offsetted tracer concentration Q, and save into DQ_W
!%%%
         !print *, "IC, K:", IC, K
      DO J = 1, JM_W
      DO I = 1, IM_W
         DQ_W(I,J,K) = Q(I+I0_W,J+J0_W,K,IC)*DELP1_W(I,J,K)
      ENDDO
      ENDDO

C E-W advective cross term 
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Pass new arguments as necessary (e.g. IM_W, JM_W, etc)
!%%%
      CALL XADV( IM_W,    JM_W,    Q(:,:,K,IC), U_W(:,:,K),
     &           JS_W(K), JN_W(K), WK_W(:,:,1), IGZD,
     &           IM,      JM,      I0_W,        J0_W,
     &           I0,      J0 )

      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         WK_W(I,J,1) = Q(I+I0_W,J+J0_W,K,IC) + 0.5*WK_W(I,J,1)
      ENDDO
      ENDDO
 
C N-S advective cross term
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits as necessary
!%%%
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         JT = XPLX(J+J0_W) - V_W(I,J,K)
         IF (JT .GT. JM-1) THEN
            JT=JM-1             
         ELSE IF (JT .LT. 1) THEN
            JT=1                
         ENDIF 
         WK_W(I,J,2) = V_W(I,J,K) * 
     &        (Q(I+I0_W,JT,K,IC) - Q(I+I0_W,JT+1,K,IC))
      ENDDO
      ENDDO
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits as necessary
!%%%
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         WK_W(I,J,2) = Q(I+I0_W,J+J0_W,K,IC) + 0.5*WK_W(I,J,2)
      ENDDO
      ENDDO

C****6***0*********0*********0*********0*********0*********0**********72
C compute flux in  E-W direction
C Return flux contribution from TPCORE in FX1_TP array (bey, 9/28/00)
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Call XTP, with window arrays etc.
!%%%
      call xtp(im_w,jm_w,igzd,JN_w(k),JS_w(k),PU_w(:,:,k),DQ_w(:,:,k),
     &       wk_w(:,:,2),CRX_w(:,:,k),fx_w(:,:,k),xmass_w(:,:,k),IORD,
     &       fx1_tp_w(:,:,k),i0_w,j0_w,i0,J0)

C compute flux in  N-S direction
C Return flux contribution from TPCORE in FY1_TP array (bey, 9/28/00)
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Call YTP with window arrays etc.
!%%%
      call ytp(IM_w,JM_w,acosp_w(:),DQ_w(:,:,k),wk_w(:,:,1),
     &         CRY_w(:,:,k),ymass_w(:,:,k),fy_w(:,:,k),JORD,
     &         fy1_tp_w(:,:,k),igzd, Jmax)
!C****6***0*********0*********0*********0*********0*********0**********72

      if(ZCROSS) then

C qz is the horizontal advection modified value for input to th=1,5
C vertical transport operator FZPPM
C Note: DQ contains only first order upwind contribution.

         DO J = 1, JM_W
         DO I = 1, IM_W
            QZ_W(I,J,K) = DQ_W(I,J,K) / DELP_W(I,J,K)
         ENDDO
         ENDDO

      ELSE

         DO J = 1, JM_W
         DO I = 1, IM_W
            QZ_W(I,J,K) = Q(I+I0_W,J+J0_W,K,IC)
         ENDDO
         ENDDO

      ENDIF

2500  continue     ! k-loop
C****6***0*********0*********0*********0*********0*********0**********72
C Compute fluxes in the vertical direction
C Return flux contribution from FZPPM in FZ1_TP for ND26 (bey, 9/28/00)
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Call FZPPM with window arrays
!%%%
      call FZPPM(qz_w,fz_w,IM_w,JM_w,NL,DQ_w,
     &               W_w,delp_w,KORD,fz1_tp_w)

!C****6***0*********0*********0*********0*********0*********0**********72

C Final update

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope
CMIC$* private(i,j,k,sum1,sum2)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, SUM1, SUM2 )
#endif
#endif 

      DO 101 K = 1, NL

      DO J = 1, JM_W
      DO I = 1, IM_W
         DQ_W(I,J,K) = DQ_W(I,J,K)   
     &               +   FX_W(I,J,K) - FX_W(I+1,J,K)
     &               + ( FY_W(I,J,K) - FY_W(I,J+1,K) ) * ACOSP_W(J)
     &               +   FZ_W(I,J,K) - FZ_W(I,J,K+1)
      ENDDO
      ENDDO

!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%      sum1 = fy(IM,j1  ,k)
!%%%      sum2 = fy(IM,J2+1,k)
!%%%
!%%%      do i=1,IM-1
!%%%         sum1 = sum1 + fy(i,j1  ,k)
!%%%         sum2 = sum2 + fy(i,J2+1,k)
!%%%      enddo
!%%% 
!%%%      DQ(1, 1,k) = DQ(1, 1,k) - sum1*RCAP + fz(1, 1,k) - fz(1, 1,k+1)
!%%%      DQ(1,JM,k) = DQ(1,JM,k) + sum2*RCAP + fz(1,JM,k) - fz(1,JM,k+1)
!%%% 
!%%%      do i=2,IM
!%%%      DQ(i, 1,k) = DQ(1, 1,k)
!%%%      DQ(i,JM,k) = DQ(1,JM,k)
!%%%      enddo
!%%%
101   continue
C****6***0*********0*********0*********0*********0*********0**********72
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Call QCKXYZ w/ window arrays, etc.
!%%%
      if(FILL) call qckxyz(DQ_w,DG3_w,IM_w,JM_w,NL,cosp_w,
     &          acosp_w,IC,NSTEP,DP,igzd, Jmax,fx_w, fy_w, fz_w )

      !=================================================================
      ! bey, 6/20/00. for mass-flux diagnostic
      ! NOTE: DIAG_FLUX is not called within a parallel loop, 
      ! so parallelization can be done within the subroutine
      !=================================================================
      CALL DIAG_FLUX( IC,   FX_w,     FX1_TP_w, FY_w,    FY1_TP_w, 
     &                FZ_w, FZ1_TP_w, NDT,      ACOSP_w, Jmax,
     &                I0_W, J0_W,     IM_W,     JM_W,    IGZD )

!C****6***0*********0*********0*********0*********0*********0**********72

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all
CMIC$* shared(q,IC,NL,j1,im,jm,jm1,DQ_w,delp2_w)
CMIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO 920 k=1,NL

         DO J = 1, JM_W
         DO I = 1, IM_W
            Q(I+I0_W,J+J0_W,K,IC) = DQ_W(I,J,K) / DELP2_W(I,J,K)
         ENDDO
         ENDDO
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%
!%%% don't need the following polar part (yxw,08/28/01)
!%%%      if(j1.ne.2) then
!%%%      DO 450 I=1,IM
!%%%      Q(I,  2,k,IC) = Q(I, 1,k,IC)
!%%%      Q(I,JM1,k,IC) = Q(I,JM,k,IC)
!%%% 450   CONTINUE
!%%%      endif
!%%%
 920  CONTINUE         
 5000 CONTINUE

      ! Return to calling program
      END SUBROUTINE TPCORE_WINDOW

!------------------------------------------------------------------------------

      subroutine cosa(cosp,cose,jm_w,j0_w,PI,DP,igzd,j0,Jmax)  
!cosa_w (yxw, 8/23/2001)
C****6***0*********0*********0*********0*********0*********0**********72
      INTEGER jm_w,j0_w,j0,Jmax,igzd,j
      TYPE (XPLEX) cosp(-10:Jmax),cose(-10:Jmax),
     &      sine(-10:jm_w+igzd+2),PI,DP,ph5
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop boundaries  
!%%%
      DO J = -IGZD, JM_W+2+IGZD    
         PH5     = -0.5*PI + (XPLX(J+J0_W+J0-1)-0.5)*DP
         SINE(J) = SIN(PH5)
      ENDDO
      
      DO J = -IGZD, JM_W+IGZD+1
         COSP(J) = (SINE(J+1)-SINE(J))/DP
      ENDDO
 
C Define cosine at edges..

      DO J = 1-IGZD, JM_W+IGZD+1
         COSE(J) = 0.5 * (COSP(J-1)+COSP(J))
      ENDDO

      ! Return to TPCORE
      END SUBROUTINE COSA

!------------------------------------------------------------------------------

      subroutine cosc(cosp,cose,jm_w,J1_w,j2_w,PI,DP,igzd,Jmax) 
C****6***0*********0*********0*********0*********0*********0**********72
      INTEGER jm_w,j2,j,j2_w,igzd,Jmax,j1_w
      TYPE (XPLEX) cosp(-10:Jmax),cose(-10:Jmax),PI,DP,phi
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop boundaries  
!%%%
      PHI = -0.5*PI+(J1_W-IGZD-3)*DP

      DO J = -IGZD, JM_W+1+IGZD
         PHI     =  PHI + DP
         COSP(J) = COS(PHI)
      ENDDO

      DO J = 1-IGZD, JM_W+IGZD+1
         COSE(J) = 0.5*(COSP(J)+COSP(J-1))
      ENDDO

      DO J = 1-IGZD, JM_W+IGZD
         COSP(J) = 0.5*(COSE(J)+COSE(J+1))
      ENDDO

      ! Return to TPCORE
      END SUBROUTINE COSC

!------------------------------------------------------------------------------

      subroutine filew(q,qtmp,IMR,JNP,ipx,tiny, fx)
      INTEGER i,j,imr,jnp,ipx
C****6***0*********0*********0*********0*********0*********0**********72
      TYPE (XPLEX) q(IMR,*),qtmp(JNP,IMR), fx(IMR+1,JNP),d0,d1,d2,
     & tiny
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C
      ipx = 0
C Copy & swap direction for vectorization.
      do i=1,imr
      do j=1,jnp
         qtmp(j,i) = q(i,j)
      ENDDO
      ENDDO

      DO I = 2, IMR-1
      DO J = 1, JNP
         IF(QTMP(J,I).LT.0.) THEN
            IPX =  1
! WEST
            D0 = MAX(0.,QTMP(J,I-1))
            D1 = MIN(-QTMP(J,I),D0)
            QTMP(J,I-1) = QTMP(J,I-1) - D1
            QTMP(J,I) = QTMP(J,I) + D1
            FX(I,J) = FX(I,J)+D1 !(yxw,02/09/2003)
! EAST
            D0 = MAX(0.,QTMP(J,I+1))
            D2 = MIN(-QTMP(J,I),D0)
            QTMP(J,I+1) = QTMP(J,I+1) - D2
            QTMP(J,I) = QTMP(J,I) + D2 + TINY
            FX(I+1,J) = FX(I+1,J) - D2 !(yxw,02/09/2003)
         ENDIF
      ENDDO
      ENDDO

      I=1
      do 65 j=1,JNP
      if(qtmp(j,i).lt.0.) then
      ipx =  1
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (bmy, 3/10/03)
!%%% Comment this out
!%%%c west
!%%%      d0 = max(0.,qtmp(j,imr))
!%%%      d1 = min(-qtmp(j,i),d0)
!%%%      qtmp(j,imr) = qtmp(j,imr) - d1
!%%%      qtmp(j,i) = qtmp(j,i) + d1
!%%%
c east
      d0 = max(0.,qtmp(j,i+1))
      d2 = min(-qtmp(j,i),d0)
      qtmp(j,i+1) = qtmp(j,i+1) - d2
      qtmp(j,i) = qtmp(j,i) + d2 + tiny
      fx(i+1, j) = fx(i+1,j) - d2   !(yxw,02/09/2003)
      endif
65    continue
      i=IMR
      do 75 j=1,JNP
      if(qtmp(j,i).lt.0.) then
      ipx =  1
c west
      d0 = max(0.,qtmp(j,i-1))
      d1 = min(-qtmp(j,i),d0)
      qtmp(j,i-1) = qtmp(j,i-1) - d1
      qtmp(j,i) = qtmp(j,i) + d1+tiny
      fx(i,j) = fx(i,j) + d1   !(yxw,02/09/2003)
c east
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (bmy, 3/10/03)
!%%% Comment this out
!%%%      d0 = max(0.,qtmp(j,1))
!%%%      d2 = min(-qtmp(j,i),d0)
!%%%      qtmp(j,1) = qtmp(j,1) - d2
!%%%
!%%%      qtmp(j,i) = qtmp(j,i) + d2 + tiny
!%%%
      endif            
75    continue
C
      if(ipx.ne.0) then
      do 85 j=1,jnp
      do 85 i=1,imr
85    q(i,j) = qtmp(j,i)
C      else
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%C Pole
!%%%      if(q(1,1).lt.0. or. q(1,JNP).lt.0.) ipx = 1
!%%%
      endif
      return
      end subroutine filew

!------------------------------------------------------------------------------

      subroutine filns(q,IMR,JNP,cosp,acosp,ipy,tiny,DP,fy,Jmax)
C****6***0*********0*********0*********0*********0*********0**********72
      INTEGER IMR,JNP,ipy,JMAX,i,j,ic
      TYPE (XPLEX) DP,dq,dn,d0,d1,ds,d2,tiny
      TYPE (XPLEX) q(IMR,*),cosp(-10:Jmax),acosp(-10:Jmax),
     & fy(IMR,JNP+1)
      LOGICAL first
      DATA first /.true./
C     SAVE cap1
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C
      if(first) then
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%      DP = 4.*ATAN(1.)/float(JNP-1)
!%%%      cap1 = IMR*(1.-COS((j1-1.5)*DP))/DP
!%%%
      first = .false.
      endif
C
      ipy = 0
      do 55 j=2,jNP-1
      DO 55 i=1,IMR
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
C North
      dn = q(i,j+1)*cosp(j+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,j+1) = (dn - d1)*acosp(j+1)
      fy(i,j+1) = fy(i,j+1)-d1  !(yxw,02/09/2003)
      dq = dq - d1
C South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      fy(i,j)= fy(i,j) + d2   !(yxw,02/09/2003)
      
      endif
55    continue
C
      
      do i=1,imr
      IF(q(i,1).LT.0.) THEN
      ipy =  1
      dq  = - q(i,1)*cosp(1)
C North
      dn = q(i,1+1)*cosp(1+1)
      d0 = max(0.,dn)
      d1 = min(dq,d0)
      q(i,1+1) = (dn - d1)*acosp(1+1)
      q(i,1) = (d1 - dq)*acosp(1) + tiny
      fy(i,1+1) = fy(i,1) - d1  !(yxw,02/09/2003)
      endif

      enddo

      j = JNP 
      do i=1,imr
      IF(q(i,j).LT.0.) THEN
      ipy =  1
      dq  = - q(i,j)*cosp(j)
C South
      ds = q(i,j-1)*cosp(j-1)
      d0 = max(0.,ds)
      d2 = min(dq,d0)
      q(i,j-1) = (ds - d2)*acosp(j-1)
      q(i,j) = (d2 - dq)*acosp(j) + tiny
      fy(i,j) = fy(i,j) + d2   !(yxw,02/09/2003)
      endif

      enddo         ! (yxw,09/26/01)

!%%%
!%%% MODIFICATIONS FOR NESTED GRID (bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%C Check Poles.
!%%%      if(q(1,1).lt.0.) then
!%%%      dq = q(1,1)*cap1/float(IMR)*acosp(j1)
!%%%      do i=1,imr
!%%%      q(i,1) = 0.
!%%%      q(i,j1) = q(i,j1) + dq
!%%%      if(q(i,j1).lt.0.) ipy = 1
!%%%      enddo
!%%%      endif
!%%%
!%%%      if(q(1,JNP).lt.0.) then
!%%%      dq = q(1,JNP)*cap1/float(IMR)*acosp(j2)
!%%%      do i=1,imr
!%%%      q(i,JNP) = 0.
!%%%      q(i,j2) = q(i,j2) + dq
!%%%      if(q(i,j2).lt.0.) ipy = 1
!%%%      enddo
!%%%      endif
!%%%
      return
      end subroutine filns

!------------------------------------------------------------------------------

      subroutine fxppm(IMR,UT,P,DC,fx1,fx2,IORD,igzd)
      INTEGER IMR,IORD,igzd,lmt,i,j,k
C****6***0*********0*********0*********0*********0*********0**********72
      TYPE (XPLEX),PARAMETER::R3 =xplex(1./3.,0d0), R23=xplex(2./3.,0d0)
      TYPE (XPLEX) UT(1-igzd:imr+igzd+1),fx1(IMR+1),P(1-igzd:imr+igzd),
     &     DC(1-igzd:imr+igzd)
      TYPE (XPLEX) AR(0:IMR+1),AL(0:IMR+1),A6(0:IMR+1),fx2(IMR+1)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C     
      LMT = IORD - 3

      DO i=0,IMR+1
         AL(i) = 0.5*(p(i-1)+p(i)) + (DC(i-1) - DC(i))*R3
      ENDDO

      DO I=0,IMR
         AR(i) = AL(i+1)
      ENDDO

      AR(IMR+1)=0.5*(p(IMR+1)+p(IMR+2))+(DC(IMR+1)-DC(IMR+2))*R3

      DO I=0,IMR+1
         A6(I) = 3.*(P(I)+P(I)  - (AL(I)+AR(I)))
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Now call LMTPM_X instead of LMTPPM
!%%%
      IF(LMT.LE.2) CALL LMTPPM_X(DC(:),A6(:),AR(:),
     &                     AL(:),P(:),IMR,IGZD,LMT)

C Abs(UT(i)) < 1
      DO I=1,IMR+1
         IF(UT(I).GT.0.) THEN
            FX1(I) = P(I-1)
            FX2(I) = AR(I-1) + 0.5*UT(I)*(AL(I-1) - AR(I-1) +
     &                         A6(I-1)*(1.-R23*UT(I)) )
         ELSE
            FX1(I) = P(I)
            FX2(I) = AL(I) - 0.5*UT(I)*(AR(I) - AL(I) +
     &                       A6(I)*(1.+R23*UT(I)))
         ENDIF
      ENDDO
C
      DO I=1,IMR+1
         FX2(I) = FX2(I) - FX1(I)
      ENDDO

      ! Return to TPCORE
      END SUBROUTINE FXPPM

!------------------------------------------------------------------------------

      SUBROUTINE fyppm(C,P,DC,fy1,fy2,IMR,JNP,A6,AR,AL,JORD,igzd)
      INTEGER IMR,JNP,JORD,igzd,IMH,JMR,LMT,i,j,k
C****6***0*********0*********0*********0*********0*********0**********72
      TYPE (XPLEX),PARAMETER::R3=xplex(1./3.,0d0), R23=xplex(2./3.,0d0 )
      TYPE (XPLEX) C(IMR,JNP+1),fy1(IMR,JNP+1),
     &     DC(IMR,-1:JNP+2) ,fy2(IMR,JNP+1), 
     &     P(IMR,-2:JNP+3)
      TYPE (XPLEX) AR(IMR,0:JNP+1), AL(IMR,0:JNP+1),A6(IMR,0:JNP+1)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================

      IMH = IMR / 2
      JMR = JNP - 1
      LMT = JORD - 3

      DO i=1,IMR*(JNP+1)
         AL(i,1) = 0.5*(p(i,0)+p(i,1)) + (DC(i,0) - DC(i,1))*R3
         AR(i,0) = AL(i,1)
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Add the following DO-loop
!%%%
      do i=1,IMR
         AL(i,0)=0.5*(p(i,-1)+p(i,0))+
     &            (DC(i,-1)-DC(i,0))*R3
         AR(i,JNP+1)=0.5*(p(i,JNP+1)+p(i,JNP+2))+
     &             (DC(i,JNP+1)-DC(i,JNP+2))*R3
      enddo
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%C Poles:
!%%%
!%%%      DO i=1,IMH
!%%%      AL(i,1) = AL(i+IMH,2)
!%%%      AL(i+IMH,1) = AL(i,2)
!%%%
!%%%      AR(i,JNP) = AR(i+IMH,JMR)
!%%%      AR(i+IMH,JNP) = AR(i,JMR)
!%%%      enddo
!%%%
!%%%
      DO I=1,IMR*(JNP+2)
         A6(I,0) = 3.*(P(I,0)+P(I,0)  - (AL(I,0)+AR(I,0)))
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Now call new routine LMTPPM_Y
!%%%
      IF(LMT.le.2) call LMTPPM_Y( DC(:,:), A6(:,:), AR(:,:),
     &                            AL(:,:), P(:,:), IMR, JNP, LMT )

      DO I = 1, IMR*(JNP+1)
         IF(C(I,1).GT.0.) THEN
            FY1(I,1) = P(I,0)
            FY2(I,1) = AR(I,0) + 0.5*C(I,1)*(AL(I,0) - AR(I,0) +
     &           A6(I,0)*(1.-R23*C(I,1)) )
         ELSE
            FY1(I,1) = P(I,1)
            FY2(I,1) = AL(I,1) - 0.5*C(I,1)*(AR(I,1) - AL(I,1) +
     &           A6(I,1)*(1.+R23*C(I,1)))
         ENDIF
      ENDDO

      DO I = 1, IMR*(JNP+1)
         FY2(I,1) = FY2(I,1) - FY1(I,1)
      ENDDO

      ! Return to TPCORE
      END SUBROUTINE FYPPM

!------------------------------------------------------------------------------

      SUBROUTINE FZPPM(P,fz,IMR,JNP,NL,DQ,WZ,delp,KORD,fz1_tp)
C****6***0*********0*********0*********0*********0*********0**********72

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h" 
      INTEGER IMR,JNP,NL,KORD,km,km1,LMT
      TYPE (XPLEX),PARAMETER:: R23=xplex(2./3.,0d0), R3=xplex(1./3.,0d0)
      TYPE (XPLEX) WZ(IMR,JNP,NL),
     &     P(IMR,JNP,NL),
     &     DQ(IMR,JNP,NL),
     &     fz(IMR,JNP,NL+1),delp(IMR,JNP,NL)
C local 2d arrays
      TYPE (XPLEX) AR(IMR,NL),AL(IMR,NL),A6(IMR,NL),delq(IMR,NL),
     &      DC(IMR,NL)

! bey, 6/20/00. for mass-flux diagnostic
      TYPE (XPLEX) fz1_tp(IMR,JNP,NL)

      TYPE (XPLEX) lac,qmp,Tmax,Tmin,Bmax,Bmin,C1,C2,Tmp,Qmax,Qmin,A1,
     & A2,D1,D2,QM,DP,C3,CMAX,CMIN,CP,CM
      INTEGER I,j,k
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Comment this out
!%%%     TYPE (XPLEX) x, y, z
!%%%     TYPE (XPLEX) median
!%%%     median(x,y,z) = min(max(x,y), max(y,z), max(z,x))
!%%%
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================ 
      km = NL
      km1 = NL-1
      LMT = max(KORD - 3, 0)
C find global min/max     
      ! VMAX1D causes bus errors on SGI.  Replace with F90 intrinsic
      ! functions "MAXVAL" and "MINVAL".  These functions produced
      ! identical results as vmax1d in testing.  "MAXVAL" and "MINVAL"
      ! should also execute more efficiently as well. (bmy, 4/24/00)
      Tmax = MAXVAL( P(:,:,1)  )
      Tmin = MINVAL( P(:,:,1)  )
      Bmax = MAXVAL( P(:,:,NL) )
      Bmin = MINVAL( P(:,:,NL) )
#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope
CMIC$* shared(LMT,Tmax,Tmin,Bmax,Bmin,JNP,IMR)
CMIC$* shared(fz,DQ,WZ,fz1_tp)
CMIC$* private(i,j,k,c1,c2,tmp,qmax,qmin,A1,A2,d1,d2,qm,dp,c3)
CMIC$* private(cmax,cmin,DC,delq,AR,AL,A6,CM,CP, qmp, lac)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K, C1, C2, TMP, QMAX, QMIN, A1, A2,
!$OMP+                     D1, D2, QM, DP, C3, CMAX, CMIN, DC, DELQ,
!$OMP+                     AR, AL, A6, CM, CP, QMP, LAC )
#endif
#endif

      do 4000 j=1,JNP

      DO K = 2, KM
      DO I = 1, IMR
         A6(I,K) = DELP(I,J,K-1) + DELP(I,J,K)
      ENDDO
      ENDDO

      DO K = 1, KM1
      DO I = 1, IMR
         DELQ(I,K) = P(I,J,K+1) - P(I,J,K)
      ENDDO
      ENDDO

      DO K = 2, KM1
      DO I = 1, IMR
         C1 = (DELP(I,J,K-1)+0.5*DELP(I,J,K))/A6(I,K+1)
         C2 = (DELP(I,J,K+1)+0.5*DELP(I,J,K))/A6(I,K)
         TMP = DELP(I,J,K)*(C1*DELQ(I,K) + C2*DELQ(I,K-1))
     &        / (A6(I,K)+DELP(I,J,K+1))
         QMAX = MAX(P(I,J,K-1),P(I,J,K),P(I,J,K+1)) - P(I,J,K)
         QMIN = P(I,J,K) - MIN(P(I,J,K-1),P(I,J,K),P(I,J,K+1))
         DC(I,K) = SIGN(MIN(ABS(TMP),QMAX,QMIN), TMP)
      ENDDO
      ENDDO

C****6***0*********0*********0*********0*********0*********0**********72
C Compute the first guess at cell interface
C First guesses are required to be continuous.
C****6***0*********0*********0*********0*********0*********0**********72
 
C Interior.
 
      DO K = 3, KM1
      DO I = 1, IMR
         C1 = DELQ(I,K-1)*DELP(I,J,K-1) / A6(I,K)
         A1 = A6(I,K-1) / (A6(I,K) + DELP(I,J,K-1))
         A2 = A6(I,K+1) / (A6(I,K) + DELP(I,J,K))
         AL(I,K) = P(I,J,K-1) + C1 + 2./(A6(I,K-1)+A6(I,K+1)) *
     &        ( DELP(I,J,K  )*(C1*(A1 - A2)+A2*DC(I,K-1)) -
     &        DELP(I,J,K-1)*A1*DC(I,K  ) )
      ENDDO
      ENDDO

C Area preserving cubic with 2nd deriv. = 0 at the boundaries
C Top
      DO 10 I=1,IMR
         D1 = DELP(I,J,1)
         D2 = DELP(I,J,2)
         QM = (D2*P(I,J,1)+D1*P(I,J,2)) / (D1+D2)
         DP = 2.*(P(I,J,2)-P(I,J,1)) / (D1+D2)
         C1 = 4.*(AL(I,3)-QM-D2*DP) / ( D2*(2.*D2*D2+D1*(D2+3.*D1)) )
         C3 = DP - 0.5*C1*(D2*(5.*D1+D2)-3.*D1**2)
         AL(I,2) = QM - 0.25*C1*D1*D2*(D2+3.*D1)
         AL(I,1) = D1*(2.*C1*D1**2-C3) + AL(I,2)
         DC(I,1) =  P(I,J,1) - AL(I,1)
C No over- and undershoot condition
         AL(I,1) = MAX(TMIN,AL(I,1))
         AL(I,1) = MIN(TMAX,AL(I,1))
         CMAX = MAX(P(I,J,1), P(I,J,2))
         CMIN = MIN(P(I,J,1), P(I,J,2))
         AL(I,2) = MAX(CMIN,AL(I,2))
         AL(I,2) = MIN(CMAX,AL(I,2))
10    CONTINUE
 
C Bottom
      DO 15 I=1,IMR
         D1 = DELP(I,J,KM )
         D2 = DELP(I,J,KM1)
         QM = (D2*P(I,J,KM)+D1*P(I,J,KM1)) / (D1+D2)
         DP = 2.*(P(I,J,KM1)-P(I,J,KM)) / (D1+D2)
         C1 = 4.*(AL(I,KM1)-QM-D2*DP) / (D2*(2.*D2*D2+D1*(D2+3.*D1)))
         C3 = DP - 0.5*C1*(D2*(5.*D1+D2)-3.*D1**2)
         AL(I,KM) = QM - 0.25*C1*D1*D2*(D2+3.*D1)
         AR(I,KM) = D1*(2.*C1*D1**2-C3) + AL(I,KM)
         DC(I,KM) = AR(I,KM) -  P(I,J,KM)
C No over- and undershoot condition
         CMAX = MAX(P(I,J,KM), P(I,J,KM1))
         CMIN = MIN(P(I,J,KM), P(I,J,KM1))
         AL(I,KM) = MAX(CMIN,AL(I,KM))
         AL(I,KM) = MIN(CMAX,AL(I,KM))
         AR(I,KM) = MAX(BMIN,AR(I,KM))
         AR(I,KM) = MIN(BMAX,AR(I,KM))
 15   CONTINUE

      DO K=1,KM1
      DO I=1,IMR
         AR(I,K) = AL(I,K+1)
      ENDDO
      ENDDO
 
C f(s) = AL + s*[(AR-AL) + A6*(1-s)]         ( 0 <= s  <= 1 )
C Top 2 layers
      DO K=1,2
         DO I=1,IMR
            A6(I,K) = 3.*(P(I,J,K)+P(I,J,K) - (AL(I,K)+AR(I,K)))
         ENDDO
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Now call new routine LMTPPM_Z 
!%%%
         CALL LMTPPM_Z(DC(1,K),A6(1,K),AR(1,K),AL(1,K),P(1,J,K),
     &                 IMR,0)
      ENDDO

C Interior.
      IF(LMT.LE.2) THEN
         DO K=3,NL-2
            DO I=1,IMR
               A6(I,K) = 3.*(P(I,J,K)+P(I,J,K) - (AL(I,K)+AR(I,K)))
            ENDDO
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Now call new routine LMTPPM_Z 
!%%%
            CALL LMTPPM_Z(DC(1,K),A6(1,K),AR(1,K),AL(1,K),P(1,J,K),
     &           IMR,LMT)
         ENDDO

      ELSEIF(LMT .EQ. 4) THEN

c****6***0*********0*********0*********0*********0*********0**********72
C Huynh's 2nd constraint
c****6***0*********0*********0*********0*********0*********0**********72

      DO K=2, NL-1
         DO I=1,IMR
            DC(I,K) = DELQ(I,K) - DELQ(I,K-1)
         ENDDO
      ENDDO

      DO  K=3, NL-2
         DO  I=1, IMR
C Right edges
         QMP   = P(I,J,K)                 + 2.0*DELQ(I,K-1)
         LAC   = P(I,J,K) + 1.5*DC(I,K-1) + 0.5*DELQ(I,K-1)
         QMIN  = MIN(P(I,J,K), QMP, LAC)
         QMAX  = MAX(P(I,J,K), QMP, LAC)
C        AR(I,K) = MEDIAN(AR(I,K), QMIN, QMAX)
         AR(I,K) = MIN(MAX(AR(I,K), QMIN), QMAX)
C Left  edges
         QMP   = P(I,J,K)                 - 2.0*DELQ(I,K)
         LAC   = P(I,J,K) + 1.5*DC(I,K+1) - 0.5*DELQ(I,K)
         QMIN  = MIN(P(I,J,K), QMP, LAC)
         QMAX  = MAX(P(I,J,K), QMP, LAC)
c        AL(i,k) = median(AL(i,k), qmin, qmax)
         AL(I,K) = MIN(MAX(AL(I,K), QMIN), QMAX)
C Recompute A6
         A6(I,K) = 3.*(2.*P(I,J,K) - (AR(I,K)+AL(I,K)))
         ENDDO
      ENDDO
      ENDIF

C Bottom 2 layers
      DO K=NL-1,NL
         DO I=1,IMR
         A6(I,K) = 3.*(P(I,J,K)+P(I,J,K) - (AL(I,K)+AR(I,K)))
         ENDDO
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Now call new routine LMTPPM_Z 
!%%%
      CALL LMTPPM_Z(DC(1,K),A6(1,K),AR(1,K),AL(1,K),P(1,J,K),
     &            IMR,0)
      enddo
 
      DO K = 2, NL
      DO I = 1, IMR
         IF(WZ(I,J,K-1).GT.0.) THEN
            CM = WZ(I,J,K-1) / DELP(I,J,K-1)
            DC(I,K) = P(I,J,K-1)
            FZ(I,J,K) = AR(I,K-1)+0.5*CM*(AL(I,K-1)-AR(I,K-1)+
     &           A6(I,K-1)*(1.-R23*CM))
         ELSE
            CP = WZ(I,J,K-1) / DELP(I,J,K)
            DC(I,K) = P(I,J,K)
            FZ(I,J,K) = AL(I,K)+0.5*CP*(AL(I,K)-AR(I,K)-
     &           A6(I,K)*(1.+R23*CP))
         ENDIF
      ENDDO
      ENDDO

      DO K = 2, NL
      DO I = 1, IMR
         FZ(I,J,K) = WZ(I,J,K-1) * (FZ(I,J,K) - DC(I,K))
         DC(I,K) = WZ(I,J,K-1) * DC(I,K)
      ENDDO
      ENDDO

      DO 350 I=1,IMR
         FZ(I,J,   1) = 0.
         FZ(I,J,NL+1) = 0.
         DQ(I,J, 1) = DQ(I,J, 1) - DC(I, 2)
         DQ(I,J,NL) = DQ(I,J,NL) + DC(I,NL)
         FZ1_TP(I,J,NL) = DC(I,NL) !(yxw, 01/21/2003)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         !%%% ERROR!  FZ1_TP is only declared with NL layers.  This line
         !%%% causes an out-of-bounds error which causes the run to die
         !%%% when running on the ALTIX platform.  Comment out. (bmy, 11/9/04)
         !%%%FZ1_TP(I,J,NL+1) = 0      !(yxw, 02/09/2003)
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         FZ1_TP(I,J,1) = 0.        !(yxw, 01/21/2003)
 350  CONTINUE
 
!-----------------------------------------------------------------------------
! bey, 6/20/00. for mass-flux diagnostic, loop had to be extended 
!      do 360 k=2,km1
!      do 360 i=1,IMR
!360   DQ(i,j,k) = DQ(i,j,k) + DC(i,k) - DC(i,k+1)
!-----------------------------------------------------------------------------
      DO K=2,KM1
      DO I=1,IMR
         DQ(I,J,K) = DQ(I,J,K) + DC(I,K) - DC(I,K+1)

         ! bey, 6/20/00. for mass-flux diagnostic
         FZ1_TP(I,J,K) = DC(I,K)
      ENDDO
      ENDDO
 
4000  CONTINUE
      
      ! Return to TPCORE
      END SUBROUTINE FZPPM

!------------------------------------------------------------------------------

      SUBROUTINE HILO(Q,IM,JM,QMAX,QMIN,BT,BD)
      INTEGER IM,JM,I,J
C****6***0*********0*********0*********0*********0*********0**********72
      TYPE (XPLEX) q(IM,JM),Qmax(IM,JM),Qmin(IM,JM),bt(IM,*),bd(IM,*)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================ 
C y-sweep
      DO J = 1, JM
      DO I = 0, IM+1
         BT(I,J) = MAX(Q(I,J-1),Q(I,J),Q(I,J+1))
         BD(I,J) = MIN(Q(I,J-1),Q(I,J),Q(I,J+1))
      ENDDO
      ENDDO
C
C x-sweep
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Comment this line out
!%%%      IM1 = IM-1
!%%%
      DO J = 1, JM
      DO I = 1, IM
         QMAX(I,J) = MAX(BT(I-1,J),BT(I,J),BT(I+1,J))
         QMIN(I,J) = MIN(BD(I-1,J),BD(I,J),BD(I+1,J))
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%
!%%%
!%%% don't need the following part (yxw,08/27/01)
!%%%
!%%%      DO j=j1,j2
!%%%     i = 1
!%%%      Qmax(1,j) = max(bt(IM,j),bt(1,j),bt(2,j))
!%%%      Qmin(1,j) = min(bd(IM,j),bd(1,j),bd(2,j))
!%%%     i = IM
!%%%      Qmax(IM,j) = max(bt(IM1,j),bt(IM,j),bt(1,j))
!%%%      Qmin(IM,j) = min(bd(IM1,j),bd(IM,j),bd(1,j))
!%%%      enddo
!%%%
!%%%C N. Pole:
!%%%      Pmax = q(1,JM)
!%%%      Pmin = q(1,JM)
!%%%      do i=1,IM
!%%%      if(q(i,j2) .gt. Pmax) then
!%%%            Pmax = q(i,j2)
!%%%       elseif(q(i,j2) .lt. Pmin) then
!%%%            Pmin = q(i,j2)
!%%%      endif
!%%%      enddo
!%%%
!%%%      do i=1,IM
!%%%      Qmax(i,JM) = Pmax
!%%%      Qmin(i,JM) = Pmin
!%%%      enddo
!%%%
!%%%C S. Pole:
!%%%      Pmax = q(1,1)
!%%%      Pmin = q(1,1)
!%%%      do i=1,IM
!%%%      if(q(i,j1) .gt. Pmax) then
!%%%            Pmax = q(i,j1)
!%%%      elseif(q(i,j1) .lt. Pmin) then
!%%%            Pmin = q(i,j1)
!%%%      endif
!%%%      enddo
!%%%
!%%%      do i=1,IM
!%%%      Qmax(i,1) = Pmax
!%%%      Qmin(i,1) = Pmin
!%%%      enddo
!%%%
!%%%      if(j1 .ne. 2) then
!%%%      JM1 = JM-1
!%%%      do i=1,IM
!%%%      Qmax(i,2) = Qmax(i,1)
!%%%      Qmin(i,2) = Qmin(i,1)
!%%%
!%%%      Qmax(i,JM1) = Qmax(i,JM)
!%%%      Qmin(i,JM1) = Qmin(i,JM)
!%%%      enddo
!%%%      endif

      ! Return to TPCORE
      END SUBROUTINE HILO

!------------------------------------------------------------------------------

      SUBROUTINE HILO3D(P,IM,JM,KM,PMAX,PMIN,QMAX,QMIN,BT,BD)
C****6***0*********0*********0*********0*********0*********0**********72

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h" 
      INTEGER I,J,K,IM,JM,KM,KM1,KM2
      TYPE (XPLEX) P(IM+2,JM+2,km),Pmax(IM,JM,km),Pmin(IM,JM,km),
     &     Qmax(IM,JM,km),Qmin(IM,JM,km),bt(im,jm),bd(im,jm)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================ 

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope
CMIC$* shared(P,Qmax,Qmin,im,jm,j1,j2)
CMIC$* private(k,bt,bd)
#else
!$OMP PARALLEL DO PRIVATE( K, BT, BD )
#endif
#endif

      DO 1000 K=1,KM
         CALL HILO(P(1,1,K),IM,JM,QMAX(1,1,K),QMIN(1,1,K),BT,BD)
1000  CONTINUE
 
      KM1 = KM-1
      KM2 = KM-2

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope
CMIC$* shared(Pmax,Pmin,Qmax,Qmin,im,jm,km,km1,km2)
CMIC$* private(i,j)
#else
!$OMP PARALLEL DO PRIVATE( I, J )
#endif
#endif
 
      DO J = 1, JM
      DO I = 1, IM
C k=1 and k=km
         PMAX(I,J, 1) = MAX(QMAX(I,J,  2),QMAX(I,J, 1))
         PMIN(I,J, 1) = MIN(QMIN(I,J,  2),QMIN(I,J, 1))
         PMAX(I,J,KM) = MAX(QMAX(I,J,KM1),QMAX(I,J,KM))
         PMIN(I,J,KM) = MIN(QMIN(I,J,KM1),QMIN(I,J,KM))
C k=2 and k=km1
         PMAX(I,J,  2) = MAX(QMAX(I,J,  3),PMAX(I,J, 1))
         PMIN(I,J,  2) = MIN(QMIN(I,J,  3),PMIN(I,J, 1))
         PMAX(I,J,KM1) = MAX(QMAX(I,J,KM2),PMAX(I,J,KM))
         PMIN(I,J,KM1) = MIN(QMIN(I,J,KM2),PMIN(I,J,KM))
      ENDDO
      ENDDO

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope
CMIC$* shared(Pmax,Pmin,Qmax,Qmin,im,jm,km,km1,km2)
CMIC$* private(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE( I, J, K )
#endif
#endif

      DO K = 3, KM2
      DO J = 1, JM
      DO I = 1, IM
         PMAX(I,J,K) = MAX(QMAX(I,J,K-1),QMAX(I,J,K),QMAX(I,J,K+1))
         PMIN(I,J,K) = MIN(QMIN(I,J,K-1),QMIN(I,J,K),QMIN(I,J,K+1))
      ENDDO
      ENDDO
      ENDDO
      
      ! Return to TPCORE
      END SUBROUTINE HILO3D

!------------------------------------------------------------------------------

      SUBROUTINE LMTPPM_X(DC,A6,AR,AL,P,IM,IGZD,LMT)
C****6***0*********0*********0*********0*********0*********0**********72
C
C A6 =  CURVATURE OF THE TEST PARABOLA
C AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
C AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
C DC =  0.5 * MISMATCH
C P  =  CELL-AVERAGED VALUE
C IM =  VECTOR LENGTH
C
C OPTIONS:
C
C LMT = 0: FULL MONOTONICITY
C LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
C LMT = 2: POSITIVE-DEFINITE CONSTRAINT
C
      INTEGER LMT,IM,IGZD,I,J,K
      TYPE (XPLEX) DA1,DA2,A6DA,FMIN
      TYPE (XPLEX),PARAMETER:: R12 = xplex(1./12.,0d0 )
      TYPE (XPLEX) A6(0:IM+1),AR(0:IM+1),AL(0:IM+1),
     &     P(1-igzd:im+igzd),DC(1-igzd:im+igzd)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C
      IF(LMT.EQ.0) THEN
C FULL CONSTRAINT
      DO 100 I=0,IM+1
         IF(DC(I).EQ.0.) THEN
            AR(I) = P(I)
            AL(I) = P(I)
            A6(I) = 0.
         ELSE
            DA1  = AR(I) - AL(I)
            DA2  = DA1**2
            A6DA = A6(I)*DA1
            IF(A6DA .LT. -DA2) THEN
               A6(I) = 3.*(AL(I)-P(I))
               AR(I) = AL(I) - A6(I)
            ELSEIF(A6DA .GT. DA2) THEN
               A6(I) = 3.*(AR(I)-P(I))
               AL(I) = AR(I) - A6(I)
            ENDIF
         ENDIF
 100  CONTINUE
      ELSEIF(LMT.EQ.1) THEN
C SEMI-MONOTONIC CONSTRAINT
         DO 150 I=0,IM+1
            IF(ABS(AR(I)-AL(I)) .GE. -A6(I)) GO TO 150
            IF(P(I).LT.AR(I) .AND. P(I).LT.AL(I)) THEN
               AR(I) = P(I)
               AL(I) = P(I)
               A6(I) = 0.
            ELSEIF(AR(I) .GT. AL(I)) THEN
               A6(I) = 3.*(AL(I)-P(I))
               AR(I) = AL(I) - A6(I)
            ELSE
               A6(I) = 3.*(AR(I)-P(I))
            AL(I) = AR(I) - A6(I)
         ENDIF
 150  CONTINUE
      ELSEIF(LMT.EQ.2) THEN
         DO 250 I=0,IM+1
            IF(ABS(AR(I)-AL(I)) .GE. -A6(I)) GO TO 250
            FMIN = P(I) + 0.25*(AR(I)-AL(I))**2/A6(I) + A6(I)*R12
            IF(FMIN.GE.0.) GO TO 250
            IF(P(I).LT.AR(I) .AND. P(I).LT.AL(I)) THEN
               AR(I) = P(I)
               AL(I) = P(I)
               A6(I) = 0.
            ELSEIF(AR(I) .GT. AL(I)) THEN
               A6(I) = 3.*(AL(I)-P(I))
               AR(I) = AL(I) - A6(I)
            ELSE
               A6(I) = 3.*(AR(I)-P(I))
               AL(I) = AR(I) - A6(I)
            ENDIF
 250     CONTINUE
      ENDIF

      ! Return to TPCORE
      END SUBROUTINE LMTPPM_X

!------------------------------------------------------------------------------

      SUBROUTINE LMTPPM_Y(DC,A6,AR,AL,P,IM,JNP,LMT)
C****6***0*********0*********0*********0*********0*********0**********72
C
C A6 =  CURVATURE OF THE TEST PARABOLA
C AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
C AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
C DC =  0.5 * MISMATCH
C P  =  CELL-AVERAGED VALUE
C IM =  VECTOR LENGTH
C
C OPTIONS:
C
C LMT = 0: FULL MONOTONICITY
C LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
C LMT = 2: POSITIVE-DEFINITE CONSTRAINT
C
      INTEGER LMT,IM,IGZD,I,J,K,JNP
      TYPE (XPLEX) DA1,DA2,A6DA,FMIN
      TYPE (XPLEX),PARAMETER:: R12 = xplex(1./12.,0d0 )
      TYPE (XPLEX) A6(IM,0:JNP+1),AR(IM,0:JNP+1),AL(IM,0:JNP+1),
     &     P(IM,-2:JNP+3),DC(IM,-1:JNP+2)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C
      IF(LMT.EQ.0) THEN
C FULL CONSTRAINT
         DO 100 I=1,IM*(JNP+2)
            IF(DC(I,0).EQ.0.) THEN
               AR(I,0) = P(I,0)
               AL(I,0) = P(I,0)
               A6(I,0) = 0.
            ELSE
               DA1  = AR(I,0) - AL(I,0)
               DA2  = DA1**2
               A6DA = A6(I,0)*DA1
               IF(A6DA .LT. -DA2) THEN
                  A6(I,0) = 3.*(AL(I,0)-P(I,0))
                  AR(I,0) = AL(I,0) - A6(I,0)
               ELSEIF(A6DA .GT. DA2) THEN
                  A6(I,0) = 3.*(AR(I,0)-P(I,0))
                  AL(I,0) = AR(I,0) - A6(I,0)
               ENDIF
            ENDIF
 100     CONTINUE
      ELSEIF(LMT.EQ.1) THEN
C SEMI-MONOTONIC CONSTRAINT
      DO 150 I=1,IM*(JNP+2)
         IF(ABS(AR(I,0)-AL(I,0)) .GE. -A6(I,0)) GO TO 150
         IF(P(I,0).LT.AR(I,0) .AND. P(I,0).LT.AL(I,0)) THEN
            AR(I,0) = P(I,0)
            AL(I,0) = P(I,0)
            A6(I,0) = 0.
         ELSEIF(AR(I,0) .GT. AL(I,0)) THEN
            A6(I,0) = 3.*(AL(I,0)-P(I,0))
            AR(I,0) = AL(I,0) - A6(I,0)
         ELSE
            A6(I,0) = 3.*(AR(I,0)-P(I,0))
            AL(I,0) = AR(I,0) - A6(I,0)
         ENDIF
 150  CONTINUE
      ELSEIF(LMT.EQ.2) THEN
         DO 250 I=1,IM*(JNP+2)
            IF(ABS(AR(I,0)-AL(I,0)) .GE. -A6(I,0)) GO TO 250
            FMIN = P(I,0) + 0.25*(AR(I,0)-AL(I,0))**2/A6(I,0) + 
     &           A6(I,0)*R12
            IF(FMIN.GE.0.) GO TO 250
            IF(P(I,0).LT.AR(I,0) .AND. P(I,0).LT.AL(I,0)) THEN
               AR(I,0) = P(I,0)
               AL(I,0) = P(I,0)
               A6(I,0) = 0.
            ELSEIF(AR(I,0) .GT. AL(I,0)) THEN
               A6(I,0) = 3.*(AL(I,0)-P(I,0))
               AR(I,0) = AL(I,0) - A6(I,0)
            ELSE
               A6(I,0) = 3.*(AR(I,0)-P(I,0))
               AL(I,0) = AR(I,0) - A6(I,0)
            ENDIF
 250     CONTINUE
      ENDIF

      ! Return to TPCORE
      END SUBROUTINE LMTPPM_Y

!------------------------------------------------------------------------------

      SUBROUTINE LMTPPM_Z(DC,A6,AR,AL,P,IM,LMT)
C****6***0*********0*********0*********0*********0*********0**********72
C
C A6 =  CURVATURE OF THE TEST PARABOLA
C AR =  RIGHT EDGE VALUE OF THE TEST PARABOLA
C AL =  LEFT  EDGE VALUE OF THE TEST PARABOLA
C DC =  0.5 * MISMATCH
C P  =  CELL-AVERAGED VALUE
C IM =  VECTOR LENGTH
C
C OPTIONS:
C
C LMT = 0: FULL MONOTONICITY
C LMT = 1: SEMI-MONOTONIC CONSTRAINT (NO UNDERSHOOTS)
C LMT = 2: POSITIVE-DEFINITE CONSTRAINT
C
      INTEGER LMT,IM,IGZD,I,J,K
      TYPE (XPLEX) DA1,DA2,A6DA,FMIN
      TYPE (XPLEX),PARAMETER::  R12 = xplex(1./12.,0d0 )
      TYPE (XPLEX) A6(IM),AR(IM),AL(IM),P(IM),DC(IM)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C
      IF(LMT.EQ.0) THEN
C FULL CONSTRAINT
         DO 100 I=1,IM
            IF(DC(I).EQ.0.) THEN
               AR(I) = P(I)
               AL(I) = P(I)
               A6(I) = 0.
            ELSE
               DA1  = AR(I) - AL(I)
               DA2  = DA1**2
               A6DA = A6(I)*DA1
               IF(A6DA .LT. -DA2) THEN
                  A6(I) = 3.*(AL(I)-P(I))
                  AR(I) = AL(I) - A6(I)
               ELSEIF(A6DA .GT. DA2) THEN
                  A6(I) = 3.*(AR(I)-P(I))
                  AL(I) = AR(I) - A6(I)
               ENDIF
            ENDIF
 100     CONTINUE
      ELSEIF(LMT.EQ.1) THEN
C SEMI-MONOTONIC CONSTRAINT
         DO 150 I=1,IM
            IF(ABS(AR(I)-AL(I)) .GE. -A6(I)) GO TO 150
            IF(P(I).LT.AR(I) .AND. P(I).LT.AL(I)) THEN
               AR(I) = P(I)
               AL(I) = P(I)
               A6(I) = 0.
            ELSEIF(AR(I) .GT. AL(I)) THEN
               A6(I) = 3.*(AL(I)-P(I))
               AR(I) = AL(I) - A6(I)
            ELSE
               A6(I) = 3.*(AR(I)-P(I))
               AL(I) = AR(I) - A6(I)
            ENDIF
 150     CONTINUE
      ELSEIF(LMT.EQ.2) THEN
         DO 250 I=1,IM
            IF(ABS(AR(I)-AL(I)) .GE. -A6(I)) GO TO 250
            FMIN = P(I) + 0.25*(AR(I)-AL(I))**2/A6(I) + A6(I)*R12
            IF(FMIN.GE.0.) GO TO 250
            IF(P(I).LT.AR(I) .AND. P(I).LT.AL(I)) THEN
               AR(I) = P(I)
               AL(I) = P(I)
               A6(I) = 0.
            ELSEIF(AR(I) .GT. AL(I)) THEN
               A6(I) = 3.*(AL(I)-P(I))
               AR(I) = AL(I) - A6(I)
            ELSE
               A6(I) = 3.*(AR(I)-P(I))
               AL(I) = AR(I) - A6(I)
            ENDIF
 250     CONTINUE
      ENDIF

      ! Return to TPCORE
      END SUBROUTINE LMTPPM_Z

!------------------------------------------------------------------------------

      SUBROUTINE QCKXYZ(Q,QTMP,IMR,JNP,NLAY,COSP,ACOSP,IC,NSTEP,DP,
     &           IGZD, JMAX,FX, FY, FZ)
C****6***0*********0*********0*********0*********0*********0**********72

      ! Added to pass C-preprocessor switches (bmy, 3/9/01)
#     include "define.h" 
      INTEGER IMR,JNP,NLAY,IC,NSTEP,IGZD,JMAX,NLM1,I,J,K,L 
      TYPE (XPLEX) ,PARAMETER::  tiny = xplex(1.D-30,0d0 )
      INTEGER , PARAMETER ::  kmax = 200 
      TYPE (XPLEX) Q(IMR,JNP,NLAY),qtmp(JNP, IMR),cosp(-10:Jmax),
     &       acosp(-10:Jmax),DP,QLY,QUP,DUP
      TYPE (XPLEX) fx(IMR+1, JNP, NLAY), fy(IMR, JNP+1, NLAY), 
     &       fz(IMR, JNP, NLAY+1)          !(yxw, 02/09/2003)
      integer IP(kmax),ipz,lpz,iu,ru
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================ 
      NLM1 = NLAY-1

C Do horizontal filling.

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope
CMIC$* private(i,j,L,qtmp)
#else
!$OMP PARALLEL DO PRIVATE( I, J, L, QTMP )
#endif
#endif

      DO 1000 L=1,NLAY
         CALL FILNS(Q(1,1,L),IMR,JNP,COSP,ACOSP,IP(L),TINY,DP, 
     &              FY(1,1,L), JMAX)
!%%% 
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Call FILEW 
         IF(IP(L).NE.0) 
     &        CALL FILEW(Q(1,1,L),QTMP,IMR,JNP,IP(L),TINY, FX(1,1,L))

1000  CONTINUE

      IPZ = 0
      DO L=1,NLAY
         IF(IP(L) .NE. 0) THEN
            IPZ = L
            GO TO 111
         ENDIF
      ENDDO
      RETURN

 111  CONTINUE

      IF(IPZ .EQ. 0) RETURN

      IF(IPZ .EQ. 1) THEN
         LPZ = 2
      ELSE
         LPZ = IPZ
      ENDIF

C Do vertical filling.

#if   defined( multitask  )
#if   defined( CRAY       )
CMIC$ do all autoscope
CMIC$* private(i,j,L,qup,qly,dup)
#else
!$OMP PARALLEL DO PRIVATE( I, J, L, QUP, QLY, DUP )
#endif
#endif

      DO 2000 J=1,JNP

         IF(IPZ .EQ. 1) THEN
C Top layer
            DO I=1,IMR
               IF(Q(I,J,1).LT.0.) THEN
                  Q(I,J,2) = Q(I,J,2) + Q(I,J,1)
                  FZ(I,J,2) = FZ(I,J,2) + Q(I,J,1) !(yxw,02/09/2003)
                  Q(I,J,1) = 0.
               ENDIF
            ENDDO
         ENDIF
 
         DO 225 L = LPZ,NLM1
            DO I=1,IMR
               IF( Q(I,J,L).LT.0.) THEN
C From above
                  QUP =  Q(I,J,L-1)
                  QLY = -Q(I,J,L)
                  DUP  = MIN(QLY,QUP)
                  Q(I,J,L-1) = QUP - DUP
                  Q(I,J,L  ) = DUP-QLY
                  FZ(I,J,L)= FZ(I,J,L-1) + DUP
C Below
                  Q(I,J,L+1) = Q(I,J,L+1) + Q(I,J,L) 
                  FZ(I,J,L+1)= FZ(I,J,L+1)+ Q(I,J,L) !(yxw, 02/09/2003)
                  Q(I,J,L)   = 0.
               ENDIF
            ENDDO
 225     CONTINUE
 
C BOTTOM LAYER
         L = NLAY
         DO I=1,IMR
            IF( Q(I,J,L).LT.0.) THEN
C From above 
               QUP = Q(I,J,NLM1)
               QLY = -Q(I,J,L)
               DUP = MIN(QLY,QUP)
               Q(I,J,NLM1) = QUP - DUP
               FZ(I,J,L) = FZ(I,J,L) + DUP !(yxw,02/09/2003)
C From "below" the surface.
               Q(I,J,L) = 0.
               FZ(I,J,L+1) = FZ(I,J,L+1) + DUP - QLY
            ENDIF
         ENDDO
 2000 CONTINUE

      ! Return to TPCORE
      END SUBROUTINE qckxyz

!------------------------------------------------------------------------------

      SUBROUTINE XADV( IM_W, JM_W, P,   UA,   JS,   JN, ADX_W, 
     &                 IGZD, IMR,  JNP, I0_W, J0_W, I0, J0 ) !( yxw,08/23/01) 
      INTEGER im_w,jm_w,js,jn,igzd,imr,jnp,i0_w,j0_w,i0,j0,i,j,k,iu,ru,
     & iiu,ip
C****6***0*********0*********0*********0*********0*********0**********72
      TYPE (XPLEX) p(IMR,JNP),adx_w(1-igzd:im_w+igzd,1-igzd:jm_w+igzd),
     & qtmp(1-igzd:im_w+igzd),UA(1-igzd:im_w+igzd,1-igzd:jm_w+igzd)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C     
      do 1309 j=1-igzd,jm_w+igzd
         if((J+j0_w+j0) .GT.JS  .and.(J+j0_w+j0).LT.JN) GO TO 1309

         do i=1-igzd,IM_w+igzd
            qtmp(i) = p(i+i0_w,j+j0_w)
         enddo

!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Comment this out      
!%%%      do i=-IML,0
!%%%      qtmp(i)       = p(IMR+i,j)
!%%%      qtmp(IMR+1-i) = p(1-i,j)
!%%%     enddo
!%%%
!%%%
         DO I=1-IGZD,IM_W+IGZD
            IU = UA(I,J)
            RU = UA(I,J) - IU
            IIU = I+I0_W-IU
            IF(UA(I,J).GE.0.) THEN
               IF (IIU .LT. 2) THEN 
                  IIU=2
               ELSEIF (IIU .GT. IMR) THEN
                  IIU=IMR
               ENDIF
               ADX_W(I,J) = P(IIU,J+J0_W)+RU*(P(IIU-1,J+J0_W)
     &              -P(IIU,J+J0_W))
            ELSE
               IF (IIU .LT. 1) THEN 
                  IIU=1
               ELSEIF (IIU .GT. IMR-1) THEN
                  IIU=IMR-1
               ENDIF
               ADX_W(I,J) = P(IIU,J+J0_W)+RU*(P(IIU,J+J0_W)-
     &              P(IIU+1,J+J0_W))
            ENDIF
         ENDDO
 
         DO I=1-IGZD,IM_W+IGZD
            ADX_W(I,J) = ADX_W(I,J) - P(I+I0_W,J+J0_W)
         ENDDO
 1309 CONTINUE
 
C Eulerian upwind
 
      DO J=JS+1-J0_W-J0,JN-1-J0_W-J0
C
      
         DO I=1-IGZD,IM_W+IGZD
            QTMP(I) = P(I+I0_W,J+J0_W)
         ENDDO
         
         DO I=1-IGZD,IM_W+IGZD
            IP = XPLX(I+I0_W) - UA(I,J)
            IF (IP .GT. IMR-1) THEN
               IP=IMR-1
            ELSE IF(IP .LT. 1) THEN
               IP=1
            ENDIF
            ADX_W(I,J) = UA(I,J)*(P(IP,J+J0_W)-P(IP+1,J+J0_W))
         ENDDO
         
      ENDDO
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/030
!%%% Polar cap stuff, useless for nested simulation
!%%%C don't need the following polar part (yxw, 08/23/01)
!%%%      if(j1.ne.2) then
!%%%      do i=1,IMR
!%%%      adx(i,  2) = 0.
!%%%      adx(i,JMR) = 0.
!%%%      enddo
!%%%      endif
!%%%
!%%%C set cross term due to x-adv at the poles to zero.
!%%%      do i=1,IMR
!%%%      adx(i,  1) = 0.
!%%%      adx(i,JNP) = 0.
!%%%      enddo
!%%%

      ! Return to TPCORE
      END SUBROUTINE XADV

!------------------------------------------------------------------------------

      SUBROUTINE XMIST(IMR,P,DC,igzd)
      INTEGER IMR,IGZD,i
      TYPE (XPLEX) TMP,PMAX,PMIN
C****6***0*********0*********0*********0*********0*********0**********72
      TYPE (XPLEX) P(1-igzd:IMR+igzd),DC(1-igzd:IMR+igzd)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C
C 2nd order version.
C
      DC(:)=0d0      ! initialization (yxw,09/28/01)

      DO I=2-IGZD,IMR+IGZD-1
         TMP = 0.25*(P(I+1) - P(I-1))
         PMAX = MAX(P(I-1), P(I), P(I+1)) - P(I)
         PMIN = P(I) - MIN(P(I-1), P(I), P(I+1))
         DC(I) = SIGN(MIN(ABS(TMP),PMAX,PMIN), TMP)
      ENDDO

      ! Return to TPCORE
      END SUBROUTINE XMIST

!------------------------------------------------------------------------------

      SUBROUTINE XTP(IM_W,JM_W,IGZD,JN,JS,PU,DQ,Q,C,FX2,XMASS,IORD,
     &                 FX1_TP,I0_W,J0_W,I0,J0)
C****6***0*********0*********0*********0*********0*********0**********72
      INTEGER im_w,jm_w,igzd,jn,js,i0_w,j0_w,imp,i,j,k,iord,i0,j0
      TYPE (XPLEX) C(1-IGZD:IM_W+1+IGZD,1-IGZD:JM_W+IGZD),
     &     FX1(IM_W+1), DC(1-IGZD:IM_W+IGZD),
     &     DQ(IM_W,JM_W),QTMP(1-IGZD:IM_W+IGZD),XMASS(IM_W+1,JM_W)
      TYPE (XPLEX) PU(IM_W+1,JM_W),Q(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD)
      TYPE (XPLEX) FX2(IM_W+1,JM_W)
      INTEGER ISAVE(IM_W),itmp,iu,ist
! bey, 6/20/00. for mass-flux diagnostic
      TYPE (XPLEX) FX1_TP(IM_W,JM_W),rut

C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
      IMP = im_w + 1
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%C van Leer at high latitudes
!%%%      jvan = max(1,jm/20)
!%%%      j1vl = j1+jvan
!%%%      j2vl = j2-jvan
!%%%
      do 1310 j=1,jm_w
C
         do i=1-igzd,im_w+igzd
            qtmp(i) = q(i,j)
         enddo
C
         if(j .ge.(JN-j0_w-J0) .or. j .le. (JS-j0_w-J0)) goto 2222
C****6***0*********0*********0*********0*********0*********0**********72
C *** Eulerian ***
C****6***0*********0*********0*********0*********0*********0**********72

!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%C don't need the following. (yxw,08/23/01)
!%%%      qtmp(0)     = q(im,J)
!%%%      qtmp(-1)    = q(im-1,J)
!%%%      qtmp(IMP)   = q(1,J)
!%%%      qtmp(IMP+1) = q(2,J)
!%%%
         IF(IORD.EQ.1) THEN
            DO I=1,IM_W+1
               IU = XPLX(I) - C(I,J)
               FX1(I) = QTMP(IU)
            ENDDO
 
C Zero high order contribution
            DO I=1,IM_W+1
               FX2(I,J) = 0.
            ENDDO
         ELSE
            CALL XMIST(IM_W,QTMP,DC,IGZD) 
            DC(1-IGZD)=DC(2-IGZD)
            DC(IM_W+IGZD)=DC(IM_W+IGZD-1)
C     
            IF(IORD.EQ.2) THEN
               DO I=1,IM_W+1
                  IU = XPLX(I) - C(I,J)
                  FX1(I  ) = QTMP(IU)
                  FX2(I,J) = DC(IU)*(SIGN(1.d0,C(I,J))-C(I,J))
               ENDDO
            ELSE
        CALL FXPPM(IM_W,C(:,J),QTMP(:),DC(:),FX1(:),FX2(:,J),IORD,IGZD)
      ENDIF
C
      ENDIF
C
      DO I=1,IM_W+1
         FX1(I  ) = FX1(I  )*XMASS(I,J)
         FX2(I,J) = FX2(I,J)*XMASS(I,J)
      ENDDO
C
      GOTO 1309
C
C****6***0*********0*********0*********0*********0*********0**********72
C *** Conservative (flux-form) Semi-Lagrangian transport ***
C****6***0*********0*********0*********0*********0*********0**********72
 
 2222 CONTINUE
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Don't need ghost zones
!%%%C ghost zone for the western edge:
!%%%      iuw =  -c(1,j)
!%%%      iuw = min(0, iuw)
!%%%
!%%%      do i=iuw, 0
!%%%         qtmp(i) = q(im+i,j)
!%%%      enddo
!%%%
!%%%C ghost zone for the eastern edge:
!%%%      iue = imp - c(im,j)
!%%%      iue = max(imp, iue)
!%%%
!%%%      do i=imp, iue
!%%%        qtmp(i) = q(i-im,j)
!%%%      enddo
!%%%
      IF(IORD.EQ.1) THEN
         DO I=1,IM_W+1
            IU = C(I,J)
            IF(C(I,J) .LE. 0.) THEN
               ITMP = I+I0_W - IU
               ISAVE(I) = ITMP - 1
            ELSE
               ITMP = I+I0_W - IU - 1
               ISAVE(I) = ITMP + 1
            ENDIF
            IF (ITMP .GT. IM_W+IGZD+I0_W) THEN
               ITMP=IM_W+IGZD+I0_W
            ELSE IF (ITMP .LT. 1-IGZD+I0_W) THEN
               ITMP=1-IGZD+I0_W
            ENDIF 
            FX1(I) = (C(I,J)-IU) * QTMP(ITMP-I0_W)
         ENDDO

C Zero high order contribution
         DO I=1,IM_W+1
            FX2(I,J) = 0.
         ENDDO
         
      ELSE
         CALL XMIST(IM_W,QTMP,DC,IGZD)
         DC(1-IGZD) = DC(2-IGZD)
         DC(IM_W+IGZD)=DC(IM_W+IGZD-1)
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Don't need ghost zones
!%%%      do i=iuw, 0
!%%%         dc(i) = dc(im+i)
!%%%      enddo
!%%%
!%%%      do i=imp, iue
!%%%         dc(i) = dc(i-im)
!%%%      enddo
!%%%
         DO I=1,IM_W+1
            IU  = C(I,J)
            RUT = C(I,J) - IU
            IF(C(I,J) .LE. 0.) THEN
               ITMP = I+I0_W+I0 - IU
               ISAVE(I) = ITMP - 1
               IF (ITMP .GT. IM_W+IGZD+I0_W+I0) THEN
                  ITMP=IM_W+IGZD+I0_W+I0
               ELSE IF (ITMP .LT. 1-IGZD+I0_W+I0) THEN
                  ITMP=1-IGZD+I0_W+I0
               ENDIF
               ITMP=ITMP-I0_W-I0
               FX2(I,J) = -RUT*DC(ITMP)*(1.+RUT)
            ELSE
               ITMP = I+I0_W+I0 - IU - 1
               ISAVE(I) = ITMP + 1
               IF (ITMP .GT. IM_W+IGZD+I0_W+I0) THEN
                  ITMP=IM_W+IGZD+I0_W+I0
               ELSE IF (ITMP .LT. 1-IGZD+I0_W+I0) THEN
                  ITMP=1-IGZD+I0_W+I0
               ENDIF
               ITMP=ITMP-I0_W-I0
               FX2(I,J) = RUT*DC(ITMP)*(1.-RUT)
            ENDIF
            FX1(I) = RUT*QTMP(ITMP)
         ENDDO
         
      ENDIF
      
      DO I=1,IM_W+1
         IF (ISAVE(I) .GT. IM_W+IGZD+I0_W+I0) THEN
            ISAVE(I)=IM_W+IGZD+I0_W+I0
         ELSE IF (ISAVE(I) .LT. 1-IGZD+I0_W+I0) THEN
            ISAVE(I)=1-IGZD+I0_W+I0
         ENDIF
         ISAVE(I)=ISAVE(I)-I0-I0_W
         IF(C(I,J).GT.1.) THEN
CDIR$ NOVECTOR
            DO IST =ISAVE(I),I-1
               FX1(I) = FX1(I) + QTMP(IST)
            ENDDO
         ELSEIF(C(I,J).LT.-1.) THEN
CDIR$ NOVECTOR
            DO IST = I,ISAVE(I)
               FX1(I) = FX1(I) - QTMP(IST)
            ENDDO
         ENDIF
      ENDDO
CDIR$ VECTOR
      DO I=1,IM_W+1
         FX1(I)   = PU(I,J)*FX1(I)
         FX2(I,J) = PU(I,J)*FX2(I,J)
      ENDDO

C use extrapolation to calculate fx1 and fx2 at grid IMP (yxw, 08/24/01)

 1309 CONTINUE 

C Update using low order fluxes.
      DO I=1,IM_W

         DQ(I,J) =  DQ(I,J) + FX1(I)-FX1(I+1)

! bey, 6/20/00. for mass-flux diagnostic
         FX1_TP(I,J) = FX1(I)
      ENDDO
 
 1310 CONTINUE

      ! Return to TPCORE
      END SUBROUTINE XTP

!------------------------------------------------------------------------------

      SUBROUTINE  YMIST(IMR,JNP,P,DC,IGZD)
      INTEGER IMR,JNP,IGZD,I
      TYPE (XPLEX) TMP,PMAX,PMIN
C****6***0*********0*********0*********0*********0*********0**********72
      TYPE (XPLEX) ,PARAMETER :: R24 = xplex(1./24.,0d0 )
      TYPE (XPLEX) P(IMR,-2:JNP+3),DC(IMR,-1:JNP+2)
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================
C
C 2nd order version for scalars
C
      DO I=1,IMR*(JNP+4)
         TMP = 0.25*(P(I,0) - P(I,-2))
         PMAX = MAX(P(I,0),P(I,-1),P(I,-2))-P(I,-1)
         PMIN = P(I,-1) - MIN(P(I,-1),P(I,-2),P(I,0))
         DC(I,-1) = SIGN(MIN(ABS(TMP),PMIN,PMAX),TMP)
      ENDDO
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%C Poles:
!%%%
!%%%      if(j1.ne.2) then
!%%%      do i=1,IMR
!%%%      DC(i,1) = 0.
!%%%      DC(i,JNP) = 0.
!%%%      enddo
!%%%      else
!%%%C Determine slopes in polar caps for scalars!
!%%%
!%%%      do 20 i=1,IMH
!%%%C South
!%%%      tmp = 0.25*(p(i,2) - p(i+imh,2))
!%%%      Pmax = max(p(i,2),p(i,1), p(i+imh,2)) - p(i,1)
!%%%      Pmin = p(i,1) - min(p(i,2),p(i,1), p(i+imh,2))
!%%%      DC(i,1)=sign(min(abs(tmp),Pmax,Pmin),tmp)
!%%%C North.
!%%%      tmp = 0.25*(p(i+imh,JMR) - p(i,JMR))
!%%%      Pmax = max(p(i+imh,JMR),p(i,jnp), p(i,JMR)) - p(i,JNP)
!%%%      Pmin = p(i,JNP) - min(p(i+imh,JMR),p(i,jnp), p(i,JMR))
!%%%      DC(i,JNP) = sign(min(abs(tmp),Pmax,pmin),tmp)
!%%% 20    continue
!%%%
!%%%C Scalars:
!%%%      do 25 i=imh+1,IMR
!%%%      DC(i,  1) =  - DC(i-imh,  1)
!%%%      DC(i,JNP) =  - DC(i-imh,JNP)
!%%% 25    continue
!%%%      endif
!%%% 
      ! Return to TPCORE
      END SUBROUTINE YMIST

!------------------------------------------------------------------------------

      SUBROUTINE YTP(IMR,JNP,ACOSP,DQ,Q,CRY,YMASS,FY2,JORD,FY1_TP,IGZD,
     &               JMAX)
C****6***0*********0*********0*********0*********0*********0**********72
      INTEGER IMR,JNP,JORD,JMAX,i,j,k,len,jt,IGZD
      TYPE (XPLEX) Q(1-IGZD:IMR+IGZD,1-IGZD:JNP+IGZD),
     &     CRY(1-IGZD:IMR+IGZD,1-IGZD:JNP+IGZD+1),
     &     YMASS(IMR,JNP+1), FY2(IMR,JNP+1),
     &     ACOSP(-10:JMAX),DQ(IMR,JNP)

! bey, 6/20/00. for mass-flux diagnostic
      TYPE (XPLEX) fy1_tp(IMR,JNP)
      
C============================================================================
C Cray NOBOUNDS directive will turn off all subscript bounds checking.
C This directive is also legal on SGI compilers (bmy, 4/24/00)
CDIR$ NOBOUNDS
C============================================================================

C Work array
      TYPE (XPLEX) fy1(IMR,JNP+1),P(IMR,-2:JNP+3),C(IMR,JNP+1), 
     &     AR(IMR,0:JNP+1),
     &     AL(IMR,0:JNP+1),
     &     A6(IMR,0:JNP+1),DC2(IMR,-1:JNP+2)
C
      LEN = IMR*(JNP+1)
      
      P(:,:)=0D0

      DO J = -2, JNP+3
      DO I =  1, IMR
         P(I,J)=Q(I,J)
      ENDDO
      ENDDO
      
      DO J = 1, JNP+1
      DO I = 1, IMR
         C(I,J)=CRY(I,J)
      ENDDO
      ENDDO

      IF(JORD.EQ.1) THEN

         DO I=1,LEN
            JT = 1. - C(I,1)
            FY1(I,1) = P(I,JT)
         ENDDO

         DO I=1,LEN
            FY2(I,1) = 0.
         ENDDO

      ELSE
         CALL YMIST(IMR,JNP,P(:,:),DC2(:,:),IGZD)

         IF(JORD.LE.0 .OR. JORD.GE.3) THEN
            CALL FYPPM(C(:,:),P(:,:),DC2(:,:),FY1(:,:),FY2(:,:),
     &           IMR,JNP,A6(:,:),AR(:,:),AL(:,:),JORD,IGZD)
         ELSE
            DO I=1,LEN
               JT = XPLX(1) - C(I,1)
               FY1(I,1) = P(I,JT)
               FY2(I,1) = (SIGN(1d0,C(I,1))-C(I,1))*DC2(I,JT)
            ENDDO
         ENDIF
      ENDIF
C
      DO I=1,LEN
         FY1(I,1) = FY1(I,1)*YMASS(I,1)
         FY2(I,1) = FY2(I,1)*YMASS(I,1)
      ENDDO
C
!=============================================================================
! This loop had to be extended for the mass-flux diagnostics (bmy, 4/26/00)
!      DO 1400 j=j1,j2
!      DO 1400 i=1,IMR
!1400  DQ(i,j) = DQ(i,j) + (fy1(i,j) - fy1(i,j+1)) * acosp(j)
!=============================================================================
!%%% 
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% use extrapolation to calculate fy1 and fy2 at grid JNP+1
!%%%      do i=1,IMR
!%%%        fy1(i,JNP+1)=2*fy1(i,JNP)-fy1(i,JNP-1)
!%%%         fy2(i,JNP+1)=2*fy2(i,JNP)-fy2(i,JNP-1)
!%%%      enddo
!%%%
      DO J = 1, JNP
      DO I = 1, IMR
         DQ(I,J) = DQ(I,J) + (FY1(I,J) - FY1(I,J+1)) * ACOSP(J)

         ! bey, 6/20/00. for mass-flux diagnostic
         FY1_TP(I,J) = FY1(I,J) 
      ENDDO
      ENDDO
!%%% 
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested simulation
!%%%C Poles
!%%%      sum1 = fy1(IMR,j1  )
!%%%      sum2 = fy1(IMR,J2+1)
!%%%      do i=1,IMR-1
!%%%      sum1 = sum1 + fy1(i,j1  )
!%%%      sum2 = sum2 + fy1(i,J2+1)
!%%%      enddo
!%%%
!%%%      sum1 = DQ(1,  1) - sum1 * RCAP
!%%%      sum2 = DQ(1,JNP) + sum2 * RCAP
!%%%      do i=1,IMR
!%%%      DQ(i,  1) = sum1
!%%%      DQ(i,JNP) = sum2
!%%%      enddo
!%%%
!%%%      if(j1.ne.2) then
!%%%      do i=1,IMR
!%%%      DQ(i,  2) = sum1
!%%%      DQ(i,JMR) = sum2
!%%%      enddo
!%%%      endif
!%%%
      ! Return to TPCORE
      END SUBROUTINE YTP

!------------------------------------------------------------------------------

      SUBROUTINE PRESS_FIX( FX,    FY,   NDT,  ACOSP, Jmax,
     &                      I0_W,  J0_W, IM_W, JM_W,  IGZD )
!
!******************************************************************************
!  Subroutine PRESS_FIX is a wrapper for the Pressure fixer DYN0.  PRESS_FIX
!  takes the mass fluxes in pressure units and converts them to [kg air/s]
!  using the correct geometry for TPCORE. (bdf, bmy, 3/10/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FX     (TYPE (XPLEX) ) : E-W flux passed from TPCORE [mb/timestep]
!  (2 ) FY     (TYPE (XPLEX) ) : N-S flux passed from TPCORE [mb/timestep]
!  (3 ) NDT    (INTEGER) : Dynamic timestep for TPCORE [s]
!  (4 ) ACOSP  (TYPE (XPLEX) ) : Array of inverse cosines    [unitless]
!  (5 ) J1     (INTEGER) : TPCORE polar cap extent     [# of boxes]
!  (6 ) I0_W   (INTEGER) : TPCORE REGION longitude offset [# boxes]
!  (7 ) J0_W   (INTEGER) : TPCORE REGION latitude  offset [# boxes]
!  (8 ) IM_W   (INTEGER) : TPCORE REGION longitude extent [# boxes]
!  (9 ) JM_W   (INTEGER) : TPCORE REGION latitude  extent [# boxes]
!  (10) IGZD   (INTEGER) : Variable equal to 1-I0_W or 1-J0_W
!
!  NOTES:
!  (1 ) Differences from "tpcore_mod" denoted by !%%% (bmy, 3/10/03)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_M2
      USE TIME_MOD, ONLY : GET_TS_DYN

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! Diagnostic switches
#     include "CMN_GCTM"  ! g0_100

      ! Arguments
!%%%
!%%% MODIFICATIONS FOR NESTED GRID
!%%% Modify array boundarys accordingly 
!%%%
      INTEGER, INTENT(IN)    :: NDT, Jmax, I0_W, J0_W, IM_W, JM_W, IGZD
      TYPE (XPLEX),  INTENT(IN)    :: ACOSP(-10:Jmax)

      TYPE (XPLEX),  INTENT(INOUT) :: FX(1-IGZD:IM_W+IGZD+1,
     &                             1-IGZD:JM_W+IGZD,LLPAR)

      TYPE (XPLEX),  INTENT(INOUT) :: FY(1-IGZD:IM_W+IGZD,
     &                             1-IGZD:JM_W+IGZD+1,LLPAR)

      ! Local variables
      INTEGER                :: I, J, K, K2, L
      TYPE (XPLEX)                 :: DTC, DTDYN, NSDYN, SUM1, SUM2
!%%%
!%%% MODIFICATION FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%      TYPE (XPLEX)                 :: NP_FLUX(IIPAR,LLPAR)
!%%%      TYPE (XPLEX)                 :: SP_FLUX(IIPAR,LLPAR)
!%%%
      TYPE (XPLEX)::ALFA( 1-IGZD:IM_W+IGZD+1,1-IGZD:JM_W+IGZD,  LLPAR  )
      TYPE (XPLEX)::BETA( 1-IGZD:IM_W+IGZD,  1-IGZD:JM_W+IGZD+1,LLPAR  )
      TYPE (XPLEX)::GAMA( 1-IGZD:IM_W+IGZD,  1-IGZD:JM_W+IGZD,  LLPAR+1)
      TYPE (XPLEX)::UMFLX(1-IGZD:IM_W+IGZD+1,1-IGZD:JM_W+IGZD,  LLPAR  )
      TYPE (XPLEX)::VMFLX(1-IGZD:IM_W+IGZD,  1-IGZD:JM_W+IGZD+1,LLPAR  )

      ! Local SAVEd variables
      LOGICAL, SAVE :: FIRST = .TRUE.
      TYPE (XPLEX),  SAVE :: DXYP(JJPAR)

      !=================================================================
      ! PRESS_FIX begins here!
      !
      ! K  is the vertical index down from the atmosphere top downwards
      ! K2 is the vertical index up from the surface
      !=================================================================

      ! Initialize arrays
      ALFA  = 0d0
      BETA  = 0d0
      GAMA  = 0d0
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Now use function GET_TS_DYN to get the dynamic timestep 
!%%%      ! NSDYN is the dynamic time step in seconds
!%%%      NSDYN = NDYN * 60d0
      ! Dynamic timestep [s]
      NSDYN = GET_TS_DYN() * 60d0
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%     ! J2 is the south polar edge
!%%%     J2    = JJPAR - J1 + 1
!%%%
!%%%
      ! DTDYN = TYPE (XPLEX) value for NDT, the dynamic timestep
      DTDYN = XPLX( NDT )

      ! Save grid box surface areas [m2] into the local DXYP array
      IF ( FIRST ) THEN
         DO J = 1, JJPAR
            DXYP(J) = GET_AREA_M2( J )
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! FX is the E-W mass flux from TPCORE in [mb/timestep].
      ! UMFLX is the mass flux in [kg air/s], which is what DYN0 needs.
      !
      ! FY is the E-W mass flux from TPCORE in [mb/timestep].
      ! VMFLX is the mass flux in [kg air/s], which is what DYN0 needs.
      !
      ! The unit conversion from [mb/timestep] to [kg air/s] is:
      !
      !   mb  | 100 Pa | 1 kg air |  s^2  |  step  | DXYP m^2     kg air
      ! ------+--------+----------+-------+--------+---------- = -------
      !  step |   mb   | Pa m s^2 | 9.8 m | DTDYN s|    s           s
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, K2, DTC )
#endif
      DO K = 1, LLPAR
         K2 = LLPAR - K + 1

         ! Compute UMFLX from FX
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits 
!%%%         DO J = 1, JJPAR
!%%%         DO I = 1, IIPAR
!%%%
         DO J = 1-IGZD, JM_W+IGZD
         DO I = 1-IGZD, IM_W+IGZD+1
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Now use DXYP(J+J0_w) instead of DXYP(J) 
!%%%         UMFLX(I,J,K2) = FX(I,J,K) * ( G0_100 * DXYP(J) ) / DTDYN
!%%%
            UMFLX(I,J,K2) = FX(I,J,K) * ( G0_100 * DXYP(J+J0_w))/ DTDYN
         ENDDO
         ENDDO

         ! Compute VMFLX from FY
         DO I = 1-IGZD, IM_W+IGZD
         DO J = 1-IGZD, JM_W+IGZD+1
            IF ( FY(I,J,K) .GE. 0 ) THEN
               DTC = FY(I,J,K) * G0_100 * ACOSP(J)  * DXYP(J+J0_w)/DTDYN
            ELSE
            DTC = FY(I,J,K) * G0_100 * ACOSP(J-1)* DXYP(J-1+J0_w)/DTDYN
            ENDIF
            
            VMFLX(I,J,K2) = DTC
         ENDDO
         ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         !=================================================================
!%%%         ! TREATMENT OF THE POLES: 1
!%%%         ! copy ymass values strait into vmflx at poles for pressure fixer
!%%%         !=================================================================
!%%%         DO I = 1, IIPAR
!%%%            VMFLX(I,1,K2)     = FY(I,1,K)
!%%%            VMFLX(I,J1-1,K2)  = FY(I,J1-1,K)
!%%%            VMFLX(I,JJPAR,K2) = FY(I,JJPAR,K)
!%%%         ENDDO
!%%%
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%      DO K = 1, LLPAR
!%%%
!%%%         !=================================================================
!%%%         ! TREATMENT OF THE POLES: 2
!%%%         ! North polar cap: J=1
!%%%         !=================================================================
!%%%         SUM1 = FY(IIPAR,J1,K)
!%%%         DO I = 1, IIPAR-1
!%%%            SUM1 = SUM1 + FY(I,J1,K)
!%%%         ENDDO
!%%%
!%%%         ! NORTH POLE FLUX IN KG.
!%%%         DO I = 1, IIPAR
!%%%            NP_FLUX(I,K) = SUM1 * G0_100 * ACOSP(1) * DXYP(1)
!%%%         ENDDO
!%%%
!%%%         !=================================================================
!%%%         ! TREATMENT OF THE POLES: 3
!%%%         ! South polar cap: J=JJPAR
!%%%         !=================================================================
!%%%         SUM2 = FY(IIPAR,J2+1,K)
!%%%         DO I = 1, IIPAR-1
!%%%            SUM2 = SUM2 + FY(I,J2+1,K)
!%%%         ENDDO
!%%%
!%%%         DO I = 1, IIPAR
!%%%            SP_FLUX(I,K) = SUM2 * G0_100 * ACOSP(JJPAR) * DXYP(JJPAR)
!%%%         ENDDO
!%%%      ENDDO

      !=================================================================
      ! Call DYN0 to fix the pressures
      !=================================================================
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Remove J1, NP_FLUX, SP_FLUX from call to DYN0
!%%%
      CALL DYN0( NSDYN, UMFLX, VMFLX, ALFA, BETA, GAMA, 
     &           I0_W,  J0_W,  IM_W,  JM_W, IGZD )

      !=================================================================
      ! ALFA is the E-W mass flux adjusted by DYN0 in [kg air/s]
      ! FX   is the E-W mass flux for TPCORE       in [mb/timestep].
      !
      ! BETA is the N-S mass flux adjusted by DYN0 in [kg air/s]
      ! FY   is the E-W mass flux for TPCORE       in [mb/timestep].
      !
      ! The unit conversion from to [kg air/s] to [mb/timestep] is:
      !
      !  kg air | Pa m s^2 | 9.8 m |    1     | DTDYN s |  mb       mb
      ! --------+----------+-------+----------+---------+------- = ----
      !    s    | 1 kg air |  s^2  | DXYP m^2 |  step   | 100 Pa   step
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, K2 )
#endif
      DO K = 1, LLPAR
         K2 = LLPAR - K + 1
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits for nested grid
!%%% Now use DXYP(J+J0_W) instead of DXYP(J) in several DO-loops below
!%%%

         ! Update FX from ALFA
         DO J = 1-IGZD, JM_W+IGZD
         DO I = 1-IGZD, IM_W+1+IGZD
            FX(I,J,K) = ALFA(I,J,K2) * DTDYN /(G0_100 * DXYP(J+J0_W) )
         ENDDO
         ENDDO

         ! Update FY from BETA
         DO I = 1-IGZD, IM_W+IGZD
         DO J = 1-IGZD, JM_W+1+IGZD
            IF ( BETA(I,J,K) .GE. 0 ) THEN
               FY(I,J,K) = BETA(I,J,K2) * DTDYN /
     &                     ( G0_100 * ACOSP(J) * DXYP(J+J0_w) )
            ELSE
               FY(I,J,K) = BETA(I,J,K2) * DTDYN /
     &                     ( G0_100 * ACOSP(J-1) * DXYP(J-1+J0_w) )
            ENDIF
         ENDDO
         ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         ! Special treatment of BETA at the poles
!%%%         DO I = 1, IIPAR
!%%%            FY(I,1,K)     = BETA(I,1,K2)
!%%%            FY(I,J1-1,K)  = BETA(I,J1-1,K2)
!%%%            FY(I,JJPAR,K) = BETA(I,JJPAR,K2)
!%%%         ENDDO
!%%%
      ENDDO

      ! Return to calling program
      END SUBROUTINE PRESS_FIX

!------------------------------------------------------------------------------

      SUBROUTINE DYN0( DTWIND, UMFLX, VMFLX, ALFA, BETA, GAMA,
     &                 I0_W,   J0_W,  IM_W,  JM_W, IGZD )
!
!******************************************************************************
!  Subroutine DYN0 is the pressure fixer for TPCORE.  DYN0 readjusts the
!  mass fluxes ALFA, BETA, GAMA, so that they are consistent with the
!  met fields. (bdf, bmy, 3/10/03, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DTWIND  (TYPE (XPLEX) ) : Time step between wind intervals [s]
!  (2 ) UMFLX   (TYPE (XPLEX) ) : Wet air mass flux in E-W     direction [kg air/s]
!  (3 ) VMFLX   (TYPE (XPLEX) ) : Wet air mass flux in N-S     direction [kg air/s]
!  (4 ) ALFA    (TYPE (XPLEX) ) : Dry air mass flux in E-W     direction [kg air/s]
!  (5 ) BETA    (TYPE (XPLEX) ) : Dry air mass flux in N-S     direction [kg air/s]
!  (6 ) GAMA    (TYPE (XPLEX) ) : Dry air mass flux in up/down direction [kg air/s]
!  (7 ) I0_W   (INTEGER) : TPCORE REGION longitude offset [# boxes]
!  (8 ) J0_W   (INTEGER) : TPCORE REGION latitude  offset [# boxes]
!  (9 ) IM_W   (INTEGER) : TPCORE REGION longitude extent [# boxes]
!  (10) JM_W   (INTEGER) : TPCORE REGION latitude  extent [# boxes]
!  (11) IGZD   (INTEGER) : Variable equal to 1-I0_W or 1-J0_W
!
!  Arguments as Output:
!  ============================================================================
!  (8 ) ALFA    (TYPE (XPLEX) ) : ALFA air mass, after pressure fix is applied
!  (9 ) BETA    (TYPE (XPLEX) ) : BETA air mass, after pressure fix is applied
!  (10) GAMA    (TYPE (XPLEX) ) : GAMA air mass, after pressure fix is applied
!
!  NOTES:
!  (1 ) Differences from "tpcore_mod.f" denoted by !%%% (bmy, 3/10/03)
!  (2 ) Removed reference to CMN (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : SPHU, PSC2, AIRDEN, AIRVOL
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PRESSURE_MOD, ONLY : GET_BP

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%      INTEGER, INTENT(IN)    :: J1
!%%%      TYPE (XPLEX),  INTENT(IN)    :: NP_FLUX(IIPAR,LLPAR)
!%%%      TYPE (XPLEX),  INTENT(IN)    :: SP_FLUX(IIPAR,LLPAR)
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Change array limits accordingly
!%%%
      INTEGER, INTENT(IN)    :: I0_W, J0_W, IM_W, JM_W, IGZD
      TYPE (XPLEX),  INTENT(IN)    :: DTWIND

      TYPE (XPLEX),  INTENT(IN)    :: UMFLX(1-IGZD:IM_W+IGZD+1,
     &                                1-IGZD:JM_W+IGZD,LLPAR)

      TYPE (XPLEX),  INTENT(IN)    :: VMFLX(1-IGZD:IM_W+IGZD,
     &                               1-IGZD:JM_W+IGZD+1,LLPAR)

      TYPE (XPLEX),  INTENT(INOUT) :: ALFA(1-IGZD:IM_W+IGZD+1,
     &                               1-IGZD:JM_W+IGZD,LLPAR)

      TYPE (XPLEX),  INTENT(INOUT) :: BETA(1-IGZD:IM_W+IGZD,
     &                               1-IGZD:JM_W+IGZD+1,LLPAR)

      TYPE (XPLEX),  INTENT(INOUT) :: GAMA(1-IGZD:IM_W+IGZD,
     &                               1-IGZD:JM_W+IGZD,LLPAR+1)

      ! Local variables
      LOGICAL       :: LSP, LNP, LEW
      INTEGER       :: IIX, JJX, KM, JB, JE, IEPZ, IMZ
      INTEGER       :: I, J, J2, K, L
      TYPE (XPLEX)        :: ALFAX, UFILT, VFILT, PCTM8, G100
      TYPE (XPLEX)        :: AIRQAV, AWE, SUMAD0, SUMAW0, AIRWET
      TYPE (XPLEX)        :: AIRH2O, AIRQKG ,SUM1, SUMA, SUMP, SUMQ
      TYPE (XPLEX)        :: SUMU, SUMV, SUMW, ZIMZ, ZDTW, G0
      TYPE (XPLEX)        :: AD_L(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)     :: AIRD(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD,LLPAR)
      TYPE (XPLEX)    :: AIRNEW(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD,LLPAR)
      TYPE (XPLEX)     :: AIRX(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD,LLPAR)
      TYPE (XPLEX)        :: AX(1-IGZD:IM_W+IGZD+1,1-IGZD:JM_W+IGZD)
      TYPE (XPLEX)        :: BX(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD+1)
      TYPE (XPLEX)        :: MERR(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD)
      TYPE (XPLEX)        :: PCTM(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD)
      TYPE (XPLEX)        :: PERR(1-IGZD:IM_W+IGZD,1-IGZD:JM_W+IGZD)
      TYPE (XPLEX)        :: SPHU_KG(IIPAR,JJPAR,LLPAR)
      TYPE (XPLEX)       :: SUMAQ(1-IGZD:IM_W+IGZD+1,1-IGZD:JM_W+IGZD+1)
      TYPE (XPLEX)        :: XYB(IIPAR,JJPAR)
      TYPE (XPLEX)        :: XYZB(IIPAR,JJPAR,LLPAR)

      ! Local saved variables
      LOGICAL, SAVE :: FIRST = .TRUE.
      TYPE (XPLEX),  SAVE :: DXYP(JJPAR)
      TYPE (XPLEX),  SAVE :: DSIG(LLPAR)

      !=================================================================
      ! DYN0 begins here!
      !
      ! UNITS OF AIR MASS AND TRACER = (kg)
      !
      ! Air mass (kg) is given by:
      !    area (m^2) * pressure thickness (Pa) / g0
      !
      ! DXYP(I,J)      = area of [I,J] [m^2]
      !
      ! PSC2(I,J)        = surf pressure [Pa] averaged in extended zone.
      !
      ! SPHU_KG(I,J,K) = specific humidity of grid box
      !                  [kg H2O/kg wet air] averaged in extended zone.
      !
      ! AIRQKG(I,J)    = Mass of H2O [kg] at each level
      !                = PSC2(I,J)) * SPHU_KG(I,J,K)
      !
      ! AIRD(I,J,K)    = dry-air mass [kg] in each box as calculated
      !                  in CTM at the beginning of each time step,
      !                  updated at end of DYN0.
      !
      ! PCTM(I,J)      = inferred wet-air (total) surf press [Pa] calc.
      !                  in CTM (using SUMAQ & AIRD-X-NEW)
      !
      ! AIRNEW(I,J,K)  = new dry-air mass in each CTM box after
      !                  horizontal divergence (ALFA+BETA) over time
      !                  step DTWIND (sec)
      !
      ! AIRX(I,J,K)    = expected dry-air mass in each CTM box after
      !                  calculating the vertical divergence (GAMA)
      !                  (also used for GCM dry mass)
      !                = XYZA(I,J,K) + XYZB(I,J,K)*PCTM(I,J) - AIRQKG
      !
      ! DTWIND         = time step [s] that applies to the averaged
      !                  wind fields (i.e., the time between successive
      !                  pressures.
      !
      !-----------------------------------------------------------------
      !
      ! Assume that we have "wet-air" mass fluxes across each boundary
      !
      !   UMFLX(I,J,K)   ==> [I,J,K] ==>  UMFLX(I+1,J,K)   [kg air/s]
      !   VMFLX(I,J,K)   ==> [I,J,K] ==>  VMFLX(I,J+1,K)   [kg air/s]
      !
      ! Convert to "dry-air" mass flux in/out of box using
      ! average Q at boundary
      !
      !   ALFA(I,J,K)   ==> [I,J,K] ==>  ALFA(I+1,J,K)     [kg air/s]
      !   BETA(I,J,K)   ==> [I,J,K] ==>  BETA(I,J+1,K)     [kg air/s]
      !
      ! Calculate convergence in each layer of dry air, compare with
      ! expected dry air mass (AIRX) and then calculate vertical
      ! dry-mass fluxes
      !
      !   GAMA(I,J,K)   ==> [I,J,K] ==>  GAMA(I,J,K+1)     [kg air/s]
      !
      ! Horizontal pressure filter adjusts UMFLX & VMFLX to reduce
      ! error in [PCTM - PSC2]
      !
      !     UMFLX + pressure filter ==> UMFLX#,
      !     VMFLX + filter ==> VMFLX# (temporary)
      !
      ! The pressure filter does nearest neighbor flux
      ! (adjusting ALFA/BETA)
      !
      !-----------------------------------------------------------------
      !
      ! Note that K->K+1 is downward (increasing pressure) and
      ! that boundaries:
      !   GAMA(I,J,1) = GAMA(I,J,KM+1)    = 0   no flux across
      !                                         upper/lower boundaries
      !
      !   BETA(I,1,K) = BETA(I,JJPAR+1,K) = 0   no flux at S & N poles
      !
      !   ALFA(1,J,K) = ALFA(IIPAR+1,J,K)       is NOT ZERO, but cyclic
      !
      ! Dimensions for ALFA, BETA, GAMA are extended by +1 beyond grid
      ! to allow simple formulation of fluxes in/out of final grid box.
      !
      ! GCM input UMFLX,VMFLX,PSG is ALWAYS of GLOBAL dimensions
      ! (IIPAR x JJPAR x LLPAR)
      !
      ! Indices of ALFA, BETA, GAMA, SPHU_KG & PS are always LOCAL
      ! (IIPAR x JJPAR x KM): FOR GEOS-CHEM, KM = LLPAR (bmy
      !
      ! Indices of tracer (STT), and diagnostics are local
      ! (w.r.t. WINDOW.  WINDOW calculations are defined by an
      ! offset and size
      !
      !         I0 .ge.0 and IIPAR+I0 .le. IIPAR
      !         J0 .ge.0 and JJPAR+J0 .le. JJPAR
      !         K0 .ge.0 and KM+K0 .le. LLPAR
      !
      ! The WINDOW calculation must allow for a boundary layer
      ! of grid boxes:
      !
      !         IG(abs. coords) = IW(in window) + I0
      !         JG(abs. coords) = JW(in window) + J0
      !         KG(abs. coords) = KW(in window) + K0
      !
      ! vertical window (NEW) allows for an upper boundary with flow
      ! across it and specified mixing ratio b.c.'s at KG = K0
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         
         ! Surface area [m2]
         DO J = 1, JJPAR
            DXYP(J) = GET_AREA_M2( J )
         ENDDO

         ! Sigma-level thickness [unitless]
         ! Assumes we are using a pure-sigma grid
         DO L = 1, LLPAR        
            DSIG(L) = GET_BP(L) - GET_BP(L+1)
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%! for tpcore poles.
!%%%      J2 = JJPAR - J1 + 1
!%%%
      ! geos code
      G0 = 9.8d0

      !=================================================================
      ! XYZB is the factor needed to get mass in kg of gridbox
      ! mass (kg) = XYZB (kg/mb) * P (mb)
      !
      ! AD_L is the dry air mass in the grid box
      !
      ! SPHU_KG is the water vapor [kg H2O/kg air]
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
#endif
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         XYZB(I,J,L)    = DSIG(L) * DXYP(J) * 1.d2 / G0
         AD_L(I,J,L)    = AIRDEN(L,I,J) * AIRVOL(I,J,L)
         SPHU_KG(I,J,L) = SPHU(I,J,L) / 1000d0
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! XYB is the factor needed to get mass in kg of column
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
#endif
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         XYB(I,J) = SUM( XYZB(I,J,1:LLPAR) )
      ENDDO
      ENDDO

      !=================================================================
      ! Define other variables
      !=================================================================
      G100   = 100.D0 / G0
      ZDTW   =   1.D0 / DTWIND
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%      LSP    = (       J0 .EQ. 0     )
!%%%      LNP    = ( JJPAR+J0 .EQ. JJPAR )
!%%%      LEW    = (    IIPAR .EQ. IIPAR )
!%%%
      !=================================================================
      ! Initialize ALFA with UMFLX and BETA with VMFLX
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
#endif
      DO L = 1, LLPAR
         DO J = 1-IGZD, IGZD+JM_W
            DO I = 1-IGZD,IM_W+IGZD+1
               ALFA(I,J,L) = UMFLX(I,J,L)
            ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% don't need easternmost edge
!%%%            ALFA(IIPAR+1,J,L) = ALFA(1,J,L)
         ENDDO

!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%         DO J = 2, JJPAR
!%%%         DO I = 1, IIPAR
!%%%
          DO J = 1-IGZD, JM_W+IGZD+1
             DO I = 1-IGZD, IM_W+IGZD
                BETA(I,J,L) = VMFLX(I,J,L)
             ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% don't need southernmost edge
!%%%            DO I = 1, IIPAR
!%%%               BETA(I,1,L)       = 0.D0
!%%%               BETA(I,JJPAR+1,L) = 0.D0
!%%%            ENDDO
!%%%
         ENDDO
      ENDDO

      !=================================================================
      ! SUMAQ(I,J): column integral of water (kg)
      ! Check on air mass
      !=================================================================
      SUMAD0 = 0.D0
      SUMAW0 = 0.D0
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits for nested-grid
!%%% Also use I+I0_W, J+J0_W to reference several arrays below
!%%%
      DO J = 1-IGZD, JM_W+IGZD+1
      DO I = 1-IGZD, IM_W+IGZD+1
         SUMAQ(I,J) = 0.D0

         DO K = 1, LLPAR
            AIRWET     = PSC2(I+I0_w,J+J0_w) * XYZB(I+I0_w,J+J0_w,K)
            AIRH2O     = SPHU_KG(I+I0_w,J+J0_w,K) * AIRWET
            SUMAQ(I,J) = SUMAQ(I,J)     + AIRH2O
            SUMAD0     = SUMAD0         + AIRWET
            SUMAW0     = SUMAW0         + AIRH2O
         ENDDO
      ENDDO
      ENDDO

      SUMAD0 = SUMAD0 - SUMAW0

      !=================================================================
      ! Initialize AIRD, the dry-air mass [kg] in each box as calculated
      ! in CTM at the start of each time step, updated at end of DYN0.
      !
      ! Compute AIRNEW, the new dry-air mass in each CTM box after
      ! horizontal divergence (ALFA+BETA) over time step DTWIND (sec)
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K )
#endif
      DO K =  1, LLPAR
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits and
!%%% Also reference AD with (I+I0_W, J+J0_W, K) instead of (I,J,K)
!%%%      DO J = J1, J2
!%%%
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         AIRD(I,J,K)   = AD_L(I+I0_W,J+J0_W,K)
         AIRNEW(I,J,K) = AIRD(I,J,K) + DTWIND *
     &                                 ( ALFA(I,J,K) - ALFA(I+1,J,K) +
     &                                   BETA(I,J,K) - BETA(I,J+1,K) )
      ENDDO
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/11/03)
!%%% Polar-cap stuff, useless for nested grid
!%%%
!%%%      !=================================================================
!%%%      ! treatment of the poles for tpcore.
!%%%      ! j=2 and j=jjpar-1 don't have any airmass change.
!%%%      !=================================================================
!%%%!$OMP PARALLEL DO
!%%%!$OMP+DEFAULT( SHARED )
!%%%!$OMP+PRIVATE( I, K )
!%%%      DO K = 1, LLPAR
!%%%
!%%%         ! J=1
!%%%         DO I = 1, IIPAR
!%%%            AIRNEW(I,1,K) = AD_L(I,1,K) - NP_FLUX(I,K)
!%%%         ENDDO
!%%%
!%%%         ! J=JJPAR
!%%%         DO I = 1, IIPAR
!%%%            AIRNEW(I,JJPAR,K) = AD_L(I,JJPAR,K) + SP_FLUX(I,K)
!%%%         ENDDO
!%%%      ENDDO
!%%%!$OMP END PARALLEL DO
!%%%
!%%%      !=================================================================
!%%%      ! Average AIRNEW at the South pole
!%%%      !=================================================================
!%%%      ZIMZ = 1.D0 / DBLE( IIPAR )
!%%%
!%%%      IF ( LSP ) THEN
!%%%        JB = 2
!%%%
!%%%!$OMP PARALLEL DO
!%%%!$OMP+DEFAULT( SHARED )
!%%%!$OMP+PRIVATE( I, K, SUMA )
!%%%        DO K = 1, LLPAR
!%%%           SUMA = SUM( AIRNEW(1:IIPAR,1,K) ) * ZIMZ
!%%%
!%%%           DO I = 1, IIPAR
!%%%              AIRNEW(I,1,K) = SUMA
!%%%           ENDDO
!%%%        ENDDO
!%%%!$OMP END PARALLEL DO
!%%%      ELSE
!%%%         JB = 1
!%%%      ENDIF
!%%%      !=================================================================
!%%%      ! Average AIRNEW at the North pole
!%%%      !=================================================================
!%%%      IF ( LNP ) THEN
!%%%        JE = JJPAR - 1
!%%%
!%%%        ! poles, just average AIRNEW
!%%%!$OMP PARALLEL DO
!%%%!$OMP+DEFAULT( SHARED )
!%%%!$OMP+PRIVATE( I, K, SUMA )
!%%%        DO K = 1, LLPAR
!%%%           SUMA = SUM( AIRNEW(1:IIPAR,JJPAR,K) ) * ZIMZ
!%%%
!%%%           DO I = 1, IIPAR
!%%%              AIRNEW(I,JJPAR,K) = SUMA
!%%%           ENDDO
!%%%        ENDDO
!%%%!$OMP END PARALLEL DO
!%%%      ELSE
!%%%         JE = JJPAR
!%%%      ENDIF
!%%%
      !================================================================
      ! BEGIN FILTER of PRESSURE ERRORS
      !
      ! Define the error in surface pressure PERR expected at end of
      ! time step filter by error in adjacent boxes, weight by areas,
      ! adjust ALFA & BETA
      !
      ! PCTM(I,J) = new CTM wet-air column based on
      !             dry-air convergence (Pascals)
      ! PERR(I,J) = pressure-error between CTM-GCM at new time
      !             (before filter)
      !================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
#endif
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IGZD+IM_W
         PCTM(I,J) = SUM( AIRNEW(I,J,:) ) / XYB(I+I0_w,J+J0_w)
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         ! special case for j=2, jjpar-1 for tpcore pole configuration.
!%%%         IF ( J .eq. 2 .OR. J .eq. JJPAR-1 ) THEN
!%%%            PCTM(I,J) = PSC2(I,J)
!%%%         ENDIF
!%%%
         PERR(I,J) = PCTM(I,J) - PSC2(I+I0_w,J+J0_w)
         MERR(I,J) = PERR(I,J) * DXYP(J+J0_w) * G100
      ENDDO
      ENDDO
      
      !### Debug
      !write(6,*) 'before PFILTR'
      !write(6,*) 'sum MERR is ', sum(MERR)
      !write(6,*) 'sum AX is ', sum(AX)
      !write(6,*) 'sum BX is ', sum(BX)

      ! Call pressure filter
      CALL PFILTR( MERR,  AX,  BX, DXYP, IIPAR, JJPAR,
     &             IM_W, JM_W, 1,  JGLOB, J0_w, IGZD )

      !### Debug
      !write(6,*) 'after PFILTR'
      !write(6,*) 'sum MERR is ', sum(MERR)
      !write(6,*) 'sum AX is ', sum(AX)
      !write(6,*) 'sum BX is ', sum(BX)

      !=================================================================
      ! Calculate corrections to ALFA from the filtered AX
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, IIX, J, K, UFILT )
#endif
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+1+IGZD
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% We don't have to worry about wrapping around the world 
!%%%         IIX   = MIN(I+I0_w,IM_W+I0_w)
!%%%
         IIX   = I+I0_w    
         UFILT = AX(I,J) / ( XYB(IIX,J+J0_w) * DTWIND )

         DO K = 1, LLPAR
            ALFA(I,J,K) = ALFA(I,J,K) + UFILT * XYZB(IIX,J+J0_w,K)
         ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! Calculate corrections to BETA from the filtered BX
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JJX, K, VFILT )
#endif
      DO J = 1-IGZD, JM_W+IGZD+1
         JJX = J+J0_W
         IF ( J+J+79+79 .GT. 181 )  JJX =J+J0_W- 1         !(YXW_1X1)

         DO I = 1-IGZD, IM_W+IGZD
            VFILT = BX(I,J) / ( XYB(I+I0_W,JJX) * DTWIND )

            DO K = 1, LLPAR
               BETA(I,J,K) = BETA(I,J,K) + VFILT * XYZB(I+I0_W,JJX,K)
            ENDDO
         ENDDO
      ENDDO

      !=================================================================
      ! Calculate the corrected AIRNEW's & PCTM after P-filter:
      ! has changed ALFA+BETAs and ctm surface pressure (PCTM)
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K )
#endif
      DO K = 1, LLPAR
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         AIRNEW(I,J,K) = AIRD(I,J,K) + DTWIND *
     &                                 ( ALFA(I,J,K) - ALFA(I+1,J,K) +
     &                                   BETA(I,J,K) - BETA(I,J+1,K) )
      ENDDO
      ENDDO
      ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/11/03)
!%%% Polar-cap stuff, useless for nested grid
!%%%      !=================================================================
!%%%      ! Average the adjusted AIRNEW at the South pole
!%%%      !=================================================================
!%%%      ZIMZ = 1.D0 / DBLE( IIPAR )
!%%%
!%%%      IF ( LSP ) THEN
!%%%         JB = 2
!%%%
!%%%!$OMP PARALLEL DO
!%%%!$OMP+DEFAULT( SHARED )
!%%%!$OMP+PRIVATE( I, K, SUMA )
!%%%         DO K = 1, LLPAR
!%%%            SUMA = SUM( AIRNEW(1:IIPAR,1,K ) ) * ZIMZ
!%%%
!%%%            DO I = 1, IIPAR
!%%%               AIRNEW(I,1,K) = SUMA
!%%%            ENDDO
!%%%         ENDDO
!%%%!$OMP END PARALLEL DO
!%%%      ELSE
!%%%         JB = 1
!%%%      ENDIF
!%%%
!%%%      !=================================================================
!%%%      ! Average the adjusted AIRNEW at the North pole
!%%%      !=================================================================
!%%%      IF ( LNP ) THEN
!%%%         JE = JJPAR -1
!%%%
!%%%!$OMP PARALLEL DO
!%%%!$OMP+DEFAULT( SHARED )
!%%%!$OMP+PRIVATE( I, K, SUMA )
!%%%         DO K = 1, LLPAR
!%%%            SUMA = SUM( AIRNEW(1:IIPAR,JJPAR,K) ) * ZIMZ
!%%%
!%%%            DO I = 1,IIPAR
!%%%               AIRNEW(I,JJPAR,K) = SUMA
!%%%            ENDDO
!%%%         ENDDO
!%%%!$OMP END PARALLEL DO
!%%%      ELSE
!%%%         JE = JJPAR
!%%%      ENDIF
!%%%
      !=================================================================
      ! END OF PRESSURE FILTER
      !
      ! GAMA:  redistribute the new dry-air mass consistent with the
      ! new CTM surface pressure, rigid upper b.c., no change in PCTM
      !
      ! AIRX(I,J,K) = dry-air mass expected, based on PCTM
      !               PCTM(I,J) & PERR(I,J)
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K, PCTM8, AIRQKG )
#endif
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         PCTM8     = ( SUM( AIRNEW(I,J,:) ) + SUMAQ(I,J) ) /
     &                               XYB(I+I0_W,J+J0_W)
         PCTM(I,J) = PCTM8
         PERR(I,J) = PCTM8 - PSC2(I+I0_W,J+J0_W)

         DO K = 1, LLPAR
            AIRQKG      = SPHU_KG(I+I0_W,J+J0_W,K) *
     &             ( XYZB(I+I0_W,J+J0_W,K) * PSC2(I+I0_W,J+J0_W) )
            AIRX(I,J,K) = PCTM8 * XYZB(I+I0_W,J+J0_W,K) - AIRQKG
         ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! GAMA from top down to be consistent with AIRX, AIRNEW not reset!
      !=================================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, K )
#endif
      DO J = 1-IGZD, JM_W+IGZD
      DO I = 1-IGZD, IM_W+IGZD
         GAMA(I,J,LLPAR+1) = 0.D0

         DO K = LLPAR, 2, -1
            GAMA(I,J,K) = GAMA(I,J,K+1) - (AIRNEW(I,J,K) - AIRX(I,J,K))
         ENDDO

         ! GAMA(I,J,1) will not be exactly ZERO, but it must be set so!
         GAMA(I,J,1) = 0.D0

         DO K = 2, LLPAR
            GAMA(I,J,K) = GAMA(I,J,K) * ZDTW
         ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE DYN0

!------------------------------------------------------------------------------

      SUBROUTINE PFILTR( MERR, ALFAX, BETAX, AXY,   ID,   JD,
     &                   IM,   JM,    NITR,  JGLOB, J0_W, IGZD )
!
!******************************************************************************
!  Subroutine PFILTR applies the pressure-filter, the pressure
!  between predicted Ps(CTM) and Ps(GCM). (bdf, yxw, bmy, 3/10/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) MERR(IM,JM)    (TYPE (XPLEX) ) : mass error
!  (2  ) ALFAX(ID+1,JD) (TYPE (XPLEX) ) : perturbed ALFA by MERR
!  (3  ) BETAX(ID,JD+1) (TYPE (XPLEX) ) : perturbed BETA by MERR
!  (4  ) AXY(ID,JD)     (TYPE (XPLEX) ) : area of grid box (I,J) in [m^2]
!  (5-6) ID, JD         (INTEGER) : "Global" array dimensions for lon, lat
!  (7-8) IM, JM         (INTEGER) : "Window" array dimensions for lon, lat
!  (9  ) NITR           (INTEGER) : number of iterations (NITR .LE. 4)
!  (10 ) JGLOB          (INTEGER) : GLOBAL REGION latitude extent [# boxes]
!  (11 ) J0_W           (INTEGER) : TPCORE REGION latitude offset [# boxes]
!  (12 ) IGZD           (INTEGER) : Variable equal to 1-I0_W or 1-J0_W
!
!  Arguments as Output:
!  ============================================================================
!  (1  ) MERR(ID,JD)    (TYPE (XPLEX) ) : adjusted mass error
!  (2  ) ALFAX(ID+1,JD) (TYPE (XPLEX) ) : adjusted ALFAX
!  (3  ) BETAX(ID,JD+1) (TYPE (XPLEX) ) : adjusted BETAX
!
!  NOTES:
!  (1 ) Differences from "tpcore_mod.f" denoted by !%%% (bmy, 3/10/03)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      ! Arguments
!%%%
!%%% MODIFICATIONS FOR NESTED GRID
!%%% Remove LSP, LNP, LEW from arg list (bmy, 3/10/03)
!%%%      LOGICAL, INTENT(IN)    :: LSP,LNP,LEW
!%%%
      INTEGER,INTENT(IN):: ID, JD, IM, JM, NITR, JGLOB,J0_w,IGzd
      TYPE (XPLEX),  INTENT(IN)    :: AXY(JGLOB)
      TYPE (XPLEX),INTENT(INOUT) :: MERR(1-igzd:IM+igzd,1-igzd:JM+igzd)
      TYPE (XPLEX),INTENT(INOUT)::ALFAX(1-igzd:IM+igzd+1,1-igzd:JM+igzd)
      TYPE (XPLEX),INTENT(INOUT)::BETAX(1-igzd:IM+igzd,1-igzd:JM+igzd+1)

      ! Local variables
      LOGICAL                :: LPOLE
      INTEGER                :: I, J, K
      TYPE (XPLEX)                 :: X0(1-igzd:IM+igzd,1-igzd:JM+igzd)

      !=================================================================
      ! PFILTR begins here!
      !=================================================================

      ! LPOLE is true if J=1 is the SOUTH POLE and J=JM is the NORTH POLE
      ! (this is the way GEOS-CHEM is set up, so LPOLE should be TRUE!)
!%%%
!%%% MODIFICATIONS FOR NESTED GRID
!%%% Polar cap stuff, useless for nested-grid simulation
!%%%      LPOLE = ( LSP .AND. LNP )
!%%%
      LPOLE = .FALSE.

      ! Zero ALFAX, BETAX, save MERR in X0
      DO J = 1-IGZD, JM+IGZD
         DO I = 1-IGZD, IM+IGZD
            ALFAX(I,J) = 0.D0
            BETAX(I,J) = 0.D0
            X0(I,J)    = MERR(I,J)
         ENDDO

         ALFAX(IM+IGZD+1,J) = 0.D0
      ENDDO

      DO I = 1-IGZD, IM+IGZD
         BETAX(I,JM+IGZD+1) = 0.D0
      ENDDO

      !=================================================================
      ! Call LOCFLT to do the local filtering
      !=================================================================
      !write(6,*) 'before LOCFLT'
      !write(6,*) 'sum of MERR ', sum(MERR)
      !write(6,*) 'sum of ALFAX ', sum(ALFAX)
      !write(6,*) 'sum of BETAX ', sum(BETAX)
      
      CALL LOCFLT( MERR, ALFAX, BETAX, AXY, ID,  JD,
     &             IM,   JM,    5,   JGLOB, J0_w,Igzd )

      !write(6,*) 'After LOCFLT'
      !write(6,*) 'sum of MERR ', sum(MERR)
      !write(6,*) 'sum of ALFAX ', sum(ALFAX)
      !write(6,*) 'sum of BETAX ', sum(BETAX)

      !=================================================================
      ! Call POLFLT to do the pole filtering (if necessary)
      !=================================================================
      IF ( LPOLE ) THEN
         CALL POLFLT( MERR, BETAX, AXY, xplx(1.D0), ID, JD, IM, JM )
      ENDIF

      !=================================================================
      ! Compute mass error MERR and return
      ! MERR, ALFAX, and BETAX are now adjusted
      !=================================================================
      DO J = 1-IGZD, JM+IGZD
      DO I = 1-IGZD, IM+IGZD
         MERR(I,J) = X0(I,J) + ALFAX(I,J) - ALFAX(I+1,J)
     &                       + BETAX(I,J) - BETAX(I,J+1)
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE PFILTR

!------------------------------------------------------------------------------

      SUBROUTINE LOCFLT( XERR, AX, BX,   AXY,   ID,   JD,
     &                   IM,   JM, NITR, JGLOB, J0_W, IGZD )
!
!******************************************************************************
!  Subroutine LOCFLT applies the pressure-filter to non-polar boxes.
!  LOCFLT is called from subroutine PFILTR above (bdf, bmy, 3/10/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) XERR(ID,JD) (TYPE (XPLEX) ) : mass error
!  (2  ) AX(ID+1,JD) (TYPE (XPLEX) ) : perturbed ALFA by XERR
!  (3  ) BX(ID,JD+1) (TYPE (XPLEX) ) : perturbed BETA by XERR
!  (4  ) AXY(ID,JD)  (TYPE (XPLEX) ) : area of grid box (I,J) in [m^2]
!  (5-6) ID, JD      (INTEGER) : "Global" array dimensions for lon, lat
!  (7-8) IM, JM      (INTEGER) : "Window" array dimensions for lon, lat
!  (9  ) NITR        (INTEGER) : number of iterations (NITR .LE. 4)
!  (10 ) JGLOB       (INTEGER) : GLOBAL REGION latitude extent [# boxes]
!  (11 ) J0_W        (INTEGER) : TPCORE REGION latitude offset [# boxes]
!  (12 ) IGZD        (INTEGER) : Variable equal to 1-I0_W or 1-J0_W
!
!  Arguments as Output:
!  ============================================================================
!  (1  ) XERR(ID,JD) (TYPE (XPLEX) ) : adjusted mass error
!  (2  ) AX(ID+1,JD) (TYPE (XPLEX) ) : adjusted AX
!  (3  ) BX(ID,JD+1) (TYPE (XPLEX) ) : adjusted BX
!
!  NOTES:
!  (1 ) Differences from "tpcore_mod.f" denoted by "!%%%" (bmy, 3/10/03)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
!%%%
!%%% MODIFICATIONS FOR NESTED GRID
!%%% Remove LSP, LNP, LEW from arg list (bmy, 3/10/03)
!%%%      LOGICAL, INTENT(IN)    :: LSP, LNP, LEW
!%%%
      ! Arguments
      INTEGER, INTENT(IN)    :: IGZD, ID, JD, IM, JM, NITR, JGLOB, J0_W
      TYPE (XPLEX),  INTENT(IN)    :: AXY(JGLOB)
      TYPE (XPLEX),INTENT(INOUT) ::XERR(1-IGZD:IM+IGZD,1-IGZD:JM+IGZD )
      TYPE (XPLEX),INTENT(INOUT)::AX(1-IGZD:IM+IGZD+1,1-IGZD:JM+IGZD  )
      TYPE (XPLEX),INTENT(INOUT)::BX(1-IGZD:IM+IGZD,  1-IGZD:JM+IGZD+1)

      ! Local variables
      INTEGER                :: I, IA, NAZ, J, J1, J2, NFLTR
      TYPE (XPLEX)                 :: SUMA, FNAZ8
      TYPE (XPLEX)              :: X0( 1-IGZD:IM+IGZD, 1-IGZD:JM+IGZD ) 

      !=================================================================
      ! LOCFLT begins here!
      !
      ! Initialize corrective column mass flows (kg):  AX->alfa, BX->beta
      !=================================================================
      DO J = 1-IGZD, JM+IGZD
      DO I = 1-IGZD, IM+IGZD
         X0(I,J) = XERR(I,J)
      ENDDO
      ENDDO

      !=================================================================
      ! Iterate over mass-error filter
      ! accumulate corrections in AX & BX
      !=================================================================
      DO NFLTR = 1, NITR

         !==============================================================
         ! calculate AX = E-W filter
         !==============================================================

         ! Compute polar box limits
         J1 = 1  - IGZD
         J2 = JM + IGZD
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         IF ( LSP ) J1 = 2
!%%%         IF ( LNP ) J2 = JM - 1

         ! Loop over non-polar latitudes
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FNAZ8 )
#endif
         DO J = J1, J2

            ! Calculate pressure-filter E-W wind between boxes [I-1] & [I].
            ! Enhance filtered wind by size of EPZ, will redistribute
            ! later within
            FNAZ8 = 0.125d0

            DO I = 2-IGZD, IM+IGZD
               AX(I,J) = AX(I,J) + FNAZ8 *(XERR(I-1,J) - XERR(I,J))
            ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% We don't need to worry about wrapping around the date line
!%%%        ! calculate pressure-filter E-W wind at edges I=1 & I=IM+1
!%%%            IF ( LEW )  THEN
!%%%               AX(IM+1,J) = AX(IM+1,J) + FNAZ8 * (XERR(IM,J) -XERR(1,J))
!%%%               AX(1,J)    = AX(1,J)    + FNAZ8 * (XERR(IM,J) -XERR(1,J))
!%%%            ELSE
!%%%
            ! WINDOW, assume zero error outside window
            AX(1-IGZD,J)   =  AX(1-IGZD,J)    - FNAZ8 * XERR(1-IGZD,J)
            AX(IM+IGZD+1,J) = AX(IM+IGZD+1,J) + FNAZ8 * XERR(IM+IGZD,J)
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% We don't need to worry about wrapping around the date line
!%%%            ENDIF
!%%%
         ENDDO

         !==============================================================
         ! calculate BX = N-S filter, N-S wind between boxes [J-1] & [J]
         !==============================================================
#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FNAZ8 )
#endif
         DO J = 2-IGZD, JM+IGZD
            FNAZ8 = 0.25D0 * AXY(J+J0_W) / ( AXY(J+J0_W-1)
     &                                     + AXY(J+J0_W) )

            DO I = 1-IGZD, IM+IGZD
               BX(I,J) = BX(I,J) + FNAZ8 * ( XERR(I,J-1) - XERR(I,J) )
            ENDDO
         ENDDO

         ! enhance the filtering by factor of 2 ONLY into/out-of polar caps
         FNAZ8 = 0.5D0 * AXY(1-IGZD+J0_W) / ( AXY(IGZD+J0_W) +
     &                                        AXY(1-IGZD+J0_W) )
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         ! When LSP=TRUE then J=1 is SOUTH POLE
!%%%         IF ( LSP ) THEN
!%%%            DO I = 1, IM
!%%%               BX(I,2) = BX(I,2) + FNAZ8 * (XERR(I,1) -XERR(I,2))
!%%%            ENDDO
!%%%         ELSE
!%%%
            DO I = 1-IGZD, IM+IGZD
               BX(I,1-IGZD)= BX(I,1-IGZD) -0.5D0 *FNAZ8 * XERR(I,1-IGZD)
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%               BX(I,2)= BX(I,2) +0.5D0 *FNAZ8 * (XERR(I,1) - XERR(I,2))
!%%%
            ENDDO
!%%% 
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         !ENDIF
!%%%
         FNAZ8 = 0.5D0 * AXY(JM+IGZD+J0_W+1) / ( AXY(JM+IGZD+J0_W) 
     &                + AXY(JM+IGZD+J0_W+1) )

!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         ! When LNP=TRUE, then J=JM is NORTH POLE
!%%%         IF ( LNP )  THEN
!%%%            DO I = 1, IM
!%%%               BX(I,JM) = BX(I,JM) +FNAZ8 *(XERR(I,JM-1) -XERR(I,JM))
!%%%            ENDDO
!%%%         ELSE
!%%%
            DO I = 1-IGZD,IM+IGZD
            BX(I,JM+IGZD+1)= BX(I,JM+IGZD+1)+0.5D0*FNAZ8*XERR(I,JM+IGZD)
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%               BX(I,JM)  = BX(I,JM)   + 0.5D0 *FNAZ8 *
!%%%      &                                  (XERR(I,JM-1) -XERR(I,JM))
!%%%
            ENDDO
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Polar cap stuff, useless for nested grid
!%%%         ENDIF
!%%%
         !==============================================================
         ! need N-S flux across boundaries if window calculation
         ! (assume XERR=0 outside)
         !
         ! JM for optimal matrix/looping, it would be best to
         ! define XERR=0 for an oversized array XERR(0:IM+1,0:JM+1)
         ! Update the mass error (XERR)
         !==============================================================
         DO J = 1-IGZD, JM+IGZD
         DO I = 1-IGZD, IM+IGZD
            XERR(I,J) = X0(I,J) + AX(I,J) - AX(I+1,J)
     &                          + BX(I,J) - BX(I,J+1)
         ENDDO
         ENDDO

      ENDDO  ! NFLTR

      ! Return to calling program
      END SUBROUTINE LOCFLT

!------------------------------------------------------------------------------

      SUBROUTINE POLFLT( XERR, BX, AXY, COEF, ID, JD, IM, JM )
!
!******************************************************************************
!  Subroutine POLFLT applies the pressure-filter to polar boxes.
!  POLFLT is called from subroutine PFILTR above (bdf, bmy, 3/10/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) XERR(ID,JD) (TYPE (XPLEX) ) : mass error
!  (2  ) BX(ID,JD+1) (TYPE (XPLEX) ) : perturbed BETA by XERR
!  (3  ) AXY(ID,JD)  (TYPE (XPLEX) ) : area of grid box (I,J) in [m^2]
!  (4  ) COEF        (TYPE (XPLEX) ) : Multiplicative coefficient ?????
!  (5-6) ID, JD      (INTEGER) : "Window" array dimensions for lon, lat
!  (7-8) IM, JM      (INTEGER) : "Global" array dimensions for lon, lat
!
!  Arguments as Output:
!  ============================================================================
!  (1  ) XERR(ID,JD) (TYPE (XPLEX) ) : adjusted mass error
!  (2  ) BX(ID,JD+1) (TYPE (XPLEX) ) : adjusted BX
!
!  NOTES:
!  (1  ) Differences from "tpcore_mod" denoted by !%%% (bmy, 3/10/03)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)    :: ID, JD, IM, JM
      TYPE (XPLEX),  INTENT(IN)    :: AXY(JD)
      TYPE (XPLEX),  INTENT(IN)    :: COEF
      TYPE (XPLEX),  INTENT(INOUT) :: XERR(ID,JD)
      TYPE (XPLEX),  INTENT(INOUT) :: BX(ID,JD+1)

      ! Local variables
      INTEGER                :: I, J
      TYPE (XPLEX)                 :: ERAV, BXJ(JD+1), TOTAL

      !=================================================================
      ! POLFLT begins here!
      !
      ! Initialize corrective column mass flows (kg):  BXJ->beta
      !=================================================================
      DO I = 1, IM

         ! Initialize
         ERAV  = 0.D0
         TOTAL = 0.D0

         ! Sum XERR in ERAV and sum AXY in TOTAL
         DO J = 1, JM
            ERAV  = ERAV  + XERR(I,J)
            TOTAL = TOTAL + AXY(J)
         ENDDO

         ! Compute area-weighted mass error total
         ERAV = ERAV / TOTAL

         ! mass-error filter, make corrections in BX
         BXJ(1) = 0.D0

         DO J = 2, JM
            BXJ(J) = BXJ(J-1) + XERR(I,J-1) - AXY(J-1) * ERAV
         ENDDO

         DO J = 2, JM
            BX(I,J) = BX(I,J) + COEF * BXJ(J)
         ENDDO

      ENDDO  ! I

      ! Update XERR
      DO J = 1, JM
      DO I = 1, IM
         XERR(I,J) = XERR(I,J) + BX(I,J) - BX(I,J+1)
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE POLFLT

!------------------------------------------------------------------------------

      SUBROUTINE DIAG_FLUX( IC,   FX,     FX1_TP, FY,    FY1_TP, 
     &                      FZ,   FZ1_TP, NDT,    ACOSP, Jmax, 
     &                      I0_W, J0_W,   IM_W,   JM_W,  IGZD )
!
!******************************************************************************
!  Subroutine DIAG_FLUX archives the mass fluxes in TPCORE version 7.1.
!  (bey, bmy, 9/20/00, 6/20/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1  ) IC        (INTEGER) : Current tracer # 
!  (2,3) FX,FX1_TP (TYPE (XPLEX) ) : Flux into the west side of grid box (I,J,K) 
!  (4,5) FY,FY1_TP (TYPE (XPLEX) ) : Flux into the south side of grid box (I,J,K)
!  (6,7) FZ,FZ1_TP (TYPE (XPLEX) ) : Flux into top of grid box (I,J,K) 
!  (8  ) NDT       (INTEGER) : Dynamic timestep in seconds
!  (9  ) ACOSP     (INTEGER) : Inverse cosine at latitude (J)
!  (10 ) Jmax      (INTEGER) : Max dimension for TPCORE internal arrays
!  (11 ) I0_W      (INTEGER) : TPCORE REGION longitude offset [# boxes]
!  (12 ) J0_W      (INTEGER) : TPCORE REGION latitude  offset [# boxes]
!  (13 ) IM_W      (INTEGER) : TPCORE REGION longitude extent [# boxes]
!  (14 ) JM_W      (INTEGER) : TPCORE REGION latitude  extent [# boxes]
!  (15 ) IGZD      (INTEGER) : Variable equal to 1-I0_W or 1-J0_W
!
!  Diagnostics archived:
!  ============================================================================
!  (1  ) ND24 : Eastward flux of tracer in kg/s
!  (2  ) ND25 : Westward flux of tracer in kg/s
!  (3  ) ND26 : Upward   flux of tracer in kg/s
!
!  NOTES:
!  (1 ) Differences from "tpcore_mod.f" denoted by !%%% (bmy, 3/10/03)
!  (2 ) Now references TCVV & ITS_A_CH4_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (3 ) Remove code for the CO-OH simulation (bmy, 6/24/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,       ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE GLOBAL_CH4_MOD, ONLY : XNUMOL_CH4, TCH4
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE TRACER_MOD,     ONLY : TCVV, ITS_A_CH4_SIM

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! Diagnostic switches
#     include "CMN_GCTM"      ! g0_100

      ! Arguments
      INTEGER, INTENT(IN) :: IC, NDT, Jmax
      INTEGER, INTENT(IN) :: I0_W, J0_W, IM_W, JM_W, IGZD
      TYPE (XPLEX),  INTENT(IN) :: FX(     IM_W+1, JM_W,   LLPAR   )  
      TYPE (XPLEX),  INTENT(IN) :: FX1_TP( IM_W,   JM_W,   LLPAR   )
      TYPE (XPLEX),  INTENT(IN) :: FY(     IM_W,   JM_W+1, LLPAR   )  
      TYPE (XPLEX),  INTENT(IN) :: FY1_TP( IM_W,   JM_W,   LLPAR   )
      TYPE (XPLEX),  INTENT(IN) :: FZ(     IM_W,   JM_W,   LLPAR+1 )   
      TYPE (XPLEX),  INTENT(IN) :: FZ1_TP( IM_W,   JM_W,   LLPAR   )
      TYPE (XPLEX),  INTENT(IN) :: ACOSP(-10:Jmax)
       
      ! Local variables
      INTEGER             :: I, J, JREF, K, K2
      TYPE (XPLEX)              :: DTC, DTDYN

      ! Local SAVEd variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      TYPE (XPLEX),  SAVE       :: DXYP(JJPAR)

      !=================================================================
      ! DIAG_FLUX begins here!
      ! 
      ! FX, FX1_TP, FY, FY1_TP, FZ, FZ1_TP have units of [mb/timestep].
      ! 
      ! To get tracer fluxes in kg/s :
      ! * (100./9.8)                 => kg/m2
      ! * DXYP(J)/(DTDYN * TCVV(IC)) => kg/s
      ! 
      ! Direction of the fluxes : 
      ! ----------------------------------------------------------------
      ! FX(I,J,K) => flux coming into the west edge of the box I 
      !              (from I-1 to I).
      !           => a positive flux goes from west to east.
      ! 
      ! FY(I,J,K) => flux coming into the south edge of the box J 
      !              (from J to J-1).
      !           => a positive flux goes from south to north     
      !              (from J-1 to J)
      ! 
      ! FZ(I,J,K) => flux coming down into the box k.
      !           => a positive flux goes down.
      !=================================================================

      ! DTDYN = TYPE (XPLEX) value for NDT, the dynamic timestep
      DTDYN = XPLX( NDT )

      ! First-time initialization
      IF ( FIRST ) THEN
         
         ! Surface area [m2]
         DO J = 1, JJPAR
            DXYP(J) = GET_AREA_M2( J )
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! ND24 Diagnostic: Eastward flux of tracer in [kg/s]
      !=================================================================
      IF ( ND24 > 0 ) THEN

#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JREF, K, K2, DTC )
#endif
         DO K = 1, LLPAR

            ! K  is the vertical index down from the atmosphere top downwards
            ! K2 is the vertical index up from the surface
            K2 = LLPAR - K + 1
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%
            DO J = 1,JM_W 
               JREF = J +J0_W
               DO I = 1, IM_W 
                  DTC = ( FX(I,J,K) + FX1_TP(I,J,K) ) *
     &                  ( G0_100    * DXYP(JREF)    ) /
     &                  ( TCVV(IC)  * DTDYN         )

                  MASSFLEW(I+I0_W,J+J0_W,K2,IC) = 
     &                MASSFLEW(I+I0_W,J+J0_W,K2,IC) + DTC
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! ND25 Diagnostic: Northward flux of tracer in [kg/s]
      !=================================================================
      IF ( ND25 > 0 ) THEN

#if   defined( multitask ) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JREF, K, K2, DTC )
#endif
         DO K = 1, LLPAR

            ! K  is the vertical index down from the atmosphere top downwards
            ! K2 is the vertical index up from the surface
            K2 = LLPAR - K + 1
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%
             DO J = 1, JM_W
               JREF = J +J0_W
               DO I = 1, IM_W
                  DTC = ( FY(I,J,K) + FY1_TP(I,J,K)       ) *
     &                  ( ACOSP(J)  * G0_100 * DXYP(JREF) ) /
     &                  ( TCVV(IC)  * DTDYN               )

                  ! Contribution for CH4 run (bmy, 1/17/01)
                  IF ( ITS_A_CH4_SIM() ) THEN
                  TCH4(I+I0_W,J+J0_W,K,10)=TCH4(I+I0_W,J+J0_W,K,10)+ 
     &                                ( DTC * DTDYN * XNUMOL_CH4 )
                  ENDIF

                  MASSFLNS(I+I0_W,J+J0_W,K2,IC) =
     &                MASSFLNS(I+I0_W,J+J0_W,K2,IC) + DTC
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! ND26 Diagnostic : Upward flux of tracer in [kg/s] 
      !=================================================================
      IF ( ND26 > 0 ) THEN

#if   defined( multitask )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, JREF, K, K2, DTC )
#endif
         DO K = 1, LLPAR

            ! K  is the vertical index down from the atmosphere top downwards
            ! K2 is the vertical index up from the surface
            ! flux through top of the atmosphere is always zero, so don't need
            ! to archive it. While flux through the surface is not zero,
            ! so move the array index downward a layer
            ! DTC is the flux through the bottom of the box (yxw,02/09/2003)
            K2 = LLPAR - K + 1
!%%%
!%%% MODIFICATIONS FOR NESTED GRID (yxw, bmy, 3/10/03)
!%%% Rewrite DO-loop limits
!%%%
            DO J = 1, JM_W
               JREF = J +J0_W
               DO I = 1, IM_W
                  DTC = ( FZ(I,J,K) + FZ1_TP(I,J,K) ) *
     &                  ( G0_100    * DXYP(JREF)    ) /
     &                  ( TCVV(IC)  * DTDYN         )

!#IF   DEFINED( LGEOSCO )
! Not really the cross-tropopause flux - there is some horizontal
! transport across tropopause as well.
!                  IF ( K2 == LPAUSE(I,J) - 1 ) THEN
!                     TCO(I,J,1,10) = TCO(I,J,1,10) + 
!     &                               DTC * DTDYN * XNUMOL_CO
!                  ENDIF
!#endif
                  MASSFLUP(I+I0_W,J+J0_W,K2,IC) = 
     &               MASSFLUP(I+I0_W,J+J0_W,K2,IC) + DTC
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE DIAG_FLUX

!------------------------------------------------------------------------------

      SUBROUTINE POSITION_WINDOW( JS, JN, JL, JH, OUT, J1_IN, J2_IN )

C Determine the relationship of (JS0,JN0) and (J1_W, J2_W) (yxw, 08/23/01)
C (ji_in,j2_in) are the part of the window inside the region (JS0,JN0)

      INTEGER JS, JN,JL, JH, J1_IN, J2_IN
      LOGICAL OUT

      IF (JH .GT. JS) THEN
         IF (JL .GT. JN) THEN
            OUT=.TRUE.
         ELSE
            OUT=.FALSE.
            IF ((JH-JS) .GE. (JN-JS)) THEN
                J2_IN=JN
            ELSE
                J2_IN=JH
            ENDIF
            IF (JL .GE. JS) THEN
                J1_IN=JL
            ELSE
                J1_IN=JS
            ENDIF
         ENDIF
       ELSE
         OUT=.TRUE.
       ENDIF

       ! Return to TPCORE
       END SUBROUTINE POSITION_WINDOW

!------------------------------------------------------------------------------
      END MODULE TPCORE_WINDOW_MOD
