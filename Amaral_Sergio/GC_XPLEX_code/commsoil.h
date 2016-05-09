! $Id: commsoil.h,v 1.2 2012/03/01 22:00:26 daven Exp $
!
!**********************************************************************
!                                                                     *
!  HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP                       *
!  MODULE FOR SOIL NOx EMISSIONS                                      *
!  by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob               *
!  (Release V2.1)                                                     *
!                                                                     *
!  Contact person: Bob Yantosca (bmy@io.harvard.edu)                  *
!                                                                     *
!**********************************************************************
! NOTES:
! (1 ) Be sure to force TYPE (XPLEX) with the DBLE function            
!       and the "D" exponent, wherever necessary (bmy, 10/6/99)             
! (2 ) Changed RCS ID tag comment character from "C" to "!" to allow 
!       freeform compilation.  Also added & continuation characters in 
!       column 73 to allow header files to be included in F90 freeform 
!       files. Updated comments, cosmetic changes. (bmy, 6/25/02)
! (3 ) Now use cpp switches to define 1x1 parameters.  Also added
!       space in the #ifdef block for the 1x125 grid (bmy, 12/1/04)
! (4 ) Bug fix: 2681 should be 2861 in NLAND (bmy, 9/22/06)
! (5 ) Set # of land boxes for GEOS-5 nested grids (yxw, dan, bmy, 11/6/08)
!**********************************************************************
!
! header file for soil NOx emissions
      
      ! The defined soil types
      INTEGER, PARAMETER :: NSOIL = 11

      ! Number of soil pulsing types
      INTEGER, PARAMETER :: NPULSE = 3  

#if   defined( GRID4x5  )

      ! There are 1118 land boxes for the 4 x 5 GLOBAL GRID
      INTEGER, PARAMETER :: NLAND = 1118

#elif defined( GRID2x25 )

      ! There are 3920 land boxes for the 2 x 2.5 GLOBAL GRID
      INTEGER, PARAMETER :: NLAND = 3920 

#elif defined( GRID1x125 )

      !%%% NOTE: still to be determined
      INTEGER, PARAMETER :: NLAND = 9999 

#elif defined( GRID1x1 ) && defined( NESTED_CH )

      ! There are 2861 land points for the 1x1 CHINA nested grid
      INTEGER, PARAMETER :: NLAND = 2861

#elif defined( GRID1x1 ) && defined( NESTED_NA )

      ! There are 2118  land points for the 1x1 N. AMERICA nested grid
      INTEGER, PARAMETER :: NLAND = 2118

#elif defined( GRID1x1 )

      ! There are 17174 land points for the 1x1 GLOBAL grid
      INTEGER, PARAMETER :: NLAND=17174

#elif defined( GRID05x0666 ) && defined( NESTED_CH ) && !defined( NESTED_SD )

      ! There are 8261 land points for the 0.5 x 0.666 CHINA nested grid
      INTEGER, PARAMETER :: NLAND = 8261
      
#elif defined( GRID05x0666 ) && defined( NESTED_NA ) && !defined( NESTED_SD )

      INTEGER, PARAMETER :: NLAND = 8568
      
#elif defined( GRID05x0666 ) && (defined( NESTED_NA ) || defined( NESTED_CH )) && defined( NESTED_SD )
       
      ! Parameters for smaller domain
      INTEGER, PARAMETER :: NLAND = 5153  

#endif


! water/desert/ice//Trop. Rain. Forst.//conifers//dry deciduous//
! other deciduous//woodland//grassland//agriculture (other than rice)
! rice paddies//wetland/tundra
      INTEGER INDEXSOIL(2,NLAND) !i,j of the grid
      TYPE (XPLEX) SOILPULS(NPULSE+1,NLAND) 
	 !tracking of wet/dry & three types of pulsing (Y&L, 94)
      TYPE (XPLEX) SOILPREP(2,NLAND)  !two month observed precip
      TYPE (XPLEX) SOILFERT(NLAND)    !ferterlizers
      TYPE (XPLEX) PULSFACT(NPULSE)   !pulsing factors
      TYPE (XPLEX) PULSDECAY(NPULSE)  !pulsing decay per timestep
      TYPE (XPLEX) SOILNOX(IGLOB,JGLOB) !stores output

      INTEGER NCONSOIL(NVEGTYPE) !olson->soil type,nvegtype in commbio.h
      TYPE (XPLEX) CANOPYNOX(MAXIJ,NTYPE) !track NOx within canopy dry dep.
      TYPE (XPLEX) SOILTA(NSOIL),SOILTB(NSOIL),SOILAW(NSOIL),
     & SOILAD(NSOIL)
      TYPE (XPLEX) SOILEXC(NSOIL)     !canopy wind extinction coeff.

      COMMON /SOIL/ SOILNOX,  INDEXSOIL, NCONSOIL, SOILPULS,            &
     &              SOILPREP, SOILFERT,  CANOPYNOX

      ! The correct sequence of PULSFACT is 5, 10, 15 
      DATA PULSFACT%r  / 5.D0,    10.D0,   15.D0   /

      ! PULSDECAY now contains the correct decay factors from Yienger & Levy
      DATA PULSDECAY%r / 0.805D0, 0.384D0, 0.208D0 / 

      ! SOILTA = Coefficient used to convert from surface temperture to  
      !          soil temperature     
      DATA SOILTA%r /0.D0,   0.84D0, 0.84D0, 0.84D0, 0.84D0,              &
     &             0.66D0, 0.66D0, 1.03D0, 1.03D0, 0.92D0,              &
     &             0.66D0/

      ! SOILTB = Coefficient used to convert from surface temperture to  
      !          soil temperature   
      DATA SOILTB%r /0.D0,   3.6D0,  3.6D0,  3.6D0,  3.6D0,               &
     &             8.8D0,  8.8D0,  2.9D0,  2.9D0,  4.4D0,               &
     &             8.8D0/
   
      ! SOILAW = Wet biome coefficient   
      DATA SOILAW%r /0.D0,   2.6D0,  0.03D0, 0.06D0, 0.03D0,              &
     &             0.17D0, 0.36D0, 0.36D0, 0.36D0, 0.003D0,             &
     &             0.05D0/

      ! SOILAD = Dry biome coefficient  
      DATA SOILAD%r /0.D0,   8.6D0,  0.22D0, 0.40D0, 0.22D0,              &
     &             1.44D0, 2.65D0, 2.65D0, 2.65D0, 0.003D0,             &
     &             0.37D0/

      ! SOILEXC = Canopy wind extinction coeff.  
      DATA SOILEXC%r /0.1D0, 4.D0,   4.D0,   4.D0,   4.D0,                &
     &              2.D0,  1.D0,   2.D0,   2.D0,   0.5D0,               &
     &              0.1D0/

