      INTEGER NTRAC             ! Number of tracers
      INTEGER NX                ! State space
      INTEGER NY                ! Measurment space
      INTEGER NYY	
      INTEGER IIGR              ! Gridsize of data
      INTEGER JJGR              ! Gridsize of data
      INTEGER LLGR              ! Gridsize of data
      INTEGER II1               ! start index for i-grid
      INTEGER II2               ! end index for i-grid
      INTEGER JJ1               ! start index for j-grid
      INTEGER JJ2               ! end index for j-grid
      INTEGER NAS               ! number of pts in nmc asia grid 

      TYPE (XPLEX) LON_lim1		! min longitude box
      TYPE (XPLEX) LON_lim2		! max longitude box	
      TYPE (XPLEX) LAT_lim1		! min latitude box
      TYPE (XPLEX) LAT_lim2 		! max latitude box


      PARAMETER ( NTRAC = 1   )  ! Tagged tracers (+1 for total CO)
!      PARAMETER ( NX    = 8    )
      PARAMETER ( NX    = 14    )
      PARAMETER ( NY    = 2500 )  ! MAX due to svdcmp.f
      PARAMETER ( NYY   = 5000 )  ! MAX due to svdcmp.f

! Window #1 : 
! lon [180W ; 15W]
! lat [ 0N ; 76N]
!      PARAMETER ( IIGR  = 67 ) ! 67 = [180W - 15W]
!      PARAMETER ( JJGR  = 39 ) ! 34 = [0N - 76N]
!      PARAMETER ( LLGR  = 1 )
!      PARAMETER ( II1 = 1 )    !  0 = 180W
!      PARAMETER ( II2 = 67 )   ! 66 = 15W
!      PARAMETER ( JJ1 = 47 )   ! 50 = 0N
!      PARAMETER ( JJ2 = 85 )   ! 83 = 76N
!      PARAMETER ( NAS = 2613 )	! IIGR*JJGR

! Window #2 : [reduced window]
! lon [167.5W ; 45W]
! lat [ 16N ; 70N]
!      PARAMETER ( IIGR  = 50 ) 
!      PARAMETER ( JJGR  = 28 ) 
!      PARAMETER ( LLGR  = 1 )
!      PARAMETER ( II1 = 6 )
!      PARAMETER ( II2 = 55 )
!      PARAMETER ( JJ1 = 54 )
!      PARAMETER ( JJ2 = 81 )
!      PARAMETER ( NAS = 1400 )	! IIGR*JJGR

! Window #3 : [reduced window for Alaska and Canada]
! lon [167.5W ; 45W]
! lat [ 30N ; 70N]
      PARAMETER ( IIGR  = 50 )
      PARAMETER ( JJGR  = 21 ) 
      PARAMETER ( LLGR  = 1 )
      PARAMETER ( II1 = 6 )
      PARAMETER ( II2 = 55 )
      PARAMETER ( JJ1 = 61 )
      PARAMETER ( JJ2 = 81 )
      PARAMETER ( NAS = 1050 )	! IIGR*JJGR
	

! Read MOPITT data and model outputs directly: (timedep inversion)
	PARAMETER (LAT_Lim1 = 30.)	
	PARAMETER (LAT_Lim2 = 75.)
	PARAMETER (LON_Lim1 = -167.5)	
	PARAMETER (LON_Lim2 = -30.)

