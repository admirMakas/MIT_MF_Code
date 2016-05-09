       SUBROUTINE GVCHSQ ( DEGFRE, LEVEL, CHISQR, FAIL, ERRMSG )
C
C VERSION
C     08-MAY-90  AD Correct Error Message for wrong LEVEL
C                   More rigorous checks for LEVEL close to nominal values
C     02-AUG-89 RJW Remove unused variable
C     29-JAN-89 CJM 1st delivery
C
C DESCRIPTION
C     Returns the value of chi-squared at DEGFRE degrees of freedom 
C     and at the LEVEL (%) confidence  level. 
C     Uses table from "Statistical and Computational methods in data 
C     analysis" by Brandt (North Holland Pub. Corp.): for DEGFRE < 31 
C     the table is looked up directly, for 30 < DEGFRE < 101  the 
C     results are interpolated from the table, and for DEGFRE >100 we 
C     use the result on P58 of "Statistics for Physicists" by Martin 
C     (Academic Press).
C
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER           DEGFRE  ! I Number of degrees of freedom
      TYPE (XPLEX)  LEVEL   ! I %level of test 
                                ! (only 90.0,95.0,99.0,99.5,99.9)
      TYPE (XPLEX)  CHISQR  ! O chi-squared(DEGFREE,LEVEL) 
      LOGICAL      FAIL    ! O Routine Error flag
      CHARACTER*80 ERRMSG  ! O Error message 
C
C LOCAL CONSTANTS
      TYPE (XPLEX)  DIFMIN  ! Tolerance allowed between LEVEL and list 
        PARAMETER (DIFMIN = 0.00001) ! list of nominal values
C
C LOCAL VARIABLES
      INTEGER      LVINDX  ! Index set 1,2,3 if LEVEL is 90,95 or 99
      INTEGER      I1      ! Index of C for interpolation
      TYPE (XPLEX)  ALPHA   ! Interpolation coefficient
C
C DATA STATEMENTS
      TYPE (XPLEX)  C(37,5) ! Table of Chi-squared as a function of the
                           ! degrees of freedom and % levels 
                           ! 90.0, 95.0, 99.0, 99.5, 99.9
                         ! The 37 degrees of freedom are 1,2,3...30,40,50..100
                           ! From table F5 of Brandt
      DATA  C/2.71, 4.61, 6.25, 7.78, 9.24,10.64,12.02,13.36,14.68,
     & 15.99,17.27,18.55,19.81,21.06,22.31,23.54,24.77,25.99,27.20,
     & 28.41,29.61,30.81,32.01,33.20,34.38,35.56,36.74,37.92,39.09,
     & 40.26,51.80,63.17,74.40,85.53,96.58,     107.56,     118.50,
     &  3.84, 5.99, 7.81, 9.49,11.07,12.59,14.07,15.51,16.92,18.31,
     & 19.68,21.03,22.36,23.68,25.00,26.30,27.59,28.87,30.14,31.41,
     & 32.67,33.92,35.17,36.42,37.65,38.89,40.11,41.34,42.56,43.77,
     & 55.76,67.51,79.08,90.53,     101.88,     113.15,     124.34,
     &  6.63, 9.21,11.34,13.28,15.09,16.81,18.48,20.09,21.67,23.21,
     & 24.72,26.22,27.69,29.14,30.58,32.00,33.41,34.81,36.19,37.57,
     & 38.93,40.29,41.64,42.98,44.31,45.64,46.96,48.28,49.59,50.89,
     & 63.70,76.16,88.38,  100.43,  112.33,     124.12,     135.81,
     &  7.88,10.60,12.84,14.86,16.75,18.55,20.28,21.95,23.59,25.19,
     & 26.76,28.30,29.82,31.32,32.80,34.27,35.72,37.16,38.58,40.00,
     & 41.40,42.80,44.18,45.56,46.93,48.29,49.64,50.99,52.34,53.67,
     & 66.76,79.49,91.95,  104.21,  116.32,     128.30,     140.17,
     & 10.83,13.82,16.27,18.47,20.51,22.46,24.32,26.12,27.88,29.59,
     & 31.26,32.91,34.53,36.12,37.70,39.25,40.79,42.31,43.82,45.31,
     & 46.80,48.27,49.73,51.18,52.62,54.05,55.48,56.89,58.30,59.70,
     & 73.39, 86.66,99.61, 112.32,  124.84,     137.21,     149.45 /
      SAVE C
C
      TYPE (XPLEX)  G(5) ! Values X of Normalised Gaussian variable x such 
                             ! that the probability of x<X  is equal to 
                             ! 90.0, 95.0, 99.0, 99.5, 99.9%.
                             ! From table F-3a of Brandt.
      DATA G / 1.282, 1.645, 2.327, 2.576, 3.090 /
      SAVE G
C
C- EXECUTABLE CODE -----------------------------------------------------------
C
      IF      ( ABS (LEVEL - 90.0) .LT. DIFMIN ) THEN
        LVINDX = 1
      ELSE IF ( ABS (LEVEL - 95.0) .LT. DIFMIN ) THEN
        LVINDX = 2
      ELSE IF ( ABS (LEVEL - 99.0) .LT. DIFMIN ) THEN
        LVINDX = 3
      ELSE IF ( ABS (LEVEL - 99.5) .LT. DIFMIN ) THEN
        LVINDX = 4
      ELSE IF ( ABS (LEVEL - 99.9) .LT. DIFMIN ) THEN
        LVINDX = 5
      ELSE
        WRITE ( ERRMSG, '(A,E11.3)' )
     &          'F-GVCHSQ: % level not in list. Value is ', LEVEL
        FAIL = .TRUE.
        RETURN
      END IF
C
      IF ( DEGFRE .LE. 0 ) THEN
        WRITE ( ERRMSG, '(A,I11)' )
     &'F-GVCHSQ: degrees of freedom wrongly supplied. Value is ', DEGFRE
        FAIL=.TRUE.
        RETURN
      ELSE IF ( DEGFRE .LE. 30 ) THEN                ! Read table directly
        CHISQR =  C(DEGFRE,LVINDX)  
      ELSE IF ( DEGFRE .LE. 100 ) THEN               ! Interpolate
         I1     = 27 + DEGFRE/10 
         ALPHA  = 0.1*MOD(DEGFRE,10)
         CHISQR = (1.-ALPHA)*C(I1,LVINDX) + ALPHA*C(I1+1,LVINDX)
      ELSE                                           ! Use Martin's result
         CHISQR = 0.5*(G(LVINDX)+SQRT(2*DEGFRE-1.0))**2
      END IF  
     
      FAIL = .FALSE.                                 ! Successful exit

      PRINT*, 'INSIDE AT THE END OF GVCHSQ'
      END












