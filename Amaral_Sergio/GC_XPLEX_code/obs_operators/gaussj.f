      SUBROUTINE gaussj( IA, NA, NB, A, B, SINGLR, FAIL, ERRMSG )
C
C VERSION
C     10-JAN-94  RJW  Remove non-ANSI common block WORKER
C     18-JAN-90  RJW  Increased IAMAX to 200
C     08-JUN-89  AD
C
C DESCRIPTION
C     General Purpose Matrix routine.
C     Invert Matrix A:  B = A**-1 
C     Using A as the output matrix is also possible: A = A**-1
C     This routine uses Gauss-Jordan Elimination with full pivoting, based on 
C     the routine GAUSSJ given on pp 28-29 of Numerical Recipes. 
C     This routine is no quicker than using MTXLEQ but saves having to construct
C     an Identity matrix.
C     There is probably scope for some more optmisation in this routine.
C     
C                                      Ref: Numerical Recipes (Press et al.)
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE
C
C ARGUMENTS
      INTEGER      IA       !  I  order of matrix A 
      INTEGER      NA       !  I  1st index of A as declared externally
      INTEGER      NB       !  I  1st index of B as declared externally
      TYPE (XPLEX)         A(NA,IA) !  I  
      TYPE (XPLEX)         B(NB,IA) !  I    Note: NA,NB >= IA is required
      LOGICAL      SINGLR   !  O  Set TRUE if A is singular 
      LOGICAL      FAIL     !  O  Set TRUE if a Fatal Error is detected
      CHARACTER*80 ERRMSG   !  O  Error message written if FAIL is TRUE
C
C LOCAL CONSTANTS
      INTEGER     IAMAX         ! Largest expected value of IA
      PARAMETER ( IAMAX = 3000 )
C
C LOCAL VARIABLES
      INTEGER  IPIV(IAMAX)      ! Workspace
      INTEGER  INDXR(IAMAX)
      INTEGER  INDXC(IAMAX)
      INTEGER  I, J, K          ! counters
      INTEGER  IROW, ICOL       ! store pivot row/column
      TYPE (XPLEX)   BIG, DUM, PIVINV
C
 1000 FORMAT ( A, I11, A, I11 )
C
C- EXECUTABLE CODE ------------------------------------------------------------
C

C Initialise B
      DO I = 1, NB
         DO J = 1, IA
            B(I,J) = 0.0
         ENDDO
      ENDDO 

C Check parameters
      IF ( IA .GT. IAMAX ) THEN                       ! Matrix A too large
        WRITE ( ERRMSG, 1000 )
     &          'F-MTXINV: IA =', IA, ' exceeds IAMAX =', IAMAX
        FAIL = .TRUE.
        RETURN                                        ! Exit with Fatal Error
      END IF
C
      IF ( IA .GT. NA ) THEN                          ! Inconsistent parameters
        WRITE ( ERRMSG, 1000 )
     &          'F-MTXINV: IA =', IA, ' exceeds NA =', NA 
        FAIL = .TRUE.
        RETURN                                         ! Exit with Fatal Error
      END IF
C
      FAIL   = .FALSE.        ! no further Fatal errors, only singular matrices
      SINGLR = .FALSE.        ! set TRUE if a Singular Matrix is detected
C
      DO I = 1, IA
        DO J = 1, IA
          B(J,I) = A(J,I)     ! copy matrix A to B (may not be necessary here?)
        END DO
        IPIV (I) = 0   
        IF (B(I,I) .EQ. 0.0) THEN 
           PRINT *, 'B EQ TO ZERO IN INPUT!!',I,B(I,I)
!           STOP 'B ZERO'
        ENDIF
      END DO

C
      DO I = 1, IA              ! Main loop over columns to be reduced
         BIG = 0.0
         DO J = 1, IA           ! Outer loop of search for a pivot element
            IF ( IPIV(J) .NE. 1 ) THEN
               DO K = 1, IA
                  IF ( IPIV(K) .EQ. 0 ) THEN
                     IF ( ABS(B(J,K)) .GE. BIG ) THEN
                        BIG  = ABS(B(J,K))
                        IROW = J
                        ICOL = K
                     END IF
                  ELSE IF ( IPIV(K) .GT. 1 ) THEN
                     print*,'SING 1!'
                     SINGLR = .TRUE.
                     RETURN     ! Exit for Singular Matrix
                  END IF
               END DO
            END IF
         END DO
         IPIV(ICOL) = IPIV(ICOL) + 1
C 
         IF ( IROW .NE. ICOL ) THEN 
C         got pivot element so interchange rows if reqd to put pivot on diagonal
            DO J = 1, IA
               DUM       = B(IROW,J)
               B(IROW,J) = B(ICOL,J)
               B(ICOL,J) = DUM
c     print *, 'irow,icol, j, b(irow,j),b(icol,j)',irow,icol, j,
c     &           b(irow,j),b(icol,j)
            END DO
         END IF
         INDXR(I) = IROW
         INDXC(I) = ICOL    
         IF ( B(ICOL,ICOL) .EQ. 0.0 ) THEN
            print*,'SING 2!'
            SINGLR = .TRUE.
            RETURN              ! Exit for Singular Matrix
         END IF


         PIVINV       = 1.0/B(ICOL,ICOL) ! Now divide pivot row by pivot element
         B(ICOL,ICOL) = 1.0
         DO J = 1, IA
            B(ICOL,J) = B(ICOL,J)*PIVINV
         END DO
         DO J = 1, IA           ! Reduce Rows except for pivot row
            IF ( J .NE. ICOL ) THEN
               DUM       = B(J,ICOL)
               B(J,ICOL) = 0.0
               DO K = 1, IA
                  B(J,K) = B(J,K)-B(ICOL,K)*DUM
               END DO
            END IF
         END DO
      END DO                    ! end of loop over column reduction
C     
      DO J = IA, 1, -1          ! unscramble column changes
         IF ( INDXR(J) .NE. INDXC(J) ) THEN
            DO K = 1, IA
               DUM           = B(K,INDXR(J))
               B(K,INDXR(J)) = B(K,INDXC(J))
               B(K,INDXC(J)) = DUM
            END DO
         END IF
      END DO
C     
      END

