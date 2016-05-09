!Subroutine to find the inverse of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Reference : Algorithm has been well explained in:
!http://math.uww.edu/~mcfarlat/inverse.htm           
!http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
      SUBROUTINE FINDInv(matrix, inverse, n, errorflag)

      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

      !Declarations
      INTEGER, INTENT(IN)  :: n
      INTEGER, INTENT(OUT) :: errorflag !Return error status. 
                                !-1 for error, 0 for normal
      TYPE (XPLEX), INTENT(IN), DIMENSION(n,n)  :: matrix !Input matrix
      TYPE (XPLEX), INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
      
      LOGICAL :: FLAGr, FLAGc
      INTEGER :: i, j, k, l
      TYPE (XPLEX) :: m
      TYPE (XPLEX), DIMENSION(n,2*n) :: augmatrix !augmented matrix
	
      !Augment input matrix with an identity matrix
      DO i = 1, n
         DO j = 1, 2*n
            IF (j <= n ) THEN
               augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
               augmatrix(i,j) = 1
            ELSE
               augmatrix(i,j) = 0
            ENDIF
         END DO
      END DO	
      !Ensure diagonal elements are non-zero
      DO k = 1, n-1
         DO j = k+1,n
            IF (augmatrix(k,k) == 0) THEN
               DO i = k+1, n
                  IF (augmatrix(i,k) /= 0) THEN
                     DO  l = 1, 2* n
                        augmatrix(k,l) = augmatrix(k,l)+augmatrix(i,l)
                     END DO
                  ENDIF
               END DO
            ENDIF
         END DO
      END DO	
	
      !Augment input matrix with an identity matrix
      DO i = 1, n
         DO j = 1, 2*n
            IF (j <= n ) THEN
               augmatrix(i,j) = matrix(i,j)
            ELSE IF ((i+n) == j) THEN
               augmatrix(i,j) = 1
            ELSE
               augmatrix(i,j) = 0
            ENDIF
         END DO
      END DO	
      !Ensure diagonal elements are non-zero

      !Reduce augmented matrix to upper traingular form
      DO k =1, n-1
         DO j = k+1, n			
            m = augmatrix(j,k)/augmatrix(k,k)
            DO i = k, 2*n
               augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
            END DO
         END DO
      END DO
	
      !Test for invertibility
      DO i = 1, n
         IF (augmatrix(i,i) == 0) THEN
!!!            PRINT*, "Matrix is non - invertible"
            inverse = 0
            errorflag = -1
            return
         ENDIF
      END DO
	
      !Make diagonal elements as 1
      DO i = 1 , n
         m = augmatrix(i,i)
         DO j = i , (2 * n)				
            augmatrix(i,j) = (augmatrix(i,j) / m)
         END DO
      END DO
	
      !Reduced right side half of augmented matrix to identity matrix
      DO k = n-1, 1, -1
         DO i =1, k
            m = augmatrix(i,k+1)
            DO j = k, (2*n)
               augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
            END DO
         END DO
      END DO				
	
      !store answer
      DO i =1, n
         DO j = 1, n
            inverse(i,j) = augmatrix(i,j+n)
         END DO
      END DO
      errorflag = 0

      END SUBROUTINE FINDinv
