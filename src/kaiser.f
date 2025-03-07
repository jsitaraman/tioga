c$$$
c$$$ Copyright TIOGA Developers. See COPYRIGHT file for details.
c$$$
c$$$ SPDX-License-Identifier: (BSD 3-Clause)


      subroutine kaiser_wrap(a,nrows,n,eigenv,trace,sume,ier)
      integer, intent(inout) :: nrows,n
      real*8, intent(inout):: a(nrows*n)
      real*8, intent(inout) :: eigenv(n)
      real*8, intent(inout) :: trace,sume
      integer, intent(inout) :: ier
      call kaiser(a,nrows,n,eigenv,trace,sume,ierr)
      end subroutine kaiser_wrap

      SUBROUTINE kaiser(a, nrows, n, eigenv, trace, sume, ier)
      
      !  EIGENVALUES AND VECTORS OF A SYMMETRIC +VE DEFINITE MATRIX,
      !  USING KAISER'S METHOD.
      !  REFERENCE: KAISER,H.F. 'THE JK METHOD: A PROCEDURE FOR FINDING THE
      !  EIGENVALUES OF A REAL*8 SYMMETRIC MATRIX', COMPUT.J., VOL.15, 271-273, 1972.
      
      !  ARGUMENTS:-
      !  A       = INPUT, AN ARRAY CONTAINING THE MATRIX
      !            OUTPUT, THE COLUMNS OF A CONTAIN THE NORMALIZED EIGENVECTORS
      !            OF A.   N.B. A IS OVERWRITTEN !
      !  NROWS   = INPUT, THE FIRST DIMENSION OF A IN THE CALLING PROGRAM.
      !  N       = INPUT, THE ORDER OF A, I.E. NO. OF COLUMNS.
      !            N MUST BE <= NROWS.
      !  EIGENV()= OUTPUT, A VECTOR CONTAINING THE ORDERED EIGENVALUES.
      !  TRACE   = OUTPUT, THE TRACE OF THE INPUT MATRIX.
      !  SUME    = OUTPUT, THE SUM OF THE EIGENVALUES COMPUTED.
      !            N.B. ANY SYMMETRIC MATRIX MAY BE INPUT, BUT IF IT IS NOT +VE
      !            DEFINITE, THE ABSOLUTE VALUES OF THE EIGENVALUES WILL BE FOUND.
      !            IF TRACE = SUME, THEN ALL OF THE EIGENVALUES ARE POSITIVE
      !            OR ZERO.   IF SUME > TRACE, THE DIFFERENCE IS TWICE THE SUM OF
      !            THE EIGENVALUES WHICH HAVE BEEN GIVEN THE WRONG SIGNS !
      !  IER     = OUTPUT, ERROR INDICATOR
      !             = 0 NO ERROR
      !             = 1 N > NROWS OR N < 1
      !             = 2 FAILED TO CONVERGE IN 10 ITERATIONS
      
      !  LATEST REVISION - 6 September 1990
      !  Fortran 90 version - 20 November 1998
      
      !*************************************************************************
      
      IMPLICIT NONE
      
      integer nrows,n
      real*8 a(nrows,n)
      real*8 eigenv(n)
      real*8 trace,sume
      integer ier
      
      
      ! Local variables
      
      REAL*8 , PARAMETER :: small = 1.0e-12, zero = 0.0
      REAL*8 , PARAMETER :: half = 0.5, one = 1.0
      INTEGER   :: i, iter, j, k, ncount, nn
      REAL*8  :: absp, absq, COS, ctn, eps, halfp, p, q, SIN, ss
      REAL*8  :: TAN, temp, xj, xk
      
      !   CALCULATE CONVERGENCE TOLERANCE, EPS.
      !   CALCULATE TRACE.   INITIAL SETTINGS.
      
      ier = 1
      IF(n < 1 .OR. n > nrows) RETURN
      ier = 0
      iter = 0
      trace = zero
      ss = zero
      DO j = 1,n
        trace = trace + a(j,j)
        DO i = 1,n
          ss = ss + a(i,j)**2
        END DO
      END DO
      sume = zero
      eps = small*ss/n
      nn = n*(n-1)/2
      ncount = nn
      
      !   ORTHOGONALIZE PAIRS OF COLUMNS J & K, K > J.
      
20    DO j = 1,n-1
        DO k = j+1,n
          
      !   CALCULATE PLANAR ROTATION REQUIRED
          
          halfp = zero
          q = zero
          DO i = 1,n
            xj = a(i,j)
            xk = a(i,k)
            halfp = halfp + xj*xk
            q = q + (xj+xk) * (xj-xk)
          END DO
          p = halfp + halfp
          absp = ABS(p)
          
      !   If P is very small, the vectors are almost orthogonal.
      !   Skip the rotation if Q >= 0 (correct ordering).
          
          IF (absp < eps .AND. q >= zero) THEN
            ncount = ncount - 1
            IF (ncount <= 0) GO TO 160
            CYCLE
          END IF
          
      !   Rotation needed.
          
          absq = ABS(q)
          IF(absp <= absq) THEN
            TAN = absp/absq
            COS = one/SQRT(one + TAN*TAN)
            SIN = TAN*COS
          ELSE
            ctn = absq/absp
            SIN = one/SQRT(one + ctn*ctn)
            COS = ctn*SIN
          END IF
          COS = SQRT((one + COS)*half)
          SIN = SIN/(COS + COS)
          IF(q < zero) THEN
            temp = COS
            COS = SIN
            SIN = temp
          END IF
          IF(p < zero) SIN = -SIN
          
      !   PERFORM ROTATION
          
          DO i = 1,n
            temp = a(i,j)
            a(i,j) = temp*COS + a(i,k)*SIN
            a(i,k) = -temp*SIN + a(i,k)*COS
          END DO
        END DO
      END DO
      ncount = nn
      iter = iter + 1
      IF(iter < 10) GO TO 20
      ier = 2
      
      !   CONVERGED, OR GAVE UP AFTER 10 ITERATIONS
      
160   DO j = 1,n
        temp = SUM( a(1:n,j)**2 )
        eigenv(j) = SQRT(temp)
        sume = sume + eigenv(j)
      END DO
      
      !   SCALE COLUMNS TO HAVE UNIT LENGTH
      
      DO j = 1,n
        IF (eigenv(j) > zero) THEN
          temp = one/eigenv(j)
        ELSE
          temp = zero
        END IF
        a(1:n,j) = a(1:n,j)*temp
      END DO
      
      RETURN
      END SUBROUTINE kaiser
     
     




























































































































































