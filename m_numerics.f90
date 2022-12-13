module m_numerics
   use iso_fortran_env

contains
   SUBROUTINE LUDCMP(A, N, NP, INDX)
      IMPLICIT REAL*8(A - H, O - Z)
      PARAMETER(NMAX=100, TINY=1.D-10)
      DIMENSION A(NP, NP), INDX(N), VV(NMAX)
      D = 1.D0
      DO I = 1, N
         AAMAX = 0.D0
         DO J = 1, N
            IF (abs(A(I, J)) .GT. AAMAX) AAMAX = ABS(A(I, J))
         end do
         IF (AAMAX .EQ. 0.D0) then
            PRINT *, 'SINGULAR MATRIX.'
            stop
         end if
         VV(I) = 1.D0/AAMAX
      end do
      DO J = 1, N
         IF (J .GT. 1) THEN
            DO I = 1, J - 1
               SUM = A(I, J)
               IF (I .GT. 1) THEN
                  DO K = 1, I - 1
                     SUM = SUM - A(I, K)*A(K, J)
                  end do
                  A(I, J) = SUM
               END IF
            end do
         END IF
         AAMAX = 0.D0
         DO I = J, N
            SUM = A(I, J)
            IF (J .GT. 1) THEN
               DO K = 1, J - 1
                  SUM = SUM - A(I, K)*A(K, J)
               end do
               A(I, J) = SUM
            END IF
            DUM = VV(I)*ABS(SUM)
            IF (DUM .GE. AAMAX) THEN
               IMAX = I
               AAMAX = DUM
            END IF
         end do
         IF (J .NE. IMAX) THEN
            DO K = 1, N
               DUM = A(IMAX, K)
               A(IMAX, K) = A(J, K)
               A(J, K) = DUM
            end do
            D = -D
            VV(IMAX) = VV(J)
         END IF
         INDX(J) = IMAX
         IF (J .NE. N) THEN
            IF (A(J, J) .EQ. 0.D0) A(J, J) = TINY
            DUM = 1./A(J, J)
            DO I = J + 1, N
               A(I, J) = A(I, J)*DUM
            end do
         END IF
      end do
      IF (A(N, N) .EQ. 0.D0) A(N, N) = TINY
      RETURN
   END SUBROUTINE LUDCMP

   SUBROUTINE LUBKSB(A, N, NP, INDX, B)
      IMPLICIT REAL*8(A - H, O - Z)
      DIMENSION A(NP, NP), INDX(N), B(N)
      II = 0
      DO I = 1, N
         LL = INDX(I)
         SUM = B(LL)
         B(LL) = B(I)
         IF (II .NE. 0) THEN
            DO J = II, I - 1
               SUM = SUM - A(I, J)*B(J)
            end do
         ELSE IF (SUM .NE. 0.D0) THEN
            II = I
         END IF
         B(I) = SUM
      end do
      DO I = N, 1, -1
         SUM = B(I)
         IF (I .LT. N) THEN
            DO J = I + 1, N
               SUM = SUM - A(I, J)*B(J)
            end do
         END IF
         B(I) = SUM/A(I, I)
      end do
      RETURN
   END

end module m_numerics
