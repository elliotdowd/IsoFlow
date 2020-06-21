!-*- f90 -*- -
! function matmul22(A, B)

!     implicit none

!     real(kind=8) :: res(2, 2)
!     real(kind=8) :: A(2, 2), B(2, 2)

!     res = matmul(A, B)

!     return res

!     end function



subroutine M22INV(A)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(2,2), INTENT(INOUT)  :: A
    !DOUBLE PRECISION, DIMENSION(2,2)  :: M22INV

    DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
    DOUBLE PRECISION :: DET
    DOUBLE PRECISION, DIMENSION(2,2) :: COFACTOR

    DET = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    COFACTOR(1,1) = +A(2,2)
    COFACTOR(1,2) = -A(2,1)
    COFACTOR(2,1) = -A(1,2)
    COFACTOR(2,2) = +A(1,1)

    A = TRANSPOSE(COFACTOR) / DET

    END SUBROUTINE