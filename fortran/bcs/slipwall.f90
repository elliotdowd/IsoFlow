
!-*- f90 -*- -
subroutine slip(u0, v0, u1, v1, s_proj, M)
    implicit none
    integer, parameter:: dp=kind(0.d0)

    integer, intent(in) :: M
    integer :: i

    real(kind=8), intent(inout) :: u0(M)
    real(kind=8), intent(inout) :: v0(M)
    real(kind=8), intent(in) :: u1(M)
    real(kind=8), intent(in) :: v1(M)
    real(kind=8), intent(in) :: s_proj(M, 2, 6)

    real(kind=8) :: matL(2, 2), matLinv(2, 2), matR(2, 2)
    real(kind=8) :: det(2, 2)
    real(kind=8) :: u0v0(2, M)
    real(kind=8) :: u1v1(2)

    do i = 1, M

        matL = reshape( (/ s_proj(i,1,4), -s_proj(i,1,3), &
                & s_proj(i,1,3), s_proj(i,1,4) /), (/2, 2/) )
        matR = reshape( (/ s_proj(i,2,4), -s_proj(i,2,3), &
                & -s_proj(i,2,3), -s_proj(i,2,4) /), (/2, 2/) )

        u1v1 = (/ u1(i), v1(i) /)

        det = 1.0_dp / (matL(1,1)*matL(2,2) - matL(1,2)*matL(2,1))
        matLinv = det * reshape( (/matL(2,2), -matL(1,2), matL(2,1), matL(1,1) /), (/2, 2/) )

        u0v0(:, i) = matmul( matmul(matLinv, matR), u1v1 )

        ! return velocity at halo cells 
        u0(i) = u0v0(1, i)
        v0(i) = u0v0(2, i)

    end do


    contains
        ! 2x2 matrix inversion
        FUNCTION M22INV(A)

        IMPLICIT NONE

        DOUBLE PRECISION, DIMENSION(2,2), INTENT(IN)  :: A
        DOUBLE PRECISION, DIMENSION(2,2)  :: M22INV

        DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
        DOUBLE PRECISION :: DET
        DOUBLE PRECISION, DIMENSION(2,2) :: COFACTOR

        DET =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

        COFACTOR(1,1) = +A(2,2)
        COFACTOR(1,2) = -A(2,1)
        COFACTOR(2,1) = -A(1,2)
        COFACTOR(2,2) = +A(1,1)

        M22INV = TRANSPOSE(COFACTOR) / DET

        END FUNCTION


end subroutine