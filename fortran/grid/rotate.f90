subroutine rotate(xx, yy, theta, M, N)

    implicit none
    integer, parameter :: dp=kind(0.d0)

    integer :: i, j
    real, intent(in) :: theta
    real(kind=8) :: ROT(2, 2)

    integer, intent(in) :: M, N

    real(kind=8), intent(inout) :: xx(M, N)
    real(kind=8), intent(inout) :: yy(M, N)
    real(kind=8) :: xy(2)

    ROT = reshape( (/ cos(theta), -sin(theta), sin(theta), cos(theta) /), (/2, 2/) )

    do i = 1, M

        do j = 1, N

            xy = matmul( ROT, (/ xx(i,j), yy(i,j) /) )
            xx(i,j) = xy(1)
            yy(i,j) = xy(2)

        end do

    end do

end subroutine