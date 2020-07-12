!-*- f90 -*- -
function dstack( A, M, N, R )

    implicit none

    integer :: i, j, k
    integer, intent(inout) :: M, N, R
    real(kind=8), intent(inout) :: A(M,N)
    real(kind=8), dimension(M,N,R) :: dstack

    do i = 1, M

        do j = 1, N

            do k = 1, R

                dstack(i,j,k) = A(i,j)

            end do

        end do

    end do

end function dstack