!-*- f90 -*- -
function product3( A, Q, M, N, R )

    implicit none

    integer :: i, j, k
    integer, intent(in) :: M, N, R
    real, intent(in) :: A(M,N)
    real, intent(in) :: Q(M,N,R)
    real :: product3(M,N,R)

    do i = 1, M

        do j = 1, N

            do k = 1, R

                product3(i,j,k) = Q(i,j,k) * A(i,j)

            end do

        end do

    end do

    end function