!-*- f90 -*- -
function Mm( M, p, q )

    implicit none
    integer, parameter :: dp=kind(0.d0)

    integer :: i, j
    integer, intent(in) :: p, q
    real(kind=8), intent(in) :: M(p,q)
    real(kind=8) :: Mm(p,q)

    do i = 1, p
        do j = 1, q
            Mm(i,j) = 0.5_dp * ( M(i,j)-abs(M(i,j)) )
        end do
    end do

end function


!-*- f90 -*- -
function Mp( M, p, q )

    implicit none
    integer, parameter :: dp=kind(0.d0)

    integer :: i, j
    integer, intent(in) :: p, q
    real(kind=8), intent(in) :: M(p,q)
    real(kind=8) :: Mp(p,q)

    do i = 1, p
        do j = 1, q
            Mp(i,j) = 0.5_dp * ( M(i,j)+abs(M(i,j)) )
        end do
    end do

end function