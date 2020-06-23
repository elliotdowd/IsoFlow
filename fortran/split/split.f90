!-*- f90 -*- -
function Mm( M )

    implicit none
    integer, parameter :: dp=kind(0.d0)

    real(kind=8), intent(in) :: M(:,:)
    real(kind=8) :: Mm(:,:)

    Mm = 0.5_dp * ( M-abs(M) )

end function


!-*- f90 -*- -
function Mp( M )

    implicit none
    integer, parameter :: dp=kind(0.d0)

    real(kind=8), intent(in) :: M(:,:)
    real(kind=8) :: Mm(:,:)

    Mm = 0.5_dp * ( M+abs(M) )

end function