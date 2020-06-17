
!-*- f90 -*- -
subroutine calc_covariant(s_proj, u, v, U, V)
! calculate covariant velocity for a given cell
! input projected cell face areas [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y, S_zeta, S_eta] and velocity components

    implicit none
    integer, parameter:: dp=kind(0.d0)

    real(kind=8), intent(in) :: s_proj(:,:,6)
    real(kind=8), intent(in) :: u(:,:), v(:,:)
    real(kind=8), intent(inout) :: U(:,:)
    real(kind=8), intent(inout) :: V(:,:)

    U(:,:) = (1.0_dp/s_proj(:,:,5)) * &
           &    (u.*s_proj(:,:,1) + v*cells.projFaceArea(:,:,2));
       
    V(:,:) = (1.0_dp/s_proj(:,:,6)) * &
           &    (u.*s_proj(:,:,3) + v*s_proj(:,:,4));

end