!-*- f90 -*- -
subroutine residual(res, dt, E_left, E_right, F_bot, F_top, s_proj, M, N)
! calculate state vector residual from flux components on each cell face

    implicit none

    integer i, j, k
    integer, parameter :: R = 4

    integer, intent(in) :: M, N

    real(kind=8), intent(inout) :: res(M,N,4)
    real(kind=8), intent(in) :: dt(M,N)
    real(kind=8), intent(in) :: E_left(M,N,R), E_right(M,N,R), F_bot(M,N,R), F_top(M,N,R)
    real(kind=8), intent(in) :: s_proj(M,N,6)

    do i = 1, M

        do j = 1, N

            do k = 1, R

                res(i,j,k) = -dt(i,j) * ( (E_right(i,j,k)*s_proj(i,j,5) - E_left(i,j,k) * s_proj(i,j,5)) &
                                        & + (F_top(i,j,k)*s_proj(i,j,6) - F_bot(i,j,k)*s_proj(i,j,6)) )
                
            end do

        end do

    end do

end subroutine


!-*- f90 -*- -
subroutine face_flux( mdot_half_zeta, mdot_half_eta, phi, p_zeta, p_eta, E_hat_left, E_hat_right, F_hat_bot, F_hat_top, M, N )
! calculate cell face fluxes from mass flow rate, phi, and p state vectors

    implicit none
    integer, parameter:: dp=kind(0.d0)

    integer :: i, j, k
    integer, parameter :: R = 4
    integer, intent(in) :: M, N

    real(kind=8), intent(in) :: mdot_half_zeta(M+1, N+2), p_zeta(M+1, N+2, 4)
    real(kind=8), intent(in) :: mdot_half_eta(M+2, N+1), p_eta(M+2, N+1, 4)
    real(kind=8), intent(in) :: phi(M+2,N+2,4)
    real(kind=8), intent(inout) :: E_hat_left(M,N,4), E_hat_right(M,N,4)
    real(kind=8), intent(inout) :: F_hat_bot(M,N,4), F_hat_top(M,N,4)

    ! loop through zeta computational direction
    do i = 1, M

        do j = 2, N+1

            do k = 1, R

                E_hat_left(i,j-1,k) = (0.5_dp) * mdot_half_zeta(i,j) * ( Phi(i,j,k) + Phi(i+1,j,k) ) &
                                    & -(0.5_dp) * abs(mdot_half_zeta(i,j)) * ( Phi(i+1,j,k) - Phi(i,j,k) ) &
                                    & + P_zeta(i,j,k)
                E_hat_right(i,j-1,k) = (0.5_dp) * mdot_half_zeta(i+1,j) * ( Phi(i+1,j,k) + Phi(i+2,j,k) ) &
                                    & -(0.5_dp) * abs(mdot_half_zeta(i+1,j)) * ( Phi(i+2,j,k) - Phi(i+1,j,k) ) &
                                    & + P_zeta(i+1,j,k)

            end do

        end do

    end do

    ! loop through eta computational direction
    do i = 2, M

        do j = 1, N

            do k = 1, R

                F_hat_bot(i-1,j,k) = (0.5_dp) * mdot_half_eta(i,j) * ( Phi(i,j,k) + Phi(i,j+1,k) ) &
                                  & -(0.5_dp) * abs(mdot_half_eta(i,j)) * ( Phi(i,j+1,k) - Phi(i,j,k) ) &
                                  & + P_eta(i,j,k)
                F_hat_top(i-1,j,k) = (0.5_dp) * mdot_half_eta(i,j+1) * ( Phi(i,j+1,k) + Phi(i,j+2,k) ) &
                                  & -(0.5_dp) * abs(mdot_half_eta(i,j+1)) * ( Phi(i,j+2,k) - Phi(i,j+1,k) ) &
                                  & + P_eta(i,j+1,k)

            end do

        end do

    end do

end subroutine


!-*- f90 -*- -
subroutine face_flux_muscl( mdot_half_zeta, mdot_half_eta, phiL, phiR, phiD, phiU, p_zeta, p_eta, E_hat_left, E_hat_right, &
                        &   F_hat_bot, F_hat_top, M, N )
! calculate cell face fluxes from mass flow rate, phi, and p state vectors using MUSCL interpolation

    implicit none
    integer, parameter:: dp=kind(0.d0)

    integer :: i, j, k
    integer, parameter :: R = 4
    integer, intent(in) :: M, N

    real(kind=8), intent(in) :: mdot_half_zeta(M+1, N+2), p_zeta(M+1, N+2, 4)
    real(kind=8), intent(in) :: mdot_half_eta(M+2, N+1), p_eta(M+2, N+1, 4)
    real(kind=8), intent(in) :: phiL(M+1,N+2,4), phiR(M+1,N+2,4)
    real(kind=8), intent(in) :: phiD(M+2,N+1,4), phiU(M+2,N+1,4)
    real(kind=8), intent(inout) :: E_hat_left(M,N,4), E_hat_right(M,N,4)
    real(kind=8), intent(inout) :: F_hat_bot(M,N,4), F_hat_top(M,N,4)

    ! loop through zeta computational direction
    do i = 1, M

        do j = 2, N+1

            do k = 1, R

                E_hat_left(i,j-1,k) = (0.5_dp) * mdot_half_zeta(i,j) * ( PhiL(i,j,k) + PhiR(i,j,k) ) &
                                    & -(0.5_dp) * abs(mdot_half_zeta(i,j)) * ( PhiR(i,j,k) - PhiL(i,j,k) ) &
                                    & + P_zeta(i,j,k)
                E_hat_right(i,j-1,k) = (0.5_dp) * mdot_half_zeta(i+1,j) * ( PhiL(i+1,j,k) + PhiR(i+1,j,k) ) &
                                    & -(0.5_dp) * abs(mdot_half_zeta(i+1,j)) * ( PhiR(i+1,j,k) - PhiL(i+1,j,k) ) &
                                    & + P_zeta(i+1,j,k)

            end do

        end do

    end do

    ! loop through eta computational direction
    do i = 2, M

        do j = 1, N

            do k = 1, R

                F_hat_bot(i-1,j,k) = (0.5_dp) * mdot_half_eta(i,j) * ( PhiD(i,j,k) + PhiU(i,j,k) ) &
                                  & -(0.5_dp) * abs(mdot_half_eta(i,j)) * ( PhiU(i,j,k) - PhiD(i,j,k) ) &
                                  & + P_eta(i,j,k)
                F_hat_top(i-1,j,k) = (0.5_dp) * mdot_half_eta(i,j+1) * ( PhiD(i,j+1,k) + PhiU(i,j+1,k) ) &
                                  & -(0.5_dp) * abs(mdot_half_eta(i,j+1)) * ( PhiU(i,j+1,k) - PhiD(i,j+1,k) ) &
                                  & + P_eta(i,j+1,k)

            end do

        end do

    end do

end