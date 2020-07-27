!-*- f90 -*- -
subroutine c_entropy_fixed( c_zeta, c_eta, ht, gam, U, V, M, N )

    implicit none
    integer, parameter:: dp=kind(0.d0)

    integer :: i, j
    integer, intent(in) :: M, N
    real(kind=8) :: c_st(M+2,N+2)
    real(kind=8) :: c_l(M+1,N+2), c_r(M+1,N+2), c_b(M+2,N+1), c_t(M+2,N+1)
    real(kind=8), intent(in) :: ht(M+2,N+2), gam(M+2,N+2), U(M+2,N+2), V(M+2,N+2)
    real(kind=8), intent(inout) :: c_zeta(M+1,N+2), c_eta(M+2,N+1)

    c_st = 2.0_dp * ht * ( (gam-1) / (gam+1) )

    c_l = c_st(1:M+1,:) / max( sqrt( c_st(1:M+1,:) ),  U(1:M+1,:) )
    c_r = c_st(2:M+2,:) / max( sqrt( c_st(2:M+2,:) ), -U(2:M+2,:) )
    c_b = c_st(:,1:N+1) / max( sqrt( c_st(:,1:N+1) ),  V(:,1:N+1) )
    c_t = c_st(:,2:N+2) / max( sqrt( c_st(:,2:N+2) ), -V(:,2:N+2) )

    c_zeta = min( c_l, c_r )
    c_eta  = min( c_b, c_t )

end subroutine


!-*- f90 -*- -
subroutine residual( res, dt, E_left, E_right, F_bot, F_top, s_proj, M, N )
! calculate state vector residual from flux components on each cell face

    implicit none

    integer i, j, k
    integer, parameter :: R = 4

    integer, intent(in) :: M, N

    real(kind=8), intent(inout) :: res(M,N,4)
    real(kind=8), intent(in) :: dt(M,N)
    real(kind=8), intent(in) :: E_left(M,N,R), E_right(M,N,R), F_bot(M,N,R), F_top(M,N,R)
    real(kind=8), intent(in) :: s_proj(M+2,N+2,6)

    do i = 1, M

        do j = 1, N

            do k = 1, R

                res(i,j,k) = -dt(i,j) * ( (E_right(i,j,k)*s_proj(i+1,j+1,5) - E_left(i,j,k) * s_proj(i,j+1,5)) &
                                      & + (F_top(i,j,k)*s_proj(i+1,j+1,6) - F_bot(i,j,k)*s_proj(i+1,j,6)) )
                
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


!-*- f90 -*- -
subroutine roe_fbar( fbar, Phi, U, normal, p, M, N )

    implicit none
    integer, parameter:: dp=kind(0.d0)

    integer :: k

    integer, intent(in) :: M, N
    real(kind=8), intent(in) :: U(M,N), p(M,N)
    real(kind=8), intent(in) :: Phi(M,N,4), normal(M,N,4)
    real(kind=8), intent(inout) :: fbar(M,N,4)

    do k = 1, 4

        fbar(:,:,k) = Phi(:,:,k) * U + normal(:,:,k) * p

    end do

end


!-*- f90 -*- -
! calculate Fd term for Roe FV formulation
subroutine roe_dissipation( Ed_half, Fd_half, p, phi, Psi_zeta, Psi_eta, d_zeta, d_eta, Pu_zeta, Pp_zeta, Uu_zeta, Up_zeta, &
                                                                &    Pu_eta, Pp_eta, Uu_eta, Up_eta, M, N )

    implicit none
    integer, parameter:: dp=kind(0.d0)

    integer :: k

    integer, intent(in) :: M, N
    real(kind=8), intent(in) :: d_zeta(M+2,N+2), Pu_zeta(M+1,N+2), Pp_zeta(M+1,N+2), &
                                            &    Uu_zeta(M+1,N+2), Up_zeta(M+1,N+2)
    real(kind=8), intent(in) :: d_eta(M+2,N+2), Pu_eta(M+2,N+1), Pp_eta(M+2,N+1), &
                                            &   Uu_eta(M+2,N+1), Up_eta(M+2,N+1)
    real(kind=8), intent(in) :: Psi_zeta(M+1,N+2,4), Psi_eta(M+2,N+1,4)
    real(kind=8), intent(in) :: p(M+2,N+2), phi(M+2,N+2,4)
    real(kind=8), intent(inout) :: Ed_half(M+1,N+2,4), Fd_half(M+2,N+1,4)

    ! numerical dissipation terms
    do k = 1, 4

        if (k .NE. 4) then
            Ed_half(:,:,k) = -0.5_dp * ( D_zeta(2:M+1,:) * (Phi(2:M+2,:,k)-Phi(1:M+1,:,k)) + (Pu_zeta+Pp_zeta)*Psi_zeta(:,:,k) + &
                                                                                        &    (Uu_zeta+Up_zeta)*Phi(2:M+2,:,k) )
            Fd_half(:,:,k) = -0.5_dp * ( D_eta(:,2:N+1) * (Phi(:,2:N+2,k)-Phi(:,1:N+1,k)) +  (Pu_eta+Pp_eta)*Psi_eta(:,:,k) + &
                                                                                        &    (Uu_eta+Up_eta)*Phi(:,2:N+2,k) )
        else
            Ed_half(:,:,k) = -0.5_dp * ( D_zeta(2:M+1,:) * ((Phi(2:M+2,:,k)-p(2:M+2,:))-(Phi(1:M+1,:,k)-p(1:M+1,:))) + &
                                                    &       (Pu_zeta+Pp_zeta)*Psi_zeta(:,:,k) + &
                                                    &       (Uu_zeta+Up_zeta)*Phi(2:M+2,:,k) )
            Fd_half(:,:,k) = -0.5_dp * ( D_eta(:,2:N+1) * ((Phi(:,2:N+2,k)-p(:,2:N+2))-(Phi(:,1:N+1,k)-p(:,1:N+1))) + &
                                                    &      (Pu_eta+Pp_eta)*Psi_eta(:,:,k) + &
                                                    &      (Uu_eta+Up_eta)*Phi(:,2:N+2,k) )

        end if

    end do

end