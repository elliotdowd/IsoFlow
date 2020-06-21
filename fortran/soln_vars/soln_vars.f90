!-*- f90 -*- -
subroutine init_Q(Q, p, T, M_in, R, gam, M, N)

    implicit none
    integer, parameter :: dp=kind(0.d0)

    integer :: i, j
    
    integer, intent(in) :: M, N

    real(kind=8), intent(in) :: p
    real(kind=8), intent(in) :: T
    real(kind=8), intent(in) :: M_in
    real(kind=8), intent(in) :: R
    real(kind=8), intent(in) :: gam
    real(kind=8), intent(inout) :: Q(M, N, 4)

    do i = 1, M

       do j = 1, N

           Q(i,j,1) = p / (R * T)
           Q(i,j,2) = Q(i,j,1) * M_in * sqrt(gam*p/Q(i,j,1))
           Q(i,j,3) = 0.0_dp
           Q(i,j,4) = rho_et(p, Q(i,j,1), Q(i,j,2)/Q(i,j,1), Q(i,j,3)/Q(i,j,1), gam)    

       end do 

    end do  

       contains
              function rho_et(p, rho, u, v, gam)
                  real(kind=8) :: rho_et, p, rho, u, v, gam
                  integer, parameter :: dp=kind(0.d0)

                  rho_et = (p/(gam-1)) + 0.5_dp*rho*(u**2 + v**2)

              end function rho_et

end subroutine


!-*- f90 -*- -
subroutine calc_covariant(s_proj, u, v, uu, vv, M, N)
! calculate covariant velocity for a given cell (uu, vv)
! input projected cell face areas [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y, S_zeta, S_eta] and velocity components

    implicit none
    integer, parameter:: dp=kind(0.d0)

    integer :: i, j
    
    integer, intent(in) :: M, N

    real(kind=8), intent(in) :: s_proj(M,N,6)
    real(kind=8), intent(in) :: u(M,N)
    real(kind=8), intent(in) :: v(M,N)
    real(kind=8), intent(inout) :: uu(M,N)
    real(kind=8), intent(inout) :: vv(M,N)

    do i = 1, M
       
       do j = 1, N

           uu(:,:) = (1.0_dp/s_proj(i,j,5)) * &
               &    (u*s_proj(i,j,1) + v*s_proj(i,j,2))
       
           vv(:,:) = (1.0_dp/s_proj(i,j,6)) * &
               &    (u*s_proj(i,j,3) + v*s_proj(i,j,4))

       end do

    end do

end


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
subroutine face_flux( mdot_zeta, mdot_eta, phi, p_zeta, p_eta, )