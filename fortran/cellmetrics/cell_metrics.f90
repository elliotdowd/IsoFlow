!-*- f90 -*- -
subroutine inv_metrics(xx, xy_inv, M, N)
! calculates inverse grid metrics
! xy_inv = [y_eta(i+.5,j) x_eta(i+.5,j) y_zeta(i,j+.5) x_zeta(i,j+.5)]

implicit none

integer :: i, j
integer, intent(in) :: M
integer, intent(in) :: N
real(kind=8), dimension(M+1, N+1, 2), intent(in) :: xx
real(kind=8), dimension(M, N, 4), intent(inout) :: xy_inv

do i = 1, M
    
    do j = 1, N
        
        xy_inv(i,j,1) = xx(i+1,j+1,2) - xx(i+1,j,2)
        xy_inv(i,j,2) = xx(i+1,j+1,1) - xx(i+1,j,1)
        xy_inv(i,j,3) = xx(i+1,j+1,2) - xx(i,j+1,2)
        xy_inv(i,j,4) = xx(i+1,j+1,1) - xx(i,j+1,1)

    end do

end do

end subroutine



!-*- f90 -*- -
subroutine face_areas(xy_inv, s_proj, M, N)

! projected cell face areas
! S_proj = [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y S_zeta, S_eta]

implicit none

integer :: i, j
integer, intent(in) :: M
integer, intent(in) :: N

real(kind=8), dimension(M, N, 4), intent(in) :: xy_inv
real(kind=8), dimension(M, N, 6), intent(inout) :: s_proj

do i = 1, M
    
    do j = 1, N
        
        s_proj(i,j,1) = xy_inv(i,j,1)
        s_proj(i,j,2) = -xy_inv(i,j,2)
        s_proj(i,j,3) = -xy_inv(i,j,3)
        s_proj(i,j,4) = xy_inv(i,j,4)
        
        s_proj(i,j,5) = sqrt(s_proj(i,j,1)**2 + s_proj(i,j,2)**2)
        s_proj(i,j,6) = sqrt(s_proj(i,j,3)**2 + s_proj(i,j,4)**2)
        
        !normal(i,j,1) = S_proj(i,j,1) / S_face(i,j,1)
        !normal(i,j,2) = S_proj(i,j,2) / S_face(i,j,1)

        
    end do

end do

end subroutine


!-*- f90 -*- -
subroutine calc_cellCentroids(xx, yy, ccx, ccy area, M, N)
    ! calculate polygon centroid and cell area given vectors of x and y points, cell area
    ! returns cx, cy in vector
    implicit none

    integer, parameter :: sides = 4
    integer :: i, j, k

    real, dimension(sides+1) :: x, y
    real :: A
    real, dimension(2) :: temp

    integer, intent(in) :: M, N

    real(kind=8), dimension(M, N), intent(inout) :: ccx, ccy
    real(kind=8), dimension(M, N), intent(in) :: area
    real(kind=8), dimension(M+1, N+1), intent(in) :: xx, yy

    do i = 1, M

        do j = 1, N

            x = (/ xx(i, j), xx(i+1, j), xx(i+1, j+1), xx(i, j+1), xx(i, j) /)
            y = (/ yy(i, j), yy(i+1, j), yy(i+1, j+1), yy(i, j+1), yy(i, j) /)


            do k = 1, sides
                ! cell area summation
                ccx = ccx(i,j) + (1/(6*area(i,j))) * ((x(k)+x(k+1)) * (x(k)*y(k+1)-x(k+1)*y(k)))
                ccy = ccy(i,j) + (1/(6*area(i,j))) * ((y(k)+y(k+1)) * (x(k)*y(k+1)-x(k+1)*y(k)))

            end do

        end do
                    
    end do        

end subroutine



function calc_cellArea(N, xx)
        ! calculate polygon area given polygon side numbers N, vectors of x and y points

    integer :: i, j, k
    integer, intent(in) :: N
        
    real, dimension(N), intent(in) :: x, y
    real :: calc_cellArea

    do i = 1, M

        do j = 1, N

            x = (/ xx(i, j, 1), xx(i+1, j, 1), xx(i+1, j+1, 1), xx(i, j+1, 1), xx(i, j, 1) /)
            y = (/ xx(i, j, 2), xx(i+1, j, 2), xx(i+1, j+1, 2), xx(i, j+1, 2), xx(i, j, 2) /)

            calc_cellArea = 0
            do i = 1, N
                ! cell area summation
                calc_cellArea = calc_cellArea + 0.5 * (x(i)*y(i+1) - x(i+1)*y(i))

            end do

        end do
    
    end do

end function calc_cellArea