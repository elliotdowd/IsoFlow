!-*- f90 -*- -
subroutine inv_metrics(xx, xy_inv, M, N)
! calculates inverse grid metrics
! xy_inv = [y_eta(i+.5,j) x_eta(i+.5,j) y_zeta(i,j+.5) x_zeta(i,j+.5)]
    implicit none
    real :: start, finish

    integer :: i, j
    integer, intent(in) :: M
    integer, intent(in) :: N
    real(kind=8), intent(in) :: xx(M+1, N+1, 2)
    real(kind=8), intent(inout) :: xy_inv(M, N, 4)

    call cpu_time(start)

    do i = 1, M
        
        do j = 1, N
            
            xy_inv(i,j,1) = xx(i+1,j+1,2) - xx(i+1,j,2)
            xy_inv(i,j,2) = xx(i+1,j+1,1) - xx(i+1,j,1)
            xy_inv(i,j,3) = xx(i+1,j+1,2) - xx(i,j+1,2)
            xy_inv(i,j,4) = xx(i+1,j+1,1) - xx(i,j+1,1)

        end do

    end do

    call cpu_time(finish)
    print '("inv_metrics time = ",f12.9," seconds.")',finish-start

end subroutine


!-*- f90 -*- -
subroutine face_areas(xy_inv, s_proj, M, N)
! projected cell face areas
! S_proj = [S_zeta.x, S_zeta.y, S_eta.x, S_eta.y S_zeta, S_eta]
    implicit none
    real :: start, finish

    integer :: i, j
    integer, intent(in) :: M
    integer, intent(in) :: N

    real(kind=8), intent(in) :: xy_inv(M, N, 4)
    real(kind=8), intent(inout) :: s_proj(M, N, 6)

    call cpu_time(start)

    do i = 1, M
        
        do j = 1, N
            
            s_proj(i,j,1) = xy_inv(i,j,1)
            s_proj(i,j,2) = -xy_inv(i,j,2)
            s_proj(i,j,3) = -xy_inv(i,j,3)
            s_proj(i,j,4) = xy_inv(i,j,4)
            
            s_proj(i,j,5) = sqrt(s_proj(i,j,1)**2 + s_proj(i,j,2)**2)
            s_proj(i,j,6) = sqrt(s_proj(i,j,3)**2 + s_proj(i,j,4)**2)
            
        end do

    end do

    call cpu_time(finish)
    print '("face_areas time = ",f12.9," seconds.")',finish-start

end subroutine


!-*- f90 -*- -
subroutine calc_cellCentroids(xx, yy, ccx, ccy, area, M, N)
    ! calculate polygon centroid and cell area given vectors of x and y points, cell area
    ! returns ccx, ccy centroid locations in arrays
    implicit none

    real :: start, finish

    integer, parameter:: dp=kind(0.d0)
    integer, parameter :: sides = 4
    integer :: i, j, k

    real :: x(sides+1), y(sides+1)

    integer, intent(in) :: M, N

    real(kind=8), intent(in) :: area(M, N)
    real(kind=8), intent(in) :: xx(M+1, N+1)
    real(kind=8), intent(in) :: yy(M+1, N+1)
    real(kind=8), intent(inout) :: ccx(M, N)
    real(kind=8), intent(inout) :: ccy(M, N)

    call cpu_time(start)

    do i = 1, M

        do j = 1, N

            x = (/ xx(i, j), xx(i+1, j), xx(i+1, j+1), xx(i, j+1), xx(i, j) /)
            y = (/ yy(i, j), yy(i+1, j), yy(i+1, j+1), yy(i, j+1), yy(i, j) /)


            do k = 1, sides
                ! cell area summation
                ccx(i,j) = ccx(i,j) + (1/(6.0_dp*area(i,j))) * ((x(k)+x(k+1)) * (x(k)*y(k+1)-x(k+1)*y(k)))
                ccy(i,j) = ccy(i,j) + (1/(6.0_dp*area(i,j))) * ((y(k)+y(k+1)) * (x(k)*y(k+1)-x(k+1)*y(k)))

            end do

        end do
                    
    end do        

    call cpu_time(finish)
    print '("calc_cellcentroids time = ",f12.9," seconds.")',finish-start

end subroutine


!-*- f90 -*- -
subroutine calc_cellArea(sides, xx, yy, area, M, N)
    ! calculate polygon area given polygon side numbers N, vectors of x and y points

    implicit none
    integer, parameter:: dp=kind(0.d0)
    real :: start, finish

    integer :: i, j, k
    integer, intent(in) :: sides
    integer, intent(in) :: M, N

    real(kind=8) :: x(sides), y(sides)
        
    real(kind=8), intent(in) :: xx(M+1, N+1)
    real(kind=8), intent(in) :: yy(M+1, N+1)
    real(kind=8), intent(inout) :: area(M, N)

    call cpu_time(start)

    do i = 1, M

        do j = 1, N

            x = (/ xx(i, j), xx(i+1, j), xx(i+1, j+1), xx(i, j+1), xx(i, j) /)
            y = (/ yy(i, j), yy(i+1, j), yy(i+1, j+1), yy(i, j+1), yy(i, j) /)

            area(i,j) = 0
            do k = 1, sides
                ! cell area summation
                area(i,j) = area(i,j) + 0.5_dp * (x(k)*y(k+1) - x(k+1)*y(k))

            end do

        end do
    
    end do

    call cpu_time(finish)
    print '("calc_cellarea time = ",f12.9," seconds.")',finish-start

end subroutine