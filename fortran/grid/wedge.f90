!-*- f90 -*- -
subroutine mod2wedge(xx, yy, height, theta, wedge_start, M, N)

    implicit none
    real :: start, finish

    integer :: i, j
    integer, intent(in) :: M, N

    real(kind=8) :: tantheta

    real(kind=8), intent(in) :: height
    real(kind=8), intent(in) :: theta
    real(kind=8), intent(in) :: wedge_start
    real(kind=8), intent(in) :: xx(M+3, N+3)
    real(kind=8), intent(inout) :: yy(M+3, N+3)

    call cpu_time(start)

    tantheta = tan(theta)

    do i = 0, M+2

        do j = 0, N+2

            ! scale y-coordinates for wedge section
            if (xx(i+1, j+1) >= wedge_start) then
                ! scale for wedge shape

                yy(i+1, j+1) = yy(i+1, j+1) - height * (real(j-1)/real(N+2)) * (xx(i+1, j+1) - wedge_start) * tantheta

            end if

            ! flip geometry about x-axis

            yy(i+1, j+1) = height - yy(i+1, j+1) 

        end do

    end do

    call cpu_time(finish)
    print '("mod2wedge time = ",f12.9," seconds.")', finish-start

end subroutine