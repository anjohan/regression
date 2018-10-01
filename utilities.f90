module mod_utilities
    use iso_fortran_env, only: dp => real64

    contains
        subroutine add_noise(y, sigma)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in) :: sigma

            integer :: N, seed(4)
            real(dp), allocatable :: noise(:)

            N = size(y)
            allocate(noise(N))

            call random_number(noise(1:4))
            seed = nint(4095*noise(1:4))
            if (mod(seed(4),2) == 0) then
                if (seed(4) > 0) then
                    seed(4) = seed(4) - 1
                else
                    seed(4) = 1
                end if
            end if

            call dlarnv(3, seed, N, noise)

            y(:) = y(:) + sigma*noise
        end subroutine

        function random_meshgrid(N) result(x)
            integer, intent(in) :: N

            real(dp), allocatable :: x1(:), x2(:), x(:,:)
            integer :: info, i, j, idx

            allocate(x1(N), x2(N), x(N**2,2))

            call random_number(x1)
            call random_number(x2)

            call dlasrt('I', N, x1, info)
            call check_info(info, "dlasrt(x1)")
            call dlasrt('I', N, x2, info)
            call check_info(info, "dlasrt(x2)")

            idx = 0
            do i = 1, N
                do j = 1, N
                    idx = idx + 1
                    x(idx, :) = [x1(i), x2(j)]
                end do
            end do
        end function

        subroutine check_info(info, routine_name)
            integer, intent(in) :: info
            character(len=*) :: routine_name

            if (info /= 0) then
                write(*, "(a, i0)") "ERROR in " // routine_name // ":, info = ", info
                error stop
            end if
        end subroutine
end module
