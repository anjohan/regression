module mod_utilities
    use iso_fortran_env, only: dp => real64

    interface shuffle
        procedure :: shuffle_realreal, shuffle_realint, shuffle_intint
    end interface

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

        function random_meshgrid(N, filebase) result(x)
            integer, intent(in) :: N
            character(len=*), optional :: filebase

            real(dp), allocatable :: x1(:), x2(:), x(:,:)
            integer :: info, i, j, idx, u1, u2

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

            if(.not. present(filebase)) return

            open(newunit=u1, file=filebase // "_x1.dat", status="replace")
            write(u1, "(*(f0.6,:,','))") x1
            close(u1)
            open(newunit=u2, file=filebase // "_x2.dat", status="replace")
            write(u2, "(*(f0.6,:,','))") x2
            close(u2)

        end function

        subroutine shuffle_realreal(x, y, column_wise)
            real(dp), intent(inout) :: x(:,:), y(:)
            logical, intent(in), optional :: column_wise

            integer :: N, i, j
            real(dp) :: tmp_real, tmp_y
            real(dp), allocatable :: tmp_x(:)

            N = size(y)

            do i = N, 1, -1
                call random_number(tmp_real)
                j = nint((N-1)*tmp_real) + 1
                tmp_y = y(i)
                y(i) = y(j)
                y(j) = tmp_y
                if (present(column_wise)) then
                    tmp_x = x(:,i)
                    x(:,i) = x(:,j)
                    x(:,j) = tmp_x(:)
                else
                    tmp_x = x(i,:)
                    x(i,:) = x(j,:)
                    x(j,:) = tmp_x(:)
                end if
            end do
        end subroutine

        subroutine shuffle_intint(x, y, column_wise)
            integer, intent(inout) :: x(:,:)
            integer, intent(inout) :: y(:)
            logical, intent(in), optional :: column_wise

            integer :: N, i, j
            real(dp) :: tmp_real
            integer :: tmp_y
            integer, allocatable :: tmp_x(:)

            N = size(y)

            do i = N, 1, -1
                call random_number(tmp_real)
                j = nint((N-1)*tmp_real) + 1
                tmp_y = y(i)
                y(i) = y(j)
                y(j) = tmp_y
                if (present(column_wise)) then
                    tmp_x = x(:,i)
                    x(:,i) = x(:,j)
                    x(:,j) = tmp_x(:)
                else
                    tmp_x = x(i,:)
                    x(i,:) = x(j,:)
                    x(j,:) = tmp_x(:)
                end if
            end do
        end subroutine

        subroutine shuffle_realint(x, y, column_wise)
            real(dp), intent(inout) :: x(:,:)
            integer, intent(inout) :: y(:)
            logical, intent(in), optional :: column_wise

            integer :: N, i, j
            real(dp) :: tmp_real
            integer :: tmp_y
            real(dp), allocatable :: tmp_x(:)

            N = size(y)

            do i = N, 1, -1
                call random_number(tmp_real)
                j = nint((N-1)*tmp_real) + 1
                tmp_y = y(i)
                y(i) = y(j)
                y(j) = tmp_y
                if (present(column_wise)) then
                    tmp_x = x(:,i)
                    x(:,i) = x(:,j)
                    x(:,j) = tmp_x(:)
                else
                    tmp_x = x(i,:)
                    x(i,:) = x(j,:)
                    x(j,:) = tmp_x(:)
                end if
            end do
        end subroutine

        subroutine check_info(info, routine_name)
            integer, intent(in) :: info
            character(len=*) :: routine_name

            if (info /= 0) then
                write(*, "(a, i0)") "ERROR in " // routine_name // ":, info = ", info
                error stop
            end if
        end subroutine

        function franke(x) result(y)
            real(dp), intent(in) :: x(:,:)
            real(dp), allocatable :: y(:)

            associate(x1 => x(:,1), x2 => x(:,2))
                y = 0.75 * exp(- (9*x1-2)**2/4  - (9*x2-2)**2/4) &
                  + 0.75 * exp(- (9*x1+1)**2/49 - (9*x2+1)**2/10) &
                  + 0.50 * exp(- (9*x1-7)**2/4  - (9*x2-3)**2/4) &
                  - 0.20 * exp(- (9*x1-4)**2    - (9*x2-7)**2)
            end associate
        end function
end module
