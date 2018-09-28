module mod_polyfitter2d
    use iso_fortran_env, only: dp => real64
    implicit none

    interface polyfitter2d
        procedure :: init_polyfitter2d
    end interface

    type, public :: polyfitter2d
        real(dp), allocatable :: beta(:), X(:,:)
        integer :: d, p, N

        contains
            procedure :: create_X, fit, predict
    end type polyfitter2d

    contains
        function init_polyfitter2d(d) result(polyfitter)
            integer, intent(in) :: d
            type(polyfitter2d) :: polyfitter

            polyfitter%d = d
            polyfitter%p = (d+1)*(d+2)/2
        end function

        subroutine create_X(self, x_values, y_values)
            class(polyfitter2d), intent(inout) :: self
            real(dp), intent(in) :: x_values(:,:), y_values(:)

            integer :: i, j, idx

            self%N = size(y_values)

            if (allocated(self%X)) deallocate(self%X)
            allocate(self%X(self%N, self%p))

            associate(d => self%d, X => self%X)
                idx = 0
                do j = 0, d
                    do i = 0, d-j
                        idx = idx + 1
                        X(:,idx) = x_values(:,1)**i * x_values(:,2)**j
                    end do
                end do
            end associate
        end subroutine

        subroutine fit(self, x_values, y_values)
            class(polyfitter2d), intent(inout) :: self
            real(dp), intent(in) :: x_values(:,:), y_values(:)

            integer :: lwork, info, N, p
            real(dp), allocatable :: work(:), X(:,:), y(:)

            call self%create_X(x_values, y_values)

            lwork = -1
            allocate(work(1))

            N = self%N
            p = self%p

            ! copies, dgels modifies
            X = self%x
            y = y_values

            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work)
            allocate(work(lwork))
            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)

            if (info /= 0) then
                write(*,"(a,i0)") "DGELS problem: info = ", info
                error stop
            end if

            self%beta = y(1:p)
        end subroutine

        subroutine predict(self, x, y, y_exact, mse, r2)
            class(polyfitter2d), intent(in) :: self
            real(dp), intent(in) :: x(:,:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(out), optional :: y_exact(:), mse, r2

            integer :: N, d, p, idx, i, j

            if (.not. allocated(self%beta)) then
                write(*,*) "ERROR: Attempting to predict without fitting"
                error stop
            end if

            N = size(x,1)
            d = self%d
            p = self%p

            y(:) = 0

            associate(beta => self%beta)
                idx = 0
                do j = 0, d
                    do i = 0, d-j
                        idx = idx + 1
                        y(:) = y(:) + beta(idx) * x(:,1)**i * x(:,2)**j
                    end do
                end do
            end associate

            if (present(y_exact)) then
                if (present(mse)) mse = sum((y - y_exact)**2)/N
                if (present(r2)) r2   = 1 - sum((y - y_exact)**2) &
                                            / sum((y_exact - sum(y_exact)/N)**2)
            end if

        end subroutine
end module
