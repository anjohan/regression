module mod_regressor
    use iso_fortran_env, only: dp => real64
    use mod_basis
    use mod_utilities, only: check_info
    implicit none

    type, public, abstract :: regressor
        class(basisfunction), allocatable :: basis(:)
        real(dp), allocatable :: beta(:), X(:,:)

        contains
            procedure :: predict, create_X
            procedure(fit_procedure), deferred :: fit
    end type

    abstract interface
        subroutine fit_procedure(self, x_values, y_values)
            import regressor, dp
            class(regressor), intent(inout) :: self
            real(dp), intent(in), optional :: x_values(:,:)
            real(dp), intent(in) :: y_values(:)
        end subroutine
    end interface

    contains
        subroutine predict(self, x, y, y_exact, mse, r2)
            class(regressor), intent(in) :: self
            real(dp), intent(in) :: x(:,:)
            real(dp), intent(inout) :: y(:)
            real(dp), intent(in), optional :: y_exact(:)
            real(dp), intent(out), optional :: mse, r2

            integer :: N, i, j, p
            p = size(self%basis)
            N = size(y)
            y(:) = 0

            associate(beta => self%beta, &
                      basis => self%basis)
                do i = 1, N
                    do j = 1, p
                        y(i) = y(i) + beta(j)*basis(j)%eval(x(i,:))
                    end do
                end do
            end associate

            if (present(y_exact)) then
                if (present(mse)) mse = sum((y - y_exact)**2)/N
                if (present(r2))  r2  = 1 - sum((y - y_exact)**2) &
                                            / sum((y_exact - sum(y_exact)/N)**2)
            end if
        end subroutine

        subroutine create_X(self, x_values)
            class(regressor), intent(inout) :: self
            real(dp), intent(in) :: x_values(:,:)

            integer :: i, j, N, p

            N = size(x_values,1)
            p = size(self%basis)

            if (allocated(self%X)) deallocate(self%X)
            allocate(self%X(N, p))

            associate(X => self%X, &
                      basis => self%basis)
                do j = 1, p
                    do i = 1, N
                        X(i,j) = basis(j)%eval(x_values(i,:))
                    end do
                end do
            end associate
        end subroutine
end module
