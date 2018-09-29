module mod_ridge2d
    use iso_fortran_env, only: dp => real64
    use mod_polyfitter2d
    implicit none

    interface ridge2d
        procedure :: init_ridge2d
    end interface

    type, public, extends(polyfitter2d) :: ridge2d
        real(dp) :: lambda

        contains
            procedure :: fit => fit_ridge
    end type ridge2d

    contains
        function init_ridge2d(d, lambda) result(ridge)
            integer, intent(in) :: d
            real(dp), intent(in) :: lambda
            type(ridge2d) :: ridge

            ridge%d = d
            ridge%p = (d+1)*(d+2)/2
            ridge%lambda = lambda
        end function

        subroutine fit_ridge(self, x_values, y_values)
            class(ridge2d), intent(inout) :: self
            real(dp), intent(in) :: x_values(:,:), y_values(:)

            integer :: N, p, info, i
            real(dp), allocatable :: X_T(:,:), A(:,:), b(:)
            real(dp) :: lambda

            call self%create_X(x_values)

            N = self%N
            p = self%p
            lambda = self%lambda

            associate(X => self%X)
                X_T = transpose(X)
                A = matmul(X_T, X)
            end associate

            do i = 1, p
                A(i,i) = A(i,i) + lambda
            end do

            b = matmul(X_T, y_values)

            call dposv('L', p, 1, A, p, b, p, info)

            if (info /= 0) then
                write(*, "(a,i0)") "ERROR in dposv: info = ", info
            end if

            self%beta = b
        end subroutine
end module
