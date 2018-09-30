module mod_lasso2d
    use iso_fortran_env, only: dp => real64
    use mod_polyfitter2d
    implicit none

    interface lasso2d
        procedure :: init_lasso2d
    end interface

    type, public, extends(polyfitter2d) :: lasso2d
        real(dp) :: lambda, tolerance

        contains
            procedure :: fit => fit_lasso
    end type lasso2d

    contains
        function init_lasso2d(d, lambda, tolerance) result(lasso)
            integer, intent(in) :: d
            real(dp), intent(in) :: lambda
            real(dp), intent(in), optional :: tolerance
            type(lasso2d) :: lasso

            lasso%d = d
            lasso%p = (d+1)*(d+2)/2
            lasso%lambda = lambda

            if (present(tolerance)) then
                lasso%tolerance = tolerance
            else
                lasso%tolerance = 1e-5
            end if
        end function

        subroutine fit_lasso(self, x_values, y_values, x_is_matrix)
            class(lasso2d), intent(inout) :: self
            real(dp), intent(in) :: x_values(:,:), y_values(:)
            logical, intent(in), optional :: x_is_matrix

            integer :: N, p, info
            real(dp), allocatable :: X_T(:,:), X_T_X(:,:), H(:,:), H_cholesky(:), &
                                     X_T_y(:), two_X_T_y(:), &
                                     beta(:), grad(:)
            real(dp) :: lambda, grad_norm, tolerance

            if (present(x_is_matrix)) then
                call self%init_X(x_values, x_is_matrix)
            else
                call self%init_X(x_values)
            end if

            N = self%N
            p = self%p
            lambda = self%lambda
            tolerance = self%tolerance

            associate(X => self%X)
                X_T = transpose(X)
                X_T_X = matmul(X_T, X)
            end associate

            H = 2*X_T_X
            X_T_y = matmul(X_T, y_values)
            two_X_T_y = 2*X_T_y

            ! create triag version for cholesky
            allocate(H_cholesky(p*(p+1)/2))
            call dtrttp('L', p, H, p, H_cholesky, info)
            call check_info(info, "dtrttp")

            call dpptrf('L', p, H_cholesky, info)
            call check_info(info, "dpptrf")

            allocate(beta(p), grad(p))
            grad_norm = 2*tolerance

            !/lassostart/!
            beta(:) = 1
            minimisation: do while (grad_norm > tolerance)
                ! calculate right-hand side
                beta(:) = two_X_T_y - lambda*sgn(beta)

                ! solve H beta_new = 2 X^T y - lambda sgn(beta_old)
                call dpptrs('L', p, 1, H_cholesky, beta, p, info)
                call check_info(info, "dpptrs")

                ! calculate gradient and its norm
                grad(:) = - two_X_T_y + matmul(H, beta) + lambda*sgn(beta)
                grad_norm = norm2(grad)
            end do minimisation
            !/lassoend/!

            self%beta = beta
        end subroutine

        elemental function sgn(x)
            real(dp), intent(in) :: x
            integer :: sgn

            if (x > 0) then
                sgn = 1
            else if (x < 0) then
                sgn = -1
            else
                sgn = 0
            end if
        end function
end module
