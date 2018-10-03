module mod_lasso
    use iso_fortran_env, only: dp => real64
    use mod_regressor
    implicit none

    interface lasso
        procedure :: init_lasso
    end interface

    type, public, extends(regressor) :: lasso
        real(dp) :: lambda, tolerance

        contains
            procedure :: fit
    end type lasso

    contains
        function init_lasso(lambda, basis, X, tolerance) result(self)
            real(dp), intent(in) :: lambda
            class(basisfunction), intent(in), optional :: basis(:)
            real(dp), intent(in), optional :: X(:,:), tolerance

            type(lasso) :: self

            self%lambda = lambda

            if (present(basis)) allocate(self%basis, source=basis)
            if (present(X)) self%X = X

            if (present(tolerance)) then
                self%tolerance = tolerance
            else
                self%tolerance = 1e-5
            end if

            self%method = "LASSO"
        end function

        subroutine fit(self, x_values, y_values)
            class(lasso), intent(inout) :: self
            real(dp), intent(in), optional :: x_values(:,:)
            real(dp), intent(in) :: y_values(:)

            integer :: N, p, info
            real(dp), allocatable :: X_T(:,:), X_T_X(:,:), H(:,:), H_cholesky(:), &
                                     X_T_y(:), two_X_T_y(:), &
                                     beta(:), grad(:)
            real(dp) :: lambda, grad_norm, tolerance

            if (present(x_values)) call self%create_X(x_values)

            N = size(y_values)
            p = size(self%basis)
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
