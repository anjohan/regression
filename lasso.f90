module mod_lasso
    use iso_fortran_env, only: dp => real64
    use mod_regressor
    implicit none

    interface lasso
        procedure :: init_lasso
    end interface

    type, public, extends(regressor) :: lasso
        real(dp) :: lambda, grad_tolerance, alpha
        integer :: max_iterations

        contains
            procedure :: fit
    end type lasso

    contains
        function init_lasso(lambda, basis, X, grad_tolerance, max_iterations, alpha) result(self)
            real(dp), intent(in) :: lambda
            class(basisfunction), intent(in), optional :: basis(:)
            real(dp), intent(in), optional :: X(:,:), grad_tolerance, alpha
            integer, intent(in), optional :: max_iterations

            type(lasso) :: self

            self%lambda = lambda

            if (present(basis)) allocate(self%basis, source=basis)
            if (present(X)) self%X = X

            if (present(grad_tolerance)) then
                self%grad_tolerance = grad_tolerance
            else
                self%grad_tolerance = 1e-5
            end if

            if (present(max_iterations)) then
                self%max_iterations = max_iterations
            else
                self%max_iterations = 100
            end if
            if (present(alpha)) then
                self%alpha = alpha
            else
                self%alpha = 1.0d0
            end if

            self%method = "LASSO"
        end function

        subroutine fit(self, x_values, y_values)
            class(lasso), intent(inout) :: self
            real(dp), intent(in), optional :: x_values(:,:)
            real(dp), intent(in) :: y_values(:)

            integer :: N, p, info, iteration, max_iterations
            real(dp), allocatable :: X_T(:,:), X_T_X(:,:), H(:,:), H_cholesky(:), &
                                     X_T_y(:), two_X_T_y(:), &
                                     beta(:), grad(:), new_beta(:)
            real(dp) :: lambda, grad_norm, grad_tolerance, min_grad_norm, alpha

            if (present(x_values)) call self%create_X(x_values)

            N = size(y_values)
            p = size(self%basis)
            lambda = self%lambda
            grad_tolerance = self%grad_tolerance
            max_iterations = self%max_iterations
            alpha = self%alpha

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
            grad_norm = 2*grad_tolerance

            !/lassostart/!
            iteration = 0
            beta(:) = 1
            minimisation: do while (grad_norm > grad_tolerance &
                                    .and. iteration < max_iterations)
                iteration = iteration + 1
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
!
!             write(*,*) iteration, max_iterations
!             if (iteration < max_iterations) return
!             min_grad_norm = grad_norm
!             iteration = 0
!             alpha = 1.0d0
!             allocate(new_beta(p))
!             do while (grad_norm > grad_tolerance &
!                       .and. iteration < 1*max_iterations)
!                 iteration = iteration + 1
!                 reduce_alpha: do
!                     new_beta(:) = beta(:) - alpha*grad(:)
!                     grad(:) = - two_X_T_y + matmul(H, new_beta) + lambda*sgn(new_beta)
!                     grad_norm = norm2(grad)
!                     if (grad_norm < min_grad_norm) then
!                         min_grad_norm = grad_norm
!                         exit reduce_alpha
!                     end if
!                     alpha = 0.9d0*alpha
!                 end do reduce_alpha
!                 write(*,*) grad_norm, alpha
!                 beta(:) = new_beta(:)
!             end do
!             self%beta = beta
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
