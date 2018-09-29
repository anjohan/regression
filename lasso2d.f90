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

        subroutine fit_lasso(self, x_values, y_values)
            class(lasso2d), intent(inout) :: self
            real(dp), intent(in) :: x_values(:,:), y_values(:)

            integer :: N, p, info, i, j
            real(dp), allocatable :: X_T(:,:), X_T_X(:,:), H(:,:), H_cholesky(:), &
                                     X_T_y(:), two_X_T_y(:), &
                                     beta(:), grad(:)
            real(dp) :: lambda, norm

            call self%create_X(x_values)

            N = self%N
            p = self%p
            lambda = self%lambda

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

            !do i = 1, p
            !    write(*,*) (H(i,j), j = 1, p)
            !end do

            call dpptrf('L', p, H_cholesky, info)
            call check_info(info, "dpptrf")

            allocate(beta(p), grad(p))
            beta(:) = 1.0d0

            do
!                 grad = matmul(H,beta)
!                 !call dsymv('L', p, 1.0d0, H, p, beta, 1, 0, grad, 1)
!
!                 grad(:) = grad(:) - two_X_T_y(:) + lambda*sgn(beta(:))
!                 norm = norm2(grad)
!
!                 call dpptrs('L', p, 1, H_cholesky, grad, p, info)
!
!                 call check_info(info, "dpptrs")
!
!                 beta(:) = beta(:) - grad(:)
!
                grad(:) = -2*matmul(X_T, y_values - matmul(self%X,beta)) + lambda * sgn(beta)
                norm = norm2(grad)

                beta(:) = beta - 0.0001 * grad
                write(*,*) "Norm of gradient:", norm
                if (norm < self%tolerance) exit
            end do

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
