module mod_binary_logreg
    use iso_fortran_env, only: dp => real64
    implicit none

    interface binary_logreg
        procedure :: init_binary_logreg
    end interface

    type, public :: binary_logreg
        real(dp), allocatable :: beta(:)
        real(dp) :: lambda, learning_rate, momentum, tolerance
        integer :: batch_size, max_iterations
        contains
            procedure :: fit
            procedure :: predict
    end type

    contains
        function init_binary_logreg(learning_rate, lambda, &
                                    momentum, batch_size, tolerance, max_iterations) &
                 result(self)
            real(dp), intent(in), optional :: lambda, momentum, tolerance
            integer, intent(in), optional :: batch_size, max_iterations
            real(dp), intent(in) :: learning_rate
            type(binary_logreg) :: self

            self%learning_rate = learning_rate

            if (present(lambda)) then
                self%lambda = lambda
            else
                self%lambda = 0
            end if
            if (present(momentum)) then
                self%momentum = momentum
            else
                self%momentum = 0
            end if
            if (present(batch_size)) then
                self%batch_size = batch_size
            else
                self%batch_size = 20
            end if
            if (present(tolerance)) then
                self%tolerance = tolerance
            else
                self%tolerance = 1.0d-4
            end if
            if (present(max_iterations)) then
                self%max_iterations = max_iterations
            else
                self%max_iterations = 100
            end if

        end function

        subroutine predict(self, X, y)
            class(binary_logreg), intent(in) :: self
            integer, intent(inout) :: y(:)
            real(dp), intent(in) :: X(:,:)

            real(dp), allocatable :: p(:)
            integer :: N

            N = size(y)
            allocate(p(N))

            if (.not. allocated(self%beta)) error stop

            p(:) = matmul(X, self%beta)
            p(:) = 1/(1 + exp(-p(:)))

            where (p > 0.5d0)
                y = 1
            else where
                y = 0
            end where
        end subroutine

        subroutine fit(self, X, y)
            class(binary_logreg), intent(inout) :: self
            integer, intent(in) :: y(:)
            real(dp), intent(in)  :: X(:,:)

            integer :: N, p, batch_size, num_batches, i, j, idx, k
            real(dp), allocatable :: X_T(:,:), tmp_real(:), &
                                     grad(:), beta(:), prev_grad(:), &
                                     y_pred(:)
            integer, allocatable :: indices(:)
            real(dp) :: grad_norm, prediction

            N = size(y)
            p = size(X, 2)

            batch_size = self%batch_size
            num_batches = N/batch_size

            X_T = transpose(X)
            allocate(tmp_real(batch_size), &
                     grad(p), beta(p), prev_grad(p))

            call random_number(beta)
            beta(:) = beta(:) - 0.5d0

            allocate(y_pred(batch_size), indices(batch_size))

            steps: do k = 1, self%max_iterations
                batches: do i = 1, num_batches
                    grad(:) = 0
                    call random_number(tmp_real)
                    indices(:) = floor(N*tmp_real(:)) + 1

                    do j = 1, batch_size
                        idx = indices(j)
                        prediction = 1/(1 + exp(-dot_product(X_T(:,idx), beta)))
                        grad(:) = grad(:) + (prediction - y(idx))*X_T(:, idx)
                    end do
!
                    grad(:) = grad(:) / batch_size
                    grad_norm = maxval(abs(grad))
!
                    if (self%lambda /= 0) grad(:) = grad(:) + self%lambda*beta(:)
                    if (k+i > 2 .and. self%momentum /= 0) &
                        grad(:) = self%momentum*prev_grad(:) + grad(:)

                    beta(:) = beta(:) - self%learning_rate/batch_size*grad(:)

                    if (grad_norm < self%tolerance) exit steps
                    prev_grad(:) = grad(:)
                end do batches
            end do steps

            self%beta = beta
        end subroutine
end module
