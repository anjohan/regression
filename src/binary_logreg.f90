module mod_binary_logreg
    use mod_regressor
    use iso_fortran_env, only: dp => real64
    implicit none

    interface binary_logreg
        procedure :: init_binary_logreg
    end interface

    type, public, extends(regressor) :: binary_logreg
        real(dp) :: lambda, learning_rate, momentum, tolerance
        integer :: batch_size, max_iterations
        contains
            procedure :: fit
            !procedure :: predict
    end type

    contains
        function init_binary_logreg(basis, learning_rate, X, lambda, &
                                    momentum, batch_size, tolerance, max_iterations) &
                 result(self)
            class(basisfunction), intent(in), optional :: basis(:)
            real(dp), intent(in), optional :: X(:,:), lambda, momentum, tolerance
            integer, intent(in), optional :: batch_size, max_iterations
            real(dp), intent(in) :: learning_rate
            type(binary_logreg) :: self

            if (present(basis)) allocate(self%basis, source=basis)
            if (present(X)) self%X = X

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

            self%method = "binary_logreg"
        end function

        subroutine fit(self, x_values, y_values)
            class(binary_logreg), intent(inout) :: self
            real(dp), intent(in) :: y_values(:)
            real(dp), intent(in), optional :: x_values(:,:)

            integer :: N, p, batch_size, num_batches, i, j, idx, k
            real(dp), allocatable :: X_T(:,:), tmp_real(:), &
                                     grad(:), beta(:), prev_grad(:)
            real(dp) :: grad_norm, prediction

            if (present(x_values)) call self%create_X(x_values)

            N = size(y_values)
            p = size(self%basis)

            batch_size = self%batch_size
            num_batches = N/batch_size

            X_T = transpose(self%X)
            allocate(tmp_real(batch_size), &
                     grad(p), beta(p), prev_grad(p))

            call random_number(beta)

            steps: do k = 1, self%max_iterations
                batches: do i = 1, num_batches
                    grad(:) = 0
                    call random_number(tmp_real)

                    do j = 1, batch_size
                        idx = floor(N*tmp_real(j)) + 1
                        prediction = 1/(1 + exp(-dot_product(X_T(:,idx), beta)))
                        grad(:) = grad(:) + (prediction - y_values(idx))*X_T(:, idx)
                    end do

                    if (self%lambda /= 0) grad(:) = grad(:) + self%lambda*beta(:)
                    if (k > 1 .and. self%momentum /= 0) &
                        grad(:) = self%momentum*prev_grad(:) + grad(:)

                    beta(:) = beta(:) - self%learning_rate/batch_size*grad(:)
                    grad_norm = maxval(abs(grad))

                    if (grad_norm < self%tolerance) exit steps
                    prev_grad(:) = grad(:)
                end do batches
            end do steps

            self%beta = beta
        end subroutine
end module
