module mod_bootstrap
    use iso_fortran_env, only: dp => real64
    use mod_utilities, only: shuffle
    use mod_regressor
    implicit none

    interface bootstrapper
        procedure :: init_bootstrapper
    end interface

    type, public :: bootstrapper
        class(regressor), allocatable :: fitter
        real(dp), allocatable :: x(:,:), x_test(:,:), y(:), y_test(:), &
                                 y_predictions(:,:), &
                                 betas(:,:), &
                                 R2s(:), MSEs(:), mean_beta(:), &
                                 beta_variance(:), final_y(:)
        real(dp) :: bias, variance, final_MSE, final_R2, mean_MSE


        contains
            procedure :: bootstrap, post_analysis
    end type bootstrapper

    contains
        function init_bootstrapper(fitter) result(bs)
            class(regressor), intent(in) :: fitter

            type(bootstrapper) :: bs

            bs%fitter = fitter
        end function

        subroutine bootstrap(self, x, y, num_bootstraps, test_fraction)
            class(bootstrapper), intent(inout) :: self
            integer, intent(in) :: num_bootstraps
            real(dp), intent(in) :: x(:,:), y(:), test_fraction

            class(regressor), allocatable :: fitter
            integer :: N, p, i, N_test, N_train
            real(dp), allocatable :: tmp_real(:), y_selection(:), &
                                     X_original(:,:), X_selection(:,:), &
                                     x_test(:,:), y_test(:), x_train(:,:), y_train(:), &
                                     betas(:,:), y_predictions(:,:), R2s(:), MSEs(:)
            integer, allocatable :: indices(:)

            N = size(y)
            p = size(self%fitter%basis)

            self%x = x; self%y = y

            N_test = nint(test_fraction * N)
            N_train = N - N_test

            call shuffle(self%x, self%y)

            x_test = self%x(1:N_test,:)
            y_test = self%y(1:N_test)
            x_train = x(N_test+1:,:)
            y_train = y(N_test+1:)

            allocate(betas(p,num_bootstraps), &
                     y_predictions(N_test, num_bootstraps), &
                     R2s(num_bootstraps), &
                     MSEs(num_bootstraps), &
                     tmp_real(N_train), &
                     indices(N_train), &
                     X_selection(N_train,p), &
                     y_selection(N_train))


            call self%fitter%create_X(x_train)
            X_original = self%fitter%X

            write(*, "(a)", advance="no") "Bootstrapping " // self%fitter%method // ": ["

            !$omp parallel do private(fitter, i, tmp_real, indices, y_selection, X_selection) &
            !$omp& shared(betas,y_predictions,MSEs,R2s)
            bootstraps: do i = 1, num_bootstraps
                if (mod(i, max(1,num_bootstraps/50)) == 0) then
                    write(*, "(a)", advance="no") "="
                end if
                if (.not. allocated(fitter)) fitter = self%fitter

                call random_number(tmp_real)
                indices(:) = nint((N_train-1)*tmp_real) + 1
                !write(*,*) size(indices)

                y_selection(:) = y_train(indices)
                X_selection(:,:) = X_original(indices,:)

                fitter%X = X_selection

                call fitter%fit(y_values=y_selection)
                call fitter%predict(x_test, y_predictions(:,i), &
                                         y_test, MSEs(i), R2s(i))

                betas(:, i) = fitter%beta
                !write(*,*) "i = ", i
            end do bootstraps
            !$omp end parallel do
            write(*,"(a)") "]"

            self%x_test = x_test
            self%y_test = y_test
            self%betas = betas
            self%y_predictions = y_predictions
            self%MSEs = MSEs
            self%R2s = R2s

            call self%post_analysis()
        end subroutine

        subroutine post_analysis(self)
            class(bootstrapper), intent(inout) :: self

            integer :: num_bootstraps, p, i, N, N_test, N_train
            real(dp), allocatable :: mean_y_test(:)

            N = size(self%y)
            N_test = size(self%y_test)
            N_train = N - N_test
            p = size(self%betas, 1)
            num_bootstraps = size(self%betas,2)
            self%mean_beta = sum(self%betas, dim=2)/num_bootstraps
            self%fitter%beta = self%mean_beta

            allocate(self%final_y(N))
            call self%fitter%create_X(self%x)
            call self%fitter%predict(self%x, self%final_y, self%y, &
                                     self%final_MSE, self%final_R2)

            self%mean_MSE = sum(self%MSEs)/num_bootstraps

            mean_y_test = sum(self%y_predictions, dim=2)/num_bootstraps
            self%bias = sum((self%y_test - mean_y_test)**2)/N_test

            self%variance = 0
            do i = 1, N_test
                self%variance = self%variance + sum((self%y_predictions(i,:) - mean_y_test(i))**2)
            end do
            self%variance = self%variance/(N_test * num_bootstraps)

            allocate(self%beta_variance(p))
            do i = 1, p
                 self%beta_variance(i) = sum((self%betas(i,:) - self%mean_beta(i))**2) &
                                         / num_bootstraps
            end do
        end subroutine
end module
