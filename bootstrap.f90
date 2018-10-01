module mod_bootstrap
    use iso_fortran_env, only: dp => real64
    use mod_regressor
    implicit none

    interface bootstrapper
        procedure :: init_bootstrapper
    end interface

    type, public :: bootstrapper
        class(regressor), allocatable :: fitter
        real(dp), allocatable :: x(:,:), y(:), y_predictions(:,:), betas(:,:), &
                                 R2s(:), MSEs(:)
        integer :: N, p

        contains
            procedure :: bootstrap
    end type bootstrapper

    contains
        function init_bootstrapper(fitter) result(bs)
            class(regressor), intent(in) :: fitter

            type(bootstrapper) :: bs

            bs%fitter = fitter
        end function

        subroutine bootstrap(self, x, y, num_bootstraps)
            class(bootstrapper), intent(inout) :: self
            integer, intent(in) :: num_bootstraps
            real(dp), intent(in) :: x(:,:), y(:)

            integer :: N, p, i
            real(dp), allocatable :: tmp_real(:), y_selection(:), &
                                     X_original(:,:), X_selection(:,:)
            integer, allocatable :: indices(:)

            N = size(y)
            p = size(self%fitter%basis)

            self%x = x; self%y = y

            if (allocated(self%betas)) deallocate(self%betas)
            if (allocated(self%y_predictions)) deallocate(self%y_predictions)

            allocate(self%betas(p,num_bootstraps), self%y_predictions(N, num_bootstraps))
            allocate(self%R2s(num_bootstraps), self%MSEs(num_bootstraps))
            allocate(tmp_real(N), indices(N))
            allocate(X_selection(N,p), y_selection(N))

            associate(betas => self%betas, &
                      y_predictions => self%y_predictions, &
                      R2s => self%R2s, MSEs => self%MSEs, &
                      fitter => self%fitter)

                call fitter%create_X(x)
                X_original = fitter%X

                bootstraps: do i = 1, num_bootstraps
                    call random_number(tmp_real)
                    indices(:) = nint((N-1)*tmp_real) + 1

                    y_selection(:) = y(indices)
                    X_selection(:,:) = X_original(indices,:)

                    fitter%X = X_selection

                    call fitter%fit(y_values=y_selection)
                    call fitter%predict(x, y_predictions(:,i), y, MSEs(i), R2s(i))

                    betas(:, i) = fitter%beta
                end do bootstraps
            end associate
        end subroutine
end module
