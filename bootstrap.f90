module mod_bootstrap
    use iso_fortran_env, only: dp => real64
    use mod_polyfitter2d
    implicit none

    interface bootstrapper
        procedure :: init_bootstrapper
    end interface

    type, public :: bootstrapper
        class(polyfitter2d), allocatable :: regressor
        real(dp), allocatable :: x(:,:), y(:), y_predictions(:,:), betas(:,:), &
                                 R2s(:), MSEs(:)
        integer :: N, p

        contains
            procedure :: bootstrap
    end type bootstrapper

    contains
        function init_bootstrapper(regressor) result(bs)
            class(polyfitter2d), intent(in) :: regressor

            type(bootstrapper) :: bs

            bs%regressor = regressor
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
            p = self%regressor%p

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
                      regressor => self%regressor)

                call regressor%create_X(x)
                X_original = regressor%X

                bootstraps: do i = 1, num_bootstraps
                    call random_number(tmp_real)
                    indices(:) = nint((N-1)*tmp_real) + 1

                    y_selection(:) = y(indices)
                    X_selection(:,:) = X_original(indices,:)

                    call regressor%fit(x_selection, y_selection, .true.)
                    call regressor%predict(x, y_predictions(:,i), y, MSEs(i), R2s(i))

                    betas(:, i) = regressor%beta
                end do bootstraps
            end associate
        end subroutine
end module
