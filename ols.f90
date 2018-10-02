module mod_ols
    use mod_regressor
    use iso_fortran_env, only: dp => real64
    implicit none

    interface ols
        procedure :: init_ols
    end interface

    type, public, extends(regressor) :: ols
        contains
            procedure :: fit
    end type ols

    contains
        function init_ols(basis, X) result(self)
            class(basisfunction), intent(in), optional :: basis(:)
            real(dp), intent(in), optional :: X(:,:)
            type(ols) :: self

            if (present(basis)) allocate(self%basis, source=basis)
            if (present(X)) self%X = X

            self%method = "ols"

        end function

        subroutine fit(self, x_values, y_values)
            class(ols), intent(inout) :: self
            real(dp), intent(in) :: y_values(:)
            real(dp), intent(in), optional :: x_values(:,:)

            integer :: lwork, info, N, p
            real(dp), allocatable :: work(:), X(:,:), y(:)

            if (present(x_values)) call self%create_X(x_values)

            lwork = -1
            allocate(work(1))

            N = size(y_values)
            p = size(self%basis)

            ! copies, dgels modifies
            X = self%x
            y = y_values

            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work)
            allocate(work(lwork))
            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)

            if (info /= 0) then
                write(*,"(a,i0)") "DGELS problem: info = ", info
                error stop
            end if

            self%beta = y(1:p)
        end subroutine
end module
