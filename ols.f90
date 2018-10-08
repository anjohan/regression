module mod_ols
    use mod_regressor
    use mod_utilities, only: check_info
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

            self%method = "OLS"

        end function

        subroutine fit(self, x_values, y_values)
            class(ols), intent(inout) :: self
            real(dp), intent(in) :: y_values(:)
            real(dp), intent(in), optional :: x_values(:,:)

            integer :: lwork, info, N, p, rank
            real(dp) :: rcond
            real(dp), allocatable :: work(:), X(:,:), y(:), s(:)

            if (present(x_values)) call self%create_X(x_values)

            lwork = -1
            allocate(work(1))

            N = size(y_values)
            p = size(self%basis)

            ! copies, dgels modifies
            X = self%X
            y = y_values

            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)
            lwork = nint(work(1))
            deallocate(work)
            allocate(work(lwork))
            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)

            call check_info(info, "DGELS")

            !allocate(s(p))
            !rcond = -1.0d0

            !call dgelss(N, p, 1, X, N, y, N, s, rcond, rank, work, lwork, info)
            !lwork = nint(work(1))
            !deallocate(work)
            !allocate(work(lwork))
            !call dgelss(N, p, 1, X, N, y, N, s, rcond, rank, work, lwork, info)

            !call check_info(info, "DEGELSS")

            self%beta = y(1:p)
        end subroutine
end module
