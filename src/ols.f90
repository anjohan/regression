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

            integer :: lwork, info, N, p, rank, liwork
            real(dp) :: rcond
            real(dp), allocatable :: work(:), X(:,:), y(:), s(:)
            integer, allocatable :: iwork(:)

            if (present(x_values)) call self%create_X(x_values)

            lwork = -1
            allocate(work(1), iwork(1))

            N = size(y_values)
            p = size(self%basis)

            ! copies, dgels modifies
            X = self%X

            if (N > p) then
                y = y_values
            else
                allocate(y(p))
                y(1:N) = y_values
            end if

!            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)
!            lwork = nint(work(1))
!            deallocate(work)
!            allocate(work(lwork))
!            call dgels("N", N, p, 1, X, N, y, N, work, lwork, info)
!
!            call check_info(info, "DGELS")

            allocate(s(p))
            rcond = -1.0d0

!            call dgelss(N, p, 1, X, N, y, N, s, rcond, rank, work, lwork, info)
            call dgelsd(N, p, 1, X, N, y, max(N,p), s, rcond, rank, work, lwork, iwork, info)
            lwork = nint(work(1))
            deallocate(work)
            allocate(work(lwork))
            liwork = iwork(1)
            deallocate(iwork)
            allocate(iwork(liwork))
            call dgelsd(N, p, 1, X, N, y, max(N,p), s, rcond, rank, work, lwork, iwork, info)
!            call dgelss(N, p, 1, X, N, y, N, s, rcond, rank, work, lwork, info)

            call check_info(info, "DGELSD")

            self%beta = y(1:p)
        end subroutine
end module
