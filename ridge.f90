module mod_ridge
    use iso_fortran_env, only: dp => real64
    use mod_regressor
    implicit none

    interface ridge
        procedure :: init_ridge
    end interface

    type, public, extends(regressor) :: ridge
        real(dp) :: lambda

        contains
            procedure :: fit
    end type ridge

    contains
        function init_ridge(lambda, basis, X) result(self)
            real(dp), intent(in) :: lambda
            class(basisfunction), intent(in), optional :: basis(:)
            real(dp), intent(in), optional :: X(:,:)

            type(ridge) :: self

            self%lambda = lambda

            if (present(basis)) allocate(self%basis, source=basis)
            if (present(X)) self%X = X
        end function

        subroutine fit(self, x_values, y_values)
            class(ridge), intent(inout) :: self
            real(dp), intent(in) :: y_values(:)
            real(dp), intent(in), optional :: x_values(:,:)

            integer :: N, p, info, i
            real(dp), allocatable :: X_T(:,:), A(:,:), b(:)
            real(dp) :: lambda

            if (present(x_values)) call self%create_X(x_values)

            N = size(y_values)
            p = size(self%basis)
            lambda = self%lambda

            allocate(A(p,p))

            associate(X => self%X)
                X_T = transpose(X)
                !call dgemm('T', 'N', p, p, N, 1.0d0, X, N, X, N, 0.0d0, A, p)
                A = matmul(X_T, X)
            end associate

            do i = 1, p
                A(i,i) = A(i,i) + lambda
            end do

            b = matmul(X_T, y_values)

            call dposv('L', p, 1, A, p, b, p, info)

            if (info /= 0) then
                write(*, "(a,i0)") "ERROR in dposv: info = ", info
            end if

            self%beta = b
        end subroutine
end module
