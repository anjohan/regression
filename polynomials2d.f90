module mod_polynomials2d
    use iso_fortran_env, only: dp => real64
    use mod_basis
    implicit none

    interface polynomial2d
        procedure :: init_polynomial2d
    end interface

    type, public, extends(basis_function) :: polynomial2d
        integer :: x1_pow, x2_pow

        contains
            procedure :: eval
    end type

    contains
        function init_polynomial2d(x1_pow, x2_pow) result(polynomial)
            integer, intent(in) :: x1_pow, x2_pow
            type(polynomial2d) :: polynomial

            polynomial%x1_pow = x1_pow
            polynomial%x2_pow = x2_pow
        end function

        function eval(self, x) result(y)
            class(polynomial2d), intent(in) :: self
            real(dp), intent(in) :: x(:)
            real(dp) :: y

            y = x(1)**self%x1_pow * x(2)**self%x2_pow
        end function

        subroutine create_basis(basis, d)
            class(polynomial2d), allocatable, intent(inout) :: basis(:)
            integer, intent(in) :: d

            integer :: p, idx, i, j
            p = (d+1)*(d+2)/2

            if (allocated(basis)) deallocate(basis)
            allocate(basis(p))

            idx = 0
            do j = 0, d
                do i = 0, d-j
                    idx = idx + 1
                    basis(idx)%x1_pow = i
                    basis(idx)%x2_pow = j
                end do
            end do
        end subroutine
end module
