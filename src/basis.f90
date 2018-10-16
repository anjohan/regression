module mod_basis
    use iso_fortran_env, only: dp => real64
    implicit none

    type, public, abstract :: basisfunction
        contains
            procedure(eval_routine), deferred :: eval
    end type

    abstract interface
        function eval_routine(self, x) result(y)
            import basisfunction, dp
            class(basisfunction), intent(in) :: self
            real(dp), intent(in) :: x(:)
            real(dp) :: y
        end function
    end interface
end module
