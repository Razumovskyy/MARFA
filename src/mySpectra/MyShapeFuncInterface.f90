module MyShapeFuncInterface
    implicit none
    abstract interface
        real function shape(nu)
            integer, parameter :: DP = selected_real_kind(15,  307)
            real(kind=DP), intent(in) :: nu ! ! cm-1, (gridWV) -- spectral point in which the total contribution from lines is calculated
        end function shape
    end interface
end module MyShapeFuncInterface