module Interfaces
    implicit none
    abstract interface
        real function shape(nu)
            integer, parameter :: DP = selected_real_kind(15, 307)
            real(kind=DP), intent(in) :: nu 
        end function shape
    end interface
end module Interfaces
