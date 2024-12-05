module Interfaces
    implicit none
    
    
    abstract interface
        real function shape(nu)
            implicit none
            integer, parameter :: DP = selected_real_kind(15, 307)
            real(kind=DP), intent(in) :: nu 
        end function shape
    end interface

    
    abstract interface
        real function chifactor(nu, moleculeIntCode)
            implicit none
            integer, parameter :: DP = selected_real_kind(15, 307)
            real(kind=DP), intent(in) :: nu 
            integer, intent(in) :: moleculeIntCode
        end function chifactor
    end interface

end module Interfaces
