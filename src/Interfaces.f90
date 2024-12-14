module Interfaces
    use Constants
    implicit none
    
    
    abstract interface
        real function shape(nu)
            import :: DP
            implicit none
            real(kind=DP), intent(in) :: nu 
        end function shape
    end interface

    
    abstract interface
        real function chifactor(nu, moleculeIntCode)
            import :: DP
            implicit none
            real(kind=DP), intent(in) :: nu 
            integer, intent(in) :: moleculeIntCode
        end function chifactor
    end interface

end module Interfaces
