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
        pure function chifactor(nu, moleculeIntCode)
            import :: DP
            implicit none
            real :: chiFactor
            real(kind=DP), intent(in) :: nu 
            integer, intent(in) :: moleculeIntCode
        end function chifactor
    end interface

end module Interfaces
