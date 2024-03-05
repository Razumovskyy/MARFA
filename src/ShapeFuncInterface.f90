module ShapeFuncInterface
    implicit none
    abstract interface
        real function shape(X)
            ! X -- distance from the line center in cm**-1
            real, intent(in) :: X
        end function shape
    end interface
end module ShapeFuncInterface
