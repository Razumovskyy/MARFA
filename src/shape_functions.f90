module shape_functions
    implicit none
    abstract interface
        real function shape_func(x)
            real, intent(in) :: x
        end function shape_func
    end interface
end module shape_functions