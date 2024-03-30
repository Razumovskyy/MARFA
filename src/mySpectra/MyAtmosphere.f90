module MyAtmosphere
    implicit none
    real :: pressure ! atm
    real :: temperature
    real :: density ! [molecule/(cm^2 * km)] -- such density units are needed for having absorption coefficient in km-1
    real :: molarMass ! g/mol -- molar mass of the treated gaseous species
end module MyAtmosphere