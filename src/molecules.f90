module molecule_vars
    ! variables related to the molecules
    use globals
    implicit none
    integer :: MOTYPE ! molecule type-integer: '1' - for SO2 and H2O, '2' -- for CO2
    integer :: NLIN ! number of lines in one HITRAN file (e.g. H16.01)
    character(len=20) :: LINE_PATH ! file names with HITRAN raw data (one file per molecule)
end module molecule_vars

module van_fleck_huber_vars
    implicit none
    integer :: JM1 ! likely just a loop index
    real :: VIVI ! modified wavenumber within the loop
    real :: EVV, EVV_ ! variables related to the exponential terms in the Van Vleck-Weisskopf-Huber factor calculations
    real :: FACTV ! actual Van-Vleck factor applied
end module van_fleck_huber_vars

module molar_masses
    implicit none
    ! Molecular masses (isotopes as in GAMACHE's 2017data for 12 gases)
    ! Double-check this !
    real, parameter :: WISO_H2O(9) = [ &
        18.010565, 20.014811, 19.014780, 19.016740, 21.020985, 20.020956, 20.0, 22.0, 21.0 &
    ]
    real, parameter :: WISO_CO2(13) = [ &
        43.989830, 44.993185, 45.994076, 44.994045, 46.997431, 45.99740, 47.998322, 46.998291, &
        46.0, 49.001675, 48.0, 47.0, 44.0 &
    ]
    ! why duplicates ??
    real, parameter :: WISO_O3(18) = [ & 
        47.984745, 49.988991, 49.988991, 48.988960, 48.988960, &
        52., 52., 51., 51., 51., 50., 50., 54., 53., 53., 52., 52., 51. &
    ]
    real, parameter :: WISO_N2O(5) = [ &
        44.001062, 44.998096, 44.998096, 46.005308, 45.005278 &
    ]
    real, parameter :: WISO_CO(9) = [ &
        27.994915, 28.998270, 29.999161, 28.999130, 31.002516, 30.002485, &
        30., 32., 33. &
    ]
    real, parameter :: WISO_CH4(4) = [ &
        16.031300, 17.034655, 17.037475, 18.040830 &
    ]
    real, parameter :: WISO_O2(6) = [ &
        31.989830, 33.994076, 32.994045, 36., 35., 34. &
    ]
    real, parameter :: WISO_NO(3) = [ &
        29.997989, 30.995023, 32.002234 &
    ]
    real, parameter :: WISO_SO2(2) = [ &
        63.961901, 65.957695 &
    ]
    real, parameter :: WISO_NO2(1) = [ &
        45.992904 &    
    ]
    real, parameter :: WISO_NH3(2) = [ &
        17.026549, 18.023583 &
    ]
    real, parameter :: WISO_HNO3(1) = [ &
        62.995644 &  
    ]
    real, parameter :: WISO_OH(3) = [ &
        17.002740, 19.006986, 18.008915 &
    ]
    real, parameter :: WISO_HF(1) = [ &
        20.006229 &    
    ]
    real, parameter :: WISO_HCl(2) = [ &
        35.976678,37.973729 &
    ]
    real, parameter :: WISO_HBr(2) = [ &
        79.926160, 81.924115 &  
    ]
    real, parameter :: WISO_HI(1) = [ &
        127.912297 &  
    ]
    real, parameter :: WISO_ClO(2) = [ &
        50.963768, 52.960819 &
    ]
    real, parameter :: WISO_OCS(5) = [ &
        59.966986, 61.962780, 60.970341, 60.966371, 61.971231 &
    ]
    real, parameter :: WISO_H2CO(3) = [ &
        30.010565, 31.013920, 32.014811 &
    ]
    real, parameter :: WISO_HOCl(2) = [ &
        51.971593, 53.968644 &  
    ]
    real, parameter :: WISO_N2(1) = [ &
        28.006147 &
    ]
    real, parameter :: WISO_HCN(3) = [ &
        27.010899, 28.014254, 28.007933 &
    ]
    real, parameter :: WISO_CH3Cl(2) = [ &
        49.992328, 51.989379 &
    ]
    real, parameter :: WISO_H2O2(1) = [ &
        34.005480 &
    ]
    real, parameter :: WISO_C2H2(2) = [ &
        26.015650, 27.019005 &
    ]
    real, parameter :: WISO_C2H6(1) = [ &
        30.046950 &
    ]
    real, parameter :: WISO_PH3(1) = [ &
        33.997238 &    
    ]
    real, parameter :: WISO_COF2(1) = [ &
        65.991722 & 
    ]
    real, parameter :: WISO_SF6(1) = [ &
        145.962492 &
    ]
    real, parameter :: WISO_H2S(3) = [ &
        33.987721, 35.983515, 34.987105 &        
    ]
    real, parameter :: WISO_HCOOH(1) = [ &
        46.005480 &
    ]
    real, parameter :: WISO_HO2(1) = [ &
        32.997655 &
    ]
    real, parameter :: WISO_O(1) = [ &
        15.994915 &
    ]
    real, parameter :: WISO_ClONO2(2) = [ &
        96.956672, 98.953723 &
    ]
    real, parameter :: WISO_NO_plus(1) = [ &
        29.997989 &     
    ]
    real, parameter :: WISO_HOBr(2) = [ &
        95.921076, 97.919027 &
    ]
    real, parameter :: WISO_C2H4(2) = [ &
        28.031300, 29.034655 &
    ]
    real, parameter :: WISO_CH3OH(1) = [ &
        32.026215 &
    ]
    real, parameter :: WISO_CH3Br(2) = [ &
        93.941811, 95.939764 &
    ]
    real, parameter :: WISO_CH3CN(1) = [ &
        41.026549 &
    ]
    real, parameter :: WISO_CF4(1) = [ &
        87.993616 &
    ]
    
    ! 125 len array in the legacy code but 124 factually
    real, parameter :: WISO(124) = [ &
        WISO_H2O, WISO_CO2, WISO_O3, WISO_N2O, WISO_CO, WISO_CH4, WISO_O2, WISO_NO, WISO_SO2, WISO_NO2, WISO_NH3, &
        WISO_HNO3, WISO_OH, WISO_HCl, WISO_HBr, WISO_HI, WISO_ClO, WISO_OCS, WISO_H2CO, WISO_HOCl, WISO_N2, WISO_HCN, &
        WISO_CH3Cl, WISO_H2O2, WISO_C2H2, WISO_C2H6, WISO_PH3, WISO_COF2, WISO_SF6, WISO_H2S, WISO_HCOOH, WISO_HO2, &
        WISO_O, WISO_ClONO2, WISO_NO_plus, WISO_HOBr, WISO_C2H4, WISO_CH3OH, WISO_CH3Br, WISO_CH3CN, WISO_CF4 &
    ]
end module molar_masses
