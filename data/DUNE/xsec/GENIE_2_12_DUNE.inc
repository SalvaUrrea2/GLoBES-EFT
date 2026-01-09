
$version="3.2.13"

/* Cross sections */
cross(#CC)<
        @cross_file = "./xsec/xsec_cc_Oxygen_LQCD_LL_divided_by_16.txt"
        @eft_coeff_file = "./xsec/new_detection_coefficients/DetAllCCQEOxygenLQCDS1Correct2.txt"
        @quark_flavors = 0 : 0
>

cross(#CCtau)<
        @cross_file = "./xsec/xsec_cc_Oxygen_LQCD_LL_divided_by_16.txt" //They have set the tau cross section to 0 in the other file
        @eft_coeff_file = "./xsec/new_detection_coefficients/DetAllCCQEOxygenLQCDS1Correct2.txt"
        @quark_flavors = 0 : 0
>

cross(#NC)<
        @cross_file = "./xsec/xsec_nc.dat"
        @eft_coeff_file = "./xsec/DetAllCCQEOxygen_energy_preserved.txt"
        @quark_flavors = 0 : 0
>


