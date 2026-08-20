__all__ = ["METABOLIC_TYPES"]

METABOLIC_TYPES = {

# ------------ AEROBIC RESPIRATION ------------ 

    "glucose aerobic respiration": {
        "ox": {"C6H12O6": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "glycerol aerobic respiration": {
        "ox": {"C3H8O3": -1, "H2O": -6, "CHO3-1": 3, "H+1": 17, "e-": 14},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "lactate aerobic respiration": {
        "ox": {"C3H5O3-1": -1, "H2O": -6, "CHO3-1": 3, "H+1": 14, "e-": 12},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "pyruvate aerobic respiration": {
        "ox": {"C3H3O3-1": -1, "H2O": -6, "CHO3-1": 3, "H+1": 12, "e-": 10},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "acetate aerobic respiration": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "ethanol aerobic respiration": {
        "ox": {"C2H6O": -1, "H2O": -5, "CHO3-1": 2, "H+1": 14, "e-": 12},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "formate aerobic respiration": {
        "ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "methanol aerobic respiration": {
        "ox": {"CH4O": -1, "H2O": -2, "CHO3-1": 1, "H+1": 7, "e-": 6},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "mannitol aerobic respiration": {
        "ox": {"C6H14O6": -1, "H2O": -12, "CHO3-1": 6, "H+1": 32, "e-": 26},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "lactose aerobic respiration": {
        "ox": {"C12H22O11": -1, "H2O": -25, "CHO3-1": 12, "H+1": 60, "e-": 48},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "formaldehyde aerobic respiration": {
        "ox": {"CH2O": -1, "H2O": -2, "CHO3-1": 1, "H+1": 5, "e-": 4},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "galactose aerobic respiration": {
        "ox": {"C6H12O6 (galactose)": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "fructose aerobic respiration": {
        "ox": {"C6H12O6 (fructose)": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "mannose aerobic respiration": {
        "ox": {"C6H12O6 (mannose)": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "succinate aerobic respiration": {
        "ox": {"C4H4O4-2": -1, "H2O": -8, "CHO3-1": 4, "H+1": 16, "e-": 14},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "maltose aerobic respiration": {
        "ox": {"C12H22O11 (maltose)": -1, "H2O": -25, "CHO3-1": 12, "H+1": 60, "e-": 48},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "sucrose aerobic respiration": {
        "ox": {"C12H22O11 (sucrose)": -1, "H2O": -25, "CHO3-1": 12, "H+1": 60, "e-": 48},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "glutamate aerobic respiration": {
        "ox": {"C5H7O4N-2": -1, "H2O": -11, "CHO3-1": 5, "NH4+1": 1, "H+1": 20, "e-": 18},
        "red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

#------------ LACTATE FERMENTATION ------------ 

    "lactate fermentation": {
        "ox": {"C3H5O3-1": -1, "H2O": -2, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 5, "e-": 4},
        "ox II": {"H2": -1, "H+1": 2, "e-": 2},
        "red": {"CHO3-1": -1, "H+1": -4.5, "e-": -4, "C2H3O2-1": .5, "H2O": 2},
    },

# ------------ METHANOGENIC FERMENTATION ------------ 

    "methanol methanogenic fermentation": {
        "ox": {"CH4O": -1, "H2O": -2, "CHO3-1": 1, "H+1": 7, "e-": 6},
        "red": {"CH4O": -1, "H+1": -2, "e-": -2, "CH4": 1, "H2O": 1}
    },

    "acetate methanogenic fermentation": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "red": {"C2H3O2-1": -1, "H+1": -9, "e-": -8, "CH4": 2, "H2O": 2}
    },

# ------------ AUTOTROPHIC METHANOGENESIS & ACETOGENESIS ------------ 

    "autotrophic methanogenesis": {
        "ox": {"H2": -1, "H+1": 2, "e-": 2},
        "red": {"CHO3-1": -1, "H+1": -9, "e-": -8, "CH4": 1, "H2O": 3}
    },

    "autotrophic acetogenesis": {
        "ox": {"H2": -1, "H+1": 2, "e-": 2},
        "red": {"CHO3-1": -1, "H+1": -4.5, "e-": -4, "C2H3O2-1": .5, "H2O": 2}
    },

# ------------ DENITRIFICATION ------------ 

    "acetate nitrate denitrification": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "red": {"NO3-1": -1, "H+1": -6, "e-": -5, "N2": .5, "H2O": 3}
    },

    "acetate nitrite denitrification": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "red": {"NO2-1": -1, "H+1": -4, "e-": -3, "N2": .5, "H2O": 2}
    },

    "formate nitrate denitrification": {
        "ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2},
        "red": {"NO3-1": -1, "H+1": -6, "e-": -5, "N2": .5, "H2O": 3}
    },

    "formate nitrite denitrification": {
        "ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2},
        "red": {"NO2-1": -1, "H+1": -4, "e-": -3, "N2": .5, "H2O": 2}
    },

# ------------ AMMONIFICATION ------------ 

    "lactate nitrate ammonification": {
        "ox": {"C3H5O3-1": -1, "H2O": -2, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 5, "e-": 4}, 
        "red": {"NO3-1": -1, "H+1": -10, "e-": -8, "NH4+1": 1, "H2O": 3}
    },

    "lactate nitrite ammonification": {
        "ox": {"C3H5O3-1": -1, "H2O": -2, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 5, "e-": 4}, 
        "red": {"NO2-1": -1, "H+1": -8, "e-": -6, "NH4+1": 1, "H2O": 2}
    },

    "acetate nitrate ammonification": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8}, 
        "red": {"NO3-1": -1, "H+1": -10, "e-": -8, "NH4+1": 1, "H2O": 3}
    },

    "acetate nitrite ammonification": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8}, 
        "red": {"NO2-1": -1, "H+1": -8, "e-": -6, "NH4+1": 1, "H2O": 2}
    },

    "formate nitrate ammonification": {
        "ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2}, 
        "red": {"NO3-1": -1, "H+1": -10, "e-": -8, "NH4+1": 1, "H2O": 3}
    },

    "formate nitrite ammonification": {
        "ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2}, 
        "red": {"NO2-1": -1, "H+1": -8, "e-": -6, "NH4+1": 1, "H2O": 2}
    },

# ------------ SULFUR RESPIRATION ------------

    "acetate sulfur respiration": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "red": {"S": -1, "H+1": -1, "e-": -2, "HS-1": 1}
    },

    "acetate sulfate respiration": {
        "ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "red": {"SO4-2": -1, "H+1": -9, "e-": -8, "HS-1": 1, "H2O": 4}
    },

    "ethanol sulfate respiration": {
        "ox": {"C2H6O": -1, "H2O": -5, "CHO3-1": 2, "H+1": 14, "e-": 12},
        "red": {"SO4-2": -1, "H+1": -9, "e-": -8, "HS-1": 1, "H2O": 4}
    },

    "ethanol sulfur respiration": {
        "ox": {"C2H6O": -1, "H2O": -5, "CHO3-1": 2, "H+1": 14, "e-": 12},
        "red": {"S": -1, "H+1": -1, "e-": -2, "HS-1": 1}
    },

# ------------ GLUCOSE FERMENTATION ------------

    # "glucose fermentation (I)": {

    #     "ox": {"C6H12O6": -1, "C3H3O3-1": 2, "H+1": 6, "e-": 4},

    #     "red": {

    #         "C2H3O2-1": {"C3H3O3-1": -1, "H2O": -1, "CHO2-1": 1, "C2H3O2-1": 1, "H+1": 1},
    #         "C2H6O": {"C3H3O3-1": -1, "H2O": -1, "H+1": -2, "e-": -2, "C2H6O": 1, "CHO3-1": 1},
    #         "C3H5O3-1": {"C3H3O3-1": -1, "H+1": -2, "e-": -2, "C3H5O3-1": 1},
    #         "C3H8O3": {"C3H3O3-1": -1, "H+1": -5, "e-": -4, "C3H8O3": 1},
    #         "C4H4O4-2": {"C3H3O3-1": -1, "CHO3-1": -1, "H+1": -4, "e-": -4, "C4H4O4-2": 1, "H2O": 2}
    #     }
    # },

    # "glucose fermentation (II)": {

    #     "ox": {"C6H12O6": -1, "C21H32N7O16P3S": -2, "H2O": -2, "C23H34N7O17P3S": 2, "CHO3-1": 2, "H+1": 10, "e-": 8},

    #     "red": {
    #         "C2H3O2-1": {"C23H34N7O17P3S": -1, "H2O": -1, "C2H3O2-1": 1, "H+1": 1, "C21H32N7O16P3S": 1},
    #         "C4H7O2-1": {"C23H34N7O17P3S": -1, "H+1": -1.5, "e-": -2, "C4H7O2-1": .5, "C21H32N7O16P3S": 1},
    #         "H2": {"H+1": -1, "e-": -1, "H2": .5},
    #         "C4H10O": {"C23H34N7O17P3S": -1, "H+1": -4, "e-": -4, "C4H10O": .5, "H2O": .5, "C21H32N7O16P3S": 1},
    #         "C2H6O": {"C23H34N7O17P3S": -1, "H+1": -4, "e-": -4, "C2H6O": 1, "C21H32N7O16P3S": 1}
    #     }
    # },


    "glucose fermentation": {

        "ox I": {"C6H12O6": -1, "C3H3O3-1": 2, "H+1": 6, "e-": 4},

        "red I": {

            "CHO2-1": {"C3H3O3-1": -1, "H2O": -1, "CHO2-1": 1, "C2H3O2-1": 1, "H+1": 1},
            "C2H6O": {"C3H3O3-1": -1, "H2O": -1, "H+1": -2, "e-": -2, "C2H6O": 1, "CHO3-1": 1},
            "C3H5O3-1": {"C3H3O3-1": -1, "H+1": -2, "e-": -2, "C3H5O3-1": 1},
            "C3H8O3": {"C3H3O3-1": -1, "H+1": -5, "e-": -4, "C3H8O3": 1},
            "C4H4O4-2": {"C3H3O3-1": -1, "CHO3-1": -1, "H+1": -4, "e-": -4, "C4H4O4-2": 1, "H2O": 2}
        },

        "ox II": {"C6H12O6": -1, "C21H32N7O16P3S": -2, "H2O": -2, "C23H34N7O17P3S": 2, "CHO3-1": 2, "H+1": 10, "e-": 8},

        "red II": {
            "C2H3O2-1": {"C23H34N7O17P3S": -1, "H2O": -1, "C2H3O2-1": 1, "H+1": 1, "C21H32N7O16P3S": 1},
            "C4H7O2-1": {"C23H34N7O17P3S": -1, "H+1": -1.5, "e-": -2, "C4H7O2-1": .5, "C21H32N7O16P3S": 1},
            "H2": {"H+1": -1, "e-": -1, "H2": .5},
            "C4H10O": {"C23H34N7O17P3S": -1, "H+1": -4, "e-": -4, "C4H10O": .5, "H2O": .5, "C21H32N7O16P3S": 1},
            "C2H6O": {"C23H34N7O17P3S": -1, "H+1": -4, "e-": -4, "C2H6O": 1, "C21H32N7O16P3S": 1}
        }
    },


#------------ OVERFLOW PATHWAYS ------------

    "overflow": {

        "C6H12O6": {

            "ox": {

                "C2H3O2-1": {"C6H12O6": -1, "H2O": -4, "C2H3O2-1": 2, "CHO3-1": 2, "H+1": 12, "e-": 8},
                "C2H6O": {"C6H12O6": -1, "C3H3O3-1": 2, "H+1": 6, "e-": 4},

            },

            "red": {

                "C2H6O": {"C3H3O3-1": -1, "H2O": -1, "H+1": -2, "e-": -2, "C2H6O": 1, "CHO3-1": 1},
                "C3H8O3": {"C3H3O3-1": -1, "H+1": -5, "e-": -4, "C3H8O3": 1},
                "C4H4O4-2": {"C3H3O3-1": -1, "CHO3-1": -1, "H+1": -4, "e-": -4, "C4H4O4-2": 1, "H2O": 2}
            },     
        },

        "C6H12O6 (galactose)": {

            "ox": {
            
                "C2H3O2-1": {"C6H12O6 (galactose)": -1, "H2O": -4, "C2H3O2-1": 2, "CHO3-1": 2, "H+1": 12, "e-": 8}
            
            }
        },

        "C3H8O3": {

            "ox": {
                        
                "C2H3O2-1": {"C3H8O3":   -1, "H2O": -2, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 8, "e-": 6}
                        
            }
        },

        "C3H5O3-1": {
        
            "ox": {
                        
                "C2H3O2-1": {"C3H5O3-1": -1, "H2O": -2, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 5, "e-": 4}
                        
            }   
        },

        "C3H3O3-1": {
        
            "ox": {
                        
                "C2H3O2-1": {"C3H3O3-1": -1, "H2O": -2, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 3, "e-": 2}
                        
            }   
        },

        "C2H6O": {
        
            "ox": {
                        
                "C2H3O2-1": {"C2H6O": -1, "H2O": -1, "C2H3O2-1": 1, "H+1": 5, "e-": 4}
                        
            }   
        },       

    }

}