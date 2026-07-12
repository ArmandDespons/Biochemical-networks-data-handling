METABOLIC_TYPES = {

# ------------ AEROBIC RESPIRATION ------------ 

    "glucose aerobic respiration": {
        "Ox": {"C6H12O6": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "glucose aerobic respiration & fermentation (E. coli overflow)": {
        "Ox": {"C6H12O6": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "Ox II": {"C6H12O6": -1, "H2O": -4, "C2H3O2-1": 2, "CHO3-1": 2, "H+1": 12, "e-": 8},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2},
    },

    "glycerol aerobic respiration": {
        "Ox": {"C3H8O3": -1, "H2O": -6, "CHO3-1": 3, "H+1": 17, "e-": 14},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "lactate aerobic respiration": {
        "Ox": {"C3H5O3-1": -1, "H2O": -6, "CHO3-1": 3, "H+1": 14, "e-": 12},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "acetate aerobic respiration": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "ethanol aerobic respiration": {
        "Ox": {"C2H6O": -1, "H2O": -5, "CHO3-1": 2, "H+1": 14, "e-": 12},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "formate aerobic respiration": {
        "Ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "methanol aerobic respiration": {
        "Ox": {"CH4O": -1, "H2O": -2, "CHO3-1": 1, "H+1": 7, "e-": 6},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "mannitol aerobic respiration": {
        "Ox": {"C6H14O6": -1, "H2O": -12, "CHO3-1": 6, "H+1": 32, "e-": 26},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "lactose aerobic respiration": {
        "Ox": {"C12H22O11 (lactose)": -1, "H2O": -25, "CHO3-1": 12, "H+1": 60, "e-": 48},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "formaldehyde aerobic respiration": {
        "Ox": {"CH2O": -1, "H2O": -2, "CHO3-1": 1, "H+1": 5, "e-": 4},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "galactose aerobic respiration": {
        "Ox": {"C6H12O6 (galactose)": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "fructose aerobic respiration": {
        "Ox": {"C6H12O6 (fructose)": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "mannose aerobic respiration": {
        "Ox": {"C6H12O6 (mannose)": -1, "H2O": -12, "CHO3-1": 6, "H+1": 30, "e-": 24},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "succinate aerobic respiration": {
        "Ox": {"C4H4O4-2": -1, "H2O": -8, "CHO3-1": 4, "H+1": 16, "e-": 14},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "maltose aerobic respiration": {
        "Ox": {"C12H22O11 (maltose)": -1, "H2O": -25, "CHO3-1": 12, "H+1": 60, "e-": 48},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "sucrose aerobic respiration": {
        "Ox": {"C12H22O11 (sucrose)": -1, "H2O": -25, "CHO3-1": 12, "H+1": 60, "e-": 48},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

    "glutamate aerobic respiration": {
        "Ox": {"C5H7O4N-2": -1, "H2O": -11, "CHO3-1": 5, "NH4+1": 1, "H+1": 20, "e-": 18},
        "Red": {"O2": -1, "H+1": -4, "e-": -4, "H2O": 2}
    },

#------------ LACTATE FERMENTATION ------------ 

    "lactate fermentation": {
        "Ox": {"C3H5O3-1": -1, "H2O": -2, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 5, "e-": 4},
        "Ox II": {"H2": -1, "H+1": 2, "e-": 2},
        "Red": {"CHO3-1": -1, "H+1": -4.5, "e-": -4, "C2H3O2-1": .5, "H2O": 2},
        "Red II": {"H2": -.5, "CHO3-1": -1, "H+1": -3.5, "e-": -3, "C2H3O2-1": .5, "H2O": 2}
    },

# ------------ METHANOGENIC FERMENTATION ------------ 

    "methanol methanogenic fermentation": {
        "Ox": {"CH4O": -1, "H2O": -2, "CHO3-1": 1, "H+1": 7, "e-": 6},
        "Red": {"CH4O": -1, "H+1": -2, "e-": -2, "CH4": 1, "H2O": 1}
    },

    "acetate methanogenic fermentation": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "Red": {"C2H3O2-1": -1, "H+1": -9, "e-": -8, "CH4": 2, "H2O": 2}
    },

# ------------ AUTOTROPHIC METHANOGENESIS & ACETOGENESIS ------------ 

    "autotrophic methanogenesis": {
        "Ox": {"H2": -1, "H+1": 2, "e-": 2},
        "Red": {"CHO3-1": -1, "H+1": -9, "e-": -8, "CH4": 1, "H2O": 3}
    },

    "autotrophic acetogenesis": {
        "Ox": {"H2": -1, "H+1": 2, "e-": 2},
        "Red": {"CHO3-1": -1, "H+1": -4.5, "e-": -4, "C2H3O2-1": .5, "H2O": 2}
    },

# ------------ DENITRIFICATION ------------ 

    "acetate-nitrate denitrification": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "Red": {"NO3-1": -1, "H+1": -6, "e-": -5, "N2": .5, "H2O": 3}
    },

    "acetate-nitrite denitrification": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "Red": {"NO2-1": -1, "H+1": -4, "e-": -3, "N2": .5, "H2O": 2}
    },

    "formate-nitrate denitrification": {
        "Ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2},
        "Red": {"NO3-1": -1, "H+1": -6, "e-": -5, "N2": .5, "H2O": 3}
    },

    "formate-nitrite denitrification": {
        "Ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2},
        "Red": {"NO2-1": -1, "H+1": -4, "e-": -3, "N2": .5, "H2O": 2}
    },

# ------------ AMMONIFICATION ------------ 

    "lactate-nitrate ammonification": {
        "Ox": {"C3H5O3-1": -1, "H2O": -3, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 5, "e-": 4}, 
        "Red": {"NO3-1": -1, "H+1": -10, "e-": -8, "NH4+1": 1, "H2O": 3}
    },

    "lactate-nitrite ammonification": {
        "Ox": {"C3H5O3-1": -1, "H2O": -3, "C2H3O2-1": 1, "CHO3-1": 1, "H+1": 5, "e-": 4}, 
        "Red": {"NO2-1": -1, "H+1": -8, "e-": -6, "NH4+1": 1, "H2O": 2}
    },

    "acetate-nitrate ammonification": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8}, 
        "Red": {"NO3-1": -1, "H+1": -10, "e-": -8, "NH4+1": 1, "H2O": 3}
    },

    "acetate-nitrite ammonification": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8}, 
        "Red": {"NO2-1": -1, "H+1": -8, "e-": -6, "NH4+1": 1, "H2O": 2}
    },

    "formate-nitrate ammonification": {
        "Ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2}, 
        "Red": {"NO3-1": -1, "H+1": -10, "e-": -8, "NH4+1": 1, "H2O": 3}
    },

    "formate-nitrite ammonification": {
        "Ox": {"CHO2-1": -1, "H2O": -1, "CHO3-1": 1, "H+1": 2, "e-": 2}, 
        "Red": {"NO2-1": -1, "H+1": -8, "e-": -6, "NH4+1": 1, "H2O": 2}
    },

# ------------ SULFUR RESPIRATION ------------

    "acetate-sulfur respiration": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "Red": {"S": -1, "H+1": -1, "e-": -2, "HS-1": 1}
    },

    "acetate-sulfate respiration": {
        "Ox": {"C2H3O2-1": -1, "H2O": -4, "CHO3-1": 2, "H+1": 9, "e-": 8},
        "Red": {"SO4-2": -1, "H+1": -9, "e-": -8, "HS-1": 1, "H2O": 4}
    },

    "ethanol-sulfate respiration": {
        "Ox": {"C2H6O": -1, "H2O": -5, "CHO3-1": 2, "H+1": 14, "e-": 12},
        "Red": {"SO4-2": -1, "H+1": -9, "e-": -8, "HS-1": 1, "H2O": 4}
    },

    "ethanol-sulfur respiration": {
        "Ox": {"C2H6O": -1, "H2O": -5, "CHO3-1": 2, "H+1": 14, "e-": 12},
        "Red": {"S": -1, "H+1": -1, "e-": -2, "HS-1": 1}
    },

# ------------ GLUCOSE FERMENTATION ------------

    "glucose fermentation": {

        "acetyl-coa": {

            "Ox": {"C6H12O6": -1, "C21H32N7O16P3S": -2, "H2O": -2, "C23H34N7O17P3S": 2, "CHO3-1": 2, "H+1": 10, "e-": 8},

            "Red": {
                "C2H3O2-1": {"C23H34N7O17P3S": -1, "H2O": -1, "C2H3O2-1": 1, "H+1": 1, "C21H32N7O16P3S": 1},
                "C4H7O2-1": {"C23H34N7O17P3S": -1, "H+1": -1.5, "e-": -2, "C4H7O2-1": .5, "C21H32N7O16P3S": 1},
                "H2": {"H+1": -1, "e-": -1, "H2": .5},
                "C4H10O": {"C23H34N7O17P3S": -1, "H+1": -4, "e-": -4, "C4H10O": .5, "H2O": .5, "C21H32N7O16P3S": 1},
                "C2H6O": {"C23H34N7O17P3S": -1, "H+1": -4, "e-": -4, "C2H6O": 1, "C21H32N7O16P3S": 1}
            }
        },

        "pyruvate": {

            "Ox": {"C6H12O6": -1, "C3H3O3-1": 2, "H+1": 6, "e-": 4},

            "Red": {

                "C2H3O2-1": {"C3H3O3-1": -1, "H2O": -1, "CHO2-1": 1, "C2H3O2-1": 1, "H+1": 1},
                "C2H6O": {"C3H3O3-1": -1, "H2O": -1, "H+1": -2, "e-": -2, "C2H6O": 1, "CHO3-1": 1},
                "C3H5O3-1": {"C3H3O3-1": -1, "H+1": -2, "e-": -2, "C3H5O3-1": 1},
                "C3H8O3": {"C3H3O3-1": -1, "H+1": -5, "e-": -4, "C3H8O3": 1},
                "C4H4O4-2": {"C3H3O3-1": -1, "CHO3-1": -1, "H+1": -4, "e-": -4, "C4H4O4-2": 1, "H2O": 2}
            }
        }
    }
}