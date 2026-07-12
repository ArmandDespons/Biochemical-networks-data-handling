import re

import numpy as np 


from pandas import DataFrame
from ast import literal_eval
from collections.abc import Iterable
from numbers import Integral, Real

from .misc import _matrix_parser, _check_coeffs

from typing import List, Dict, Union, Optional

from .Species import Species

__all__ = ["Equation", "ConservationLaws"]


class Equation:

    def __init__(
        self, 
        reactants: List[Union[str, Species]], 
        products: List[Union[str, Species]], 
        known_coefficients: Optional[Dict[str, float]] = None
    ):
        
        if known_coefficients is None:
            known_coefficients = {}


        # 1. Convert to Species objects
        self.reactants = [s if isinstance(s, Species) else Species(s) for s in reactants]
        self.products = [s if isinstance(s, Species) else Species(s) for s in products]
        

        _dummy_coeff = {s.formula: 1 for s in self.reactants + self.products}
        self.stoichiometry = ConservationLaws(self, _coeffs = _dummy_coeff)


        _check_coeffs(known_coefficients, self.reactants, self.products)

        # 2. Store user-defined coefficients
        self.known_coefficients = known_coefficients
        

        # 3. Store defaults coefficients (1.0) for anything not explicitly defined
        self.default_coefficients = {}
        for s in self.reactants + self.products:
            if s.formula not in self.known_coefficients:
                self.default_coefficients[s.formula] = -1.0 if s in self.reactants else 1.0

        self.balanced_coefficients = None 
        self.conservation_laws = ConservationLaws(self)
        
        if self.is_balanced:
            self.balanced_coefficients = self.coefficients
        
        self.result_inference = None
        

    @property
    def known_coeff(self):
        """Access or assign known coefficient(s) for species in this equation."""
        return self.known_coefficients


    @known_coeff.setter
    def known_coeff(self, value):
        """Set one or more known coefficients and refresh derived equation state."""
        if isinstance(value, tuple) and len(value) == 2:
            value = {value[0]: value[1]}

        if not isinstance(value, dict):
            raise TypeError("known_coefficient must be set with a dict or a (species, coefficient) tuple")

        _check_coeffs(value, self.reactants, self.products)

        if self.balanced_coefficients: 
            self.default_coefficients = {s: 1.0 for s in self.balanced_coefficients.keys()}
            self.balanced_coefficients = None

        for species in value.keys():
            self.default_coefficients.pop(species, None)

        _to_delete = []
        for s in self.known_coefficients.keys(): 
            if not s in value.keys(): _to_delete.append(s) 

        for s in _to_delete: 
            self.known_coefficients.pop(s)
            self.default_coefficients.update({s: 1.0})

        _check_coeffs(self.default_coefficients, self.reactants, self.products)

        print(f'Warning: Deleting previously kwown coefficient for species {_to_delete}, now set to default value (+/-)1.')
        del _to_delete

        self.known_coefficients.update(value)
        self.conservation_laws = ConservationLaws(self)


    @classmethod
    def from_string(cls, equation_str: str):
        """
        Initialize an Equation instance from a string. The species involved in the reaction are
        stored as Species objects if theyr have not been previously created.

            - Reactants and products must be spearated by '='
            - Differents species must be separated with ' + ' with spaces
            - Stoichiometric coefficients are optional (set by default to 1)
        
            [Coeffs1] Reactant1 + [Coeffs2] Reactant1 + ... = [Coeffs3] Product1 + [Coeffs4] Product2 + ...

        Examples: 

            NaCl = Na+ + Cl-

            .5 C6H12O6 + 3 O2 = 3 CO2 + 3 H2O

            2 A + 2 B = 2 AB
        """

        if "=" not in equation_str:
            raise ValueError("Equation must contain '='.")
        
        left_side, right_side = equation_str.split("=")

        def parse_side(side_str: str):
            species_list = []
            known_coeffs = {}
            
            # Split by '+' only when it is surrounded by spaces 
            items = [item.strip() for item in re.split(r'\s+\+\s+', side_str.strip())]
            
            for item in items:
                if not item: continue
                
                # Regex strng parsing:
                # ^([\d\.]+)?    -> Optional leading coefficient
                # \s* -> Optional space
                # (              -> Start Formula Group
                #  [A-Za-z0-9\.]+ -> Elements, numbers, and decimals (CH1.8O.5)
                #  (?:[+-]\d*)?  -> Optional charge at the very end (+1, -, +)
                # )
                match = re.match(r'^([\d\.]+)?\s*([A-Za-z0-9\.]+(?:[+-]\d*)?)$', item)
                
                if match:
                    coeff_str = match.group(1)
                    formula = match.group(2)
                    
                    obj = Species(formula)
                    species_list.append(obj)
                    
                    if coeff_str:
                        known_coeffs[formula] = float(coeff_str)
                else:
                    # Fallback for single ions like "H+" if the regex is too strict
                    obj = Species(item)
                    species_list.append(obj)
            
            return species_list, known_coeffs

        reac_objs, reac_known = parse_side(left_side)
        prod_objs, prod_known = parse_side(right_side)
        
        return cls(reac_objs, prod_objs, {**reac_known, **prod_known})


    @classmethod
    def from_dict(cls, dict: Dict[str, float] | Dict[Species, float]):

        reactants, products = [], []
        to_delete = []
        for species, coeff in dict.items():

            if isinstance(coeff, float) or isinstance(coeff, Integral):

                if coeff>0:
                    products.append(species)

                elif coeff<0: 
                    reactants.append(species)

                elif np.isnan(coeff): 
                    to_delete.append(species)

            elif coeff=='-':
                reactants.append(species)
                to_delete.append(species)

            elif coeff=='+':
                products.append(species)
                to_delete.append(species)

            else: 
                raise ValueError(f"Unauthorized value {coeff}. Only integers or floats are allowed for kwown coefficients or '+'/'-' for products/reactants with unknown coefficients")

        for species in to_delete: dict.pop(species)

        return cls(reactants, products, known_coefficients=dict)



    @property
    def is_balanced(self):
        return self.conservation_laws.is_balanced
    

    @property
    def coefficients(self) -> Dict[str, float]:
        """Combines all coefficient sources for calculations."""

        if self.balanced_coefficients:
            return {**self.balanced_coefficients, **self.known_coefficients}

        else: 
            return {**self.default_coefficients, **self.known_coefficients}
        

    def element_recovery(self, element, in_percent=True):

        if not element in self.conservation_laws.list:
            raise ValueError(f'{element} is not listed as an element.')
        
        
        in_reactant = -np.fromiter((self.conservation_laws.get[element][s.formula] for s in self.reactants), dtype=np.float64).sum()
        in_product = np.fromiter((self.conservation_laws.get[element][s.formula] for s in self.products), dtype=np.float64).sum()

        base = 100 if in_percent else 1
        
        return in_product/in_reactant*base


    def infer_coeffs(self, known_coeff: Optional[Dict[str, float]] = {}, verbose=True):
        """
        Solves the linear system A * x = b for the unknown coefficients.
        Accepts an optional dictionary of additional known fluxes/coefficients.
        """

        if self.balanced_coefficients: 
            self.default_coefficients = {s: 1.0 for s in self.balanced_coefficients.keys()}
            self.balanced_coefficients = None

        result_inference = {}

        if known_coeff != {}: self.known_coeff = known_coeff
            

        # Refresh the conservation laws matrix to reflect the new knowns
        # self.conservation_laws = ConservationLaws(self)
            
        
        result_inference['Unknowns coefficients'] = list(self.default_coefficients.keys())
        result_inference['Knowns coefficients'] = self.known_coefficients

        cl = self.conservation_laws
        known_species = list(self.known_coefficients.keys())
        unknown_species = list(self.default_coefficients.keys())

        # Get column indices for knowns and unknowns
        known_indices = [cl.species_to_index[s] for s in known_species]
        unknown_indices = [cl.species_to_index[s] for s in unknown_species]

        L_matrix = cl.L

        if 'Charge' in cl.rows_L.keys():
            L_matrix[:-1, unknown_indices] = np.abs(L_matrix[:-1, unknown_indices])

            for r in self.reactants:

                if r.formula in unknown_species and r.charge!=0: L_matrix[-1, cl.species_to_index[r.formula]] *= -1

        else: 
            L_matrix[:, unknown_indices] = np.abs(L_matrix[:, unknown_indices])

        
        # Construct the LHS matrix by using _matrix_parser that return the row-independent L_matrix whose columns associated 
        # with the known coefficients have been removed
        LHS, rows_LHS = _matrix_parser(L_matrix[:, unknown_indices], cl.rows_L)


        result_inference['LHS'] = LHS
        result_inference['Rows LHS/RHS'] = rows_LHS

        # Construct vector RHS vector (the negative sum of the known columns)
        RHS = -np.sum(L_matrix[:, known_indices], axis=1)[[cl.rows_L[key] for key in rows_LHS.keys()]]

        result_inference['RHS'] = L_matrix[:, known_indices][[cl.rows_L[key] for key in rows_LHS.keys()], :]

        # Check if the matrix is square
        _shape = LHS.shape[0] - LHS.shape[1] 
        

        if _shape<0: 
            
            if verbose: print('Unable to infer coefficients: Under-determined system')

            result_inference['Status'] = 'Failed'
            result_inference['Message'] = f'Under-determined system: not enough known coefficients'
            result_inference['Inferred coefficients'] = None

        if _shape==0:

            sol = np.linalg.solve(LHS, RHS)

            self.balanced_coefficients = dict(zip(self.default_coefficients.keys(), sol))
            self.conservation_laws = ConservationLaws(self)

            result_inference['Status'] = 'Success'
            result_inference['Message'] = None
            result_inference['Inferred coefficients'] = self.balanced_coefficients


            _temp = self.balanced_coefficients.copy()
            _check_coeffs(_temp, self.reactants, self.products)
            for key, val in self.balanced_coefficients.items():
            
                _sign_mismatch = []
                if val!=_temp[key]:
                    _sign_mismatch.append(key)

                if _sign_mismatch!=[]: 
                    msg = f"Warning! Sign mismatch for species: {[i for i in _sign_mismatch]}."
                    result_inference['Message'] = msg
                
                    if verbose: print(msg)

        
        if _shape>0: 

            if verbose: print('Unable to infer coefficients')

            result_inference['Status'] = 'Failed'
            result_inference['Inferred coefficients'] = None


        self.result_inference = result_inference


    def reset_coeffs(self, which='inferred'):
        if not which in ['inferred', 'all']: 
            raise ValueError("'which' should be either 'inferred' to remove inferred coeffcients or 'all' to set all coefficients to their default values (+/- 1)") 

        if which=='inferred':

            if not self.result_inference is None:  
                _new_default = {s: 1 for s in self.result_inference['Inferred coefficients'].keys()}

        if which=='all': 
    
            _new_default = {s.formula: 1 for s in self.reactants+self.products}

            self.known_coefficients = {}

        self.balanced_coefficients = None

        _check_coeffs(_new_default, self.reactants, self.products)
        self.default_coefficients = _new_default

        self.conservation_laws = ConservationLaws(self)


    def list_species(self, which='all'): 

        if not which in ['all', 'reactants', 'products']:
            raise ValueError(f"'which' takes the following values: 'all', 'reactants', 'products'")
        
        if which=='all':
            return [s.formula for s in self.reactants+self.products]
        
        elif which=='reactants':
            return [s.formula for s in self.reactants]
        
        else: 
            return [s.formula for s in self.products]


    # ---- Gibbs free-energy difference of the reaction (DrG) ----

    def gibbs_fe_diff(self, 
                      activity: Optional[Dict[str, float]] = None, 
                      T: Optional[float] = 298.15, 
                      pH: Optional[float] = 7,
                      method = 'eQ pH=0', 
                      verbose = True
                      ): 
        
        if not self.is_balanced and verbose: 
            print('Warning: Computing the Gibbs free-energy difference for unbalanced equation.')
            
        return sum([self.coefficients.get(species.formula)*species.mu(activity=activity, T=T, pH=pH, method=method) for species in self.reactants+self.products])


    def DrG(self, activity: Optional[Dict[str, float]] = None, T: Optional[float] = 298.15, pH: Optional[float] = 7, method = 'eQ pH=0', verbose=True): 

        return self.gibbs_fe_diff(activity=activity, T=T, pH=pH, method=method, verbose=verbose)



    # ---- Enthalpy difference of the reaction (DrH) ---- 

    def enthalpy_diff(self): 
         
        if not self.is_balanced: 
            print('Warning: Computing the Gibbs free-energy difference for unbalanced equation.')

        return sum([self.coefficients.get(species.formula)*species.standard_enthaply for species in self.reactants+self.products])


    def DrH(self): 

        return self.enthalpy_diff()



    def __mul__(self, scalar):
        if not isinstance(scalar, Real) or isinstance(scalar, bool):
            return NotImplemented

        if scalar == 0:
            raise ValueError("Cannot multiply an equation by zero.")

        flip = scalar < 0
        new_reactants = list(self.products) if flip else list(self.reactants)
        new_products = list(self.reactants) if flip else list(self.products)

        if self.balanced_coefficients and self.is_balanced:
            source_coeffs = self.coefficients
        else:
            source_coeffs = self.known_coefficients

        new_known_coefficients = {formula: coeff * scalar for formula, coeff in source_coeffs.items()}

        return Equation(new_reactants, new_products, known_coefficients=new_known_coefficients)


    def __rmul__(self, scalar):
        return self.__mul__(scalar)


    def __add__(self, other):
        if not isinstance(other, Equation):
            return NotImplemented

        def is_fully_determined(eq):
            return bool(eq.balanced_coefficients) and eq.is_balanced

        c1 = self.coefficients if is_fully_determined(self) else self.known_coefficients
        c2 = other.coefficients if is_fully_determined(other) else other.known_coefficients

        species_registry = {s.formula: s for s in self.reactants + self.products + other.reactants + other.products}

        side1 = {s.formula: ('reactant' if s in self.reactants else 'product') for s in self.reactants + self.products}
        side2 = {s.formula: ('reactant' if s in other.reactants else 'product') for s in other.reactants + other.products}

        new_known_coefficients = {}
        new_reactants, new_products = [], []

        for formula in species_registry.keys():
            in1, in2 = formula in c1, formula in c2
            present1, present2 = formula in side1, formula in side2

            if in1 and in2:
                net = c1[formula] + c2[formula]

                if np.isclose(net, 0.0, atol=1e-6):
                    continue

                new_known_coefficients[formula] = net
                (new_reactants if net < 0 else new_products).append(species_registry[formula])
                continue

            if present1 and present2 and side1[formula] != side2[formula]:
                raise ValueError(
                    f"Cannot combine species '{formula}': it is a substrate in one equation and a product "
                    f"in the other, but its coefficient is not known in both equations."
                )

            if in1:
                new_known_coefficients[formula] = c1[formula]
                (new_reactants if c1[formula] < 0 else new_products).append(species_registry[formula])

            elif in2:
                new_known_coefficients[formula] = c2[formula]
                (new_reactants if c2[formula] < 0 else new_products).append(species_registry[formula])

            else:
                side = side1.get(formula, side2.get(formula))
                (new_reactants if side == 'reactant' else new_products).append(species_registry[formula])

        return Equation(new_reactants, new_products, known_coefficients=new_known_coefficients)


    def __repr__(self):
        coefficients = self.coefficients
        LHS, RHS = "", ""

        for formula, coeff in coefficients.items():
            if coeff<0:
                LHS += str(-coeff) + " " + formula + " + "

            else:
                RHS += str(coeff) + " " + formula + " + "

        return LHS[:-3] + " = " + RHS[:-3]


    

        



class ConservationLaws:


    def __init__(self, eq: Equation, _coeffs = None): 
        """
        Construct a nested dictionnary containing the conservation laws of all the element parsed from
        an Equation instance
        """
        _order_element = ["C", "H", "O", "N", "P", "S", "Mg", "Na", "Li", "Cl", "Mn", "Ca", "K", "Fe"]
        _order_dict = {element: index for index, element in enumerate(_order_element)}


        # 1. Identify all unique species and elements
        all_species = eq.reactants + eq.products
        all_elements = set()
        for s in all_species:
            all_elements.update(s.stoich.keys())
        
        all_elements = sorted(all_elements, key=lambda x: _order_dict.get(x, len(_order_element)))


        # 2. Initialize the nested structure
        self.get = {el: {} for el in all_elements}
        self.get['Charge'] = {}

        coeffs = eq.coefficients if _coeffs is None else _coeffs

        for element in all_elements:

            is_empty = True
            for s in all_species:
                if s.stoich.get(element, 0) and is_empty:
                    is_empty = False

                self.get[element][s.formula] = coeffs[s.formula] * s.stoich.get(element, 0) 

            if is_empty:
                self.get.pop(element)
        
        is_empty = True
        for s in all_species:
            if s.charge and is_empty:
                    is_empty = False
                    
            self.get['Charge'][s.formula] = coeffs[s.formula] * s.charge 

        if is_empty:
                self.get.pop('Charge')
        
        self.list = list(self.get.keys())

        self.element_to_index = dict(zip(self.list, [i for i in range(len(self.list))]))
        self.species_to_index = dict(zip([s.formula for s in all_species], [i for i in range(len(all_species))]))

        self.L,  self.rows_L = _matrix_parser(self.to_array, self.element_to_index)


    @property
    def to_array(self):
        out = np.zeros((len(self.list), len(self.species_to_index)), dtype=float)

        for element in self.list:
            for species in self.species_to_index.keys():
                out[self.element_to_index[element], self.species_to_index[species]] += self.get[element][species]

        return out


    @property
    def to_pandas(self):
        return DataFrame(self.to_array, index=self.list, columns=self.species_to_index.keys())


    @property
    def is_balanced(self): 
        # A perfectly balanced equation will have all sums equal to 0. 
        # Therefore, np.sum(..., axis=1) == 0, and (1 + sum) == 1. 
        # Using np.isclose is much safer for floating point math than exact == or bool()

        return bool( np.allclose(np.sum(self.to_array, axis=1), 0.0, atol=1e-6) ) 
    
    
    @property
    def rk(self):
        return np.linalg.matrix_rank(self.to_array)
    

    def __repr__(self):
        return DataFrame.__repr__(self.to_pandas)