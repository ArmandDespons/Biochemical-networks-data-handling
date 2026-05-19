import re

import numpy as np 

from pandas import DataFrame
from ast import literal_eval
from collections.abc import Iterable

from .misc import _matrix_parser, _check_coeffs

from typing import List, Dict, Union, Optional
from .Species import Species


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

            if isinstance(coeff, float) or isinstance(coeff, int):

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


    def infer_coeffs(self, known_coeff: Optional[Dict[str, float]] = None, verbose=True):
        """
        Solves the linear system A * x = b for the unknown coefficients.
        Accepts an optional dictionary of additional known fluxes/coefficients.
        """

        result_inference = {}

        # Erase the previous results
        if self.balanced_coefficients: 

            self.default_coefficients = {s: 1.0 for s in self.balanced_coefficients.keys()}
            self.balanced_coefficients = None


        # Update knowns if the user provided additional fluxes
        if known_coeff:

            _check_coeffs(known_coeff, self.reactants, self.products)

             # Remove these newly defined species from the defaults
            for species in known_coeff.keys():
                self.default_coefficients.pop(species, None)

            _to_delete = []
            for s in self.known_coefficients.keys():
                
                if not s in known_coeff.keys():
                    _to_delete.append(s)
                    self.default_coefficients[s] = 1.0
                
            for s in _to_delete:
                self.known_coefficients.pop(s)

                if verbose: print(f'Warning: Overriding previously kwown coefficient for species {s}, now set to unknown')

            del _to_delete 

            self.known_coefficients.update(known_coeff)
            
            _check_coeffs(self.default_coefficients, self.reactants, self.products)
            

        # Refresh the conservation laws matrix to reflect the new knowns
        self.conservation_laws = ConservationLaws(self)
            
        
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

        result_inference['RHS'] = RHS

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
            # result_inference['Message'] = None if all(sol >= 0) else f'Warning: linear problem yield negative coefficients'
            result_inference['Inferred coefficients'] = self.balanced_coefficients

        
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



    def __repr__(self):
        def fmt(s_list):
            return " + ".join([f"{abs(self.coefficients[s.formula])} {s.formula}" for s in s_list])
        

        return f"{fmt(self.reactants)} = {fmt(self.products)}"




class ConservationLaws:

    def __init__(self, eq: Equation): 
        """
        Construct a nested dictionnary containing the conservation laws of all the element parsed from
        an Equation instance
        """

        # 1. Identify all unique species and elements
        all_species = eq.reactants + eq.products
        all_elements = set()
        for s in all_species:
            all_elements.update(s.stoich.keys())

        # 2. Initialize the nested structure
        self.get = {el: {} for el in all_elements}
        self.get['Charge'] = {}

        coeffs = eq.coefficients

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





class MacrochemEquation(Equation):
    
    def __init__(self, 
        substrates: List[Union[str, Species]], 
        products: List[Union[str, Species]], 

        biomass_stoichiometry: Optional[Dict[str, float]] = None, 

        electron_donor: Optional[str] = None,
        electron_acceptor: Optional[str] = None,
        carbon_source: Optional[str] = None,

        yields: Optional[Dict[str, float]] = None,
        biomass_alias: Optional[str] = 'CH1.8O0.5N0.2'
    ):
        
        if yields is None:
            yields = {}
        
        
        if biomass_stoichiometry: 
            _biomass_stoichiometry = {}
            _biomass_stoichiometry.update({'name': 'biomass'})

            _biomass_elements = ['C', 'H', 'O', 'N', 'P', 'S']

            for element in _biomass_elements: 
                for key in biomass_stoichiometry: 

                    if element in key: 
                        _biomass_stoichiometry.update({element: biomass_stoichiometry.get(key)})

            if not 'C' in _biomass_stoichiometry.keys():
                _biomass_stoichiometry = {'C': 1, **_biomass_stoichiometry}


            self.biomass = Species.from_dict(_biomass_stoichiometry)

            products[0] = self.biomass

        else: 
            if isinstance(products[0], Species): 
                products[0]._update_name('biomass')
            else: 
                products[0] = Species(products[0], name='biomass')

            self.biomass = products[0]

        self.biomass_stoich = self.biomass.stoich


        if isinstance(biomass_alias, str): 
            self._biomass_alias = biomass_alias

        elif biomass_alias is None:
            self._biomass_alias = self.biomass.formula


        super().__init__(reactants = substrates, 
                         products = products, 
                         known_coefficients = self._parse_biomass_alias(yields, self._biomass_alias, self.biomass.formula)
                         )


        self.substrates = self.reactants

        self._yields = self.known_coefficients


        self._electron_donor = self._validate_species(electron_donor) if electron_donor else None
        self._electron_acceptor = self._validate_species(electron_acceptor) if electron_acceptor else None
        self._carbon_source = self._validate_species(carbon_source) if carbon_source else None


    @classmethod
    def from_string(cls, eqn_str: str, **kwargs):

        eqn = Equation.from_string(eqn_str)
     
        return cls(substrates=eqn.reactants, products=eqn.products, yields=eqn.known_coefficients, **kwargs)
    

    @classmethod
    def from_dict(cls, dict: Dict[str, float], biomass_alias='CH1.8O0.5N0.2', biomass_key=None, **kwargs):
        
        if biomass_key is None: biomass_key = biomass_alias

        if not biomass_key in dict.keys():
            raise ValueError(f"Biomass key '{biomass_key}' not found among the dictionary keys: {list(dict.keys())}")

        _element_biomass = ['C biomass', 'H biomass', 'O biomass', 'N biomass', 'P biomass', 'S biomass']
        _dict = dict.copy()
        biomass_stoich = {k: _dict.pop(k) for k in _element_biomass if k in dict.keys()}

        if biomass_stoich=={}:
            biomass_stoich = None

        eqn = Equation.from_dict({biomass_alias: _dict.pop(biomass_key), **_dict})

        return cls(substrates=eqn.reactants, 
                   products=eqn.products, 
                   biomass_stoichiometry=biomass_stoich, 
                   yields=eqn.known_coefficients, 
                   biomass_alias=biomass_alias,
                   **kwargs
                   )
    
    @classmethod
    def from_parsed_df(cls, df: DataFrame, index: Optional[int] = None, **kwargs):

        if index is None: 
            return {index: cls.from_parsed_df(df, index) for index in df.index}
            
        elif isinstance(index, Iterable):
            return {i: cls.from_parsed_df(df, i) for i in index}
        
        elif isinstance(index, int):

            _necessary_cols = ['Substrates', 'Products', 'Electron donor formula', 'Electron acceptor formula', 'Carbon source formula', 'H biomass', 'O biomass', 'N biomass']
            if not all([i in df.columns for i in _necessary_cols]):
                raise ValueError(f"The parsed DataFrame must contain columns: {_necessary_cols}.")


            substrates = df.loc[index, 'Substrates'] if isinstance(df.loc[index, 'Substrates'], list) else literal_eval(df.loc[index, 'Substrates']) 
            products = df.loc[index, 'Products'] if isinstance(df.loc[index, 'Products'], list) else literal_eval(df.loc[index, 'Products']) 

            Ed, Ea, Cs = df.loc[index, 'Electron donor formula'], df.loc[index, 'Electron acceptor formula'], df.loc[index, 'Carbon source formula']

            biomass_stoichiometry = df.loc[index, ['H biomass', 'O biomass', 'N biomass']].rename(index={'H biomass': 'H', 'O biomass': 'O', 'N biomass': 'N'}).to_dict()
                        
            return  cls(substrates=substrates, 
                        products=products, 
                        biomass_stoichiometry=biomass_stoichiometry,
                        electron_donor=Ed,
                        electron_acceptor=Ea,
                        carbon_source=Cs,
                        **kwargs
                        )
        
        else: 
            raise ValueError("'Index' should be either an iterable object or an int. ")


    def _validate_species(self, formula: str) -> Species:
        """
            Convert a species formula (str) to a Species instance, provided it is 
            among the substrates/products of the macrochemical equation. 
        """

        if formula is None: 
            return formula

        for s in self.reactants + self.products:
            if s.formula == formula:
                return s
            
        raise ValueError(f"Role Error: {formula} must be one of the reactants: {[r.formula for r in self.reactants]}")


    @staticmethod
    def _use_biomass_alias(dict, biomass_alias, biomass_formula):
        out = dict.copy()
        out.update({biomass_alias: dict.get(biomass_formula)})
        out.pop(biomass_formula, None)
    
        return out
    
    @staticmethod
    def _parse_biomass_alias(dict, biomass_alias, biomass_formula):

        if biomass_alias in dict.keys(): 
            out = dict.copy()
            out.update({biomass_formula: out.pop(biomass_alias)})

            return out
    
        else: 

            return dict


    def infer_yields(self, known_coeff: Optional[Dict[str, float]]=None, verbose=True):

        _known_coeff = None
        if known_coeff:
            _known_coeff = self._parse_biomass_alias(known_coeff, self._biomass_alias, self.biomass.formula)

        self.infer_coeffs(_known_coeff, verbose=verbose)


    
    def list_species(self, which='all', use_biomass_alias=True): 

        if not which in ['all', 'reactants', 'products']:
            raise ValueError(f"'which' takes the following values: 'all', 'reactants', 'products'")
        
        if which=='all':
            if use_biomass_alias:
               return [s.formula if s.formula!=self.biomass.formula else self._biomass_alias for s in self.reactants+self.products]
            
            else: 
                return [s.formula for s in self.reactants+self.products]

        elif which=='products':
            if use_biomass_alias:
               return [s.formula if s.formula!=self.biomass.formula else self._biomass_alias for s in self.products]
            
            else: 
                return [s.formula for s in self.products]
        
        else: 
            return [s.formula for s in self.reactants]

    # --- Electron Donor (Ed) ---

    @property
    def electron_donor(self) -> Optional[Species]:
        return self._electron_donor

    @electron_donor.setter
    def electron_donor(self, value: Union[str, Species]):
        formula = value.formula if isinstance(value, Species) else value
        self._electron_donor = self._validate_species(formula)

    @property
    def Ed(self): return self.electron_donor

    @Ed.setter
    def Ed(self, value): self.electron_donor = value


    # --- Electron Acceptor (Ea) ---

    @property
    def electron_acceptor(self) -> Optional[Species]:
        return self._electron_acceptor

    @electron_acceptor.setter
    def electron_acceptor(self, value: Union[str, Species]):
        formula = value.formula if isinstance(value, Species) else value
        self._electron_acceptor = self._validate_species(formula)

    @property
    def Ea(self): return self.electron_acceptor

    @Ed.setter
    def Ea(self, value): self.electron_acceptor = value


    # --- Carbon source (Cs) ---

    @property
    def carbon_source(self) -> Optional[Species]:
        return self._carbon_source

    @carbon_source.setter
    def carbon_source(self, value: Union[str, Species]):
        formula = value.formula if isinstance(value, Species) else value
        self._carbon_source = self._validate_species(formula)

    @property
    def Cs(self): return self.carbon_source

    @Ed.setter
    def Cs(self, value): self.carbon_source = value


    # --- Yields (y) ---

    @property
    def yields(self) -> Dict[str, np.float64]:
        return self._use_biomass_alias(self.coefficients, self._biomass_alias, self.biomass.formula)
    
    @yields.setter
    def yields(self, value: Dict[str, float]):
        """Updates the known coefficients and refreshes conservation laws."""
       
        _new_yields = self._parse_biomass_alias(value.copy(), self._biomass_alias, self.biomass.formula)
        for key in value.keys(): 
            if not key in self.list_species(): _new_yields.pop(key)

        
        _check_coeffs(_new_yields, self.reactants, self.products)

        self.known_coefficients.update(_new_yields)

        for s in value.keys():
            if s in self.default_coefficients.keys(): self.default_coefficients.pop(s)

        # We re-run the conservation law setup to account for new knowns
        self.conservation_laws = ConservationLaws(self)

    @property
    def y(self): return self.yields

    @y.setter
    def y(self, value): self.yields = value



    def __repr__(self):
        base = super().__repr__()

        Ed = self._electron_donor.formula if self._electron_donor else 'Unknown'
        Ea = self._electron_acceptor.formula if self._electron_acceptor else 'Unknown'
        Cs = self._carbon_source.formula if self._carbon_source else 'Unknown'

        return f"{base}" + "\n \n" + f"Electron donor: {Ed}" + " | " + f"Electron acceptor: {Ea}" + " | " + f"Carbon source: {Cs}"

