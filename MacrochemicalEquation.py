import re

import numpy as np 

from pandas import DataFrame
from ast import literal_eval
from collections.abc import Iterable

from .misc import _check_coeffs

from typing import List, Dict, Union, Optional

from .Species import *
from .Equation import *

__all__ = ["MacrochemEquation"]




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
            raise TypeError("'Index' should be either an iterable object or an int. ")


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

    @Ea.setter
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

    @Cs.setter
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
    