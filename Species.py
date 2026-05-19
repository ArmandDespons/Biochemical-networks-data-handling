import re
from numpy import log, isnan
import warnings
from typing import Dict, Optional
from ast import literal_eval
from .CHEMICAL_DB import CHEMICAL_DB  
from .Atom import Atom

class Species:
    _DATA = CHEMICAL_DB
    _REGISTRY: Dict[str, 'Species'] = {}

    Atom._initialize_REGISTRY()


    def __new__(cls, formula: str, name: Optional[str] = None, **kwargs):

        if formula not in cls._REGISTRY:
            cls._REGISTRY[formula] = super().__new__(cls)

        return cls._REGISTRY[formula]


    def __init__(self, formula: str, name: Optional[str] = None, **kwargs):

        if getattr(self, '_initialized', False):
            # Update attributes if already initialized, excluding properties
            for key, value in kwargs.items():
                if not isinstance(getattr(type(self), key, None), property):
                    setattr(self, key, value)
            return

        self._initialized = True
        self.formula = formula
        
        db_entry = self._DATA.get(formula, {})
        if db_entry != {}:
            self.db_entry = db_entry
            self._in_db = True
        else: 
            self._in_db = False
            self.db_entry = 'Not found in the database'

        self.name = name or db_entry.get('Name', formula)
        
        self.stoichiometry: Dict[str, float] = {}
        self.charge: int = 0
        self._infer_from_formula()
        

        
        # Internal storage for overridden mass
        self._override_molar_mass = kwargs.get('molar_mass', None)

        combined_data = {**db_entry, **kwargs}
        for key, value in combined_data.items():
            safe_key = key.lower().split(' (')[0].replace(" ", "_")

            # Avoid overwriting properties or element-specific counts
            if safe_key in ['molar_mass', 'm', 'stoich', 'g', 'h', 'c', 'o', 'h', 'n', 's']:
                continue
            setattr(self, safe_key, value)


    @classmethod
    def from_dict(cls, dic: dict):
        formula = ''
        name = None
        for element, stoich in zip(dic.keys(), dic.values()): 

            if element.lower() != 'charge' and element.lower() != 'name':
                formula += element + str(stoich)
            
            elif element.lower() == 'charge': 
                if stoich<0:
                    formula += str(stoich)
                else: 
                    formula += '+' + str(stoich)

            elif element.lower() == 'name': 
                name = stoich

        return cls(formula, name=name)


    @classmethod
    def from_name(cls, name: str): 
        _name = name.lower()

        if not _name in cls._DATA.keys():
            raise ValueError(f'{_name} not found in the dataset.')

        else:
            if isinstance(cls._DATA[_name], str):
        
                return cls(cls._DATA[_name], name=_name)


    def _infer_from_formula(self):
        """Parses formula and charge from the string."""

        working_formula = self.formula
        charge_match = re.search(r'([+-]\d*)$', working_formula)
        if charge_match:
            raw_charge = charge_match.group(1)
            self.charge = 1 if raw_charge == '+' else -1 if raw_charge == '-' else int(raw_charge)
            working_formula = working_formula[:charge_match.start()]

        elements = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', working_formula)
        if not elements and working_formula:
            self.stoichiometry[working_formula] = 1.0
        else:
            for symbol, count in elements:
                val = float(count) if (count and count.strip()) else 1.0
                self.stoichiometry[symbol] = self.stoichiometry.get(symbol, 0.0) + val


    @property
    def _made_of_atoms(self):
        return all(symbol in Atom._REGISTRY for symbol in self.stoichiometry.keys())
    

    def _update_name(self, value): 
        self.name = value

    # ---- Molar Mass (M) ----

    @property
    def molar_mass(self) -> float:
        """Returns the overridden mass if present, otherwise computes it."""
        if self._override_molar_mass is not None:
            return self._override_molar_mass
        
        # Standard dynamic calculation
        return sum(Atom.get(s).M * c for s, c in self.stoichiometry.items())

    @molar_mass.setter
    def molar_mass(self, value: float):
        """Set a custom molar mass with context-specific behavior."""

        if self._made_of_atoms:
            print(f"Warnings: Species '{self.formula}' is composed of known atoms. Overriding computed molar mass with {value}.")

        self._override_molar_mass = value

    @property
    def M(self) -> float:
        return self.molar_mass
    
    @M.setter
    def M(self, value):
        return self.molar_mass(value)


    # ---- Shortland to obtain the stoichiometry dictionary ----

    @property
    def stoich(self) -> Dict[str, float]:
        return self.stoichiometry


    # ---- THERMODYNAMIC PROPERTIES ----

    # ---- Standard enthalpy of formation (Hf0) ----

    @property
    def standard_enthaply(self): 
        if self._in_db: 
            return self.db_entry.get('Enthalpy of formation (kJ/mol)')
        
        elif self.name=='biomass': 

            A = 115 # kJ/Cmol 
            q_c = A*(4*self.stoich.get('C', 0) + self.stoich.get('H', 0) - 2*self.stoich.get('O', 0) )

            # a_O2 = 1 + biomass_stoich['H']/4 - biomass_stoich['O']/2
            # a_N2 = biomass_stoich['N'] 

            h_CO2 = -393.5
            h_H2O = -285.8

            return self.stoich.get('C', 0)*h_CO2 + self.stoich.get('H', 0)/2*h_H2O + q_c 
        
        else: 
            raise NotImplementedError('Not implemented yet')
        

    @property
    def H0f(self): 
        return self.standard_enthaply


    # ---- Standard Gibbs free-energy of formation (Gf0) ----

    def standard_free_energy(self, method: Optional[str] = 'eQ pH=0', T: Optional[float] = 298.15):
        
        if isnan(T): T = 298.15

        if self._in_db: 
            availaible_methods = self.db_entry['Free-energy of formation (kJ/mol)'].keys()
            
            if method in availaible_methods:
                return self.db_entry['Free-energy of formation (kJ/mol)'][method]
            
            else: 
                raise ValueError(f"Unknwon '{method}' method. Availaible methods for the given species: {list(availaible_methods)}")


        elif self.name=='biomass':
            S_C, S_H, S_O, S_N = 5.74e-3, 65.34e-3, 102.58e-3, 91.81e-3

            entropy_of_formation  = -0.813*(self.stoich.get('C', 0)*S_C + self.stoich.get('H', 0)*S_H + self.stoich.get('O', 0)*S_O + self.stoich.get('N', 0)*S_N)

            return self.standard_enthaply - T*entropy_of_formation
        
        else: 
            raise NotImplementedError('Not implemented yet')
        

    def G0f(self, method: Optional[str] = 'eQ pH=0', T: Optional[float] = 298.15): 
        
        return self.standard_free_energy(method=method, T=T)
    

    # ---- Chemical potential (mu) ----

    def chemical_potential(self, 
                           activity: Optional[float] = None,
                           T: Optional[float] = 298.15, 
                           pH: Optional[float] = 7,
                           method: Optional[str] = 'eQ pH=0'): 
        
        if isnan(T): T = 298.15
        if isnan(pH): pH = 7

        if isinstance(activity, str): activity = literal_eval(activity)
        elif activity is None: pass
        elif isnan(activity): activity = None

        if self._in_db:

            _activity = None

            if isinstance(activity, dict):
                if self.formula in activity.keys(): 
                    _activity = activity.get(self.formula)

            if isinstance(activity, float) or isinstance(activity, int): 
                _activity = activity

            if _activity is None: 

                if self.db_entry['Phase']=='aq' and self.formula!='H+1': 
                    _activity = 1e-3

                elif self.formula=='H+1': 
                    _activity = 10**(-pH)

                else: 
                    _activity = 1

            R = 8.31446261815324e-3
            mu0 = self.standard_free_energy(method=method)*(T/298.15) + self.standard_enthaply*(1 - T/298.15)

            return mu0 + R*T*log(_activity)
        

        elif self.name=='biomass':

            return self.standard_free_energy(method=method, T=T)*(T/298.15) + self.standard_enthaply*(1 - T/298.15)
        

        else: 
            raise NotImplementedError('Not implemented yet')


    def mu(self, activity: Optional[float] = None, T: Optional[float] = 298.15, pH: Optional[float] = 7, method: Optional[str] = 'eQ pH=0'):

        return self.chemical_potential(activity, T, pH, method)
    
    def __repr__(self):
        return f"Species(Name={self.name}, Formula={self.formula})"