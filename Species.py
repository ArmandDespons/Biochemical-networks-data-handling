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


    @classmethod
    def _resolve_formula(cls, formula: Optional[str], name: Optional[str]) -> str:
        """Determine the canonical DB formula key for a (formula, name) pair.

        This is the single source of truth for identity: __new__ uses it to key
        the registry and __init__ uses it to populate attributes, so a species
        looked up by formula or by name always maps to the same instance.
        """

        if formula is None and name is None:
            raise ValueError("Cannot initialize a species without a formula or name.")

        if formula is not None:
            entry = cls._DATA.get(formula.lower())
            # A formula can itself be an alias pointing at the canonical formula.
            return entry if isinstance(entry, str) else formula

        entry = cls._DATA.get(name.lower())
        if isinstance(entry, str):
            return entry

        raise ValueError(f"Species named {name} not found in the database")


    def __new__(cls, formula: Optional[str] = None, name: Optional[str] = None):

        key = cls._resolve_formula(formula, name)

        if key not in cls._REGISTRY:
            instance = super().__new__(cls)
            instance._initialized = False
            cls._REGISTRY[key] = instance

        return cls._REGISTRY[key]


    def __init__(self, formula: Optional[str] = None, name: Optional[str] = None):

        if self._initialized:
            return

        self.formula = self._resolve_formula(formula, name)
        entry = self._DATA.get(self.formula)
        self._override_molar_mass = None
        self._override_phase = None
        self._override_organic = None
        self.alias = name if name else None

        if isinstance(entry, dict):
            self._in_db = True
            self.db_entry = entry

            self.name = self.db_entry.get("Name")
            self.stoichiometry = dict(self.db_entry.get("Stoichiometry") or {})
            self.charge = self.db_entry.get("Charge")
            
            self._organic = self.db_entry.get("Organic")
            self._phase = self.db_entry.get("Phase")
            self._aliases = self.db_entry.get("Aliases", [])

        else:
            self._in_db = False
            self.db_entry = None
            self.stoichiometry = {}
            self.charge = 0  # overwritten by _infer_from_formula if the formula has a charge suffix
            self._organic = None
            self._phase = None
            self._aliases = []

            self._infer_from_formula()

            self.name = name

        self._initialized = True



    @classmethod
    def from_dict(cls, dic: dict):
        formula = ''
        name = None
        for element, stoich in dic.items(): 

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
    def from_list(cls, lst_species: list[str]):
        """Instantiate species from a non-homogeneous list mixing chemical formulas and/or chemical names."""

        if not isinstance(lst_species, list):
            raise TypeError("'lst_species' argument should be an instance of list[str] containing chemical formulas and/or names in the database.")

        species_list = []
        for item in lst_species:
            if not isinstance(item, str):
                raise TypeError(f"Each item of 'lst_species' must be a string (formula or name), got {type(item).__name__}.")

            # A name key maps to a formula string in the DB; a formula key maps to its info dict.
            # DB name keys are stored lower-case, so match case-insensitively before falling back to a formula lookup.
            species_list.append(cls(name=item) if isinstance(cls._DATA.get(item.lower()), str) else cls(item))

        return species_list
            


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


    def __getattr__(self, attr):
        # Guard against recursion: bail out fast for internal/dunder attributes,
        # and before `stoichiometry` exists yet (e.g. mid-__init__ or on
        # bootstrap-only lookups), since `_made_of_atoms` needs it to be set.
        if attr.startswith('_') or 'stoichiometry' not in self.__dict__:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{attr}'")

        if attr in Atom._REGISTRY and self._made_of_atoms:
            return self.stoichiometry.get(Atom.get(attr).formula, 0)

        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{attr}'")


    # ---- Degree of reductance (gamma) ----

    def gamma(self, per_carbon=False): 
        if not self._made_of_atoms: 
            raise TypeError("Cannot compute the degree of reductance of species that are not made of known Atoms.")

        DoR = -self.charge

        for el, s in self.stoich.items():
            DoR += Atom.get(el).v*s

        if per_carbon:
            DoR /= self.stoich.get('C', 1)

        return DoR


    # ---- Molar Mass (M) ----

    @property
    def molar_mass(self) -> float:
        """Returns the overridden mass if present, otherwise computes it."""
        if self._override_molar_mass is not None:
            return self._override_molar_mass
        
        if self._made_of_atoms:
            return sum(Atom.get(s).M * c for s, c in self.stoichiometry.items())
        

    @molar_mass.setter
    def molar_mass(self, value: float):
        """Set a custom molar mass with context-specific behavior."""

        if self._made_of_atoms:
            print(f"Warnings: Species '{self.formula}' is composed of known atoms. Overriding computed molar mass with {value}.")

        self._override_molar_mass = value


    # ---- Phase ----

    @property
    def phase(self):
        """Returns the overridden phase if present, otherwise the DB-defined one."""
        if self._override_phase is not None:
            return self._override_phase

        return self._phase


    @phase.setter
    def phase(self, value):
        """Set a custom phase with context-specific behavior."""

        if self._phase is not None:
            print(f"Warnings: Species '{self.formula}' already has a defined phase ('{self._phase}'). Overriding with {value}.")

        self._override_phase = value


    # ---- Organic ----

    @property
    def organic(self):
        """Returns the overridden organic flag if present, otherwise the DB-defined one."""
        if self._override_organic is not None:
            return self._override_organic

        return self._organic


    @organic.setter
    def organic(self, value):
        """Set a custom organic flag with context-specific behavior."""

        if self._organic is not None:
            print(f"Warnings: Species '{self.formula}' already has a defined organic flag ('{self._organic}'). Overriding with {value}.")

        self._override_organic = value

    @property
    def M(self) -> float:
        return self.molar_mass
    
    @M.setter
    def M(self, value):
        self.molar_mass = value


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
        elif isinstance(activity, dict): pass
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
            raise NotImplementedError(f'Species {self.formula} not found in the database. Not implemented yet')


    def mu(self, activity: Optional[float] = None, T: Optional[float] = 298.15, pH: Optional[float] = 7, method: Optional[str] = 'eQ pH=0'):

        return self.chemical_potential(activity, T, pH, method)
    
    def __repr__(self):
        return f"Species(Name={self.name}, Formula={self.formula})"
    