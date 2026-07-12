from dataclasses import dataclass
from typing import Dict, Optional

@dataclass(frozen=True)
class Atom:
    formula: str
    name: str
    molar_mass: float
    atomic_number: int
    electronegativity: float
    valence_electrons: int
    reference_valence_electrons: int

    _REGISTRY: Dict[str, 'Atom'] = None


    def __post_init__(self):
        """Runs automatically after A = Atom(...) to register the new atom."""

        # Ensure the standard table is loaded first
        if not getattr(Atom, '_initializing', False) and Atom._REGISTRY is None:
            Atom._initialize_REGISTRY()
        
        # Add this instance to the global registry
        Atom._REGISTRY[self.formula] = self

        # Also register by lowercase name for flexible lookup
        Atom._REGISTRY[self.name.lower()] = self


    @classmethod
    def _initialize_REGISTRY(cls):
        cls._initializing = True
        cls._REGISTRY = {}
        
        # Format: (Formula, Name, Mass, Atomic number, Electronegativity, Valence e-, Reference valence e-)
        atoms_data = [
            ("H",  "Hydrogen",   1.008,  1,  2.20, 1, 0),
            ("C",  "Carbon",    12.011,  6,  2.55, 4, 0),
            ("O",  "Oxygen",    15.999,  8,  3.44, 6, 8),
            ("N",  "Nitrogen",  14.007,  7,  3.04, 5, 8),
            ("P",  "Phosphorus",30.974, 15,  2.19, 5, 0),
            ("S",  "Sulfur",    32.06,  16,  2.58, 6, 0),
            ("Mg", "Magnesium", 24.305, 12,  1.31, 2, 0),
            ("Na", "Sodium",    22.990, 11,  0.93, 1, 0),
            ("Li", "Lithium",    6.94,   3,  0.98, 1, 0),
            ("Cl", "Chlorine",  35.45,  17,  3.16, 7, 8),
            ("Mn", "Manganese", 54.938, 25,  1.55, 7, 0),
            ("Ca", "Calcium",   40.078, 20,  1.00, 2, 0),
            ("K",  "Potassium", 39.098, 19,  0.82, 1, 0),
            ("Fe", "Iron",      55.845, 26,  1.83, 8, 5),
        ]

        cls._REGISTRY.update({
            symbol: cls(symbol, name, mass, z, en, valence, ref_valence) 
            for symbol, name, mass, z, en, valence, ref_valence in atoms_data
        })

        # Add lowercase names to registry for flexible lookup
        cls._REGISTRY.update({
            name.lower(): cls._REGISTRY[symbol] 
            for symbol, name, mass, z, en, valence, ref_valence in atoms_data
        })
        
        cls._initializing = False


    @classmethod
    def get(cls, key: str|list[str]) -> 'Atom':
        if cls._REGISTRY is None:
            cls._initialize_REGISTRY()
        
        if isinstance(key, str):
            result = cls._REGISTRY.get(key) or cls._REGISTRY.get(key.lower())
        
            if not result:
                raise ValueError(f"Atom '{key}' not found. Check spelling or add to registry.")
            
        elif isinstance(key, list) and all(isinstance(item, str) for item in key):

            result = tuple([cls.get(el) for el in key])

        else: 
            
            raise ValueError('Argument must be a string or list of string')

        return result
    


    @property
    def v(self) -> float:
        """Available valence electron"""

        return self.valence_electrons - self.reference_valence_electrons

    # --- Shorthand Aliases ---
    @property
    def M(self) -> float:
        """Molar mass."""
        return self.molar_mass

    @property
    def Z(self) -> int:
        """Atomic number."""
        return self.atomic_number

    @property
    def chi(self) -> Optional[float]:
        """Electronegativity (Pauling scale)."""
        return self.electronegativity
    

    def __repr__(self):
        return (f"{self.name} ({self.formula}): Z={self.atomic_number}, "
                f"Mass={self.M}, Valence={self.valence_electrons}e-, "
                f"Electronegativity={self.chi}"
                )