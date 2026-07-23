from numbers import Real
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from .Species import Species

__all__ = ["Complex"]

_ELEMENT_ORDER = ["C", "H", "O", "N", "P", "S", "Mg", "Na", "Li", "Cl", "Mn", "Ca", "K", "Fe"]


class Complex:

    def __init__(
            self,
            species: List[Species],
            coefficients: Optional[Dict[Union[str, Species], float]] = None,
            default_coeff = 1
            ):

        if not isinstance(species, list):
            raise TypeError("'species' must be a list of Species instances, chemical formulas, and/or chemical names.")

        # Set once at construction; there is no setter, so it cannot change afterwards.
        # True: coefficients keep their signed (algebraic) value. False: only their absolute value is stored.
        self._default = default_coeff

        # A complex is a *set* of species: resolve mixed items (Species instances, formulas,
        # and/or names), then de-duplicate while preserving order.
        self.species = list(dict.fromkeys(self._coerce_species(item) for item in species))

        if coefficients is None:
            coefficients = {}

        self._known_coefficients = {}
        for key, value in coefficients.items():
            s = self._resolve_species(key)
            self._known_coefficients[s.formula] = value

        self._default_coefficients = {
            s.formula: self._default for s in self.species if s.formula not in self._known_coefficients
        }


    @staticmethod
    def _coerce_species(item: Union[str, Species]) -> Species:
        """Resolve a species-list item (Species instance, chemical formula, or chemical name) to a Species."""

        if isinstance(item, Species):
            return item

        if isinstance(item, str):
            # DB name keys are stored lower-case; match case-insensitively before falling back to a formula lookup.
            if isinstance(Species._DATA.get(item.lower()), str):
                return Species(name=item)

            return Species(item)

        raise TypeError(f"Each species must be a Species instance, chemical formula, or chemical name (str), got {type(item).__name__}.")


    def _resolve_species(self, key: Union[str, Species]) -> Species:
        """Resolve a coefficient dict key to one of this complex's Species instances.

        Accepted key forms: chemical formula (str),
        a Species instance, species name (str), then species alias (str).
        """

        if isinstance(key, Species):
            if key in self.species:
                return key

            raise ValueError(f"{key!r} is not part of this complex.")

        if isinstance(key, str):
            for s in self.species:
                if s.formula == key:
                    return s

                if s.name == key.casefold():
                    return s

                if s.alias == key.lower():
                    return s
                
                if key.lower() in s._aliases:
                    return s

            raise ValueError(f"'{key}' does not match the formula, name, or alias of any species in this complex.")

        raise TypeError(f"Coefficient keys must be a formula/name/alias (str) or a Species instance, got {type(key).__name__}.")


    def pop(self, species: Union[str, Species]) -> Species:
        """Resolve `species` (formula, name, alias, or Species instance) and remove it from this complex.

        Returns the removed Species instance.
        """
        try:
            s = self._resolve_species(species)

        except ValueError:
            return None

        self.species.remove(s)
        self._known_coefficients.pop(s.formula, None)
        self._default_coefficients.pop(s.formula, None)

        return s


    @property
    def coefficients(self) -> Dict[str, float]:
        """Combines default and user-defined coefficients for every species in the complex."""
        return {**self._default_coefficients, **self._known_coefficients}


    def _update(self, new_coeff: Dict[Union[str, Species], float]):
        """Replace `_known_coefficients` with `new_coeff` and recompute `_default_coefficients` to match.
        """

        self._known_coefficients = {
            self._resolve_species(key).formula: value for key, value in new_coeff.items()
        }

        self._default_coefficients = {
            s.formula: self._default for s in self.species if s.formula not in self._known_coefficients
        }


    @property
    def elemental_composition(self) -> Tuple[np.ndarray, Dict[str, int], Dict[str, int]]:
        """Coefficient-weighted elemental composition matrix.

        Returns (matrix, species_to_col, element_to_row): rows are the union of
        elements found across every species' `stoich`, plus a 'Charge' row only
        if at least one species carries a nonzero charge; columns are the
        species (in `self.species` order), each entry weighted by the species'
        effective coefficient.
        """

        coeffs = self.coefficients

        all_elements = set()
        for s in self.species:
            all_elements.update(s.stoich.keys())

        order_index = {el: i for i, el in enumerate(_ELEMENT_ORDER)}
        elements = sorted(all_elements, key=lambda el: (order_index.get(el, len(_ELEMENT_ORDER)), el))

        has_charge = any(s.charge for s in self.species)
        if has_charge:
            elements.append("Charge")

        element_to_row = {el: i for i, el in enumerate(elements)}
        species_to_col = {s.formula: i for i, s in enumerate(self.species)}

        matrix = np.zeros((len(elements), len(self.species)))

        for s in self.species:
            col = species_to_col[s.formula]
            coeff = coeffs[s.formula]

            for el, amount in s.stoich.items():
                matrix[element_to_row[el], col] = coeff * amount

            if has_charge:
                matrix[element_to_row["Charge"], col] = coeff * s.charge

        return matrix, species_to_col, element_to_row


    # ---- Element-wise proxy for Species attributes/properties ----

    def __getattr__(self, attr):
        """Proxy a Species attribute/property element-wise across `self.species`.

        Numeric results are weighted by each species' effective coefficient;
        anything else (str, bool, None, ...) is returned as-is, coefficients
        ignored. This transparently supports Species' dynamically-resolved atom
        attributes too (e.g. `.C`, `.carbon`), since `getattr` triggers
        Species.__getattr__ the same way direct attribute access would.
        """

        if attr.startswith('_'):
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{attr}'")

        coeffs = self.coefficients
        values = [getattr(s, attr) for s in self.species]

        return [
            coeffs[s.formula] * v if isinstance(v, Real) and not isinstance(v, bool) else v
            for s, v in zip(self.species, values)
        ]


    def __add__(self, other):
        """Merge two complexes: union of species, coefficients combined (ties broken by `other`)."""

        if not isinstance(other, Complex):
            return NotImplemented

        species = list(dict.fromkeys(self.species + other.species))
        known = {**self._known_coefficients, **other._known_coefficients}
        default = {**self._default_coefficients, **other._default_coefficients}

        merged = Complex(species)
        merged._known_coefficients = known
        merged._default_coefficients = {s.formula: default[s.formula] for s in species if s.formula not in known}

        return merged


    def __len__(self):
        
        return len(self.species)


    def __repr__(self):
        coeffs = self.coefficients
        return "Complex(" + " + ".join(f"{abs(coeffs[s.formula])} {s.formula}" for s in self.species) + ")"
