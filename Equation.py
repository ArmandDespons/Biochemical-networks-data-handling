import re

import numpy as np

from numbers import Real
from typing import Dict, List, Optional, Union

from .Species import Species
from .Complex import Complex


__all__ = ["Equation"]


class Equation:
    """
    A chemical equation: a `reactants` complex balanced against a `products`
    complex, each species carrying a signed stoichiometric coefficient
    (negative for reactants, positive for products).

    A coefficient is either known (set explicitly — via `known_coefficients`
    at construction, `update_coefficients`, or the `coefficients`/`coeff`
    setter) or default (an unset placeholder of magnitude 1); `coefficients`/
    `coeff` merge the two, preferring `balanced_coefficients` once the
    equation has been solved. Construct via `__init__` (species lists, e.g.
    `Equation(['A', 'B'], ['AB'])`), `from_string` ("2 A + 2 B = 2 AB"), or
    `from_dict` ({species: signed coefficient, or a '+'/'-' sentinel for an
    unknown one}).

    `element_conservation` (an `ElementalConservation`) is rebuilt from the
    current coefficients on every change; `is_balanced` and `element_recovery`
    read off of it. `infer_coeff` solves for an independent set of species'
    coefficients from that conservation matrix — correcting any reactant/
    product role mismatch it finds along the way — and reports what it did in
    `result_inference`; `clear` undoes it.

    Supports `*` (either side) to scale by a nonzero scalar (negative flips
    reactants/products) and `+`/`-` to combine two equations; see each
    operator's docstring for the exact rules.
    """

    def __init__(
        self,
        reactants: List[Union[str, Species]],
        products: List[Union[str, Species]],
        known_coefficients: Optional[Dict[Union[str, Species], float]] = None
    ):

        if known_coefficients is None:
            known_coefficients = {}

        # Bare complexes first (every species at its default coefficient), so that
        # `_resolve_coeff` can match each known-coefficient key against the actual
        # resolved species of each side.
        self.reactants = Complex(reactants, default_coeff=-1)
        self.products = Complex(products, default_coeff=1)

        overlap = set(self.reactants.formula) & set(self.products.formula)

        if overlap:
            raise NotImplementedError(
                f"Species present in both reactants and products are not yet supported: {sorted(overlap)}"
            )

        reactant_known, product_known = self._resolve_coeff(known_coefficients)

        self.reactants = Complex(self.reactants.species, coefficients=reactant_known, default_coeff=-1)
        self.products = Complex(self.products.species, coefficients=product_known, default_coeff=1)


        from .ElementalConservation import ElementalConservation
        self.element_conservation = ElementalConservation(self)

        self.balanced_coefficients = self._default_coefficients.copy() if self.is_balanced else None

        self._overriden_complexes = False
        self._last_independent_set = set()

        self.pH = 7
        self.T = 298.15
        self.activity = None


    def _resolve_coeff(self, known_coefficients: Dict[Union[str, Species], float]):
        """Split a combined {species: coefficient} dict between reactants and products.

        Each key (formula, Species instance, name, or alias) is matched against
        whichever side it belongs to; the resolved coefficient's sign is then
        forced to the conventional one for that side (negative for reactants,
        positive for products) regardless of the sign the caller used.
        """

        reactant_known, product_known = {}, {}

        for key, value in known_coefficients.items():

            try:
                s = self.reactants._resolve_species(key)
                reactant_known[s.formula] = -abs(value)
                continue
            except ValueError:
                pass

            try:
                s = self.products._resolve_species(key)
                product_known[s.formula] = abs(value)
                continue
            except ValueError:
                pass

            raise ValueError(f"'{key}' does not match the formula, name, or alias of any species in the reactants or products.")

        return reactant_known, product_known


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
                    # Fallback for single ions like "H+" if the regex is too strictDict[Union[str, Species]
                    obj = Species(item)
                    species_list.append(obj)

            return species_list, known_coeffs

        reac_objs, reac_known = parse_side(left_side)
        prod_objs, prod_known = parse_side(right_side)

        return cls(reac_objs, prod_objs, {**reac_known, **prod_known})


    @classmethod
    def from_dict(cls, coefficients: Dict[Union[str, Species], float]):
        """
        Initialize an Equation instance from a dict mapping species (formula, Species
        instance, name, or alias) to either a known coefficient (a signed number,
        the sign is only used to assign the side, not its final value) or a bare
        '+'/'-' sentinel for a product/reactant whose coefficient is left unknown.
        """

        coefficients = dict(coefficients)

        reactants, products, to_delete = [], [], []

        for species, coeff in coefficients.items():

            if isinstance(coeff, Real) and not isinstance(coeff, bool):

                if np.isnan(coeff):
                    to_delete.append(species)

                elif coeff > 0:
                    products.append(species)

                elif coeff < 0:
                    reactants.append(species)

            elif coeff == '-':
                reactants.append(species)
                to_delete.append(species)

            elif coeff == '+':
                products.append(species)
                to_delete.append(species)

            else:
                raise ValueError(
                    f"Unauthorized value {coeff}. Only numbers are allowed for known coefficients, "
                    "or '+'/'-' for products/reactants with an unknown coefficient."
                )

        for species in to_delete:
            coefficients.pop(species)

        return cls(reactants, products, known_coefficients=coefficients)
    

    @property
    def all_species(self): 

        return self.reactants + self.products

    @property
    def is_balanced(self):

        return bool( np.allclose(np.sum(self.element_conservation.matrix, axis=1), 0.0) )
    

    def element_recovery(self, element: str, in_percent: bool = True):

        if not element in self.element_conservation.elements:
            raise ValueError(f"Element {element} not a valid element {self.element_conservation.elements}.")
        
        base = 100 if in_percent else 1
        el = self.element_conservation.matrix[self.element_conservation.rows[element], :]

        in_reactant = - el[el < 0].sum()
        in_product = el[el > 0].sum()
        
        return float( in_product / in_reactant * base )


    # ----- Thermodynamics (pH, T, activity) -----

    @property
    def pH(self):
        return self._pH

    @pH.setter
    def pH(self, value: float):
        self._pH = float(value)

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, value: float):
        self._T = float(value)

    @property
    def activity(self):
        return self._activity

    @activity.setter
    def activity(self, value: Optional[Union[float, Dict[Union[str, Species], float]]]):
        """
        None: `chemical_potential`'s own phase-based default for every species.
        A single float: applied to every species. A {species: activity} dict
        (keys resolved to formulas): applied per species, mixing with the
        phase-based default for any species left unmentioned.
        """

        if value is None:
            self._activity = None

        elif isinstance(value, dict):
            self._activity = {self.all_species._resolve_species(k).formula: v for k, v in value.items()}

        elif isinstance(value, Real) and not isinstance(value, bool):
            self._activity = float(value)

        else:
            raise TypeError(
                f"'activity' must be a float, a dict of {{species: activity}}, or None, got {type(value).__name__}."
            )


    def DrH(self, _warn: bool = True):
        """Standard reaction enthalpy: sum of each species' signed coefficient times its `standard_enthaply`."""

        if _warn and not self.is_balanced:
            print(f"Warnings: equation '{self}' is not balanced; DrH may not be physically meaningful.")

        coeffs = self.coefficients

        return sum(coeffs[s.formula] * s.standard_enthaply for s in self.all_species.species)


    def DrG(self, method: str = 'eQ pH=0', _warn: bool = True):
        """
        Reaction Gibbs free energy: sum of each species' signed coefficient
        times its `chemical_potential`, evaluated at `self.T`, `self.pH`, and
        `self.activity`.
        """

        if _warn and not self.is_balanced:
            print(f"Warnings: equation '{self}' is not balanced; DrG may not be physically meaningful.")

        coeffs = self.coefficients

        return sum(
            coeffs[s.formula] * s.chemical_potential(activity=self.activity, T=self.T, pH=self.pH, method=method)
            for s in self.all_species.species
        )


    def DrS(self, method: str = 'eQ pH=0'):

        if not self.is_balanced:
            print(f"Warnings: equation '{self}' is not balanced; DrS may not be physically meaningful.")

        return ( self.DrH(_warn=False) - self.DrG(method=method, _warn=False) ) / self.T


    # ----- Coefficients properties/setter -----

    @property
    def _known_coefficients(self): 
        
        return self.all_species._known_coefficients
    
    @property
    def _default_coefficients(self):

        return self.all_species._default_coefficients


    def update_coefficients(self, new_coeff: Dict[Union[str, Species], float], verbose: bool = True, _clear: bool = True):
        """
        Replace the known coefficients with `new_coeff` (see `_resolve_coeff`: the
        caller's sign is ignored and forced to each species' conventional role —
        negative for reactants, positive for products).

        Species previously in `_known_coefficients` but absent from `new_coeff`
        fall back to their default coefficient; if `verbose`, report which ones.

        Resets `balanced_coefficients` to None and rebuilds `element_conservation`.
        """

        if _clear: self.clear()

        reactant_known, product_known = self._resolve_coeff(new_coeff)

        reverted = set(self._known_coefficients) - set(reactant_known) - set(product_known)

        self.reactants._update(reactant_known)
        self.products._update(product_known)

        if verbose and reverted:
            print(f"Coefficients of species {sorted(reverted)} erased and set with default coefficients +/- 1")

        self.balanced_coefficients = None

        from .ElementalConservation import ElementalConservation
        self.element_conservation = ElementalConservation(self)


    @property
    def coefficients(self) -> Dict[str, float]:
        """
        All species coefficients, keyed by formula.

        Once the equation has been balanced, `balanced_coefficients` is merged
        with `known_coefficients` (the latter taking precedence); otherwise
        falls back to `default_coefficients`.
        """

        if self.balanced_coefficients:
            return {**self.balanced_coefficients, **self._known_coefficients}

        return {**self._default_coefficients, **self._known_coefficients}

    @coefficients.setter
    def coefficients(self, new_coeff: Dict[Union[str, Species], float]):
        self.update_coefficients(new_coeff)

    @property
    def coeff(self):
        return self.coefficients

    @coeff.setter
    def coeff(self, new_coeff: Dict[Union[str, Species], float]):
        self.update_coefficients(new_coeff)

    # ----- Infer missing coefficents ----

    def infer_coeff(
            self,
            known_coeff: Optional[Dict[Union[Species, str], float]] = None,
            independent_set: List[Union[str, Species]] = None,
            include: Optional[List[Union[str, Species]]] = None,
            exclude: Optional[List[Union[str, Species]]] = None,
            strict_include: bool = False,
            strict_exclude: bool = False,
            verbose: Optional[bool] = True
            ):
        """
        Solve for the coefficients of the *independent set* — the species whose
        coefficients get computed — via `LHS @ x = RHS` (`ElementalConservation
        ._linear_sys`, `np.linalg.solve`), then store the full result in
        `self.balanced_coefficients`.

        `independent_set` is used as-is if given (`include`/`exclude`/`strict_*`
        are then ignored). Otherwise it's the first result of `get_independent_sets`
        (`include`, `exclude` plus the species named in `known_coeff` — or, if
        that isn't given, in `_known_coefficients` — `strict_include`, `strict_exclude`).

        Every other species becomes known: `known_coeff`'s value where given,
        else its current effective coefficient. A species whose inferred sign
        contradicts its assumed reactant/product role is moved to the correct
        complex instead (`_resolve_complex`).

        Side effects: `self._previous_coefficients` is always set; on a role
        swap, `self._overriden_complexes` and `self._previous_reactants`/
        `_previous_products` are set. `self.result_inference` (success/failure
        details) is always set, and errors are re-raised after being recorded there.

        Raises if there's no spare species beyond the elemental rank to anchor
        the solve (every species would be needed just to reach full rank).
        """

        result = {
            'success': False,
            'independent_set': None,
            'known_set': None,
            'known_coeffs': None,
            'inferred_coeffs': None,
            'sign_mismatch': None,
            'message': None,
        }

        try:
            
            self.update_coefficients(self._known_coefficients, _clear=True) 

            n_species = len(self.element_conservation.cols)
            rk = self.element_conservation.rk

            if n_species <= rk:
                raise ValueError(
                    f"This equation has {n_species} species but its elemental matrix already has rank "
                    f"{rk}: every species is needed just to reach full rank, so the only solution is "
                    "the trivial all-zero one — there is no species left to act as a known anchor."
                )

            all_species = self.all_species

            coeff_source = known_coeff if known_coeff is not None else self._known_coefficients
            coeff_source = {all_species._resolve_species(k).formula: v for k, v in coeff_source.items()}

            if independent_set is None:
                include_formulas = {all_species._resolve_species(s).formula for s in (include or [])}

                # A species named in `include` should win over it also being in
                # `coeff_source` (known, so auto-excluded by default) — otherwise
                # the two lists would contradict each other for no reason the
                # caller asked for.
                auto_exclude = [f for f in coeff_source if f not in include_formulas]

                sets = self.element_conservation.get_independent_sets(
                    include=include,
                    exclude=list(exclude or []) + auto_exclude,
                    strict_include=strict_include,
                    strict_exclude=strict_exclude
                )

                if not sets:
                    raise ValueError("No invertible set of species could be found to infer coefficients from.")

                independent_set = sets[0]

            else:
                independent_set = [all_species._resolve_species(s).formula for s in independent_set]

            result['independent_set'] = list(independent_set)

            overridden = set(independent_set) & set(coeff_source)
            if verbose and overridden:
                print(
                    f"Species {sorted(overridden)} have a known coefficient but are part of the "
                    "independent set: treated as unknowns to infer instead."
                )

            current = self.coefficients
            new_known = {f: coeff_source.get(f, current[f]) for f in current if f not in independent_set}

            result['known_set'] = list(new_known.keys())
            result['known_coeffs'] = dict(new_known)

            self._previous_coefficients = dict(current)

            self.update_coefficients(new_known, verbose=False, _clear=False)

            LHS, RHS = self.element_conservation._linear_sys(independent_set)
            x = np.linalg.solve(LHS, RHS)
            inferred = dict(zip(independent_set, x))

            mismatched = self._resolve_complex(inferred, verbose=verbose)
            result['sign_mismatch'] = sorted(mismatched) if mismatched else None

            solved_magnitude = {f: abs(v) for f, v in inferred.items()}

            self.update_coefficients({**new_known, **solved_magnitude}, verbose=False, _clear=False)

            self.balanced_coefficients = dict(self._known_coefficients)

            result['inferred_coeffs'] = {f: self._known_coefficients[f] for f in independent_set}
            result['success'] = True

        except Exception as e:
            result['message'] = str(e)
            raise

        finally:
            self.result_inference = ResultDict(result)


    def _resolve_complex(self, coeff: Dict[str, float], verbose: bool = True):
        """
        Check `coeff` (a {formula: signed value} dict, e.g. from solving the
        linear system) for species whose sign doesn't match their current
        reactant/product role, and move any such species to the other complex.

        The pre-move `reactants`/`products` are saved to `self._previous_reactants`/
        `self._previous_products` only the first time this happens (when
        `self._overriden_complexes` is still False) — once set, they keep the
        complexes as they were at instance creation instead of being overwritten
        by an intermediate, already-moved state from a prior call.
        """

        reactant_formulas = set(self.reactants.formula)

        mismatched = [
            f for f, v in coeff.items()
            if (f in reactant_formulas and v > 0) or (f not in reactant_formulas and v < 0)
        ]

        if not mismatched:
            return mismatched

        if verbose:
            print(
                f"Sign mismatch for species {sorted(mismatched)}: the inferred coefficient has the "
                "opposite sign expected from their reactant/product role; moving them to the other side."
            )

        if not self._overriden_complexes:
            self._previous_reactants = self.reactants
            self._previous_products = self.products

        to_products = {f for f in mismatched if f in reactant_formulas}
        to_reactants = set(mismatched) - to_products

        new_reactant_species = [s for s in self.reactants.species if s.formula not in to_products] + \
                                [self.products._resolve_species(f) for f in to_reactants]
        new_product_species = [s for s in self.products.species if s.formula not in to_reactants] + \
                               [self.reactants._resolve_species(f) for f in to_products]

        self.reactants = Complex(new_reactant_species, default_coeff=-1)
        self.products = Complex(new_product_species, default_coeff=1)

        self._overriden_complexes = True

        return mismatched


    def clear(self):

        if self._overriden_complexes:
            self.reactants = self._previous_reactants
            self.products = self._previous_products
            self._overriden_complexes = False

        self.result_inference = None

        self.update_coefficients({}, verbose=False, _clear=False)
            



    def __repr__(self):
        coeffs = self.coefficients

        def side(species):
            return " + ".join(f"{abs(coeffs[s.formula]):g} {s.formula}" for s in species)

        return f"{side(self.reactants.species)} = {side(self.products.species)}"


    def __mul__(self, scalar):
        """
        Scale this equation by a non-zero scalar, returning a new `Equation`.

        Known coefficients are scaled by `abs(scalar)` (default coefficients
        are untouched); a negative scalar also flips reactants and products.
        """

        if not isinstance(scalar, Real) or isinstance(scalar, bool):
            return NotImplemented

        if scalar == 0:
            raise ValueError("Cannot multiply an Equation by a zero scalar.")

        known_coefficients = {f: abs(v) * abs(scalar) for f, v in self._known_coefficients.items()}

        reactants, products = self.reactants, self.products
        if scalar < 0:
            reactants, products = products, reactants

        return Equation(reactants.species, products.species, known_coefficients=known_coefficients)


    def __rmul__(self, scalar):
        return self.__mul__(scalar)


    def __add__(self, other):
        """
        Combine two equations into a new one.

        Species unique to either equation carry over unchanged. For a species
        shared by both:
        - both known: their signed coefficients are summed (the sum's sign
          decides the new role); a sum of exactly zero drops the species
          entirely (it cancels out between the two equations).
        - only one known: that value is used as-is, the other's default coefficient
          is discarded.
        - neither known: keeps the shared side if both agree on it; if they
          conflict, `self`'s (the first equation's) side wins.
        """

        if not isinstance(other, Equation):
            return NotImplemented

        self_known = self._known_coefficients
        other_known = other._known_coefficients

        self_reactants = set(self.reactants.formula)
        self_products = set(self.products.formula)
        other_reactants = set(other.reactants.formula)
        other_products = set(other.products.formula)

        species_map = {
            s.formula: s for s in
            list(self.reactants.species) + list(self.products.species) +
            list(other.reactants.species) + list(other.products.species)
        }

        shared = (self_reactants | self_products) & (other_reactants | other_products)

        new_known = {}
        new_reactants, new_products = set(), set()

        for f in (self_reactants | self_products) - shared:
            if f in self_known:
                new_known[f] = self_known[f]
            (new_reactants if f in self_reactants else new_products).add(f)

        for f in (other_reactants | other_products) - shared:
            if f in other_known:
                new_known[f] = other_known[f]
            (new_reactants if f in other_reactants else new_products).add(f)

        for f in shared:
            if f in self_known and f in other_known:
                combined = self_known[f] + other_known[f]

                if combined == 0:
                    continue

                new_known[f] = combined
                (new_reactants if combined < 0 else new_products).add(f)

            elif f in self_known:
                new_known[f] = self_known[f]
                (new_reactants if self_known[f] < 0 else new_products).add(f)

            elif f in other_known:
                new_known[f] = other_known[f]
                (new_reactants if other_known[f] < 0 else new_products).add(f)

            else:
                # Neither known: `self`'s side is used whether it agrees with
                # `other`'s (keeping the shared side) or not (self takes over).
                (new_reactants if f in self_reactants else new_products).add(f)

        return Equation(
            [species_map[f] for f in new_reactants],
            [species_map[f] for f in new_products],
            known_coefficients=new_known
        )


    def __sub__(self, other):
        """`self - other` is `self + (-1) * other` (see `__add__`, `__mul__`)."""

        if not isinstance(other, Equation):
            return NotImplemented

        return self + (-1) * other


    def __len__(self):

        return len(self.all_species)


    def __getattr__(self, name):

        try:
            s = self.all_species._resolve_species(name)
        except ValueError:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

        return self.coefficients[s.formula]


class ResultDict:

    def __init__(self, dict):
        self.dict = dict

    def __getattr__(self, name):

        if name in self.dict:
            return self.dict[name]

        raise AttributeError(f"Available attributes are {list(self.dict.keys())}.")