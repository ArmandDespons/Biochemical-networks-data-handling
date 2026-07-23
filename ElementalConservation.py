from itertools import combinations
from typing import Dict, List, Optional, Set, Tuple, Union

import numpy as np
from pandas import DataFrame

from .Equation import Equation
from .Species import Species

__all__ = ["ElementalConservation"]


class ElementalConservation:

    def __init__(self, eq: Equation):
        """
        Build the elemental conservation matrix of an Equation instance from the
        coefficient-weighted stoichiometry of its merged `all_species` complex.
        """

        self.matrix, self.cols, self.rows = eq.all_species.elemental_composition

        target_rank = np.linalg.matrix_rank(self.matrix)

        kept_idx = []
        independent_rows = {}
        for element, idx in self.rows.items():
            trial_idx = kept_idx + [idx]

            if np.linalg.matrix_rank(self.matrix[trial_idx, :]) > len(kept_idx):
                kept_idx = trial_idx
                independent_rows[element] = len(independent_rows)

                if len(kept_idx) == target_rank:
                    break

        # `L`: `self.matrix` reduced to its linearly independent rows only (full row rank).
        # `independent_rows`: {element: new row index in `L`}, same dict format as `self.rows`.
        self.L = self.matrix[kept_idx, :]
        self.independent_rows = independent_rows

        self.elements = list(self.rows.keys())

        # Kept to resolve `get_independent_sets`' includes/exclude (formula/Species/name/alias)
        # against the exact same species set the matrix/cols were built from.
        self._complex = eq.all_species


    @property
    def el(self):
        return self.elements


    @property
    def to_pandas(self):
        return DataFrame(self.matrix, index=self.elements, columns=self.cols.keys())


    @property
    def rk(self):
        return np.linalg.matrix_rank(self.matrix)

    def get_independent_sets(
            self,
            include: Optional[List[Union[str, Species]]] = None,
            exclude: Optional[List[Union[str, Species]]] = None,
            strict_include: bool = True,
            strict_exclude: bool = True
            ) -> List[List[str]]:
        """
        Enumerate every set of species (one per column of `L`, sized to its
        number of rows) whose columns of `L` are linearly independent.

        `include` species are forced into every returned set; `exclude` species
        are never allowed in any of them. Each returned group is a list of
        formulas ordered the same way as the keys of `self.cols`.

        By default (`strict_include`/`strict_exclude` True), an `include`/
        `exclude` that cannot be fully honored (too many forced species, or
        too few candidates left once excluded) raises. Set either to False to
        instead treat that list as a priority order and silently keep only as
        many entries from its front as can actually be honored.
        """

        include = include or []
        exclude = exclude or []

        n_rows, n_cols = self.L.shape

        if n_rows > n_cols:
            raise ValueError(
                f"The elemental matrix has more rows ({n_rows}) than columns ({n_cols}): "
                "no set of species can span every row."
            )

        include_all = list(dict.fromkeys(self._complex._resolve_species(s).formula for s in include))
        exclude_all = list(dict.fromkeys(self._complex._resolve_species(s).formula for s in exclude))

        overlap = set(include_all) & set(exclude_all)
        if overlap:
            raise ValueError(f"Species cannot be both included and excluded: {sorted(overlap)}")

        if not strict_exclude:
            exclude_all = exclude_all[:n_cols - n_rows]

        exclude_formulas = set(exclude_all)
        candidates = [f for f in self.cols if f not in exclude_formulas]

        if len(candidates) < n_rows:
            raise ValueError(
                f"Excluding {sorted(exclude_formulas)} leaves only {len(candidates)} candidate "
                f"species, fewer than the {n_rows} independent rows of the elemental matrix."
            )

        if not strict_include:
            include_all = include_all[:n_rows]

        include_formulas = set(include_all)

        if len(include_formulas) > n_rows:
            raise ValueError(
                f"'include' already specifies more species ({len(include_formulas)}) than "
                f"there are independent rows ({n_rows}) in the elemental matrix."
            )

        free_candidates = [f for f in candidates if f not in include_formulas]
        n_free = n_rows - len(include_formulas)

        independent_sets = []
        for combo in combinations(free_candidates, n_free):
            subset = include_formulas | set(combo)
            ordered_subset = [f for f in self.cols if f in subset]
            col_idx = [self.cols[f] for f in ordered_subset]
            submatrix = self.L[:, col_idx]

            is_square = submatrix.shape[0] == submatrix.shape[1]
            is_invertible = is_square and np.linalg.matrix_rank(submatrix) == submatrix.shape[0]

            if is_invertible:
                independent_sets.append(ordered_subset)

        return independent_sets


    def _linear_sys(
            self,
            independent_set: Union[List[Union[str, Species]], Set[Union[str, Species]]]
            ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Build the LHS, RHS pair of the linear system for `independent_set`.

        LHS is a rk x rk matrix (rk = number of rows of `L`) whose columns span
        `independent_set`, in the same order as `independent_set` itself when
        it is a list, or in `self.cols` order (a set has none of its own) when
        it is a set. Its rows follow `self.independent_rows`, each the absolute
        value of the matching row of `L`, except the 'Charge' row (if any)
        which holds each species' raw charge instead.

        RHS is the 1-d array obtained by summing every column of `L` not in
        `independent_set`.
        """

        rk = self.L.shape[0]

        resolved = [self._complex._resolve_species(s).formula for s in independent_set]
        formulas = [f for f in self.cols if f in resolved] if isinstance(independent_set, set) else resolved

        if len(set(formulas)) != rk:
            raise ValueError(
                f"'independent_set' must contain exactly {rk} distinct species (the rank of L), "
                f"got {len(set(formulas))}: {sorted(set(formulas))}"
            )

        in_idx = [self.cols[f] for f in formulas]
        out_idx = [i for i in self.cols.values() if i not in in_idx]

        species = [self._complex._resolve_species(f) for f in formulas]

        LHS = np.empty((rk, rk))
        for element, row in self.independent_rows.items():
            if element == "Charge":
                LHS[row, :] = [s.charge for s in species]
            else:
                LHS[row, :] = np.abs(self.L[row, in_idx])

        RHS = -self.L[:, out_idx].sum(axis=1)

        return LHS, RHS


    def __getattr__(self, name):
        if name in self.elements:
            return self.matrix[self.rows[name], :]
        
        else:
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")
        

    def __repr__(self):
        return DataFrame.__repr__(self.to_pandas)
