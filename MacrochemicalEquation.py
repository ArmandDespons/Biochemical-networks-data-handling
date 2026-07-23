from typing import Dict, List, Optional, Union
from collections.abc import Iterable
from ast import literal_eval

from pandas import DataFrame

from .Species import Species
from .Equation import Equation


__all__ = ["MacrochemEquation"]


class MacrochemEquation(Equation):
    """
    An `Equation` specialized for microbial growth-yield stoichiometry.

    The first product is distinguished as `biomass` — named `'biomass'` so
    `Species.standard_enthaply`/`chemical_potential` dispatch to its
    elemental-composition-based formulas — addressable via a friendlier
    `biomass_alias` wherever a species key is otherwise expected (`yields`,
    `infer_yields`, `list_species`), in place of its real (often unwieldy)
    formula.

    `yields`/`y` mirror `coefficients`/`coeff` (alias substituted for
    biomass); `infer_yields` mirrors `infer_coeff`. Also tracks
    `electron_donor`/`electron_acceptor`/`carbon_source` (`Ed`/`Ea`/`Cs`) and
    `metabolic_type` as validated species/metadata references.
    """

    def __init__(
        self,
        substrates: List[Union[str, Species]],
        products: List[Union[str, Species]],
        biomass_stoichiometry: Optional[Dict[str, float]] = None,
        electron_donor: Optional[Union[str, Species]] = None,
        electron_acceptor: Optional[Union[str, Species]] = None,
        carbon_source: Optional[Union[str, Species]] = None,
        metabolic_type: Optional[str] = None,
        yields: Optional[Dict[Union[str, Species], float]] = None,
        biomass_alias: Optional[str] = 'CH1.8O0.5N0.2'
    ):

        if yields is None:
            yields = {}

        # Either build biomass fresh from its elemental composition, or tag an
        # existing product as biomass in place: `.name` is set to the literal
        # 'biomass', the hook `Species`'s thermodynamic properties dispatch on.
        if biomass_stoichiometry:
            _biomass_stoichiometry = {'name': 'biomass'}

            for element in ['C', 'H', 'O', 'N', 'P', 'S']:
                for key in biomass_stoichiometry:
                    if element in key:
                        _biomass_stoichiometry[element] = biomass_stoichiometry.get(key)

            if 'C' not in _biomass_stoichiometry:
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

        self._biomass_alias = biomass_alias if isinstance(biomass_alias, str) else self.biomass.formula

        super().__init__(
            reactants=substrates,
            products=products,
            known_coefficients=self._parse_biomass_alias(yields, self._biomass_alias, self.biomass.formula)
        )

        self.substrates = self.reactants

        self._electron_donor = self._validate_species(electron_donor) if electron_donor else None
        self._electron_acceptor = self._validate_species(electron_acceptor) if electron_acceptor else None
        self._carbon_source = self._validate_species(carbon_source) if carbon_source else None

        self._metabolic_type = metabolic_type if metabolic_type else None


    # ----- Construction from string/dict/DataFrame -----

    @classmethod
    def from_string(cls, eqn_str: str, **kwargs):

        eqn = Equation.from_string(eqn_str)

        return cls(
            substrates=eqn.reactants.species,
            products=eqn.products.species,
            yields=eqn._known_coefficients,
            **kwargs
        )


    @classmethod
    def from_dict(cls, coefficients: Dict[Union[str, Species], float], biomass_alias: Optional[str] = 'CH1.8O0.5N0.2',
                  biomass_key: Optional[str] = None, **kwargs):

        if biomass_key is None:
            biomass_key = biomass_alias

        if biomass_key not in coefficients:
            raise ValueError(
                f"Biomass key '{biomass_key}' not found among the dictionary keys: {list(coefficients.keys())}"
            )

        _element_biomass = ['C biomass', 'H biomass', 'O biomass', 'N biomass', 'P biomass', 'S biomass']
        _dict = dict(coefficients)
        biomass_stoich = {k: _dict.pop(k) for k in _element_biomass if k in coefficients}

        if not biomass_stoich:
            biomass_stoich = None

        eqn = Equation.from_dict({biomass_alias: _dict.pop(biomass_key), **_dict})

        return cls(
            substrates=eqn.reactants.species,
            products=eqn.products.species,
            biomass_stoichiometry=biomass_stoich,
            yields=eqn._known_coefficients,
            biomass_alias=biomass_alias,
            **kwargs
        )


    @classmethod
    def from_parsed_df(cls, df: DataFrame, index: Optional[int] = None, **kwargs):

        if index is None:
            return {i: cls.from_parsed_df(df, i, **kwargs) for i in df.index}

        elif isinstance(index, Iterable):
            return {i: cls.from_parsed_df(df, i, **kwargs) for i in index}

        elif isinstance(index, int):

            _necessary_cols = ['Substrates', 'Products', 'Metabolic type', 'Electron donor formula',
                                'Electron acceptor formula', 'Carbon source formula', 'H biomass', 'O biomass', 'N biomass']
            if not all(col in df.columns for col in _necessary_cols):
                raise ValueError(f"The parsed DataFrame must contain columns: {_necessary_cols}.")

            substrates = df.loc[index, 'Substrates'] if isinstance(df.loc[index, 'Substrates'], list) \
                else literal_eval(df.loc[index, 'Substrates'])
            products = df.loc[index, 'Products'] if isinstance(df.loc[index, 'Products'], list) \
                else literal_eval(df.loc[index, 'Products'])

            Ed = df.loc[index, 'Electron donor formula']
            Ea = df.loc[index, 'Electron acceptor formula']
            Cs = df.loc[index, 'Carbon source formula']

            metabolic_type = df.loc[index, 'Metabolic type']

            biomass_stoichiometry = df.loc[index, ['H biomass', 'O biomass', 'N biomass']].rename(
                index={'H biomass': 'H', 'O biomass': 'O', 'N biomass': 'N'}
            ).to_dict()

            return cls(
                substrates=substrates,
                products=products,
                metabolic_type=metabolic_type,
                biomass_stoichiometry=biomass_stoichiometry,
                electron_donor=Ed,
                electron_acceptor=Ea,
                carbon_source=Cs,
                **kwargs
            )

        else:
            raise TypeError("'index' should be either an iterable object or an int.")


    # ----- Species resolution & biomass aliasing -----

    def _validate_species(self, key: Optional[Union[str, Species]]) -> Optional[Species]:
        """Resolve `key` (formula, name, alias, or Species instance) against this equation's own species."""

        if key is None:
            return None

        try:
            return self.all_species._resolve_species(key)
        except ValueError:
            raise ValueError(
                f"Role Error: {key} must be one of the reactants or products: {sorted(self.all_species.formula)}"
            )


    @staticmethod
    def _use_biomass_alias(coefficients: dict, biomass_alias: str, biomass_formula: str) -> dict:
        """Rename the real biomass-formula key to `biomass_alias`, for display."""

        out = dict(coefficients)
        out[biomass_alias] = out.pop(biomass_formula)

        return out


    @staticmethod
    def _parse_biomass_alias(coefficients: dict, biomass_alias: str, biomass_formula: str) -> dict:
        """Rename a `biomass_alias` key back to the real biomass formula, for internal use."""

        if biomass_alias in coefficients:
            out = dict(coefficients)
            out[biomass_formula] = out.pop(biomass_alias)

            return out

        return coefficients

    # ----- Infer yields -----

    def infer_yields(
            self,
            known_coeff: Optional[Dict[Union[str, Species], float]] = None,
            independent_set: Optional[List[Union[str, Species]]] = None,
            include: Optional[List[Union[str, Species]]] = None,
            exclude: Optional[List[Union[str, Species]]] = None,
            strict_include: bool = False,
            strict_exclude: bool = False,
            verbose: bool = True
            ):
        """`infer_coeff`, with `known_coeff`'s keys and every species list accepting `biomass_alias`."""

        def _swap(items):
            return [self.biomass.formula if s == self._biomass_alias else s for s in items] if items else items

        if known_coeff:
            known_coeff = self._parse_biomass_alias(dict(known_coeff), self._biomass_alias, self.biomass.formula)

        return self.infer_coeff(
            known_coeff=known_coeff,
            independent_set=_swap(independent_set),
            include=_swap(include),
            exclude=_swap(exclude),
            strict_include=strict_include,
            strict_exclude=strict_exclude,
            verbose=verbose
        )


    # ----- Metabolic type -----

    @property
    def metabolic_type(self) -> Optional[str]:
        return self._metabolic_type

    @metabolic_type.setter
    def metabolic_type(self, value: Optional[str]):
        if value is None:
            self._metabolic_type = None
            return

        if not isinstance(value, str):
            raise TypeError(f"Metabolic type must be a string. Got {type(value)} instead.")

        self._metabolic_type = value


    # ----- Electron donor (Ed) -----

    @property
    def electron_donor(self) -> Optional[Species]:
        return self._electron_donor

    @electron_donor.setter
    def electron_donor(self, value: Optional[Union[str, Species]]):
        self._electron_donor = self._validate_species(value)

    @property
    def Ed(self):
        return self.electron_donor

    @Ed.setter
    def Ed(self, value):
        self.electron_donor = value


    # ----- Electron acceptor (Ea) -----

    @property
    def electron_acceptor(self) -> Optional[Species]:
        return self._electron_acceptor

    @electron_acceptor.setter
    def electron_acceptor(self, value: Optional[Union[str, Species]]):
        self._electron_acceptor = self._validate_species(value)

    @property
    def Ea(self):
        return self.electron_acceptor

    @Ea.setter
    def Ea(self, value):
        self.electron_acceptor = value


    # ----- Carbon source (Cs) -----

    @property
    def carbon_source(self) -> Optional[Species]:
        return self._carbon_source

    @carbon_source.setter
    def carbon_source(self, value: Optional[Union[str, Species]]):
        self._carbon_source = self._validate_species(value)

    @property
    def Cs(self):
        return self.carbon_source

    @Cs.setter
    def Cs(self, value):
        self.carbon_source = value


    # ----- Yields (y) -----

    @property
    def _yields(self) -> Dict[str, float]:
        return self._known_coefficients

    @property
    def yields(self) -> Dict[str, float]:
        return self._use_biomass_alias(self.coefficients, self._biomass_alias, self.biomass.formula)

    @yields.setter
    def yields(self, value: Dict[Union[str, Species], float]):
        self.update_coefficients(self._parse_biomass_alias(dict(value), self._biomass_alias, self.biomass.formula))

    @property
    def y(self):
        return self.yields

    @y.setter
    def y(self, value):
        self.yields = value


    def __repr__(self):
        base = super().__repr__()

        Ed = self._electron_donor.formula if self._electron_donor else 'Unknown'
        Ea = self._electron_acceptor.formula if self._electron_acceptor else 'Unknown'
        Cs = self._carbon_source.formula if self._carbon_source else 'Unknown'
        MT = self._metabolic_type if self._metabolic_type else 'Unknown'

        return f"{base}" + "\n \n" + f"Metabolic type: {MT}" + " | " + f"Electron donor: {Ed}" + \
               " | " + f"Electron acceptor: {Ea}" + " | " + f"Carbon source: {Cs}"
