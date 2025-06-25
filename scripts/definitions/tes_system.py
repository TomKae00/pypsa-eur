# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from enum import Enum
from typing import Iterator, Tuple, List

from scripts.definitions.heat_system import HeatSystem
from scripts.definitions.tes_system_case import TesSystemCase


class TesSystem(Enum):
    """
    Enumeration representing different types of thermal energy storage (TES) systems.

    Attributes
    ----------
    TTES : str
        Tank Thermal Energy Storage.
    PTES : str
        Pit Thermal Energy Storage.
    ATES : str
        Aquifer Thermal Energy Storage.
    """

    TTES = "tank thermal energy storage"
    PTES = "pit thermal energy storage"
    ATES = "aquifer thermal energy storage"

    def __str__(self) -> str:
        """
        Returns the string representation of the TES system.
        """
        return self.value

    def component_name(
        self,
        boosting: TesSystemCase = TesSystemCase.BASE
    ) -> str:
        """
        Returns the TES system name, appending ' boosting' if BOOSTED.

        Parameters
        ----------
        boosting : TesSystemCase, optional
            One of TesSystemCase.BASE or TesSystemCase.BOOSTED.
            Default is BASE.

        Returns
        -------
        str
        """
        if boosting == TesSystemCase.BOOSTED:
            return f"{self.value} boosting"
        return self.value

    @staticmethod
    def boosting_cases(is_boosted: bool) -> List[TesSystemCase]:
        """
        Which boosting cases to iterate.

        - If is_boosted is True, returns [BASE, BOOSTED].
        - If is_boosted is False, returns [BASE].
        """
        if is_boosted:
            return [TesSystemCase.BASE, TesSystemCase.BOOSTED]
        return [TesSystemCase.BASE]

    @classmethod
    def _candidates_for(
        cls,
        heat_system: HeatSystem
    ) -> List["TesSystem"]:
        """
        Which TES types are valid for this heat system.

        Parameters
        ----------
        heat_system : HeatSystem

        Returns
        -------
        List[TesSystem]
            [TTES, PTES, ATES] if URBAN_CENTRAL, else [TTES].
        """
        if heat_system == HeatSystem.URBAN_CENTRAL:
            return list(cls)
        return [cls.TTES]

    @classmethod
    def _should_be_boosted(
        cls,
        heat_system: HeatSystem,
        tes_system: "TesSystem",
        supplemental_heating_storages: List[str]
    ) -> bool:
        """
        Decide whether this TES system should run with boosting.

        Parameters
        ----------
        heat_system : HeatSystem
        tes_system : TesSystem
        supplemental_heating_storages : List[str]
            From config:
            options["district_heating"]["supplemental_heating_storages"]

        Returns
        -------
        bool
            True if heat_system is URBAN_CENTRAL and
            tes_system.value is listed in supplemental_heating_storages.
        """
        return (
            heat_system == HeatSystem.URBAN_CENTRAL
            and tes_system.value in supplemental_heating_storages
        )

    @classmethod
    def variants_for(
        cls,
        heat_system: HeatSystem,
        supplemental_heating_storages: List[str],
    ) -> Iterator[Tuple["TesSystem", TesSystemCase]]:
        """
        Iterate all (TES system, boost case) pairs valid for a given
        heat system *and* your config’s supplemental_heating_storages.

        Yields
        ------
        Tuple[TesSystem, TesSystemCase]
        """
        # pick which TES types apply
        for tes_system in cls._candidates_for(heat_system):
            # decide if boosting is enabled for this pair
            boosted_enabled = cls._should_be_boosted(
                heat_system, tes_system, supplemental_heating_storages
            )
            # yield base—and boosted if allowed
            for boost_case in cls.boosting_cases(boosted_enabled):
                yield tes_system, boost_case

    @classmethod
    def names_for(
        cls,
        heat_system: HeatSystem,
        supplemental_heating_storages: List[str],
    ) -> Iterator[str]:
        """
        Iterate the exact TES‐system names you’ll use in your loops,
        automatically including both base and boosted variants when applicable.

        Yields
        ------
        str
            The component_name() of each (tes_system, boost_case) pair.
        """
        for tes_system, boost_case in cls.variants_for(
            heat_system, supplemental_heating_storages
        ):
            yield tes_system.component_name(boost_case)