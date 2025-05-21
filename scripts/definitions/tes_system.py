# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from enum import Enum

from scripts.definitions.tes_system_case import TesSystemCase # gucken wie genau bennen

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

    Methods
    -------
    __str__()
        Returns the string representation of the TES system.
    component_name(boosting)
        Returns the TES system name, with ' boosting' appended if boosting is enabled.
    """

    TTES = "tank thermal energy storage"
    PTES = "pit thermal energy storage"
    ATES = "aquifer thermal energy storage"

    def __str__(self) -> str:
        """
        Returns the string representation of the TES system.

        Returns
        -------
        str
            The string representation of the TES system.
        """
        return self.value


    def component_name(self, boosting: TesSystemCase = TesSystemCase.BASE) -> str: # mal schauen ob das case oder boosting genommen wird
        """
        Returns the name of the TES system, with ' boosting' appended if boosting is enabled.

        Parameters
        ----------
        boosting : Boosting, optional
            Whether boosting is enabled. Default is Boosting.DISABLED.

        Returns
        -------
        str
            The TES system name, with or without ' boosting' suffix.
        """
        if boosting == TesSystemCase.BOOSTED:
            return f"{self.value} boosting"
        else:
            return self.value


    @staticmethod
    def boosting_cases(is_boosted: bool):
        if is_boosted:
            return [TesSystemCase.BASE, TesSystemCase.BOOSTED]
        else:
            return [TesSystemCase.BASE]