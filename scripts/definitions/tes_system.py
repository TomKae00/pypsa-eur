# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from enum import Enum

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

    TTES = "tank_thermal_energy_storage"
    PTES = "pit_thermal_energy_storage"
    ATES = "aquifer_thermal_energy_storage"

    def __str__(self) -> str:
        """
        Returns the string representation of the TES system.
        """
        return self.value