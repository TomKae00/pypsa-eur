# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

from enum import Enum

class TesSystemCase(Enum):
    """
    Enumeration representing TES system modeling cases.

    Attributes
    ----------
    BASE : str
        TES system without boosting.
    BOOSTED : str
        TES system with boosting (higher output temperature via supplemental heating).
    """

    BASE = "base"
    BOOSTED = "boosted"

    def __str__(self) -> str:
        """
        Returns the string representation of the TES system case.

        Returns
        -------
        str
            String representation of the TES system case.
        """
        return self.value
