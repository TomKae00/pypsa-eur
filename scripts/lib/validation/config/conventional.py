# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Conventional generators configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#conventional
"""

from pydantic import ConfigDict, Field

from scripts.lib.validation.config._base import ConfigModel


class ConventionalConfig(ConfigModel):
    """Configuration for `conventional` settings."""

    model_config = ConfigDict(extra="allow")

    unit_commitment: bool = Field(
        False,
        description="Allow the overwrite of ramp_limit_up, ramp_limit_start_up, ramp_limit_shut_down, p_min_pu, min_up_time, min_down_time, and start_up_cost of conventional generators. Refer to the CSV file 'unit_commitment.csv'.",
    )
    dynamic_fuel_price: bool = Field(
        False,
        description="Use time-varying fuel prices: the native monthly history or a prepared hourly CSV supplied with fuel_price_file.",
    )
    fuel_price_file: str | None = Field(
        None,
        description="Optional hourly CSV in EUR/MWh thermal with a gas column and optional carrier columns. Requires dynamic_fuel_price; supports {horizon}. Bypasses monthly price building and smoothing.",
    )
    fuel_price_rolling_window: int = Field(
        6,
        description="Monthly rolling mean window for fossil fuel prices smoothing.",
        ge=1,
    )
    nuclear: dict[str, str | float] = Field(
        default_factory=lambda: {"p_max_pu": "data/nuclear_p_max_pu.csv"},
        description="For any carrier/technology overwrite attributes as listed below.",
    )
