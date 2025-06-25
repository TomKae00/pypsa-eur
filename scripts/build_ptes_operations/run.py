# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Approximate the top temperature of the pit thermal energy storage (PTES), ensuring that the temperature does not
exceed the operational limit.

Determine whether supplemental heating is needed. A binary indicator is generated:
    - 1: The forward temperature is less than or equal to the TES maximum; direct usage is possible.
    - 0: The forward temperature exceeds the TES maximum; supplemental heating (e.g., via a heat pump) is required.

Calculate dynamic PTES capacity profiles based on district heating forward and return flow temperatures.
The linear relation between temperature difference and capacity is taken from Sorknaes (2018).

The capacity of thermal energy storage systems varies with the temperature difference
between the forward and return flows in district heating networks assuming a direct
integration of the storage. This script calculates normalized capacity factors (e_max_pu)
for PTES systems based on these temperature differences.

Relevant Settings
-----------------
.. code:: yaml
    sector
        district_heating:
            ptes:
                dynamic_ptes_capacity:
                supplemental_heating:
                    enable:
                max_top_temperature:
                min_bottom_temperature:

Inputs
------
- `resources/<run_name>/forward_temperature.nc`
    Forward temperature profiles for the district heating networks.
- `resources/<run_name>/central_heating_return_temperature_profiles.nc`:
    Return temperature profiles for the district heating networks.

Outputs
-------
- `resources/<run_name>/ptes_top_temperature_profiles.nc`
    Clipped PTES top temperature profile (in °C).
- `resources/<run_name>/ptes_supplemental_heating_required.nc`
    Binary indicator for additional heating (1 = direct PTES use, 0 = supplemental heating required).
- `resources/<run_name>/ptes_e_max_pu_profiles.nc`
    Normalized PTES capacity profiles.
- `resources/<run_name>/ptes_temperature_boost_ratio_profiles.nc`
    Ratio of PTES charge that requires additional heating due to temperature differences.

Source
------
Sorknæs, P. 2018. "Simulation method for a pit seasonal thermal energy storage system with a heat pump in a district heating system", Energy, Volume 152, https://doi.org/10.1016/j.energy.2018.03.152.
Approximate thermal energy storage (TES) top temperature and identify need for supplemental heating.
"""

import logging

import xarray as xr
from _helpers import set_scenario_config

from scripts.build_tes_operations.tes_temperature_approximator import (
    TesTemperatureApproximator,
)
from scripts.definitions.tes_system import TesSystem

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tes_operations",
            clusters=8,
            planning_horizons="2050",
        )

    set_scenario_config(snakemake)

    # Load temperature profiles
    logger.info(
        "Loading district heating temperature profiles and constructing PTES temperature approximator"
    )
    forward_temperature = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    return_temperature = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    # Compute temperature-dependent profiles for each TES system separately
    top_temperatures = []
    direct_utilisation_profiles = []
    e_max_pu_profiles = []
    boost_ratios = []

    for tes_system in TesSystem:
        approximator = TesTemperatureApproximator(
            forward_temperature=forward_temperature,
            return_temperature=return_temperature,
            max_tes_top_temperature=snakemake.params.max_ptes_top_temperature,
            min_tes_bottom_temperature=snakemake.params.min_ptes_bottom_temperature,
        )

        top_temperatures.append(
            approximator.top_temperature.expand_dims(tes_system=[tes_system.name])
        )
        direct_utilisation_profiles.append(
            approximator.direct_utilisation_profile.expand_dims(
                tes_system=[tes_system.name]
            )
        )
        e_max_pu_profiles.append(
            approximator.e_max_pu.expand_dims(tes_system=[tes_system.name])
        )
        boost_ratios.append(
            approximator.temperature_boost_ratio.expand_dims(
                tes_system=[tes_system.name]
            )
        )

    top_temperatures = xr.concat(top_temperatures, dim="tes_system")
    direct_utilisation_profiles = xr.concat(
        direct_utilisation_profiles, dim="tes_system"
    )
    e_max_pu_profiles = xr.concat(e_max_pu_profiles, dim="tes_system")
    boost_ratios = xr.concat(boost_ratios, dim="tes_system")

    logger.info(
        f"Saving TES top temperature profile to {snakemake.output.tes_top_temperature_profiles}"
    )
    top_temperatures.to_netcdf(snakemake.output.tes_top_temperature_profiles)

    logger.info(
        f"Saving TES direct utilisation profile to {snakemake.output.tes_direct_utilisation_profiles}"
    )
    direct_utilisation_profiles.to_netcdf(
        snakemake.output.tes_direct_utilisation_profiles
    )

    logger.info("Calculating dynamic TES capacity profiles")
    logger.info(
        f"Saving TES capacity profiles to {snakemake.output.tes_e_max_pu_profiles}"
    )
    e_max_pu_profiles.to_netcdf(snakemake.output.tes_e_max_pu_profiles)

    logger.info(
        f"Saving TES reheat ratio profiles to {snakemake.output.tes_temperature_boost_ratio_profiles}"
    )
    boost_ratios.to_netcdf(snakemake.output.tes_temperature_boost_ratio_profiles)