# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Approximate the top temperature of the pit thermal energy storage (PTES), ensuring that the temperature does not
exceed the operational limit.

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
                dynamic_capacity:
                discharger_temperature_boosting_required:
                charger_temperature_boosting_required:
                max_top_temperature:
                min_bottom_temperature:

Inputs
------
- `resources/<run_name>/forward_temperature.nc`
    Forward temperature profiles for the district heating networks.
- `resources/<run_name>/central_heating_return_temperature_profiles.nc`:
    Return temperature profiles for the district heating networks.
- `resources/<run_name>/ptes_temperature_boost_ratio_profiles.nc`
    Ratio of PTES charge that requires additional heating due to temperature differences.

Outputs
-------
- `resources/<run_name>/ptes_top_temperature_profiles.nc`
    Clipped PTES top temperature profile (in °C).
- `resources/<run_name>/ptes_e_max_pu_profiles.nc`
    Normalized PTES capacity profiles.
- `ptes_temperature_boost_ratio_profiles` (netCDF):
    Charging temperature boost ratio time series.
- `ptes_forward_temperature_boost_ratio_profiles` (netCDF):
    Forward flow temperature boost ratio time series.

Source
------
Sorknæs, P. 2018. "Simulation method for a pit seasonal thermal energy storage system with a heat pump in a district heating system", Energy, Volume 152, https://doi.org/10.1016/j.energy.2018.03.152.
Approximate thermal energy storage (TES) top temperature and identify need for supplemental heating.
"""

import logging

import xarray as xr
import pandas as pd
from scripts._helpers import set_scenario_config

from scripts.build_tes_operations.tes_temperature_approximator import (
    TesTemperatureApproximator,
    TesTemperatureMode,
)

from scripts.definitions.tes_system  import TesSystem

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake(
            "build_tes_operations",
            clusters=8,
            planning_horizons="2030",
        )

    set_scenario_config(snakemake)

#    if (snakemake.params.charge_boosting_required and
#            TesTemperatureMode(snakemake.params.tes_temperature_profile) is TesTemperatureMode.DYNAMIC):
#        raise ValueError(
#            "Charger boosting cannot be used with 'dynamic' temperature profile"
#        )

    central_heating_forward_temperature = xr.open_dataarray(
        snakemake.input.central_heating_forward_temperature_profiles
    )
    central_heating_return_temperature = xr.open_dataarray(
        snakemake.input.central_heating_return_temperature_profiles
    )

    variable_to_output = {
        "top_temperature": snakemake.output.tes_top_temperature_profile,
        "e_max_pu": snakemake.output.tes_e_max_pu_profile,
        "boost_per_discharge": snakemake.output.boost_per_discharge_profile,
        "boost_per_charge": snakemake.output.boost_per_charge_profile,
    }

    profiles_all_tes_systems = {v: [] for v in variable_to_output}
    tes_system_all_names = []

    for tes_system in TesSystem:
        params_this_tes_system = snakemake.params.tes[tes_system.name]
        if not params_this_tes_system.get("enable", True):
            continue

        tes_temperature_approximator = TesTemperatureApproximator(
            forward_temperature=central_heating_forward_temperature,
            return_temperature=central_heating_return_temperature,
            max_top_temperature=params_this_tes_system["max_top_temperature"],
            min_bottom_temperature=params_this_tes_system["min_bottom_temperature"],
            temperature_profile=TesTemperatureMode(params_this_tes_system["temperature_profile"]),
            charge_boosting_required=params_this_tes_system["charge_boosting_required"],
            discharge_boosting_required=params_this_tes_system["discharge_boosting_required"],
            dynamic_capacity=params_this_tes_system["dynamic_capacity"],
        )

        for variable in profiles_all_tes_systems:
            profiles_all_tes_systems[variable].append(getattr(tes_temperature_approximator, variable))
        tes_system_all_names.append(tes_system.value)

        # hier muss das tes_system anders bennant werden, damit es richtig ausgewählt wird
    tes_system_index = pd.Index(tes_system_all_names, name="tes_system")

    for variable, dataarrays in profiles_all_tes_systems.items():
        xr.concat(dataarrays, dim=tes_system_index).to_netcdf(variable_to_output[variable])
