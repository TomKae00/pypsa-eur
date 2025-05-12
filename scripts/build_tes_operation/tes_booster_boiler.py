# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

import xarray as xr


class ComputeReheatRatio:
    """
    Computes the ratio of utilized temperature lift to remaining lift in a district
    heating system, with optional supplemental heating adjustment.

    Attributes
    ----------
    forward_temperature_celsius : xr.DataArray
        The forward temperature profile (°C) from the district heating network.
    bottom_temperature : xr.DataArray
        The return-flow temperature profile (°C).
    clipped_top_temperature : xr.DataArray
        The forward temperature clipped to the maximum allowable PTES temperature.
    supplemental_heating : bool
        If True, applies an additional multiplication of the computed ratio to
        account for supplemental reheat requirements.
    """

    def __init__(
        self,
        forward_temperature_celsius: xr.DataArray,
        bottom_temperature: xr.DataArray,
        clipped_top_temperature: xr.DataArray,
        supplemental_heating: bool = False,
    ):
        """
        Initialize the ComputeReheatRatio calculation.

        Parameters
        ----------
        forward_temperature_celsius : xr.DataArray
            The forward temperature profile (°C) from the district heating network.
        bottom_temperature : xr.DataArray
            The return-flow temperature profile (°C).
        clipped_top_temperature : xr.DataArray
            The forward temperature profile clipped at the maximum PTES temperature (°C).
        supplemental_heating : bool, optional
            Whether supplemental heating is required. If True, the computed ratio
            will be squared to reflect additional reheat effort (default: False).
        """
        self.forward_temperature = forward_temperature_celsius
        self.bottom_temperature = bottom_temperature
        self.clipped_top_temperature = clipped_top_temperature
        self.supplemental_heating = supplemental_heating

    @property
    def reheat_ratio(self) -> xr.DataArray:
        """
        Calculate the reheat ratio:

            (clipped_top_temperature - bottom_temperature)
            -----------------------------------------------
            (forward_temperature_celsius - clipped_top_temperature)

        If supplemental_heating is True, the ratio will be squared.

        Returns
        -------
        xr.DataArray
            The reheat ratio profile; squared if supplemental heating is enabled.
        """
        base_ratio = (
            self.clipped_top_temperature - self.bottom_temperature
        ) / (
            self.forward_temperature - self.clipped_top_temperature
        )
        if self.supplemental_heating:
            return base_ratio * base_ratio
        return base_ratio
