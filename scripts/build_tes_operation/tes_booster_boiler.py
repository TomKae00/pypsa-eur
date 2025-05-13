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
        forward_temperature: xr.DataArray,
        return_temperature: xr.DataArray,
        clipped_top_temperature: xr.DataArray,
        supplemental_heating_profile: xr.DataArray,
    ):
        """
        Initialize the ComputeReheatRatio calculation.

        Parameters
        ----------
        forward_temperature : xr.DataArray
            The forward temperature profile (°C) from the district heating network.
        return_temperature : xr.DataArray
            The return temperature profile (°C) from the district heating network.
        clipped_top_temperature : xr.DataArray
            The forward temperature profile clipped at the maximum PTES temperature (°C).
        supplemental_heating_profile : xr.DataArray
            Whether supplemental heating is required. If True, the computed ratio
            will be squared to reflect additional reheat effort (default: False).
        """
        self.forward_temperature = forward_temperature
        self.return_temperature = return_temperature
        self.clipped_top_temperature = clipped_top_temperature
        self.supplemental_heating = supplemental_heating_profile

    @property
    def reheat_ratio(self) -> xr.DataArray:
        # hier noch herausfinden, wo die +1 genua herkommt -> mathematische herleitung
        # überlegen, wie die capital cost berechnet werden müssen. Efficiency wird zwar mit +1 beschrieben z.b 4 MW brauhen 2,4 MW nacherhitzung mit Kessel bei VLT 120 max ptes temp 90 und RLT 40
        # , aber der Kessel hat dann eigentlich eine installierte Leistung von 2,4 MW und nicht 4 wie es jetzt im Modell der Fall wäre. Daher capital cost * (base_ratio - 1), das müsste das problem lösen. Aber nochmal überprüfen, ob das stimmt.
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
        base_ratio = 1 + (
            self.forward_temperature - self.clipped_top_temperature
        ) / (
            self.clipped_top_temperature - self.return_temperature
        )   #* self.supplemental_heating
        return base_ratio
