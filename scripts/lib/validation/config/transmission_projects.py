# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT

"""
Transmission projects configuration.

See docs in https://pypsa-eur.readthedocs.io/en/latest/configuration.html#transmission-projects
"""

from typing import Literal

from pydantic import BaseModel, Field

from scripts.lib.validation.config._base import ConfigModel


class _IncludeConfig(ConfigModel):
    """Configuration for `transmission_projects.include` settings."""

    tyndp2020: bool = Field(
        True,
        description="Whether to integrate the TYNDP 2020 dataset.",
    )
    tyndp2026: bool = Field(
        False,
        description="Include the pinned TYNDP 2026 subset for the Bornholm 2035 study; see doc/bornholm_2035_grid.md.",
    )
    nep: bool = Field(
        True,
        description="Whether to integrate the German network development plan dataset.",
    )
    manual: bool = Field(
        True,
        description="Whether to integrate the manually added transmission projects. They are taken from the previously existing links_tyndp.csv file.",
    )


class TransmissionProjectsConfig(BaseModel):
    """Configuration for `transmission_projects` settings."""

    enable: bool = Field(
        True,
        description="Whether to integrate this transmission projects or not.",
    )
    reference_year: int | None = Field(
        None,
        description="Optional static grid year. Remove not-yet-built and expired AC/DC assets after applying projects, before clustering. Intended for overnight studies; leave unset for multi-period grids.",
    )
    include: _IncludeConfig = Field(
        default_factory=_IncludeConfig,
        description="Name of the transmission projects. They should be unique and have to be provided in the `data/transmission_projects` folder.",
    )
    skip: list[str] = Field(
        default_factory=lambda: ["upgraded_lines", "upgraded_links"],
        description="Type of lines to skip from all transmission projects. Possible values are: `upgraded_lines`, `upgraded_links`, `new_lines`, `new_links`.",
    )
    status: list[str] | dict[str, list[str]] = Field(
        default_factory=lambda: ["under_construction", "in_permitting", "confirmed"],
        description="Status to include into the model as list or as dict with name of project and status to include. Possible values for status are `under_construction`, `in_permitting`, `confirmed`, `planned_not_yet_permitted`, `under_consideration`.",
    )
    new_link_capacity: Literal["zero", "keep"] = Field(
        "zero",
        description="Whether to set the new link capacity to the provided capacity or set it to zero.",
    )
