# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
# SPDX-License-Identifier: MIT
"""Native pipeline switches retain local sector coupling and existing gas assets."""

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pypsa
import pytest
import yaml

from scripts.prepare_sector_network import add_h2_gas_infrastructure


@pytest.mark.parametrize("connectivity,hydrogen", [(0, False), (1, False), (0, True)])
def test_pipeline_candidates_are_optional(tmp_path, connectivity, hydrogen):
    options = yaml.safe_load(Path("config/config.default.yaml").read_text())["sector"]
    options.update(
        gas_network=True,
        gas_network_connectivity_upgrade=connectivity,
        H2_network=hydrogen,
        H2_retrofit=False,
        hydrogen_fuel_cell=False,
        hydrogen_turbine=False,
        hydrogen_underground_storage=False,
        methanation=False,
        coal_cc=False,
        SMR=False,
        SMR_cc=False,
    )
    n = pypsa.Network()
    nodes = pd.Index(["DKE1", "DE00", "SE04"])
    n.add("Bus", nodes, carrier="AC", location=nodes, x=[12, 13, 14], y=[55, 54, 56])
    costs = pd.DataFrame(
        {"capital_cost": 1.0, "efficiency": 0.7, "lifetime": 30.0, "fuel": 20.0},
        index=[
            "electrolysis",
            "hydrogen storage tank type 1 including compressor",
            "gas storage",
            "gas",
            "CH4 (g) pipeline",
            "H2 (g) pipeline",
        ],
    )
    pipes = pd.DataFrame(
        {
            "bus0": ["DKE1"],
            "bus1": ["DE00"],
            "p_nom": [5000.0],
            "p_min_pu": [-1.0],
            "length": [200.0],
            "name": ["existing gas pipeline"],
        },
        index=["gas pipeline DK-DE"],
    )
    pipes_path, caverns_path = tmp_path / "pipes.csv", tmp_path / "caverns.csv"
    pipes.to_csv(pipes_path)
    pd.DataFrame(columns=["salt"]).to_csv(caverns_path)
    inputs = pd.DataFrame(
        {"lng": 0.0, "pipeline": 0.0, "production": 1000.0, "storage": 1.0},
        index=nodes,
    )
    spatial = SimpleNamespace(
        nodes=nodes,
        gas=SimpleNamespace(nodes=nodes + " gas", locations=nodes),
    )
    add_h2_gas_infrastructure(
        n,
        costs,
        pd.DataFrame(index=nodes),
        caverns_path,
        ["salt"],
        pipes_path,
        inputs,
        spatial,
        options,
    )

    existing = n.links.loc["gas pipeline DK-DE"]
    assert existing.p_nom == 5000
    assert not existing.p_nom_extendable
    assert existing.bus0 == "DKE1 gas" and existing.bus1 == "DE00 gas"
    assert existing.p_min_pu == -1
    new_gas = n.links[n.links.carrier.eq("gas pipeline new")]
    assert len(new_gas) == connectivity
    assert new_gas.p_nom_extendable.all()
    assert n.links.carrier.eq("H2 pipeline").any() == hydrogen
    assert not n.links.carrier.eq("H2 pipeline retrofitted").any()
    # Disabling pipelines must not remove regional hydrogen or conversion/storage.
    assert n.buses[n.buses.carrier.eq("gas")].location.tolist() == list(nodes)
    assert n.buses[n.buses.carrier.eq("H2")].location.tolist() == list(nodes)
    assert n.links.loc[n.links.carrier.eq("H2 Electrolysis"), "p_nom_extendable"].all()
    assert n.stores.loc[n.stores.carrier.eq("H2 Store"), "e_nom_extendable"].all()
    assert not n.generators.loc[
        n.generators.carrier.eq("gas"), "p_nom_extendable"
    ].any()
