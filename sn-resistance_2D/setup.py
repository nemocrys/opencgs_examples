import numpy as np
import os
import sys
import yaml

import pyelmer.elmer as elmer
from objectgmsh import (
    Model,
    Shape,
    MeshControlConstant,
    MeshControlLinear,
    MeshControlExponential,
    GeometryError
)

from opencgs.setup import ElmerSetupCz
import opencgs.control as ctrl
import opencgs.geo as geo


THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def geometry(config, sim_dir="./simdata/_test", name="cz_resistance", visualize=False, mat_config_file="config_mat.yml"):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    config_mat = ctrl.load_config(mat_config_file)
    matdata_melt = config_mat[config["melt"]["material"]]

    # initialize geometry model
    model = Model(name)

    # draw phase interfac if predefined (for transient simulation)
    if "phase_if" in config:
        phase_if_crystal = geo.line_from_points(
            model, config["phase_if"], "phase_if_crystal"
        )
        phase_if_melt = geo.line_from_points(model, config["phase_if"], "phase_if_melt")
    else:
        phase_if_crystal = None
        phase_if_melt = None
    # geometry
    crucible = geo.crucible(model, 2, **config["crucible"])
    melt = geo.melt(
        model,
        2,
        crucible,
        **config["melt"],
        crystal_radius=config["crystal"]["r"],
        phase_if=phase_if_melt,
        beta=matdata_melt["Beta"],
        gamma=matdata_melt["Surface Tension"],
        rho=matdata_melt["Density"],
    )
    crystal = geo.crystal(
        model, 2, **config["crystal"], melt=melt, phase_if=phase_if_crystal
    )
    seed = geo.seed(model, 2, **config["seed"], crystal=crystal)
    ins_crc = geo.crucible_support(
        model, 2, **config["crucible_insulation"], top_shape=crucible, name="crucible_insulation"
    )
    adp = geo.crucible_adapter(model, 2, **config["crucible_adapter"], top_shape=ins_crc)
    ax_bt = geo.crucible_support(
        model, 2, **config["axis_bt"], top_shape=adp, name="axis_bt"
    )
    vessel = geo.vessel(model, 2, **config["vessel"], adjacent_shapes=[ax_bt])
    ax_top = geo.axis_top(model, 2, **config["axis_top"], bot_shape=seed, vessel=vessel)
    heater = geo.heater(model, 2, **config["heater"], vessel=vessel)
    if config["insulation"]["bottom_t"] != 0:
        insulaton = geo.insulation(model, 2, **config["insulation"], vessel=vessel)

    model.synchronize()

    # interfaces
    if_crucible_melt = Shape(model, 1, "if_crucible_melt", crucible.get_interface(melt))
    if_melt_crystal = Shape(model, 1, "if_melt_crystal", melt.get_interface(crystal))
    if_crystal_seed = Shape(model, 1, "if_crystal_seed", crystal.get_interface(seed))
    if_seed_axtop = Shape(model, 1, "if_seed_axtop", seed.get_interface(ax_top))
    if_axtop_vessel = Shape(model, 1, "if_axtop_vessel", ax_top.get_interface(vessel))
    if_crucible_ins = Shape(model, 1, "if_crucible_ins", crucible.get_interface(ins_crc))
    if_ins_adp = Shape(model, 1, "if_ins_adp", ins_crc.get_interface(adp))
    if_adp_axbt = Shape(model, 1, "if_adp_axbt", adp.get_interface(ax_bt))
    if_axbt_vessel = Shape(model, 1, "if_axbt_vessel", ax_bt.get_interface(vessel))
    if config["insulation"]["bottom_t"] != 0:
        if_insulation_vessel = Shape(model, 1, "if_insulation_vessel", insulaton.get_interface(vessel))

    # boundaries
    bnd_symmetry_axis = Shape(model, 1, "bnd_symmetry_axis", model.symmetry_axis)

    bnd_melt = Shape(model, 1, "bnd_melt", melt.boundaries)
    bnd_melt -= if_melt_crystal
    bnd_melt -= if_crucible_melt
    bnd_melt -= bnd_symmetry_axis
    bnd_seed = Shape(model, 1, "bnd_seed", seed.boundaries)
    bnd_seed -= if_crystal_seed
    bnd_seed -= if_seed_axtop
    bnd_seed -= bnd_symmetry_axis
    # split up the boundaries of crystal, seed, ax_top for movement
    bnd_crystal_top = Shape(
        model,
        1,
        "bnd_crystal_top",
        crystal.get_boundaries_in_box(
            [seed.params.r, crystal.params.r],
            [
                crystal.params.X0[1] + crystal.params.l,
                crystal.params.X0[1] + crystal.params.l,
            ],
        ),
    )
    bnd_crystal_side = Shape(model, 1, "bnd_crystal_side", crystal.boundaries)
    bnd_crystal_side -= bnd_crystal_top
    bnd_crystal_side -= if_melt_crystal
    bnd_crystal_side -= if_crystal_seed
    bnd_crystal_side -= bnd_symmetry_axis
    bnd_axtop_bt = Shape(
        model,
        1,
        "bnd_axtop_bt",
        ax_top.get_boundaries_in_box(
            [seed.params.r, ax_top.params.r], [ax_top.params.X0[1], ax_top.params.X0[1]]
        ),
    )
    bnd_axtop_side = Shape(model, 1, "bnd_axtop_side", ax_top.boundaries)
    bnd_axtop_side -= if_seed_axtop
    bnd_axtop_side -= if_axtop_vessel
    bnd_axtop_side -= bnd_axtop_bt
    bnd_axtop_side -= bnd_symmetry_axis
    bnd_crucible_bt = Shape(
        model,
        1,
        "bnd_crucible_bt",
        crucible.get_boundaries_in_box(
            [0, ins_crc.params.r_in], [crucible.params.X0[1], crucible.params.X0[1]]
        ),
    )
    bnd_crucible_outside = Shape(model, 1, "bnd_crucible_outside", crucible.boundaries    )
    bnd_crucible_outside -= bnd_crucible_bt
    bnd_crucible_outside -= if_crucible_melt
    bnd_crucible_outside -= if_crucible_ins
    bnd_crucible_outside -= bnd_symmetry_axis

    bnd_ins_crc = Shape(model, 1, "bnd_ins_crc", ins_crc.boundaries)
    bnd_ins_crc -= if_crucible_ins
    bnd_ins_crc -= if_ins_adp
    bnd_ins_crc -= bnd_symmetry_axis
    bnd_adp = Shape(model, 1, "bnd_adp", adp.boundaries)
    bnd_adp -= if_ins_adp
    bnd_adp -= if_adp_axbt
    bnd_adp -= bnd_symmetry_axis
    bnd_axbt = Shape(model, 1, "bnd_axbt", ax_bt.boundaries)
    bnd_axbt -= if_adp_axbt
    bnd_axbt -= if_axbt_vessel
    bnd_axbt -= bnd_symmetry_axis
    if config["insulation"]["bottom_t"] != 0:
        bnd_insulation = Shape(model, 1, "bnd_insulation", insulaton.boundaries)
        bnd_insulation -= if_insulation_vessel
    bnd_heater = Shape(model, 1, "bnd_heater", heater.boundaries)
    bnd_vessel_outside = Shape(
        model,
        1,
        "bnd_vessel_outside",
        [vessel.bottom_boundary, vessel.top_boundary, vessel.right_boundary],
    )
    bnd_vessel_inside = Shape(model, 1, "bnd_vessel_inside", vessel.boundaries)
    bnd_vessel_inside -= bnd_vessel_outside
    bnd_vessel_inside -= if_axtop_vessel
    bnd_vessel_inside -= if_axbt_vessel
    if config["insulation"]["bottom_t"] != 0:
        bnd_vessel_inside -= if_insulation_vessel
    bnd_vessel_inside -= bnd_symmetry_axis

    model.make_physical()

    # mesh
    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()
    for shape in [melt, crystal, seed, ax_top, crucible, ins_crc, adp, ax_bt, vessel]:
        MeshControlLinear(model, shape, shape.mesh_size, vessel.params.r_in)
    MeshControlExponential(
        model, if_melt_crystal, crystal.params.r / 30, exp=1.6, fact=3
    )
    MeshControlExponential(model, bnd_melt, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(model, if_crucible_melt, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(
        model, bnd_crucible_outside, crucible.mesh_size / 3, exp=1.6, fact=3
    )
    MeshControlExponential(model, heater, heater.mesh_size)
    model.generate_mesh(**config["mesh"])

    if visualize:
        model.show()
    model.write_msh(f"{sim_dir}/case.msh")
    print(model)
    model.close_gmsh()
    return model


def simulation(model, config, sim_dir="./simdata/_test", mat_config={}):
    # simulation
    sim = ElmerSetupCz(
        **config["general"],
        sim_dir=sim_dir,
        probes=config["probes"],
        heating=config["heating_resistance"],
        smart_heater=config["smart-heater"],
        materials_dict=mat_config
    )
    if "solver-update" in config:
        sim.solver_update=config["solver-update"]
    
    # load additional solver to compute heat flux
    solver_flux = elmer.Solver(sim.sim, "solver_flux")
    solver_flux.data = config["FluxSolver"]
    sim._eqn_main.solvers.append(solver_flux)
    
    # functions for T-dependent matrices
    for function in mat_config["functions"]:
        sim.sim.intro_text += function

    # bodies
    sim.add_resistance_heater(model["heater"])
    sim.add_crystal(model["crystal"])
    sim.add_body(model["melt"])
    sim.add_body(model["crucible"])
    sim.add_body(model["crucible_insulation"])
    sim.add_body(model["crucible_adapter"])
    try:
        sim.add_body(model["insulation"])
    except GeometryError:
        pass
    sim.add_body(model["axis_bt"])
    sim.add_body(model["vessel"])
    sim.add_body(model["seed"])
    sim.add_body(model["axis_top"])

    # phase interface
    sim.add_phase_interface(model["if_melt_crystal"])

    # boundaries with convection
    sim.add_radiation_boundary(
        model["bnd_crucible_outside"], **config["boundaries"]["crucible_outside"]
    )
    sim.add_radiation_boundary(model["bnd_melt"], **config["boundaries"]["melt"])
    sim.add_radiation_boundary(
        model["bnd_crystal_side"],
        **config["boundaries"]["crystal"],
    )
    sim.add_radiation_boundary(
        model["bnd_crystal_top"],
        **config["boundaries"]["crystal"],
    )
    sim.add_radiation_boundary(
        model["bnd_heater"],
        **config["boundaries"]["heater"]
    )
    sim.add_radiation_boundary(model["bnd_seed"])
    sim.add_radiation_boundary(model["bnd_axtop_bt"])
    sim.add_radiation_boundary(model["bnd_axtop_side"])
    sim.add_interface(model["if_crystal_seed"])
    sim.add_interface(model["if_seed_axtop"])
    # stationary boundaries
    for bnd in [
        "bnd_crucible_bt",
        "bnd_ins_crc",
        "bnd_adp",
        "bnd_axbt",
        "bnd_vessel_inside",
    ]:
        sim.add_radiation_boundary(model[bnd])
    try:
        sim.add_radiation_boundary(model["bnd_insulation"])
    except GeometryError:
        pass
    # stationary interfaces
    for bnd in [
        "if_crucible_melt",
        "if_axtop_vessel",
        "if_crucible_ins",
        "if_ins_adp",
        "if_adp_axbt",
        "if_axbt_vessel",
    ]:
        sim.add_interface(model[bnd])
    try:
        sim.add_interface(model["if_insulation_vessel"])
    except GeometryError:
        pass
    # outside boundaries
    sim.add_temperature_boundary(
        model["bnd_vessel_outside"], **config["boundaries"]["vessel_outside"]
    )

    # symmetry axis
    sim.add_interface(model["bnd_symmetry_axis"], [0, None])

    # heat flux computation
    # doesn't work with 2D T-dependent heat conductivity of insulation
    # sim.heat_flux_computation(sim["crucible"], sim["bnd_crucible_outside"])
    # sim.heat_flux_computation(sim["crucible"], sim["bnd_crucible_bt"])
    # sim.heat_flux_computation(sim["crucible"], sim["if_crucible_melt"])
    # sim.heat_flux_computation(sim["crucible"], sim["if_crucible_ins"])

    # sim.heat_flux_computation(sim["melt"], sim["if_crucible_melt"])
    # sim.heat_flux_computation(sim["melt"], sim["if_melt_crystal"])
    # sim.heat_flux_computation(sim["melt"], sim["bnd_melt"])

    # sim.heat_flux_computation(sim["crystal"], sim["if_melt_crystal"])
    # sim.heat_flux_computation(sim["crystal"], sim["bnd_crystal_side"])
    # sim.heat_flux_computation(sim["crystal"], sim["bnd_crystal_top"])
    # sim.heat_flux_computation(sim["crystal"], sim["if_crystal_seed"])

    # sim.heat_flux_computation(sim["seed"], sim["if_crystal_seed"])
    # sim.heat_flux_computation(sim["seed"], sim["bnd_seed"])
    # sim.heat_flux_computation(sim["seed"], sim["if_seed_axtop"])

    # sim.heat_flux_computation(sim["axis_top"], sim["if_seed_axtop"])
    # sim.heat_flux_computation(sim["axis_top"], sim["bnd_axtop_bt"])
    # sim.heat_flux_computation(sim["axis_top"], sim["bnd_axtop_side"])
    # sim.heat_flux_computation(sim["axis_top"], sim["if_axtop_vessel"])

    # sim.heat_flux_computation(sim["crucible_insulation"], sim["if_crucible_ins"])
    # sim.heat_flux_computation(sim["crucible_insulation"], sim["bnd_ins_crc"])
    # sim.heat_flux_computation(sim["crucible_insulation"], sim["if_ins_adp"])

    # sim.heat_flux_computation(sim["insulation"], sim["bnd_insulation"])
    # sim.heat_flux_computation(sim["insulation"], sim["if_insulation_vessel"])
    
    # sim.heat_flux_computation(sim["heater"], sim["bnd_heater"])

    # sim.heat_flux_computation(sim["crucible_adapter"], sim["if_ins_adp"])
    # sim.heat_flux_computation(sim["crucible_adapter"], sim["bnd_adp"])
    # sim.heat_flux_computation(sim["crucible_adapter"], sim["if_adp_axbt"])

    # sim.heat_flux_computation(sim["axis_bt"], sim["if_adp_axbt"])
    # sim.heat_flux_computation(sim["axis_bt"], sim["bnd_axbt"])
    # sim.heat_flux_computation(sim["axis_bt"], sim["if_axbt_vessel"])

    # sim.heat_flux_computation(sim["vessel"], sim["if_axbt_vessel"])
    # sim.heat_flux_computation(sim["vessel"], sim["if_axtop_vessel"])
    # sim.heat_flux_computation(sim["vessel"], sim["if_insulation_vessel"])
    # sim.heat_flux_computation(sim["vessel"], sim["bnd_vessel_inside"])
    # sim.heat_flux_computation(sim["vessel"], sim["bnd_vessel_outside"])

    # export
    sim.export()
    
    return sim


if __name__ == "__main__":
    geo_config = ctrl.load_config("./config_geo.yml")
    sim_config = ctrl.load_config("./config_sim.yml")
    model = geometry(geo_config, visualize=True)
    # simulation(model, sim_config)
