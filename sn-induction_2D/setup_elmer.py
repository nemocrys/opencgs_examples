import numpy as np
import os

from objectgmsh import (
    Model,
    Shape,
    MeshControlLinear,
    MeshControlExponential,
)
from copy import deepcopy
import opencgs.control as ctrl
from opencgs.setup import ElmerSetupCz
import opencgs.geo as geo
from pyelmer import elmer

THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def geometry(config, sim_dir="./simdata/_test", name="cz_induction", visualize=False, mat_config_file="config_mat.yml"):
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
    melt.mesh_size /= 4
    crystal = geo.crystal(
        model, 2, **config["crystal"], melt=melt, phase_if=phase_if_crystal
    )
    inductor = geo.inductor(model, 2, **config["inductor"])
    seed = geo.seed(model, 2, **config["seed"], crystal=crystal)
    ins = geo.crucible_support(
        model, 2, **config["insulation"], top_shape=crucible, name="insulation"
    )
    adp = geo.crucible_adapter(model, 2, **config["crucible_adapter"], top_shape=ins)
    ax_bt = geo.crucible_support(
        model, 2, **config["axis_bt"], top_shape=adp, name="axis_bt"
    )
    vessel = geo.vessel(model, 2, **config["vessel"], adjacent_shapes=[ax_bt])
    ax_top = geo.axis_top(model, 2, **config["axis_top"], bot_shape=seed, vessel=vessel)
    filling = geo.filling(model, 2, **config["filling"], vessel=vessel)
    # filling_in_inductor = Shape(
    #     model,
    #     2,
    #     "filling_in_inductor",
    #     filling.get_part_in_box(
    #         [
    #             inductor.params.X0[0] - inductor.params.d_in / 2,
    #             inductor.params.X0[0] + inductor.params.d_in / 2,
    #         ],
    #         [inductor.params.X0[1] - inductor.params.d_in, 1e6],
    #     ),
    # )
    # filling -= filling_in_inductor
    model.synchronize()

    # boundaries
    bnd_melt = Shape(model, 1, "bnd_melt", melt.get_interface(filling))
    bnd_seed = Shape(model, 1, "bnd_seed", seed.get_interface(filling))
    # split up the boundaries of crystal, seed, ax_top for movement
    bnd_crystal_side = Shape(
        model, 1, "bnd_crystal_side", crystal.get_interface(filling)
    )
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
    bnd_crystal_side -= bnd_crystal_top
    bnd_axtop_side = Shape(model, 1, "bnd_axtop_side", ax_top.get_interface(filling))
    bnd_axtop_bt = Shape(
        model,
        1,
        "bnd_axtop_bt",
        ax_top.get_boundaries_in_box(
            [seed.params.r, ax_top.params.r], [ax_top.params.X0[1], ax_top.params.X0[1]]
        ),
    )
    bnd_axtop_side -= bnd_axtop_bt
    bnd_crucible_bt = Shape(
        model,
        1,
        "bnd_crucible_bt",
        crucible.get_boundaries_in_box(
            [0, ins.params.r_in], [crucible.params.X0[1], crucible.params.X0[1]]
        ),
    )
    bnd_crucible_outside = Shape(
        model, 1, "bnd_crucible_outside", crucible.get_interface(filling)
    )
    bnd_crucible_outside -= bnd_crucible_bt

    bnd_ins = Shape(model, 1, "bnd_ins", ins.get_interface(filling))
    bnd_adp = Shape(model, 1, "bnd_adp", adp.get_interface(filling))
    bnd_axbt = Shape(model, 1, "bnd_axbt", ax_bt.get_interface(filling))
    bnd_vessel_inside = Shape(
        model,
        1,
        "bnd_vessel_inside",
        [
            vessel.get_boundaries_in_box(
                [ax_bt.params.r_out, vessel.params.r_in],
                [ax_bt.params.X0[1], ax_bt.params.X0[1]],
                one_only=True,
            ),  # bottom
            vessel.get_boundaries_in_box(
                [ax_top.params.r, vessel.params.r_in],
                [
                    ax_bt.params.X0[1] + vessel.params.h_in,
                    ax_bt.params.X0[1] + vessel.params.h_in,
                ],
                one_only=True,
            ),  # top
            vessel.get_boundaries_in_box(
                [vessel.params.r_in, vessel.params.r_in],
                [ax_bt.params.X0[1], ax_bt.params.X0[1] + vessel.params.h_in],
                one_only=True,
            ),  # wall
        ],
    )
    bnd_vessel_outside = Shape(
        model,
        1,
        "bnd_vessel_outside",
        [vessel.bottom_boundary, vessel.top_boundary, vessel.right_boundary],
    )
    bnd_inductor = Shape(
        model, 1, "bnd_inductor", inductor.get_interface(filling)
    )
    bnd_symmetry_axis = Shape(model, 1, "bnd_symmetry_axis", model.symmetry_axis)
    
    # interfaces
    if_crucible_melt = Shape(model, 1, "if_crucible_melt", crucible.get_interface(melt))
    if_melt_crystal = Shape(model, 1, "if_melt_crystal", melt.get_interface(crystal))
    if_crystal_seed = Shape(model, 1, "if_crystal_seed", crystal.get_interface(seed))
    if_seed_axtop = Shape(model, 1, "if_seed_axtop", seed.get_interface(ax_top))
    if_axtop_vessel = Shape(model, 1, "if_axtop_vessel", ax_top.get_interface(vessel))
    if_crucible_ins = Shape(model, 1, "if_crucible_ins", crucible.get_interface(ins))
    if_ins_adp = Shape(model, 1, "if_ins_adp", ins.get_interface(adp))
    if_adp_axbt = Shape(model, 1, "if_adp_axbt", adp.get_interface(ax_bt))
    if_axbt_vessel = Shape(model, 1, "if_axbt_vessel", ax_bt.get_interface(vessel))

    model.make_physical()

    # mesh
    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()
    for shape in [melt, crystal, seed, ax_top, crucible, ins, adp, ax_bt, vessel]:
        MeshControlLinear(model, shape, shape.mesh_size, filling.mesh_size)
    MeshControlExponential(
        model, if_melt_crystal, crystal.params.r / 30, exp=1.6, fact=3
    )
    MeshControlExponential(model, bnd_melt, melt.mesh_size / 5, exp=1.8, fact=3)
    MeshControlExponential(model, if_crucible_melt, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(
        model, bnd_crucible_outside, crucible.mesh_size / 3, exp=1.6, fact=3
    )
    MeshControlExponential(model, inductor, inductor.mesh_size)
    model.generate_mesh(**config["mesh"])

    if visualize:
        model.show()
    model.write_msh(f"{sim_dir}/case.msh")
    print(model)
    model.close_gmsh()
    return model


def simulation(
    model,
    config,
    sim_dir="./simdata/_test",
    mat_config={},
    with_flow=False,
    u_scaling=1,
    flow_variable="UMean_avg300s",
    flow_variable_old="UMean_avg300s_old",
    flow_relaxation_factor=1,
    heat_conductivity_scaling=1,
):
    # simulation
    sim = ElmerSetupCz(
        **config["general"],
        sim_dir=sim_dir,
        probes=config["probes"],
        heating=config["heating_induction"],
        smart_heater=config["smart-heater"],
        materials_dict=mat_config
    )
    if "solver-update" in config:
        sim.solver_update=config["solver-update"]
    sim.sim.solvers["ResultOutputSolver"].data.update(
        {
            "Vtu Part collection": True,
            "Scalar Field 1": "temperature",
            "Scalar Field 2": "potential re",
            "Scalar Field 3": "potential im",
            "Scalar Field 4": "heat conductivity",
            "Scalar Field 5": "temperature loads",
            "Scalar Field 6": "current density re e",
            "Scalar Field 7": "current density im e",
            "Scalar Field 8": "joule heating e",
            "Vector Field 1": "temperature grad",
            "Vector Field 2": "temperature flux",
            "Vector Field 3": "electric field re e",
            "Vector Field 4": "electric field im e",
            "Vector Field 5": "magnetic flux density re e",
            "Vector Field 6": "magnetic flux density im e",
            "Vector Field 7": "jxb re e",
            "Vector Field 8": "jxb im e",
        }
    )
    if with_flow:
        sim.sim.solvers["ResultOutputSolver"].data.update(
            {"Vector Field 9": flow_variable}
        )

    eqn_filling = elmer.Equation(
        sim.sim,
        "eqn_filling",
        [x for i, x in enumerate(sim._eqn_main.solvers) if i != 2],
    )  # all solvers but HeatSolver (which is at position 2)

    if with_flow:
        solver_gmsh_input = elmer.Solver(sim.sim, "gmsh_input", config["GmshInput"])
        sim.solver_heat.data.update(
            {"Exported Variable 1": flow_variable, "Exported Variable 1 Dofs": 3}
        )
        if flow_relaxation_factor < 1:
            sim.solver_heat.data.update(
            {"Exported Variable 2": flow_variable_old, "Exported Variable 2 Dofs": 3}
        )

    # load additional solver to save melt interfaces
    solver_flux = elmer.Solver(sim.sim, "solver_flux")
    solver_flux.data = config["FluxSolver"]
    eqn_melt = elmer.Equation(
        sim.sim, "eqn_melt", sim._eqn_main.solvers + [solver_flux]
    )
    if with_flow:
        eqn_melt.data.update({"Convection": "Constant"})
    # melt boundary export
    save_line_crc_melt = elmer.Solver(sim.sim, "save_line_crc_melt")
    save_line_crc_melt.data = deepcopy(config["save_line_solver"])
    save_line_crc_melt.data.update(
        {
            "Equation": '"SaveLineCrcMelt"',
            "Save Mask": 'String "Save Line Crc Melt"',
            "Filename": '"save_line_crc_melt.dat"',
        }
    )
    save_line_melt_surf = elmer.Solver(sim.sim, "save_line_melt_surf")
    save_line_melt_surf.data = deepcopy(config["save_line_solver"])
    save_line_melt_surf.data.update(
        {
            "Equation": '"SaveLineMeltSurf"',
            "Save Mask": 'String "Save Line Melt Surf"',
            "Filename": '"save_line_melt_surf.dat"',
        }
    )
    save_line_melt_crys = elmer.Solver(sim.sim, "save_line_melt_crys")
    save_line_melt_crys.data = deepcopy(config["save_line_solver"])
    save_line_melt_crys.data.update(
        {
            "Equation": '"SaveLineMeltCrys"',
            "Save Mask": 'String "Save Line Melt Crys"',
            "Filename": '"save_line_melt_crys.dat"',
        }
    )
    # functions for T-dependent matrices
    for function in mat_config["functions"]:
        sim.sim.intro_text += function

    # forces
    joule_heat = sim.joule_heat
    # bodies
    sim.add_inductor(model["inductor"])
    sim.add_crystal(model["crystal"], force=joule_heat)
    melt = sim.add_body(model["melt"], force=joule_heat)
    melt.data.update({"name": "melt"})
    if with_flow:
        if flow_relaxation_factor == 1:
            melt.material.data.update(
                {
                    "Convection velocity 1": f'Variable {flow_variable} 1\n    real lua "{u_scaling}*tx[0]"',
                    "Convection velocity 2": f'Variable {flow_variable} 2\n    real lua "{u_scaling}*tx[0]"',
                }
            )
        else:
            melt.material.data.update(
                {
                    "Convection velocity 1": f'Variable {flow_variable} 1, {flow_variable_old} 1\n    real lua "{u_scaling}*({flow_relaxation_factor}*tx[0] + {1 - flow_relaxation_factor}*tx[1])"',
                    "Convection velocity 2": f'Variable {flow_variable} 2, {flow_variable_old} 2\n    real lua "{u_scaling}*({flow_relaxation_factor}*tx[0] + {1 - flow_relaxation_factor}*tx[1])"',
                }
            )
    melt.equation = eqn_melt
    sim.add_body(model["crucible"], force=joule_heat)
    sim.add_body(model["insulation"], force=joule_heat)
    sim.add_body(model["crucible_adapter"], force=joule_heat)
    sim.add_body(model["axis_bt"], force=joule_heat)
    sim.add_body(model["vessel"], force=joule_heat)
    sim.add_body(model["seed"], force=joule_heat)
    sim.add_body(model["axis_top"], force=joule_heat)
    filling = sim.add_body(model["filling"], force=joule_heat)
    filling.equation = eqn_filling
    filling.initial_condition.data.pop(
        "Temperature"
    )  # it's not required and makes the simulation diverging

    # phase interface
    sim.add_phase_interface(model["if_melt_crystal"])

    # boundaries with convection (+ movement)
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
    # moving boundaries
    sim.add_radiation_boundary(model["bnd_seed"])
    sim.add_radiation_boundary(model["bnd_axtop_bt"])
    sim.add_radiation_boundary(model["bnd_axtop_side"])
    # moving interfaces
    sim.add_interface(model["if_crystal_seed"])
    sim.add_interface(model["if_seed_axtop"])
    # stationary boundaries
    for bnd in [
        "bnd_crucible_bt",
        "bnd_ins",
        "bnd_adp",
        "bnd_axbt",
        "bnd_vessel_inside",
        "bnd_inductor",
    ]:
        sim.add_radiation_boundary(model[bnd])
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
    # outside boundaries
    sim.add_temperature_boundary(
        model["bnd_inductor"], **config["boundaries"]["inductor_inside"]
    )
    bc = sim.add_temperature_boundary(
        model["bnd_vessel_outside"], **config["boundaries"]["vessel_outside"]
    )
    bc.zero_potential = True

    # symmetry axis
    sim.add_interface(model["bnd_symmetry_axis"], movement=[0, None])

    # set flags for boundary data export for coupling
    sim.sim.boundaries["if_crucible_melt"].data.update(
        {"Save Line Crc Melt": "Logical True"}
    )
    sim.sim.boundaries["bnd_melt"].data.update({"Save Line Melt Surf": "Logical True"})
    sim.sim.boundaries["if_melt_crystal"].data.update(
        {"Save Line Melt Crys": "Logical True"}
    )

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

    # sim.heat_flux_computation(sim["insulation"], sim["if_crucible_ins"])
    # sim.heat_flux_computation(sim["insulation"], sim["bnd_ins"])
    # sim.heat_flux_computation(sim["insulation"], sim["if_ins_adp"])

    # sim.heat_flux_computation(sim["crucible_adapter"], sim["if_ins_adp"])
    # sim.heat_flux_computation(sim["crucible_adapter"], sim["bnd_adp"])
    # sim.heat_flux_computation(sim["crucible_adapter"], sim["if_adp_axbt"])

    # sim.heat_flux_computation(sim["axis_bt"], sim["if_adp_axbt"])
    # sim.heat_flux_computation(sim["axis_bt"], sim["bnd_axbt"])
    # sim.heat_flux_computation(sim["axis_bt"], sim["if_axbt_vessel"])

    # sim.heat_flux_computation(sim["vessel"], sim["if_axbt_vessel"])
    # sim.heat_flux_computation(sim["vessel"], sim["if_axtop_vessel"])
    # sim.heat_flux_computation(sim["vessel"], sim["bnd_vessel_inside"])
    # sim.heat_flux_computation(sim["vessel"], sim["bnd_vessel_outside"])

    # export
    sim.export()

    return sim


if __name__ == "__main__":
    geo_config = ctrl.load_config("./config_geo.yml")
    sim_config = ctrl.load_config("./config_sim.yml")
    model1 = geometry(geo_config, visualize=True)
