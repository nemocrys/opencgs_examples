"""
Geometry and simulation setup (in two separate functions).
"""

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
)
import opencgs.control as ctrl
import opencgs.geo as geo
from opencgs.setup import ElmerSetupCz
from copy import deepcopy


THIS_DIR = os.path.dirname(os.path.realpath(__file__))


def geometry(
    config,
    sim_dir="./simdata/_test",
    name="cz_induction",
    visualize=False,
    with_afterheater=False,
    mat_config_file="config_mat.yml",
):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    config_mat = ctrl.load_config(mat_config_file)

    # initialize geometry model
    model = Model(name)

    # draw phase interface if predefined (for transient simulation)
    if "phase_if" in config:
        phase_if_crystal = geo.line_from_points(
            model, config["phase_if"], "phase_if_crystal"
        )
        phase_if_melt = geo.line_from_points(model, config["phase_if"], "phase_if_melt")
    else:
        phase_if_crystal = None
        phase_if_melt = None

    # compute melt height
    matdata_melt = config_mat[config["melt"]["material"]]
    matdata_crystal = config_mat[config["crystal"]["material"]]
    density_crystal = matdata_crystal["Density"]
    density_melt = matdata_melt["Density"]
    crystal_surface, crystal_volume, crystal_angle = geo.compute_current_crystal_surface(
        **config["crystal"]
    )
    crystal_radius = crystal_surface[-1, 0]
    meniscus_x, meniscus_y, _, _, _ = geo.menicus_hurle(
        crystal_radius,
        crystal_angle,
        matdata_melt["Beta"],
        matdata_melt["Surface Tension"],
        matdata_melt["Density"],
        res=100,
    )
    meniscus_volume = geo.compute_meniscus_volume(meniscus_x, meniscus_y)
    if "total_mass" in config:  # melt height as result from mass balance
        melt_mass = (
            config["total_mass"]
            - meniscus_volume * density_melt
            - crystal_volume * density_crystal
        )
        config["melt"]["h"] = float(
            melt_mass / (density_melt * np.pi * config["crucible"]["r_in"] ** 2)
        )

    # geometry
    crucible = geo.crucible(model, 2, **config["crucible"])
    melt = geo.melt(
        model,
        2,
        crucible,
        **config["melt"],
        crystal_radius=crystal_radius,
        alpha_crys=crystal_angle,
        phase_if=phase_if_melt,
        beta=matdata_melt["Beta"],
        gamma=matdata_melt["Surface Tension"],
        rho=matdata_melt["Density"],
    )
    crystal = geo.crystal_from_points(
        model, 2, **config["crystal"], melt=melt, phase_if=phase_if_crystal
    )

    inductor = geo.inductor(model, 2, **config["inductor"])
    # inductor = geo.inductor(model, 2, **config["inductor_bot"], name="inductor")
    # inductor2 = geo.inductor(model, 2, **config["inductor_top"], name="inductor2")
    # inductor.geo_ids += inductor2.geo_ids
    # inductor.params.area += inductor2.params.area
    # model._shapes.remove(inductor2)

    seedholder = geo.cylinder_shape(
        model, 2, **config["seedholder"], bot_shape=crystal, name="seedholder"
    )
    shaft_al2o3 = geo.cylinder_shape(
        model, 2, **config["shaft_al2o3"], bot_shape=seedholder, name="shaft_al2o3"
    )
    axtop_adapter = geo.cylinder_shape(
        model, 2, **config["axtop_adapter"], bot_shape=shaft_al2o3, name="axtop_adapter"
    )
    if with_afterheater:
        afterheater_bot = geo.cylinder_shape(
            model,
            2,
            bot_shape=crucible,
            name="afterheater_bot",
            **config["afterheater_bot"],
        )
        afterheater_bot.set_interface(crucible)
        afterheater_top = geo.cylinder_shape(
            model,
            2,
            bot_shape=afterheater_bot,
            name="afterheater_top",
            **config["afterheater_top"],
        )
        afterheater_cover = geo.cylinder_shape(
            model,
            2,
            bot_shape=afterheater_top,
            name="afterheater_cover",
            **config["afterheater_cover"],
        )
    ins_bot = geo.stacked_cylinder_shape(
        model,
        2,
        **config["insulation_bot"],
        top_shape=crucible,
        name="insulation_bot",
    )
    if with_afterheater:
        ins_side = geo.cylinder_shape(
            model,
            2,
            bot_shape=ins_bot,
            name="insulation_side",
            **config["insulation_side"],
        )
    else:
        ins_side = geo.cylinder_shape(
            model,
            2,
            bot_shape=ins_bot,
            name="insulation_side",
            **config["insulation_side_two_rings"],
        )
    ins_side.set_interface(ins_bot)
    if with_afterheater:
        ins_top = geo.cylinder_shape(
            model,
            2,
            bot_shape=ins_side,
            name="insulation_top",
            **config["insulation_top"],
        )
        ins_top.set_interface(ins_side)
    adp = geo.stacked_cylinder_shape(
        model, 2, **config["axbot_adapter"], top_shape=ins_bot, name="axbot_adapter"
    )
    ax_bt = geo.cylinder_shape(
        model, 2, **config["axis_bt"], top_shape=adp, name="axis_bt"
    )
    vessel = geo.vessel(model, 2, **config["vessel"], adjacent_shapes=[ax_bt])
    ax_top = geo.axis_top(
        model, 2, **config["axis_top"], bot_shape=axtop_adapter, vessel=vessel
    )
    ax_top.set_interfaces([axtop_adapter, shaft_al2o3])
    exhaust_hood = geo.exhaust_hood(
        model, 2, **config["exhaust_hood"], y_bot=vessel.params.X0[1] + vessel.params.t
    )
    filling = geo.filling(model, 2, **config["filling"], vessel=vessel)
    model.synchronize()

    # boundaries
    if with_afterheater:
        bnd_afterheater_bot = Shape(
            model, 1, "bnd_afterheater_bot", afterheater_bot.get_interface(filling)
        )
        bnd_afterheater_top = Shape(
            model, 1, "bnd_afterheater_top", afterheater_top.get_interface(filling)
        )
        bnd_afterheater_cover = Shape(
            model, 1, "bnd_afterheater_cover", afterheater_cover.get_interface(filling)
        )
    bnd_melt = Shape(model, 1, "bnd_melt", melt.get_interface(filling))
    bnd_seedholder = Shape(
        model, 1, "bnd_seedholder", seedholder.get_interface(filling)
    )
    bnd_shaft_al2o3 = Shape(
        model, 1, "bnd_shaft_al2o3", shaft_al2o3.get_interface(filling)
    )
    bnd_axtop_adapter = Shape(
        model, 1, "bnd_axtop_adapter", axtop_adapter.get_interface(filling)
    )
    bnd_crystal = Shape(model, 1, "bnd_crystal", crystal.get_interface(filling))
    bnd_axtop = Shape(model, 1, "bnd_axtop", ax_top.get_interface(filling))
    bnd_crucible_outside = Shape(
        model, 1, "bnd_crucible_outside", crucible.get_interface(filling)
    )
    bnd_ins_bot = Shape(model, 1, "bnd_ins_bot", ins_bot.get_interface(filling))
    bnd_ins_side = Shape(model, 1, "bnd_ins_side", ins_side.get_interface(filling))
    if with_afterheater:
        bnd_ins_top = Shape(model, 1, "bnd_ins_top", ins_top.get_interface(filling))
    bnd_adp = Shape(model, 1, "bnd_adp", adp.get_interface(filling))
    bnd_axbt = Shape(model, 1, "bnd_axbt", ax_bt.get_interface(filling))
    bnd_exhaust_hood = Shape(
        model, 1, "bnd_exhaust_hood", exhaust_hood.get_interface(filling)
    )
    bnd_vessel_inside = Shape(
        model,
        1,
        "bnd_vessel_inside",
        [
            vessel.get_boundaries_in_box(
                [0, ax_bt.params.r_in],
                [ax_bt.params.X0[1], ax_bt.params.X0[1]],
                one_only=True,
            ),  # bottom inside
            vessel.get_boundaries_in_box(
                [ax_bt.params.r_out, vessel.params.r_in],
                [ax_bt.params.X0[1], ax_bt.params.X0[1]],
                one_only=True,
            ),  # bottom outside
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
    bnd_inductor = Shape(model, 1, "bnd_inductor", inductor.get_interface(filling))
    bnd_symmetry_axis = Shape(model, 1, "bnd_symmetry_axis", model.symmetry_axis)

    # interfaces
    if_crucible_melt = Shape(model, 1, "if_crucible_melt", crucible.get_interface(melt))
    if_melt_crystal = Shape(model, 1, "if_melt_crystal", melt.get_interface(crystal))
    if_crystal_seedholder = Shape(
        model, 1, "if_crystal_seedholder", crystal.get_interface(seedholder)
    )
    if_axtop_vessel = Shape(model, 1, "if_axtop_vessel", ax_top.get_interface(vessel))
    if_crucible_ins_bot = Shape(
        model, 1, "if_crucible_ins_bot", crucible.get_interface(ins_bot)
    )
    if_ins_bot_adp = Shape(model, 1, "if_ins_bot_adp", ins_bot.get_interface(adp))
    if_ins_bot_side = Shape(
        model, 1, "if_ins_bot_side", ins_bot.get_interface(ins_side)
    )
    if with_afterheater:
        if_ins_side_top = Shape(
            model, 1, "if_ins_side_top", ins_side.get_interface(ins_top)
        )
        if_crucible_afterheater = Shape(
            model,
            1,
            "if_crucible_afterheater",
            crucible.get_interface(afterheater_bot),
        )
    if_adp_axbt = Shape(model, 1, "if_adp_axbt", adp.get_interface(ax_bt))
    if_axbt_vessel = Shape(model, 1, "if_axbt_vessel", ax_bt.get_interface(vessel))

    model.make_physical()

    # mesh
    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()
    shapes = [
        melt,
        crystal,
        seedholder,
        ax_top,
        crucible,
        ins_bot,
        ins_side,
        adp,
        ax_bt,
        exhaust_hood,
        vessel,
    ]
    if with_afterheater:
        shapes += [afterheater_bot, afterheater_top, afterheater_cover, ins_top]
    for shape in shapes:
        MeshControlLinear(model, shape, shape.mesh_size, filling.mesh_size)
    MeshControlExponential(model, if_melt_crystal, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(model, bnd_melt, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(model, if_crucible_melt, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(
        model, bnd_crucible_outside, melt.mesh_size / 2, exp=1.6, fact=3
    )
    bnds = [bnd_crucible_outside]
    if with_afterheater:
        bnds += [
            bnd_afterheater_bot,
            bnd_afterheater_top,
            bnd_afterheater_cover,
        ]

    for bnd in bnds:
        MeshControlExponential(
            model,
            bnd,
            crucible.mesh_size / 3,
            exp=1.6,
            fact=3,
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
    transparent_crys_melt=True,
    u_scaling=1,
    flow_variable="UMean_avg300s",
    flow_variable_old="UMean_avg300s_old",
    flow_relaxation_factor=1,
    heat_conductivity_scaling=1,
    with_afterheater=False,
):
    # simulation
    sim = ElmerSetupCz(
        **config["general"],
        sim_dir=sim_dir,
        probes=config["probes"],
        heating=config["heating_induction"],
        smart_heater=config["smart-heater"],
        materials_dict=deepcopy(mat_config),
    )
    if "solver-update" in config:
        sim.solver_update = config["solver-update"]
    sim.solver_heat.data.update({"Update Gebhart Factors": True})
    sim.sim.solvers["ResultOutputSolver"].data.update(
        {"Vtu Part collection": True}
    )  # , "Save Bulk Only": True})
    eqn_filling = elmer.Equation(
        sim.sim,
        "eqn_filling",
        [x for i, x in enumerate(sim._eqn_main.solvers) if i != 1],
    )  # all solvers but HeatSolver (which is at position 1)

    if with_flow:
        solver_gmsh_input = elmer.Solver(sim.sim, "gmsh_input", config["GmshInput"])
        sim.solver_heat.data.update(
            {"Exported Variable 1": flow_variable, "Exported Variable 1 Dofs": 3}
        )
        if flow_relaxation_factor < 1:
            sim.solver_heat.data.update(
            {"Exported Variable 2": flow_variable_old, "Exported Variable 2 Dofs": 3}
        )

    # set correct growth velocity
    area_melt = np.pi * model["crucible"].params.r_in**2
    area_crystal = np.pi * model["crystal"].params.r**2
    dty_crystal = mat_config[model["crystal"].params.material]["Density"]
    dty_melt = mat_config[model["melt"].params.material]["Density"]
    print(f"Pulling velocity: {sim.v_pull:.3f} mm/min")
    print(f"Area ratio: {area_crystal / area_melt:.4f}")
    print(f"Density ratio: {dty_crystal / dty_melt}")
    sim.v_pull += sim.v_pull*  area_crystal / area_melt * dty_crystal / dty_melt
    print(f"Growth velocity: {sim.v_pull:.4f} mm/min")

    # load additional solver to save melt interfaces
    solver_flux = elmer.Solver(sim.sim, "solver_flux")
    solver_flux.data = config["FluxSolver"]
    eqn_melt = elmer.Equation(
        sim.sim, "eqn_melt", sim._eqn_main.solvers + [solver_flux]
    )
    if with_flow:
        eqn_melt.data.update({"Convection": "Constant"})
    eqn_crys = elmer.Equation(
        sim.sim, "eqn_crys", sim._eqn_main.solvers + [solver_flux]
    )

    # boundary export for separate gas flow simulation
    save_line_gasflow = elmer.Solver(sim.sim, "save_line_gasflow")
    save_line_gasflow.data = deepcopy(config["save_line_solver"])
    save_line_gasflow.data.update(
        {
            "Equation": '"SaveLineGasflow"',
            "Save Mask": 'String "Save Line Gasflow"',
            "Filename": '"save_line_gasflow.dat"',
        }
    )
    save_line_gasflow.data.pop("Variable 2")  # don't export fluxes (they are not computed everywhere)
    save_line_gasflow.data.pop("Variable 3")  # don't export fluxes (they are not computed everywhere)

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
    # crystal boundary export
    save_line_crys_surf = elmer.Solver(sim.sim, "save_line_crys_surf")
    save_line_crys_surf.data = deepcopy(config["save_line_solver"])
    save_line_crys_surf.data.update(
        {
            "Equation": '"SaveLineCrysSurf"',
            "Save Mask": 'String "Save Line Crys Surf"',
            "Filename": '"save_line_crys_surf.dat"',
        }
    )
    # functions for T-dependent matrices
    for function in mat_config["functions"]:
        sim.sim.intro_text += function
    # forces
    joule_heat = sim.joule_heat
    # bodies
    sim.add_inductor(model["inductor"])
    crystal = sim.add_crystal(model["crystal"], force=joule_heat)
    crystal.data.update({"name": "crystal"})
    crystal.equation = eqn_crys
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
    melt.material.data["Heat Conductivity"] *= heat_conductivity_scaling
    for body in [
        "crucible",
        "insulation_bot",
        "insulation_side",
        "axbot_adapter",
        "axis_bt",
        "vessel",
        "seedholder",
        "shaft_al2o3",
        "axtop_adapter",
        "axis_top",
        "exhaust_hood",
    ]:
        bdy = sim.add_body(model[body], force=joule_heat)
        bdy.data.update({"name": body})

    # we don't want to have filling in our thermal simulation
    filling = sim.add_body(model["filling"])
    filling.data.update({"name": "filling"})
    filling.equation = eqn_filling
    filling.initial_condition.data.pop(
        "Temperature"
    )  # it's not required and makes the simulation diverging

    if with_afterheater:
        for body in [
            "insulation_top",
            "afterheater_bot",
            "afterheater_top",
            "afterheater_cover",
        ]:
            bdy = sim.add_body(model[body], force=joule_heat)
            bdy.data.update({"name": body})

    # phase interface
    sim.add_phase_interface(model["if_melt_crystal"])

    # boundaries with convection
    sim.add_radiation_boundary(
        model["bnd_crucible_outside"], **config["boundaries"]["crucible_outside"]
    )
    if transparent_crys_melt:
        bc = sim.add_interface(model["bnd_melt"])
        bc.heat_transfer_coefficient = config["boundaries"]["melt"]["htc"]
        bc.T_ext = config["boundaries"]["melt"]["T_ext"]

        bc = sim.add_interface(model["bnd_crystal"])
        bc.heat_transfer_coefficient = config["boundaries"]["crystal"]["htc"]
        bc.T_ext = config["boundaries"]["crystal"]["T_ext"]

        # raise ValueError("Not maintained, some features may be missing.")
    else:
        sim.add_radiation_boundary(model["bnd_melt"], **config["boundaries"]["melt"])
        sim.add_radiation_boundary(
            model["bnd_crystal"],
            **config["boundaries"]["crystal"],
        )
    if with_afterheater:
        sim.add_radiation_boundary(
            model["bnd_afterheater_bot"], **config["boundaries"]["afterheater"]
        )
        sim.add_radiation_boundary(
            model["bnd_afterheater_top"], **config["boundaries"]["afterheater"]
        )
        sim.add_radiation_boundary(
            model["bnd_afterheater_cover"], **config["boundaries"]["afterheater"]
        )
    # boundaries without convection
    bnds = [
        "bnd_seedholder",
        "bnd_shaft_al2o3",
        "bnd_axtop_adapter",
        "bnd_axtop",
        "bnd_ins_bot",
        "bnd_ins_side",
        "bnd_adp",
        "bnd_axbt",
        "bnd_exhaust_hood",
        "bnd_vessel_inside",
    ]
    if with_afterheater:
        bnds += [
            "bnd_ins_top",
        ]
    for bnd in bnds:
        sim.add_radiation_boundary(model[bnd])
    bc = sim.add_radiation_boundary(model["bnd_inductor"])
    bc.fixed_temperature = config["boundaries"]["inductor"]["T"]
    # interfaces
    for bnd in [
        "if_axtop_vessel",
        "if_crucible_ins_bot",
        "if_ins_bot_adp",
        "if_ins_bot_side",
        "if_adp_axbt",
        "if_axbt_vessel",
    ]:
        sim.add_interface(model[bnd])
    if with_afterheater:
        sim.add_interface(model["if_ins_side_top"])
        sim.add_interface(model["if_crucible_afterheater"])
    if transparent_crys_melt:
        bc = sim.add_radiation_boundary(model["if_crystal_seedholder"])
        bc.data.update({"Radiation Target Body": elmer.StringFromList([crystal])})
        bc = sim.add_radiation_boundary(model["if_crucible_melt"])
        bc.data.update({"Radiation Target Body": elmer.StringFromList([melt])})
    else:
        sim.add_interface(model["if_crystal_seedholder"])
        sim.add_interface(model["if_crucible_melt"])
    # outside boundaries
    sim.add_temperature_boundary(
        model["bnd_vessel_outside"], **config["boundaries"]["vessel_outside"]
    )

    # symmetry axis
    sim.add_interface(model["bnd_symmetry_axis"], movement=[0, None])

    sim.sim.boundaries["if_crucible_melt"].data.update(
        {"Save Line Crc Melt": "Logical True"}
    )
    sim.sim.boundaries["bnd_melt"].data.update({"Save Line Melt Surf": "Logical True"})
    sim.sim.boundaries["if_melt_crystal"].data.update(
        {"Save Line Melt Crys": "Logical True"}
    )
    for bnd in ["if_crystal_seedholder", "bnd_crystal"]:
        sim.sim.boundaries[bnd].data.update({"Save Line Crys Surf": "Logical True"})

    for bnd in sim.sim.boundaries.values():
        bnd.data.update({"Save Line Gasflow": "Logical True"})

    # # heat flux computation
    # sim.heat_flux_computation(sim["crucible"], sim["bnd_crucible_outside"])
    # sim.heat_flux_computation(sim["crucible"], sim["if_crucible_melt"])
    # sim.heat_flux_computation(sim["crucible"], sim["if_crucible_ins_bot"])
    # if with_afterheater:
    #     sim.heat_flux_computation(sim["crucible"], sim["if_crucible_afterheater"])

    #     sim.heat_flux_computation(sim["afterheater_bot"], sim["bnd_afterheater_bot"])
    #     sim.heat_flux_computation(sim["afterheater_bot"], sim["if_crucible_afterheater"])

    # sim.heat_flux_computation(sim["melt"], sim["if_crucible_melt"])
    # sim.heat_flux_computation(sim["melt"], sim["if_melt_crystal"])
    # sim.heat_flux_computation(sim["melt"], sim["bnd_melt"])

    # sim.heat_flux_computation(sim["crystal"], sim["if_melt_crystal"])
    # sim.heat_flux_computation(sim["crystal"], sim["bnd_crystal"])
    # sim.heat_flux_computation(sim["crystal"], sim["if_crystal_seedholder"])

    # # sim.heat_flux_computation(sim["axis_top"], sim["if_seed_axtop"])
    # sim.heat_flux_computation(sim["axis_top"], sim["bnd_axtop"])
    # sim.heat_flux_computation(sim["axis_top"], sim["if_axtop_vessel"])

    # sim.heat_flux_computation(sim["insulation_bot"], sim["if_crucible_ins_bot"])
    # sim.heat_flux_computation(sim["insulation_bot"], sim["bnd_ins_bot"])
    # sim.heat_flux_computation(sim["insulation_bot"], sim["if_ins_bot_adp"])
    # sim.heat_flux_computation(sim["insulation_bot"], sim["if_ins_bot_side"])

    # sim.heat_flux_computation(sim["insulation_side"], sim["bnd_ins_side"])
    # sim.heat_flux_computation(sim["insulation_side"], sim["if_ins_bot_side"])
    # if with_afterheater:
    #     sim.heat_flux_computation(sim["insulation_side"], sim["if_ins_side_top"])

    #     sim.heat_flux_computation(sim["insulation_top"], sim["bnd_ins_top"])
    #     sim.heat_flux_computation(sim["insulation_top"], sim["if_ins_side_top"])

    # sim.heat_flux_computation(sim["axbot_adapter"], sim["if_ins_bot_adp"])
    # sim.heat_flux_computation(sim["axbot_adapter"], sim["bnd_adp"])
    # sim.heat_flux_computation(sim["axbot_adapter"], sim["if_adp_axbt"])

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
    mat_config = ctrl.load_config("./config_mat.yml")
    model = geometry(geo_config, visualize=True)
    # simulation(model, sim_config, mat_config=mat_config)
