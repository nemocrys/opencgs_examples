import numpy as np
import os
import yaml

from objectgmsh import (
    Model,
    Shape,
    MeshControlConstant,
    MeshControlLinear,
    MeshControlExponential,
    cut,
    GeometryError,
)
import opencgs.control as ctrl
from opencgs.setup import ElmerSetupCz
from opencgs.geo import line_from_points, create_2d_shape, create_melt_shape, create_crystal_shape, create_axis_shape, create_heater_shape
import gmsh
from copy import deepcopy
from pyelmer import elmer


occ = gmsh.model.occ
      

def geometry(config, sim_dir="./simdata_elmer/_test", name="cz_resistance", visualize=False):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)
    # initialize geometry model
    model = Model(name)

    crucible_position = config["process_condition"]["crucible_position"]
    coord_update_bot = [0, crucible_position - config["default"]["base_crucible_position"]]

    vessel = create_2d_shape(model, "vessel", **config["vessel"])
    vessel_inside = create_2d_shape(model, "vessel_inside", **config["vessel_inside"])
    occ.cut(vessel.dimtags, vessel_inside.dimtags, removeTool=False)
    axis_bot_steel = create_axis_shape(model, "axis_bot_steel", y_top_update=coord_update_bot[1], **config["axis_bot_steel"])
    axis_bot_steel.set_interface(vessel)
    axis_bot_adapter = create_2d_shape(model, "axis_bot_adapter", coordinate_update=coord_update_bot, **config["axis_bot_adapter"])
    axis_bot_adapter.set_interface(axis_bot_steel)
    axis_bot_graphite = create_2d_shape(model, "axis_bot_graphite", coordinate_update=coord_update_bot, **config["axis_bot_graphite"])
    axis_bot_graphite.set_interface(axis_bot_adapter)
    crucible_support = create_2d_shape(model, "crucible_support", coordinate_update=coord_update_bot, **config["crucible_support"])
    crucible_support.set_interface(axis_bot_graphite)
    crucible = create_2d_shape(model, "crucible", coordinate_update=coord_update_bot, **config["crucible"])
    crucible.set_interface(crucible_support)
    if "phase_if" in config:
        # Fix coordinates, this is required because mgdyn solver slightly shifts them
        config["phase_if"][-1][0] = 0
        config["phase_if"][0][0] = config["process_condition"]["crystal_radius"]

        phase_if_crystal = line_from_points(
            model, config["phase_if"], "phase_if_crystal"
        )
        phase_if_melt = line_from_points(model, config["phase_if"], "phase_if_melt")
    else:
        phase_if_crystal = None
        phase_if_melt = None
    melt, y_meniscus = create_melt_shape(
        model,
        "melt",
        coordinate_update=coord_update_bot,
        bottom_y_coordinate=crucible_position,
        melt_height=config["process_condition"]["melt_level"],
        crystal_radius=config["process_condition"]["crystal_radius"],
        phase_if=phase_if_melt,
        **config["melt"]
    )
    melt.set_interface(crucible)
    crystal, y_crystal_top = create_crystal_shape(
        model,
        "crystal",
        crystal_length=config["process_condition"]["crystal_length"],
        base_crystal_length=config["default"]["base_crystal_length"],
        meniscus_position=y_meniscus,
        phase_if=phase_if_crystal,
        **config["crystal"]
    )
    crystal.set_interface(melt)
    seed_holder = create_2d_shape(model, "seed_holder", coordinate_update=[0, y_crystal_top], **config["seed_holder"])
    seed_holder.set_interface(crystal)
    axis_top_pt1 = create_2d_shape(model, "axis_top_pt1", coordinate_update=[0, y_crystal_top], **config["axis_top_pt1"])
    axis_top_pt1.set_interface(seed_holder)
    axis_top_pt2 = create_2d_shape(model, "axis_top_pt2", coordinate_update=[0, y_crystal_top], **config["axis_top_pt2"])
    axis_top_pt2.set_interface(axis_top_pt1)
    axis_top_pt3 = create_axis_shape(model, "axis_top_pt3", y_bot_update=y_crystal_top, **config["axis_top_pt3"])
    axis_top_pt3.set_interface(axis_top_pt2)
    axis_top_pt3.set_interface(vessel)

    side_heater_bot = create_heater_shape(model, "side_heater_bot", **config["side_heater_bot"])
    side_heater_mid = create_heater_shape(model, "side_heater_mid", **config["side_heater_mid"])
    side_heater_top = create_heater_shape(model, "side_heater_top", **config["side_heater_top"])
    bot_heater = create_heater_shape(model, "bot_heater", **config["bot_heater"])

    insulation_bot_1 = create_2d_shape(model, "insulation_bot_1", **config["insulation_bot_1"])
    insulation_bot_1.set_interface(vessel)
    insulation_bot_2 = create_2d_shape(model, "insulation_bot_2", **config["insulation_bot_2"])
    insulation_bot_2.set_interface(insulation_bot_1)
    insulation_bot_3 = create_2d_shape(model, "insulation_bot_3", **config["insulation_bot_3"])
    insulation_bot_3.set_interface(insulation_bot_2)
    graphite_bot_inside_1 = create_2d_shape(model, "graphite_bot_inside_1", **config["graphite_bot_inside_1"])
    graphite_bot_inside_1.set_interfaces([vessel, insulation_bot_1, insulation_bot_2, insulation_bot_3])
    graphite_bot_inside_2 = create_2d_shape(model, "graphite_bot_inside_2", **config["graphite_bot_inside_2"])
    graphite_bot_inside_2.set_interfaces([graphite_bot_inside_1, insulation_bot_3])
    graphite_bot_outside_1 = create_2d_shape(model, "graphite_bot_outside_1", **config["graphite_bot_outside_1"])
    graphite_bot_outside_1.set_interfaces([vessel, insulation_bot_1, insulation_bot_2, insulation_bot_3])
    graphite_bot_outside_2 = create_2d_shape(model, "graphite_bot_outside_2", **config["graphite_bot_outside_2"])
    graphite_bot_outside_2.set_interfaces([graphite_bot_outside_1, insulation_bot_3])
    insulation_outside = create_2d_shape(model, "insulation_outside", **config["insulation_outside"])
    insulation_outside.set_interface(vessel)
    heater_support_ring = create_2d_shape(model, "heater_support_ring", **config["heater_support_ring"])
    heater_support_ring.set_interface(insulation_outside)
    heatshield_support = create_2d_shape(model, "heatshield_support", **config["heatshield_support"])
    heatshield_support.set_interface(insulation_outside)
    insulation_heatshield_support = create_2d_shape(model, "insulation_heatshield_support", **config["insulation_heatshield_support"])
    insulation_heatshield_support.set_interface(heatshield_support)
    heatshield_outside = create_2d_shape(model, "heatshield_outside", **config["heatshield_outside"])
    heatshield_outside.set_interfaces([heatshield_support, insulation_outside])
    heatshield_inside = create_2d_shape(model, "heatshield_inside", **config["heatshield_inside"])
    heatshield_inside.set_interface(heatshield_outside)
    heatshield_filling = create_2d_shape(model, "heatshield_filling", **config["heatshield_filling"])
    heatshield_filling.set_interfaces([heatshield_inside, heatshield_outside])

    # optional insulation ring on top of heat shield
    if "insulation_ring_top" in config:
        insulation_ring_top = create_2d_shape(model, "insulation_ring_top", **config["insulation_ring_top"])
        insulation_ring_top.set_interface(insulation_outside)
    
    dimtags = []
    for shape in model.get_shapes(dimension=2):
        if shape != vessel_inside: # and shape != vessel:
            dimtags += shape.dimtags
    vessel_inside.geo_ids = cut(vessel_inside.dimtags, dimtags, remove_tool=False)
    # fragment with all shapes at once required, probably related to
    # gmsh issue https://gitlab.onelab.info/gmsh/gmsh/-/issues/2100
    occ.fragment(vessel_inside.dimtags, dimtags)
    model.synchronize()

    # detect interface boundaries between bodies
    if_vessel_axbotsteel = Shape(model, 1, "if_vessel_axbotsteel", vessel.get_interface(axis_bot_steel))
    if_axbotsteel_adapter = Shape(model, 1, "if_axbotsteel_adapter", axis_bot_steel.get_interface(axis_bot_adapter))
    if_axbotgraph_adapter = Shape(model, 1, "if_axbotgraph_adapter", axis_bot_graphite.get_interface(axis_bot_adapter))
    if_axbotgraph_crcsup = Shape(model, 1, "if_axbotgraph_crcsup", axis_bot_graphite.get_interface(crucible_support))
    if_crucsup_cruc = Shape(model, 1, "if_crucsup_cruc", crucible_support.get_interface(crucible))
    if_crucible_melt = Shape(model, 1, "if_crucible_melt", crucible.get_interface(melt))
    if_melt_crystal = Shape(model, 1, "if_melt_crystal", melt.get_interface(crystal))
    if_crystal_seedholder = Shape(model, 1, "if_crystal_seedholder", crystal.get_interface(seed_holder))
    if_seedholder_axtop1 = Shape(model, 1, "if_seedholder_axtop1", seed_holder.get_interface(axis_top_pt1))
    if_axtop1_2 = Shape(model, 1, "if_axtop1_2", axis_top_pt1.get_interface(axis_top_pt2))
    if_axtop2_3 = Shape(model, 1, "if_axtop2_3", axis_top_pt2.get_interface(axis_top_pt3))
    if_axtop3_vessel = Shape(model, 1, "if_axtop3_vessel", axis_top_pt3.get_interface(vessel))
    if_insbot1_vessel = Shape(model, 1, "if_insbot1_vessel", insulation_bot_1.get_interface(vessel))
    if_insbot1_insbot2 = Shape(model, 1, "if_insbot1_insbot2", insulation_bot_1.get_interface(insulation_bot_2))
    if_insbot2_insbot3 = Shape(model, 1, "if_insbot2_insbot3", insulation_bot_2.get_interface(insulation_bot_3))
    if_graphbotin1_vessel = Shape(model, 1, "if_graphbotin1_vessel", graphite_bot_inside_1.get_interface(vessel))
    if_graphbotin1_2 = Shape(model, 1, "if_graphbotin1_2", graphite_bot_inside_1.get_interface(graphite_bot_inside_2))
    if_graphbotin1_insbot1 = Shape(model, 1, "if_graphbotin1_insbot1", graphite_bot_inside_1.get_interface(insulation_bot_1))
    if_graphbotin1_insbot2 = Shape(model, 1, "if_graphbotin1_insbot2", graphite_bot_inside_1.get_interface(insulation_bot_2))
    if_graphbotin1_insbot3 = Shape(model, 1, "if_graphbotin1_insbot3", graphite_bot_inside_1.get_interface(insulation_bot_3))
    if_graphbotin2_insbot3 = Shape(model, 1, "if_graphbotin2_insbot3", graphite_bot_inside_2.get_interface(insulation_bot_3))
    if_graphbotout1_vessel = Shape(model, 1, "if_graphbotout1_vessel", graphite_bot_outside_1.get_interface(vessel))
    if_graphbotout1_2 = Shape(model, 1, "if_graphbotout1_2", graphite_bot_outside_1.get_interface(graphite_bot_outside_2))
    if_graphbotout1_insbot1 = Shape(model, 1, "if_graphbotout1_insbot1", graphite_bot_outside_1.get_interface(insulation_bot_1))
    if_graphbotout1_insbot2 = Shape(model, 1, "if_graphbotout1_insbot2", graphite_bot_outside_1.get_interface(insulation_bot_2))
    if_graphbotout1_insbot3 = Shape(model, 1, "if_graphbotout1_insbot3", graphite_bot_outside_1.get_interface(insulation_bot_3))
    if_graphbotout2_insbot3 = Shape(model, 1, "if_graphbotout2_insbot3", graphite_bot_outside_2.get_interface(insulation_bot_3))
    if_insout_vessel = Shape(model, 1, "if_insout_vessel", insulation_outside.get_interface(vessel))
    if_insout_graphbotout1 = Shape(model, 1, "if_insout_graphbotout1", insulation_outside.get_interface(graphite_bot_outside_1))
    if_insout_graphbotout2 = Shape(model, 1, "if_insout_graphbotout2", insulation_outside.get_interface(graphite_bot_outside_2))
    if_insout_heaterring = Shape(model, 1, "if_insout_heaterring", insulation_outside.get_interface(heater_support_ring))
    if_insout_heatshildsup = Shape(model, 1, "if_insout_heatshildsup", insulation_outside.get_interface(heatshield_support))
    if_insout_heatshildout = Shape(model, 1, "if_insout_heatshildout", insulation_outside.get_interface(heatshield_outside))
    if_heatshildsup_insheatshieldsup = Shape(model, 1, "if_heatshildsup_insheatshieldsup", heatshield_support.get_interface(insulation_heatshield_support))
    if_heatshildsup_heatshieldout = Shape(model, 1, "if_heatshildsup_heatshildout", heatshield_support.get_interface(heatshield_outside))
    if_heatshildout_heatshieldin = Shape(model, 1, "if_heatshildout_heatshieldin", heatshield_outside.get_interface(heatshield_inside))
    if_heatshieldout_filling = Shape(model, 1, "if_heatshieldout_filling", heatshield_outside.get_interface(heatshield_filling))
    if_heatshildin_filling = Shape(model, 1, "if_heatshildin_filling", heatshield_inside.get_interface(heatshield_filling))

    if "insulation_ring_top" in config:
        if_insringtop_insout = Shape(model, 1, "if_insringtop_insout", insulation_outside.get_interface(insulation_ring_top))

    # detect surface boundaries
    bnd_symmetry_axis = Shape(model, 1, "bnd_symmetry_axis", model.symmetry_axis)
    surf_vesselinside = Shape(model, 1, "surf_vesselinside", vessel.get_interface(vessel_inside))
    surf_axbotsteel = Shape(model, 1, "surf_axbotsteel", axis_bot_steel.get_interface(vessel_inside))
    surf_axbotadpter = Shape(model, 1, "surf_axbotadpter", axis_bot_adapter.get_interface(vessel_inside))
    surf_axbotgraph = Shape(model, 1, "surf_axbotgraph", axis_bot_graphite.get_interface(vessel_inside))
    surf_crucsup = Shape(model, 1, "surf_crucsup", crucible_support.get_interface(vessel_inside))
    surf_cruc = Shape(model, 1, "surf_cruc", crucible.get_interface(vessel_inside))
    surf_melt = Shape(model, 1, "surf_melt", melt.get_interface(vessel_inside))
    surf_crys = Shape(model, 1, "surf_crys", crystal.get_interface(vessel_inside))
    surf_seedhold = Shape(model, 1, "surf_seedhold", seed_holder.get_interface(vessel_inside))
    surf_axtop1 = Shape(model, 1, "surf_axtop1", axis_top_pt1.get_interface(vessel_inside))
    surf_axtop2 = Shape(model, 1, "surf_axtop2", axis_top_pt2.get_interface(vessel_inside))
    surf_axtop3 = Shape(model, 1, "surf_axtop3", axis_top_pt3.get_interface(vessel_inside))
    surf_botheater = Shape(model, 1, "surf_botheater", bot_heater.get_interface(vessel_inside))
    surf_sideheaterbot = Shape(model, 1, "surf_sideheaterbot", side_heater_bot.get_interface(vessel_inside))
    surf_sideheatermid = Shape(model, 1, "surf_sideheatermid", side_heater_mid.get_interface(vessel_inside))
    surf_sideheatertop = Shape(model, 1, "surf_sideheatertop", side_heater_top.get_interface(vessel_inside))
    surf_insbot3 = Shape(model, 1, "surf_insbot3", insulation_bot_3.get_interface(vessel_inside))
    surf_graphbotin1 = Shape(model, 1, "surf_graphbotin1", graphite_bot_inside_1.get_interface(vessel_inside))
    surf_graphbotin2 = Shape(model, 1, "surf_graphbotin2", graphite_bot_inside_2.get_interface(vessel_inside))
    surf_graphbotout1 = Shape(model, 1, "surf_graphbotout1", graphite_bot_outside_1.get_interface(vessel_inside))
    surf_graphbotout2 = Shape(model, 1, "surf_graphbotout2", graphite_bot_outside_2.get_interface(vessel_inside))
    surf_insout = Shape(model, 1, "surf_insout", insulation_outside.get_interface(vessel_inside))
    surf_heaterring = Shape(model, 1, "surf_heaterring", heater_support_ring.get_interface(vessel_inside))
    surf_heatshieldsup = Shape(model, 1, "surf_heatshieldsup", heatshield_support.get_interface(vessel_inside))
    surf_insheatshieldsup = Shape(model, 1, "surf_insheatshieldsup", insulation_heatshield_support.get_interface(vessel_inside))
    surf_heatshieldout = Shape(model, 1, "surf_heatshieldout", heatshield_outside.get_interface(vessel_inside))
    surf_heatshieldin = Shape(model, 1, "surf_heatshieldin", heatshield_inside.get_interface(vessel_inside))
    bnd_vessel_outside = Shape(
        model,
        1,
        "bnd_vessel_outside",
        [
            x for x in vessel.boundaries if x not in 
                model.symmetry_axis
                + if_vessel_axbotsteel.geo_ids
                + if_axtop3_vessel.geo_ids
                + if_insout_vessel.geo_ids
                + if_insbot1_vessel.geo_ids
                + if_graphbotin1_vessel.geo_ids
                + if_graphbotout1_vessel.geo_ids
                + surf_vesselinside.geo_ids
        ]
    )
    # model.remove_shape(vessel_inside)
    if "insulation_ring_top" in config:
        surf_insringtop = Shape(model, 1, "surf_insringtop", insulation_ring_top.get_interface(vessel_inside))

    model.make_physical()
    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()
    
    max_meshsize = 20.e-3
    for shape in model.get_shapes(3):
        MeshControlLinear(model, shape, shape.mesh_size, max_meshsize)
    MeshControlExponential(model, if_melt_crystal, melt.mesh_size / 3, exp=1.6, fact=3)
    MeshControlExponential(model, surf_crys, crystal.mesh_size / 3, exp=1.6, fact=3)
    MeshControlExponential(model, if_crucible_melt, melt.mesh_size / 3, exp=1.6, fact=3)
    MeshControlExponential(model, surf_melt, melt.mesh_size / 3, exp=1.6, fact=3)


    model.generate_mesh(**config["mesh"])

    if visualize:
        model.show()
    model.write_msh(f"{sim_dir}/case.msh")
    print(model)
    model.close_gmsh()
    return model

def simulation(model,
               config,
               sim_dir="./simdata_elmer/_test",
               mat_config={},
               with_flow=False,
               u_scaling=1,
               flow_variable="UMean_avg300s",
               flow_variable_old="UMean_avg300s_old",
               flow_relaxation_factor=1,
               heat_conductivity_scaling=1,  # not implemented
               ):
    sim = ElmerSetupCz(
        **config["general"],
        sim_dir=sim_dir,
        probes=config["probes"],
        smart_heater=config["smart-heater"],
        materials_dict=mat_config
    )
    if "solver-update" in config:
        sim.solver_update=config["solver-update"]

    sim.solver_heat.data.update({"Update Gebhardt Factors": True})
    sim.sim.solvers["ResultOutputSolver"].data.update({"Vtu Part collection": True}) #, "Save Bulk Only": True})
    # sim.sim.solvers["ResultOutputSolver"].data.update({"Discontinuous Bodies": True})  # TODO FIXME "Discontinous Bodies = True" breaks lorentz_force_to_txt

    solver_mgdyn = elmer.Solver(sim.sim, "mgdyn", config["MagnetoDynamics2DHarmonic"])
    solver_mgdyn_calc = elmer.Solver(sim.sim, "mgdyncalc", config["MagnetoDynamicsCalcFields"])
    omega = 2*np.pi* config["kristmag_f1"]["frequency"]
    solver_mgdyn.data.update({"Angular Frequency": omega})
    solver_mgdyn_calc.data.update({"Angular Frequency": omega})
    sim._eqn_main.solvers += [solver_mgdyn, solver_mgdyn_calc]
    eqn_vacuum = elmer.Equation(sim.sim, "eqn_vac", sim._eqn_main.solvers[1:])  # all solvers but HeatSolver (which is at position 0)

    if with_flow:
        solver_gmsh_input = elmer.Solver(
            sim.sim, 
            "gmsh_input",
            config["GmshInput"]
        )
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
    # sim._eqn_main.solvers += [solver_flux]  # this causes wrong T-gradient at melt boundary. It may be fixed using a different order when adding the bodies.
    eqn_melt = elmer.Equation(sim.sim, "eqn_melt", sim._eqn_main.solvers + [solver_flux])
    if with_flow:
        eqn_melt.data.update({"Convection": "Constant"})

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

    # bodies
    bot_heater = sim.add_resistance_heater(model["bot_heater"], config["heating"]["power_bot"])
    bot_heater.data.update({"name": "bot_heater"})
    side_heater_bot = sim.add_resistance_heater(model["side_heater_bot"], config["heating"]["bower_side_bot"])
    side_heater_bot.data.update({"name": "side_heater_bot"})
    side_heater_mid = sim.add_resistance_heater(model["side_heater_mid"], config["heating"]["power_side_mid"])
    side_heater_mid.data.update({"name": "side_heater_mid"})
    side_heater_top = sim.add_resistance_heater(model["side_heater_top"], config["heating"]["power_side_top"])
    side_heater_top.data.update({"name": "side_heater_top"})

    mat = elmer.Material(sim.sim,
                         "mat_heater_zero-electric-conductivity",
                         deepcopy(bot_heater.material.data)
                         )
    mat.data.update({"Electric Conductivity": 0})

    for bdy, name in zip(
        [bot_heater, side_heater_bot, side_heater_mid, side_heater_top],
        ["bot", "side_bot", "side_mid", "side_top"],
    ):
        current_dty = config["kristmag_f1"][name]["current_dty"]
        phase = config["kristmag_f1"][name]["phase"]
        if current_dty != 0:
            bdy.material = mat
            current_dty_re = np.cos(2 * np.pi * phase / 360) * current_dty
            current_dty_im = np.sin(2 * np.pi * phase / 360) * current_dty
            bdy.body_force.data.update({
                "Current density": f"real {current_dty_re}",
                "Current density Im": f"real {current_dty_im}",
            })

    crystal = sim.add_crystal(model["crystal"])
    crystal.data.update({"name": "crystal"})
    for body in [
        "vessel",
        "axis_bot_steel",
        "axis_bot_adapter",
        "axis_bot_graphite",
        "crucible_support",
        "crucible",
        "melt",
        "seed_holder",
        "axis_top_pt1",
        "axis_top_pt2",
        "axis_top_pt3",
        "insulation_bot_1",
        "insulation_bot_2",
        "insulation_bot_3",
        "graphite_bot_inside_1",
        "graphite_bot_inside_2",
        "graphite_bot_outside_1",
        "graphite_bot_outside_2",
        "insulation_outside",
        "heater_support_ring",
        "heatshield_support",
        "insulation_heatshield_support",
        "heatshield_outside",
        "heatshield_inside",
        "heatshield_filling",
    ]:
        body_ = sim.add_body(model[body])
        body_.data.update({"name": body})
    
    try:
        body_ = sim.add_body(model["insulation_ring_top"])
        body_.data.update({"name": "insulation_ring_top"})
    except GeometryError:
        pass  # simulation without this additional ring

    vacuum = sim.add_body(model["vessel_inside"])
    vacuum.data.update({"name": "vessel_inside"})
    vacuum.equation = eqn_vacuum
    vacuum.initial_condition.data.pop("Temperature")  # it's not required and makes the simulation diverging
    melt = sim.sim.bodies["melt"]
    melt.equation = eqn_melt
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
    # phase interface
    sim.add_phase_interface(model["if_melt_crystal"])

    # other interfaces
    for interface in [
        'if_vessel_axbotsteel',
        'if_axbotsteel_adapter',
        'if_axbotgraph_adapter',
        'if_axbotgraph_crcsup',
        'if_crucsup_cruc',
        'if_crucible_melt',
        'if_crystal_seedholder',
        'if_seedholder_axtop1',
        'if_axtop1_2',
        'if_axtop2_3',
        'if_axtop3_vessel',
        'if_insbot1_vessel',
        'if_insbot1_insbot2',
        'if_insbot2_insbot3',
        'if_graphbotin1_vessel',
        'if_graphbotin1_2',
        'if_graphbotin1_insbot1',
        'if_graphbotin1_insbot2',
        'if_graphbotin1_insbot3',
        'if_graphbotin2_insbot3',
        'if_graphbotout1_vessel',
        'if_graphbotout1_2',
        'if_graphbotout1_insbot1',
        'if_graphbotout1_insbot2',
        'if_graphbotout1_insbot3',
        'if_graphbotout2_insbot3',
        'if_insout_vessel',
        'if_insout_graphbotout1',
        'if_insout_graphbotout2',
        'if_insout_heaterring',
        'if_insout_heatshildsup',
        'if_insout_heatshildout',
        'if_heatshildsup_insheatshieldsup',
        'if_heatshildsup_heatshildout',
        'if_heatshildout_heatshieldin',
        'if_heatshieldout_filling',
        'if_heatshildin_filling',
    ]:
        bnd = sim.add_interface(model[interface])
        # bnd.save_scalars = True
        # bnd.save_line = True
    # sim.add_temperature_boundary(model["if_crucible_melt"], 1700)
    
    try:
        sim.add_interface(model["if_insringtop_insout"])
    except GeometryError:
        pass  # simulation without this additional ring

    # surfaces with radiation
    for surface in [
        'surf_vesselinside',
        'surf_axbotsteel',
        'surf_axbotadpter',
        'surf_axbotgraph',
        'surf_crucsup',
        'surf_cruc',
        'surf_melt',
        'surf_crys',
        'surf_seedhold',
        'surf_axtop1',
        'surf_axtop2',
        'surf_axtop3',
        'surf_botheater',
        'surf_sideheaterbot',
        'surf_sideheatermid',
        'surf_sideheatertop',
        'surf_insbot3',
        'surf_graphbotin1',
        'surf_graphbotin2',
        'surf_graphbotout1',
        'surf_graphbotout2',
        'surf_insout',
        'surf_heaterring',
        'surf_heatshieldsup',
        'surf_insheatshieldsup',
        'surf_heatshieldout',
        'surf_heatshieldin',
    ]:
        bnd = sim.add_radiation_boundary(model[surface])
        # bnd.save_scalars = True
        # bnd.save_line = True
    try:
        sim.add_radiation_boundary(model["surf_insringtop"])
    except GeometryError:
        pass  # simulation without this additional ring

    # bc = sim.add_radiation_boundary(model["surf_heater_bot"])
    # bc.data.update({"Temperature": 500})
    # bc = sim.add_radiation_boundary(model["surf_heater_side"])
    # bc.data.update({"Temperature": 500})

    bnd_outside = sim.add_temperature_boundary(
        model["bnd_vessel_outside"],
        **config["boundaries"]["vessel_outside"]
    )
    bnd_outside.zero_potential = True
    # bnd_outside.save_line = True
    # bnd_outside.save_scalars = True
    # symmetry axis
    sim.add_interface(model["bnd_symmetry_axis"], movement=[0, None])

    sim.sim.boundaries["if_crucible_melt"].data.update({"Save Line Crc Melt": "Logical True"})
    sim.sim.boundaries["surf_melt"].data.update({"Save Line Melt Surf": "Logical True"})
    sim.sim.boundaries["if_melt_crystal"].data.update({"Save Line Melt Crys": "Logical True"})

    # export
    sim.export()
    return sim

if __name__ == "__main__":
    geo_config = ctrl.load_config("./config_geo.yml")
    sim_config = ctrl.load_config("./config_sim.yml")
    mat_config = ctrl.load_config("./config_mat.yml")
    model = geometry(geo_config, visualize=True)
    simulation(model, sim_config, mat_config=mat_config)
