"""
Geometry and simulation setup (in two separate functions).
"""

import numpy as np
import os
import sys
import yaml
import gmsh
import pickle

import pyelmer.elmer as elmer
from pyelmer.execute import run_elmer_grid, run_elmer_solver
from pyelmer.post import scan_logfile
from objectgmsh import (
    Model,
    Shape,
    MeshControlConstant,
    MeshControlLinear,
    MeshControlExponential,
    cut
)
import opencgs.control as ctrl
from opencgs.setup import ElmerSetupCz
import opencgs.geo as geo
from copy import deepcopy

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
factory = gmsh.model.occ

def geometry(
    config,
    sim_dir="./simdata/_test",
    name="cz_induction",
    visualize=False,
    with_afterheater=False,
    mat_config_file="config_mat.yml",
    dim=3,
):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    config_mat = ctrl.load_config(mat_config_file)

    # initialize geometry model
    model = Model(name)
    gmsh.option.setNumber("General.NumThreads", 8)  # omp parallelized meshing

    # helper shape to cut inductor ends at outside of vessel
    surrounding_temp = geo.surrounding(
        model, dim, **config["surrounding_temp"], name="surrounding_temp"
    )
    cutbox = geo.surrounding(model, dim, **config["cutbox"], name="cutbox")
    model.remove_shape(surrounding_temp)

    # draw phase interface if predefined (for transient simulation)
    if config["crystal"]["current_length"] == 0.01:
        interface = config["phase_if_length=0.01"]
    elif config["crystal"]["current_length"] == 0.087:
        interface = config["phase_if_length=0.087"]
    else:
        raise ValueError("Crystal length not supported.")
    phase_if_crystal = geo.line_from_points(model, interface, "phase_if_crystal")
    phase_if_melt = geo.line_from_points(model, interface, "phase_if_melt")

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
    crucible = geo.crucible(model, dim, **config["crucible"])
    melt = geo.melt(
        model,
        dim,
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
        model, dim, **config["crystal"], melt=melt, phase_if=phase_if_crystal
    )
    seedholder = geo.cylinder_shape(
        model, dim, **config["seedholder"], bot_shape=crystal, name="seedholder"
    )
    shaft_al2o3 = geo.cylinder_shape(
        model, dim, **config["shaft_al2o3"], bot_shape=seedholder, name="shaft_al2o3"
    )
    axtop_adapter = geo.cylinder_shape(
        model, dim, **config["axtop_adapter"], bot_shape=shaft_al2o3, name="axtop_adapter"
    )
    ins_bot = geo.stacked_cylinder_shape(
        model,
        dim,
        **config["insulation_bot"],
        top_shape=crucible,
        name="insulation_bot",
    )
    ins_side = geo.cylinder_shape(
        model,
        dim,
        bot_shape=ins_bot,
        name="insulation_side",
        **config["insulation_side_two_rings"],
    )
    ins_side.set_interface(ins_bot)
    adp = geo.stacked_cylinder_shape(
        model, dim, **config["axbot_adapter"], top_shape=ins_bot, name="axbot_adapter"
    )
    ax_bt = geo.cylinder_shape(
        model, dim, **config["axis_bt"], top_shape=adp, name="axis_bt"
    )
    inductor = geo.coil_from_points(
        model, dim, dy=ax_bt.params.X0[1], **config["inductor_measurement_calliper"]
    )
    factory.rotate(inductor.dimtags, 0, 0, 0, 0, 1, 0, np.pi/2)
    cut(inductor.dimtags, cutbox.dimtags, remove_tool=False)
    inductor.set_interface(cutbox)
    vessel = geo.vessel(model, dim, **config["vessel"], adjacent_shapes=[ax_bt])
    cut(vessel.dimtags, inductor.dimtags, False)
    inductor.set_interface(vessel)
    # lower power supply
    circle = factory.add_circle(
        0, inductor.params.supply_bot_h, 0.16, r=inductor.params.radius_pipe * 2
    )
    pipe = factory.extrude([(1, circle)], 0, 0, -0.05)[1][1]
    factory.rotate([(2, pipe)], 0, 0, 0, 0, 1, 0, np.pi/2)
    tag = factory.fragment(vessel.dimtags, [(2, pipe)])[0][1][1]
    vessel_power_supply_bottom = Shape(model, 3, "vessel_power_supply_bottom", [tag])
    vessel_power_supply_bottom.mesh_size = vessel.mesh_size
    vessel_power_supply_bottom.params.material = "insulating-steel"
    # upper power supply
    circle = factory.add_circle(
        0, inductor.params.supply_top_h, 0.16, r=inductor.params.radius_pipe * 2
    )
    pipe = factory.extrude([(1, circle)], 0, 0, -0.05)[1][1]
    factory.rotate([(2, pipe)], 0, 0, 0, 0, 1, 0, np.pi/2)
    dimtags = factory.fragment(vessel.dimtags, [(2, pipe)])
    vessel.geo_ids = [dimtags[0][0][1]]
    vessel_power_supply_top = Shape(
        model, 3, "vessel_power_supply_top", [dimtags[0][1][1]]
    )
    vessel_power_supply_top.mesh_size = vessel.mesh_size
    vessel_power_supply_top.params.material = "insulating-steel"
    factory.synchronize()

    ax_top = geo.axis_top(
        model, dim, **config["axis_top"], bot_shape=axtop_adapter, vessel=vessel
    )
    ax_top.set_interfaces([axtop_adapter, shaft_al2o3])
    exhaust_hood = geo.exhaust_hood(
        model, dim, **config["exhaust_hood"], y_bot=vessel.params.X0[1] + vessel.params.t
    )
    filling = geo.filling(model, dim, **config["filling"], vessel=vessel)

    # fix issue with duplicate surface
    shapes = model.get_shapes(3)
    shapes[0].set_interfaces(shapes[1:])

    model.synchronize()

    # boundaries
    bnd_melt = Shape(model, dim-1, "bnd_melt", melt.get_interface(filling))
    bnd_seedholder = Shape(
        model, 1, "bnd_seedholder", seedholder.get_interface(filling)
    )
    bnd_shaft_al2o3 = Shape(
        model, dim-1, "bnd_shaft_al2o3", shaft_al2o3.get_interface(filling)
    )
    bnd_axtop_adapter = Shape(
        model, dim-1, "bnd_axtop_adapter", axtop_adapter.get_interface(filling)
    )
    bnd_crystal = Shape(model, dim-1, "bnd_crystal", crystal.get_interface(filling))
    bnd_axtop = Shape(model, dim-1, "bnd_axtop", ax_top.get_interface(filling))
    bnd_crucible_outside = Shape(
        model, dim-1, "bnd_crucible_outside", crucible.get_interface(filling)
    )
    bnd_ins_bot = Shape(model, dim-1, "bnd_ins_bot", ins_bot.get_interface(filling))
    bnd_ins_side = Shape(model, dim-1, "bnd_ins_side", ins_side.get_interface(filling))
    bnd_adp = Shape(model, dim-1, "bnd_adp", adp.get_interface(filling))
    bnd_axbt = Shape(model, dim-1, "bnd_axbt", ax_bt.get_interface(filling))
    bnd_exhaust_hood = Shape(
        model, dim-1, "bnd_exhaust_hood", exhaust_hood.get_interface(filling)
    )
    bnd_vessel_inside = Shape(
        model, dim - 1, "bnd_vessel_inside", vessel.get_interface(filling)
    )
    bnd_inductor = Shape(
        model, dim - 1, "bnd_inductor", inductor.get_interface(filling)
    )
    bnd_inductor_end_bottom = Shape(
        model, dim - 1, "bnd_inductor_end_bottom", [inductor.get_interface(cutbox)[0]]
    )
    bnd_inductor_end_top = Shape(
        model, dim - 1, "bnd_inductor_end_top", [inductor.get_interface(cutbox)[1]]
    )
    model.remove_shape(cutbox)  # not required any longer

    # interfaces
    if_crucible_melt = Shape(model, dim-1, "if_crucible_melt", crucible.get_interface(melt))
    if_melt_crystal = Shape(model, dim-1, "if_melt_crystal", melt.get_interface(crystal))
    if_crystal_seedholder = Shape(
        model, dim-1, "if_crystal_seedholder", crystal.get_interface(seedholder)
    )
    if_axtop_vessel = Shape(model, dim-1, "if_axtop_vessel", ax_top.get_interface(vessel))
    if_crucible_ins_bot = Shape(
        model, dim-1, "if_crucible_ins_bot", crucible.get_interface(ins_bot)
    )
    if_ins_bot_adp = Shape(model, dim-1, "if_ins_bot_adp", ins_bot.get_interface(adp))
    if_ins_bot_side = Shape(
        model, dim-1, "if_ins_bot_side", ins_bot.get_interface(ins_side)
    )
    if_adp_axbt = Shape(model, dim-1, "if_adp_axbt", adp.get_interface(ax_bt))
    if_axbt_vessel = Shape(model, dim-1, "if_axbt_vessel", ax_bt.get_interface(vessel))
    if_inductor_vessel_power_supply = Shape(
        model,
        dim - 1,
        "if_inductor_vessel_power_supply",
        inductor.get_interface(vessel_power_supply_bottom)
        + inductor.get_interface(vessel_power_supply_top),
    )
    
    # helper boundaries
    if_vessel_powersupply = Shape(
        model,
        dim - 1,
        "if_vessel_powersupply",
        vessel.get_interface(vessel_power_supply_bottom)
        + vessel.get_interface(vessel_power_supply_top),
    )
    bnd_vessel_outside = Shape(
        model,
        dim - 1,
        "bnd_vessel_outside",
        [
            x
            for x in vessel.boundaries
            if x
            not in if_axbt_vessel.geo_ids
            + if_axtop_vessel.geo_ids
            + bnd_vessel_inside.geo_ids
            + if_vessel_powersupply.geo_ids
        ],
    )
    bnd_vessel_powersupply_inside = Shape(
        model,
        dim - 1,
        "bnd_vessel_powersupply_inside",
        vessel_power_supply_top.get_interface(filling)
        + vessel_power_supply_bottom.get_interface(filling),
    )
    bnd_vessel_powersupply_outside = Shape(
        model,
        dim - 1,
        "bnd_vessel_powersupply_outside",
        [
            x
            for x in vessel_power_supply_top.boundaries
            + vessel_power_supply_bottom.boundaries
            if x
            not in bnd_vessel_powersupply_inside.geo_ids
            + if_vessel_powersupply.geo_ids
            + if_inductor_vessel_power_supply.geo_ids
        ],
    )

    model.make_physical()
    # model.show()

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
    for shape in shapes:
        MeshControlLinear(model, shape, shape.mesh_size, filling.mesh_size, dist_end=filling.mesh_size*2)

    MeshControlExponential(model, bnd_melt, melt.mesh_size / 3, exp=1.6, shapes=[melt, crucible])
    MeshControlExponential(model, if_melt_crystal, crystal.mesh_size / 2, exp=1.6, shapes=[crystal, melt])


    # MeshControlExponential(model, if_crucible_melt, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(
        model, bnd_crucible_outside, crucible.mesh_size / 3, exp=1.6, shapes=[crucible, melt, ins_bot, ins_side, adp]
    )
    MeshControlExponential(model, inductor, inductor.mesh_size, exp=1.6, fact=2)
    model.generate_mesh(**config["mesh"])

    if visualize:
        model.show()
    model.write_msh(f"{sim_dir}/case.msh")
    print(model)
    model.close_gmsh()
    with open(f"{sim_dir}/model.pickle", "wb") as f:
        pickle.dump(model, f)
    return model


def simulation(
    config,
    model=None,
    sim_dir="./simdata/_test",
):
    if model is None:
        with open(f"{sim_dir}/model.pickle", "rb") as f:
            model = pickle.load(f)

    # read parameters from config
    omega = 2 * np.pi * config["heating"]["frequency"]
    rad_s2s = config["general"]["rad_s2s"]

    # setup solvers
    sim = elmer.load_simulation("3D_steady", "config_elmer.yml")
    sim.settings.update({"Angular Frequency": omega})

    solver_MG = elmer.load_solver("MGDynamics", sim, "config_elmer.yml")
    solver_MG.data.update({"Angular Frequency": omega})
    solver_MGCalc = elmer.load_solver("MGDynamicsCalc", sim, "config_elmer.yml")
    solver_MGCalc.data.update({"Angular Frequency": omega})
    solver_heat = elmer.load_solver("HeatSolver", sim, "config_elmer.yml")
    solver_out = elmer.load_solver("ResultOutputSolver", sim, "config_elmer.yml")
    eqn_main = elmer.Equation(
        sim,
        "eqn_main",
        [
            solver_MG,
            solver_MGCalc,
            solver_heat,
        ],
        {"name": "eqn_main"},
    )
    eqn_filling = elmer.Equation(
        sim,
        "eqn_filling",
        [
            solver_MG,
            solver_MGCalc,
        ],
        {"name": "eqn_filling"},
    )

    solver_coil = elmer.load_solver("CoilSolver", sim, "config_elmer.yml")
    solver_coil.data.update({"Desired Coil Current": config["heating"]["current"]})
    eqn_coil = elmer.Equation(
        sim, "eqn_coil", [solver_coil, solver_MG, solver_MGCalc, solver_heat]
    )
    solver_MG.data.update(
        {
            "Use Elemental CoilCurrent": "Logical True",
            "Fix Input Current Density": False,
        }
    )
    current_source = elmer.BodyForce(sim, "current source", {"name": "current source"})
    joule_heat = elmer.BodyForce(
        sim,
        "joule_heat",
        {
            "Joule Heat": "Logical True",
            "Smart Heater Control": "Logical True",
            "Heat Source": 0,
        },
    )

    # output for coupling
    solver_flux = elmer.Solver(sim, "solver_flux")
    solver_flux.data = config["FluxSolver"]
    eqn_melt = elmer.Equation(sim, "eqn_melt", eqn_main.solvers + [solver_flux])
    save_line_crc_melt = elmer.Solver(sim, "save_line_crc_melt")
    save_line_crc_melt.data = deepcopy(config["save_line_solver"])
    save_line_crc_melt.data.update(
        {
            "Equation": '"SaveLineCrcMelt"',
            "Save Mask": 'String "Save Line Crc Melt"',
            "Filename": '"save_line_crc_melt.dat"',
        }
    )
    save_line_melt_surf = elmer.Solver(sim, "save_line_melt_surf")
    save_line_melt_surf.data = deepcopy(config["save_line_solver"])
    save_line_melt_surf.data.update(
        {
            "Equation": '"SaveLineMeltSurf"',
            "Save Mask": 'String "Save Line Melt Surf"',
            "Filename": '"save_line_melt_surf.dat"',
        }
    )
    save_line_melt_crys = elmer.Solver(sim, "save_line_melt_crys")
    save_line_melt_crys.data = deepcopy(config["save_line_solver"])
    save_line_melt_crys.data.update(
        {
            "Equation": '"SaveLineMeltCrys"',
            "Save Mask": 'String "Save Line Melt Crys"',
            "Filename": '"save_line_melt_crys.dat"',
        }
    )

    # compute correct growth velocity
    mat_config = ctrl.load_config("config_mat.yml")
    area_melt = np.pi * model["crucible"].params.r_in**2
    area_crystal = np.pi * model["crystal"].params.r**2
    dty_crystal = mat_config[model["crystal"].params.material]["Density"]
    dty_melt = mat_config[model["melt"].params.material]["Density"]
    v_pull = config["general"]["v_pull"]
    print(f"Pulling velocity: {v_pull:.3f} mm/min")
    print(f"Area ratio: {area_crystal / area_melt:.4f}")
    print(f"Density ratio: {dty_crystal / dty_melt}")
    v_pull += v_pull*  area_crystal / area_melt * dty_crystal / dty_melt
    print(f"Growth velocity: {v_pull:.4f} mm/min")

    # add bodies and materials
    material = elmer.load_material(
        model["inductor"].params.material, sim, "config_mat.yml"
    )
    material.data.update({"name": model["inductor"].params.material})
    ic = elmer.InitialCondition(sim, "inductor", {"Temperature": model["inductor"].params.T_init})
    inductor = elmer.Body(sim, "inductor", [model["inductor"].ph_id])
    inductor.material = material
    inductor.equation = eqn_coil
    inductor.body_force = current_source
    inductor.initial_condition = ic
    inductor.data.update({"name": "inductor"})

    # functions for T-dependent matrices
    for function in mat_config["functions"]:
        sim.intro_text += function

    for name in [
        "crystal",
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
        "vessel_power_supply_bottom",
        "vessel_power_supply_top",
    ]:
        if model[name].params.material not in sim.materials:
            material = elmer.load_material(
                model[name].params.material, sim, "config_mat.yml"
            )
            material.data.update({"name": model[name].params.material})
        else:
            material = sim.materials[model[name].params.material]
        ic = elmer.InitialCondition(sim, name, {"Temperature": model[name].params.T_init})
        body = elmer.Body(sim, name, [model[name].ph_id])
        body.material = material
        body.equation = eqn_main
        body.body_force = joule_heat
        body.initial_condition = ic
        body.data.update({"name": name})


    melt = elmer.Body(sim, "melt", [model["melt"].ph_id])
    material = elmer.load_material(
        model["melt"].params.material, sim, "config_mat.yml"
    )
    material.data.update({"name": model["melt"].params.material})
    material.data.pop("Surface Tension")  # get rid of parameters not used by elmer
    material.data.pop("Beta")  # get rid of parameters not used by elmer
    ic = elmer.InitialCondition(sim, "melt", {"Temperature": model["melt"].params.T_init})
    melt.material = material
    melt.initial_condition = ic
    melt.equation = eqn_melt
    melt.body_force = joule_heat
    melt.data.update({"name": "melt"})

    filling = elmer.Body(sim, "filling", [model["filling"].ph_id])
    material = elmer.load_material(
        model["filling"].params.material, sim, "config_mat.yml"
    )
    material.data.update({"name": model["filling"].params.material})
    filling.material = material
    filling.equation = eqn_filling
    filling.data.update({"name": "filling"})

    # boundaries without convection
    for name in [
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
        "bnd_vessel_powersupply_inside",
    ]:
        bnd = elmer.Boundary(sim, name, [model[name].ph_id])
        bnd.data = {
            "name": name,
        }
        if rad_s2s:
            bnd.data.update(
                {
                    "Radiation": "Diffuse Gray",
                }
            )
        else:
            bnd.data.update(
                {
                    "Radiation": "Idealized",
                    "External Temperature": config["boundaries"]["T_ext"],
                }
            )
    # crucible outside
    bnd = elmer.Boundary(sim, "bnd_crucible_outside", [model["bnd_crucible_outside"].ph_id])
    bnd.data = {
        "name": "bnd_crucible_outside",
        "External Temperature": config["boundaries"]["crucible_outside"]["T_ext"],
        "Heat transfer coefficient": config["boundaries"]["crucible_outside"]["htc"],
    }
    if rad_s2s:
        bnd.data.update(
            {
                "Radiation": "Diffuse Gray",
            }
        )
    else:
        bnd.data.update(
            {
                "Radiation": "Idealized",
                "External Temperature": config["boundaries"]["T_ext"],
            }
        )

    # melt surface, no radiation because melt is transparent
    bnd = elmer.Boundary(sim, "bnd_melt", [model["bnd_melt"].ph_id])
    bnd.data = {
        "name": "bnd_melt",
        "External Temperature": config["boundaries"]["melt"]["T_ext"],
        "Heat transfer coefficient": config["boundaries"]["melt"]["htc"],
        "Save Line Melt Surf": "Logical True"
    }
    
    # crystal surface, no radiation because crystal is transparent
    bnd = elmer.Boundary(sim, "bnd_crystal", [model["bnd_crystal"].ph_id])
    bnd.data = {
        "name": "bnd_crystal",
        "External Temperature": config["boundaries"]["crystal"]["T_ext"],
        "Heat transfer coefficient": config["boundaries"]["crystal"]["htc"],
    }

    # if crucible melt, radiation because melt is transparent
    bnd = elmer.Boundary(sim, "if_crucible_melt", [model["if_crucible_melt"].ph_id])
    bnd.data = {
        "name": "if_crucible_melt",
        "Save Line Crc Melt": "Logical True",
    }
    if rad_s2s:
        bnd.data.update(
            {
                "Radiation": "Diffuse Gray",
                "Radiation Target Body": elmer.StringFromList([sim.bodies["melt"]])
            }
        )
    else:
        bnd.data.update(
            {
                "Radiation": "Idealized",
                "External Temperature": config["boundaries"]["T_ext"],
            }
        )

    # if crystal-seedholder, radiation because crystal is transparent
    bnd = elmer.Boundary(sim, "if_crystal_seedholder", [model["if_crystal_seedholder"].ph_id])
    bnd.data = {
        "name": "if_crystal_seedholder",
    }
    if rad_s2s:
        bnd.data.update(
            {
                "Radiation": "Diffuse Gray",
                "Radiation Target Body": elmer.StringFromList([sim.bodies["crystal"]])
            }
        )
    else:
        bnd.data.update(
            {
                "Radiation": "Idealized",
                "External Temperature": config["boundaries"]["T_ext"],
            }
        )

    # phase interface
    mat_crystal = sim.materials[model["crystal"].params.material]
    if_melt_crystal = elmer.Boundary(
        sim,
        "if_melt_crystal",
        [model["if_melt_crystal"].ph_id],
        {
            "Smart Heater Boundary": "Logical True",
            "Smart Heater Temperature": mat_crystal.data["Melting Point"],
            "Heat flux": v_pull
            / 6e4
            * mat_crystal.data["Density"]
            * mat_crystal.data["Latent Heat"],
            "Heat Flux BC": True,
            "Save Line Melt Crys": "Logical True",
        },
    )

    # inductor
    bnd_inductor = elmer.Boundary(sim, "bnd_inductor", [model["bnd_inductor"].ph_id])
    bnd_inductor.data = {
        "Temperature": config["boundaries"]["T_ext"],
        "name": "bnd_inductor",
    }
    if rad_s2s:
        bnd_inductor.data.update(
            {
                "Radiation": "Diffuse Gray",
            }
        )
    else:
        bnd_inductor.data.update(
            {
                "Radiation": "Idealized",
                "External Temperature": config["boundaries"]["T_ext"],
            }
        )
    bnd_inductor_end_bottom = elmer.Boundary(
        sim, "inductor_end_bottom", [model["bnd_inductor_end_bottom"].ph_id]
    )
    bnd_inductor_end_bottom.data = {
        "AV re {e}": "Real 0.0",
        "AV im {e}": "Real 0.0",
        "AV re": "Real 0.0",
        "AV im": "Real 0.0",
        "name": "inductor_end_bottom",
        "Coil End": True,
    }
    bnd_inductor_end_top = elmer.Boundary(
        sim, "inductor_end_top", [model["bnd_inductor_end_top"].ph_id]
    )
    bnd_inductor_end_top.data = {
        "Coil Start": True,
        "AV re": "Real 0.0",
        "AV im": "Real 0.0",
        "AV re {e}": "Real 0.0",
        "AV im {e}": "Real 0.0",
        "name": "inductor_end_top",
    }
    if_inductor_powersupply = elmer.Boundary(
        sim, "if_inductor_powersupply", [model["if_inductor_vessel_power_supply"].ph_id]
    )
    if_inductor_powersupply.data = {
        "Temperature": config["boundaries"]["T_ext"],
        "name": "if_inductor_powersupply",
    }

    # vessel outside
    bnd_vessel_outside = elmer.Boundary(
        sim, "bnd_vessel_outside", [model["bnd_vessel_outside"].ph_id]
    )
    bnd_vessel_outside.data = {
        "AV re {e}": "Real 0.0",
        "AV im {e}": "Real 0.0",
        "Temperature": config["boundaries"]["T_ext"],
        "name": "bnd_vessel_outside",
    }
    bnd_powersupply_outside = elmer.Boundary(
        sim, "bnd_powersupply_outside", [model["bnd_vessel_powersupply_outside"].ph_id]
    )
    bnd_powersupply_outside.data = {
        "AV re {e}": "Real 0.0",
        "AV im {e}": "Real 0.0",
        "Temperature": config["boundaries"]["T_ext"],
        "name": "bnd_powersupply_outside",
    }

    sim.write_sif(sim_dir)
    print("Starting ElmerGrid...")
    run_elmer_grid(sim_dir, "case.msh")
    print("Starting ElmerSolver...")
    run_elmer_solver(sim_dir)
    err, warn, stats = scan_logfile(sim_dir)
    print("Errors:", err)
    print("Warnings:", warn)
    print("Statistics:", stats)
    


if __name__ == "__main__":
    geo_config = ctrl.load_config("./config_geo.yml")
    sim_config = ctrl.load_config("./config_sim.yml")
    mat_config = ctrl.load_config("./config_mat.yml")

    # # option 1:
    # geo_config["crystal"]["current_length"] = 0.01
    # sim_dir = "simdata/crystal_length_10mm_rotated"
    
    # option 2:
    geo_config["crystal"]["current_length"] = 0.087
    sim_dir = "simdata/crystal_length_87mm_rotated"

    model = geometry(geo_config, sim_dir)
    # simulation(sim_config, model, sim_dir)

    # simulation(sim_config, None, "simdata/crystal_length_87mm")
