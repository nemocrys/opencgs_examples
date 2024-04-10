import numpy as np
import os
import pickle
import yaml

import pyelmer.elmer as elmer
from pyelmer.execute import run_elmer_grid, run_elmer_solver
from pyelmer.post import scan_logfile
from objectgmsh import (
    Model,
    Shape,
    MeshControlLinear,
    MeshControlExponential,
    cut,
    factory,
)
import opencgs.control as ctrl
import opencgs.geo as geo
import gmsh
from copy import deepcopy

# # aluminum crucible
# sim_name = "aluminum_crucible"
# crucible_material = "aluminum"
# frequency = 16.6e3
# crystal_r = 2.45e-3
# crystal_l = 0.03
# interface = [  # interface shape from 2D iteratively coupled simulation
#     [0.00244998418487, 0.0324939282956],
#     [0.0023735401471, 0.0324851040641],
#     [0.00229707632746, 0.0324762775478],
#     [0.00229706105514, 0.0324762755932],
#     [0.00222082527826, 0.0324655591951],
#     [0.00214462074783, 0.03245484719],
#     [0.00214460550671, 0.0324548449955],
#     [0.0020684413141, 0.0324436186215],
#     [0.00199230243481, 0.0324323959783],
#     [0.00199228719996, 0.032432393737],
#     [0.0019161271738, 0.03242121068],
#     [0.00183996529591, 0.0324100273503],
#     [0.00183995005485, 0.0324100251501],
#     [0.00176375393857, 0.0323992140592],
#     [0.00168751983042, 0.0323883975772],
#     [0.00168750457604, 0.0323883954718],
#     [0.00161124681338, 0.0323781652127],
#     [0.00153491695411, 0.0323679252813],
#     [0.0015349016825, 0.032367923306],
#     [0.00145856404326, 0.0323584160869],
#     [0.00138213039282, 0.0323488969105],
#     [0.0013821151022, 0.0323488950895],
#     [0.00130568562301, 0.0323402097627],
#     [0.00122915312211, 0.0323315127291],
#     [0.00122913781231, 0.0323315110801],
#     [0.00115261070644, 0.0323237212073],
#     [0.0010759875293, 0.0323159215563],
#     [0.00107597220096, 0.0323159200921],
#     [0.00099934648165, 0.032309080705],
#     [0.00092264339822, 0.032302234414],
#     [0.000922628052819, 0.0323022331445],
#     [0.000845910599054, 0.0322963868087],
#     [0.000769137769757, 0.0322905362547],
#     [0.000769122409218, 0.0322905351875],
#     [0.000692323485353, 0.0322857153611],
#     [0.000615490832107, 0.0322808934197],
#     [0.000615475458766, 0.0322808925602],
#     [0.00053860992817, 0.0322771212845],
#     [0.000461727632084, 0.0322733491894],
#     [0.000461712248522, 0.0322733485413],
#     [0.000384796039663, 0.0322706410247],
#     [0.000307874450776, 0.0322679333193],
#     [0.000307859059783, 0.0322679328858],
#     [0.000230909892027, 0.032266306334],
#     [0.000153960591604, 0.0322646797818],
#     [0.000153945196105, 0.032264679565],
#     [7.69801286035e-05, 0.0322641385951],
#     [0, 0.0322635976308],
# ]

# graphite crucible
sim_name = "graphite_crucible"
crucible_material = "graphite-CZ3R6300"
frequency = 15.5e3
crystal_r = 6.2e-3
crystal_l = 0.03
interface = [  # interface shape from 2D iteratively coupled simulation
    [0.00619986295268, 0.0330372050182],
    [0.00600587238604, 0.0330141287808],
    [0.00581188244056, 0.0329910526164],
    [0.00581184375009, 0.0329910471691],
    [0.00561903430081, 0.032959665271],
    [0.00542621848314, 0.0329282823353],
    [0.00542617995667, 0.032928275789],
    [0.00523383167885, 0.032894213036],
    [0.00504148250764, 0.0328601501235],
    [0.00504144404024, 0.0328601432437],
    [0.00484922119852, 0.0328254266732],
    [0.00465699823557, 0.0327907100787],
    [0.00465695977132, 0.0327907031892],
    [0.00446464248074, 0.0327565431123],
    [0.00427232534073, 0.0327223830609],
    [0.00427228684436, 0.0327223763593],
    [0.00407974258107, 0.0326895388327],
    [0.00388719858026, 0.0326567013496],
    [0.00388716003, 0.0326566949687],
    [0.00369430396046, 0.032625739921],
    [0.00350144818565, 0.0325947849189],
    [0.00350140956831, 0.032594778958],
    [0.00330819464872, 0.032566140931],
    [0.00311498001505, 0.0325375029449],
    [0.0031149413238, 0.0325374974825],
    [0.00292134636969, 0.032511526272],
    [0.00272775165384, 0.0324855550938],
    [0.00272771288669, 0.0324855501931],
    [0.00253373918887, 0.0324625279644],
    [0.00233976565199, 0.0324395057559],
    [0.00233972681121, 0.032439501468],
    [0.00214539583684, 0.0324196559963],
    [0.00195106495018, 0.0323998105355],
    [0.00195102604151, 0.0323998069011],
    [0.001756373705, 0.0323833178571],
    [0.00156172139983, 0.0323668288189],
    [0.00156168243169, 0.0323668258717],
    [0.00136675694433, 0.0323538502208],
    [0.00117183143358, 0.0323408745724],
    [0.00117179241672, 0.0323408723392],
    [0.000976653040591, 0.0323315214674],
    [0.000781513636812, 0.0323221705993],
    [0.000781474583797, 0.0323221690993],
    [0.000586187688012, 0.0323165235566],
    [0.00039090076836, 0.0323108780183],
    [0.000390861693088, 0.0323108772656],
    [0.000195499792852, 0.032308997814],
    [0, 0.0323071183624],
]


def geometry(
    config,
    sim_dir="./simdata",
    name="cz_induction",
    visualize=False,
    mat_config_file="config_mat.yml",
):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    config_mat = ctrl.load_config(mat_config_file)
    matdata_melt = config_mat[config["melt"]["material"]]

    # initialize geometry model
    model = Model(name)
    gmsh.option.setNumber("General.NumThreads", 8)  # omp parallelized meshing
    dim = 3  # TODO parameterize or remove

    # geometry
    # helper shape to cut inductor ends at outside of vessel
    surrounding_temp = geo.surrounding(
        model, dim, **config["surrounding_temp"], name="surrounding_temp"
    )
    cutbox = geo.surrounding(model, dim, **config["cutbox"], name="cutbox")
    model.remove_shape(surrounding_temp)

    phase_if_crystal = geo.line_from_points(model, interface, "phase_if_crystal")
    phase_if_melt = geo.line_from_points(model, interface, "phase_if_melt")

    crucible = geo.crucible(model, dim, **config["crucible"], material=crucible_material)
    melt = geo.melt(
        model,
        3,
        crucible,
        **config["melt"],
        crystal_radius=crystal_r,
        phase_if=phase_if_melt,
        beta=matdata_melt["Beta"],
        gamma=matdata_melt["Surface Tension"],
        rho=matdata_melt["Density"],
    )
    crystal = geo.crystal(
        model,
        dim,
        **config["crystal"],
        melt=melt,
        phase_if=phase_if_crystal,
        r=crystal_r,
        l=crystal_l,
    )
    seed = geo.seed(model, dim, **config["seed"], crystal=crystal)
    ins = geo.crucible_support(
        model, dim, **config["insulation"], top_shape=crucible, name="insulation"
    )
    adp = geo.crucible_adapter(model, dim, **config["crucible_adapter"], top_shape=ins)
    ax_bt = geo.crucible_support(
        model, dim, **config["axis_bt"], top_shape=adp, name="axis_bt"
    )
    # inductor = geo.inductor_parametric(model, dim, **config["inductor_parametric"])
    inductor = geo.coil_from_points(
        model, dim, dy=ax_bt.params.X0[1], **config["inductor_measurement_calliper"]
    )
    cut(inductor.dimtags, cutbox.dimtags, remove_tool=False)
    inductor.set_interface(cutbox)
    vessel = geo.vessel(model, dim, **config["vessel"], adjacent_shapes=[ax_bt])
    cut(vessel.dimtags, inductor.dimtags, False)
    inductor.set_interface(vessel)
    # lower power supply
    # for parametric coil  # TODO parameterize geometry, use proper flange position & dimension
    # circle = factory.add_circle(
    #     0, -0.015, 0.16, r=inductor.params.radius_pipe * 2
    # )
    circle = factory.add_circle(
        0, inductor.params.supply_bot_h, 0.16, r=inductor.params.radius_pipe * 2
    )
    pipe = factory.extrude([(1, circle)], 0, 0, -0.05)[1][1]
    tag = factory.fragment(vessel.dimtags, [(2, pipe)])[0][1][1]
    vessel_power_supply_bottom = Shape(model, 3, "vessel_power_supply_bottom", [tag])
    vessel_power_supply_bottom.mesh_size = vessel.mesh_size
    vessel_power_supply_bottom.params.material = "insulating-steel"
    # upper power supply
    # for parametric coil  # TODO parameterize geometry, use proper flange position & dimension
    # circle = factory.add_circle(
    #     0, 0.035, 0.16, r=inductor.params.radius_pipe * 2
    # )
    circle = factory.add_circle(
        0, inductor.params.supply_top_h, 0.16, r=inductor.params.radius_pipe * 2
    )
    pipe = factory.extrude([(1, circle)], 0, 0, -0.05)[1][1]
    dimtags = factory.fragment(vessel.dimtags, [(2, pipe)])
    vessel.geo_ids = [dimtags[0][1][1]]
    vessel_power_supply_top = Shape(
        model, 3, "vessel_power_supply_top", [dimtags[0][0][1]]
    )
    vessel_power_supply_top.mesh_size = vessel.mesh_size
    vessel_power_supply_top.params.material = "insulating-steel"
    factory.synchronize()
    ax_top = geo.axis_top(model, dim, **config["axis_top"], bot_shape=seed, vessel=vessel)
    filling = geo.filling(model, dim, **config["filling"], vessel=vessel)

    # helper boundaries
    if_vessel_powersupply = Shape(
        model,
        dim - 1,
        "if_vessel_powersupply",
        vessel.get_interface(vessel_power_supply_bottom)
        + vessel.get_interface(vessel_power_supply_top),
    )

    # surfaces for radiation / convective cooling
    bnd_crystal = Shape(model, dim - 1, "bnd_crystal", crystal.get_interface(filling))
    bnd_melt = Shape(model, dim - 1, "bnd_melt", melt.get_interface(filling))
    bnd_crucible = Shape(
        model, dim - 1, "bnd_crucible", crucible.get_interface(filling)
    )
    bnd_crucible_bot = Shape(
        model,
        dim - 1,
        "bnd_crucible_bot",
        crucible.get_boundaries_in_box(
            [-ins.params.r_in, ins.params.r_in],
            [-crucible.params.t_bt] * 2,
            [-ins.params.r_in, ins.params.r_in],
        ),
    )
    bnd_crucible -= bnd_crucible_bot
    bnd_inductor = Shape(
        model, dim - 1, "bnd_inductor", inductor.get_interface(filling)
    )
    bnd_inductor_end_bottom = Shape(
        model, dim - 1, "bnd_inductor_end_bottom", inductor.get_interface(cutbox)[0:2]
    )
    bnd_inductor_end_top = Shape(
        model, dim - 1, "bnd_inductor_end_top", inductor.get_interface(cutbox)[2:]
    )
    model.remove_shape(cutbox)  # not required any longer
    bnd_ins = Shape(model, dim - 1, "bnd_ins", ins.get_interface(filling))
    bnd_adp = Shape(model, dim - 1, "bnd_adp", adp.get_interface(filling))
    bnd_axbt = Shape(model, dim - 1, "bnd_axbt", ax_bt.get_interface(filling))
    bnd_axtop = Shape(model, dim - 1, "bnd_axtop", ax_top.get_interface(filling))
    bnd_seed = Shape(model, dim - 1, "bnd_seed", seed.get_interface(filling))
    bnd_vessel_inside = Shape(
        model, dim - 1, "bnd_vessel_inside", vessel.get_interface(filling)
    )

    # interfaces
    if_crucible_melt = Shape(
        model, dim - 1, "if_crucible_melt", crucible.get_interface(melt)
    )
    if_melt_crystal = Shape(
        model, dim - 1, "if_melt_crystal", melt.get_interface(crystal)
    )
    if_crystal_seed = Shape(
        model, dim - 1, "if_crystal_seed", crystal.get_interface(seed)
    )
    if_seed_axtop = Shape(model, dim - 1, "if_seed_axtop", seed.get_interface(ax_top))
    if_axtop_vessel = Shape(
        model, dim - 1, "if_axtop_vessel", ax_top.get_interface(vessel)
    )
    if_crucible_ins = Shape(
        model, dim - 1, "if_crucible_ins", crucible.get_interface(ins)
    )
    if_ins_adp = Shape(model, dim - 1, "if_ins_adp", ins.get_interface(adp))
    if_adp_axbt = Shape(model, dim - 1, "if_adp_axbt", adp.get_interface(ax_bt))
    if_axbt_vessel = Shape(
        model, dim - 1, "if_axbt_vessel", ax_bt.get_interface(vessel)
    )
    if_inductor_vessel_power_supply = Shape(
        model,
        dim - 1,
        "if_inductor_vessel_power_supply",
        inductor.get_interface(vessel_power_supply_bottom)
        + inductor.get_interface(vessel_power_supply_top),
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
    # mesh
    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()

    # model.generate_mesh(2)
    # model.show()
    # exit()

    # melt
    MeshControlLinear(
        model,
        bnd_melt,
        melt.mesh_size / 3,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    MeshControlExponential(
        model, bnd_melt, melt.mesh_size / 3, exp=1.6, shapes=[melt, crucible]
    )
    # MeshControlLinear(model, if_crucible_melt, melt.mesh_size / 3, filling.mesh_size, dist_end=filling.mesh_size*2)
    MeshControlExponential(
        model, if_crucible_melt, crystal.mesh_size / 2, exp=1.6, shapes=[crucible, melt]
    )
    # crystal
    MeshControlLinear(
        model,
        bnd_crystal,
        crystal.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 3,
    )
    MeshControlLinear(
        model,
        if_melt_crystal,
        crystal.mesh_size / 2,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    MeshControlExponential(
        model, if_melt_crystal, crystal.mesh_size / 2, exp=1.6, shapes=[crystal, melt]
    )
    # seed
    MeshControlLinear(
        model,
        bnd_seed,
        seed.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 5,
    )
    MeshControlLinear(
        model,
        if_crystal_seed,
        seed.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 5,
    )
    MeshControlLinear(
        model,
        if_seed_axtop,
        seed.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 5,
    )
    # MeshControlExponential(model, seed, seed.mesh_size, exp=1.6)
    # MeshControlExponential(model, if_crystal_seed, seed.mesh_size)
    # MeshControlExponential(model, if_seed_axtop, seed.mesh_size)
    # top axis
    MeshControlLinear(
        model,
        bnd_axtop,
        ax_top.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    MeshControlLinear(
        model,
        if_axtop_vessel,
        ax_top.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    # crucible
    MeshControlLinear(
        model,
        bnd_crucible,
        crucible.mesh_size / 3,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    MeshControlExponential(
        model,
        bnd_crucible,
        crucible.mesh_size / 3,
        exp=1.6,
        shapes=[crucible, melt, ins, adp],
    )
    MeshControlLinear(
        model,
        if_crucible_ins,
        crucible.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    # insulation
    MeshControlLinear(
        model, bnd_ins, ins.mesh_size, filling.mesh_size, dist_end=filling.mesh_size * 2
    )
    MeshControlLinear(
        model,
        if_ins_adp,
        ins.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    # adapter
    MeshControlLinear(
        model, bnd_adp, adp.mesh_size, filling.mesh_size, dist_end=filling.mesh_size * 2
    )
    MeshControlLinear(
        model,
        if_adp_axbt,
        adp.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    # bottom axis
    MeshControlLinear(
        model,
        bnd_axbt,
        ax_bt.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    MeshControlLinear(
        model,
        if_axbt_vessel,
        ax_bt.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    # vessel
    MeshControlLinear(
        model,
        bnd_vessel_inside,
        vessel.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
    )
    # inductor
    MeshControlLinear(
        model,
        bnd_inductor,
        inductor.mesh_size,
        filling.mesh_size,
        dist_end=filling.mesh_size * 2,
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


def simulation(config, model=None, sim_dir="./simdata", name="cz_induction"):
    if model is None:
        with open(f"{sim_dir}/model.pickle", "rb") as f:
            model = pickle.load(f)

    # read parameters from config
    omega = 2 * np.pi * frequency
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
    # add bodies and materials
    material = elmer.load_material(
        model["inductor"].params.material, sim, "config_mat.yml"
    )
    material.data.update({"name": model["inductor"].params.material})
    inductor = elmer.Body(sim, "inductor", [model["inductor"].ph_id])
    inductor.material = material
    inductor.equation = eqn_coil
    inductor.body_force = current_source
    inductor.data.update({"name": "inductor"})

    for name in [
        "crucible",
        "crystal",
        "seed",
        "insulation",
        "crucible_adapter",
        "axis_bt",
        "vessel",
        "axis_top",
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
        body = elmer.Body(sim, name, [model[name].ph_id])
        body.material = material
        body.equation = eqn_main
        body.body_force = joule_heat
        body.data.update({"name": name})

    melt = elmer.Body(sim, "melt", [model["melt"].ph_id])
    material = elmer.load_material(
        model["melt"].params.material, sim, "config_mat.yml"
    )
    material.data.update({"name": model["melt"].params.material})
    melt.material = material
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

    # get rid of parameters not used by elmer
    sim.materials["tin-liquid"].data.pop("Surface Tension")
    sim.materials["tin-liquid"].data.pop("Beta")

    for name in [
        "bnd_ins",
        "bnd_adp",
        "bnd_axbt",
        "bnd_axtop",
        "bnd_seed",
        "bnd_vessel_inside",
        "bnd_vessel_powersupply_inside",
        "bnd_crucible_bot",
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

    for name in [
        "bnd_crystal",
        "bnd_melt",
        "bnd_crucible",
    ]:
        bnd = elmer.Boundary(sim, name, [model[name].ph_id])
        bnd.data = {
            "name": name,
            "External Temperature": config["boundaries"][name]["T_ext"],
            "Heat transfer coefficient": config["boundaries"][name]["htc"],
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
                }
            )
    sim.boundaries["bnd_melt"].data.update({"Save Line Melt Surf": "Logical True"})

    mat_crystal = sim.materials[model["crystal"].params.material]
    if_melt_crystal = elmer.Boundary(
        sim,
        "if_melt_crystal",
        [model["if_melt_crystal"].ph_id],
        {
            "Smart Heater Boundary": "Logical True",
            "Smart Heater Temperature": mat_crystal.data["Melting Point"],
            "Heat flux": config["general"]["v_pull"]
            / 6e4
            * mat_crystal.data["Density"]
            * mat_crystal.data["Latent Heat"],
            "Heat Flux BC": True,
            "Save Line Melt Crys": "Logical True",
        },
    )
    if_crucible_melt = elmer.Boundary(
        sim,
        "if_crucible_melt",
        [model["if_crucible_melt"].ph_id],
        {"Save Line Crc Melt": "Logical True"},
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
    if_inductor_powersupply = elmer.Boundary(
        sim, "if_inductor_powersupply", [model["if_inductor_vessel_power_supply"].ph_id]
    )
    if_inductor_powersupply.data = {
        "Temperature": config["boundaries"]["T_ext"],
        "name": "if_inductor_powersupply",
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
    with open("config_geo.yml") as f:
        config_geo = yaml.safe_load(f)
    with open("config_sim.yml") as f:
        config_sim = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir=f"simdata_{sim_name}")
    simulation(config_sim, model, sim_dir=f"simdata_{sim_name}")
    # simulation(config_sim, None, sim_dir=f"simdata_{sim_name}")
