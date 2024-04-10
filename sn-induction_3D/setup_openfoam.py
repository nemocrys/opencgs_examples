import shutil
import matplotlib.pyplot as plt
import numpy as np
from nemoblock import *
from pyelmer.post import dat_to_dataframe
from extract_lorentz_force import lorentz_force_to_txt
import yaml


def prepare_simulation(
    dimension=2,
    sim_dir="./simdata_openfoam",
    elmer_results="./simdata_elmer/results",
    of_template="./openfoam_template",
):
    # read elmer output
    try:
        df_crc_melt = dat_to_dataframe(elmer_results + "/save_line_crc_melt_relaxed.dat")
    except FileNotFoundError:
        df_crc_melt = dat_to_dataframe(elmer_results + "/save_line_crc_melt.dat")
    df_crc_melt_finalstep = df_crc_melt.loc[
        df_crc_melt["Call count"] == df_crc_melt["Call count"][df_crc_melt.index[-1]]
    ].drop_duplicates(subset=["Node index"])
    try:
        df_melt_surf = dat_to_dataframe(elmer_results + "/save_line_melt_surf_relaxed.dat")
    except FileNotFoundError:
        df_melt_surf = dat_to_dataframe(elmer_results + "/save_line_melt_surf.dat")
    df_melt_surf_finalstep = df_melt_surf.loc[
        df_melt_surf["Call count"] == df_melt_surf["Call count"][df_melt_surf.index[-1]]
    ].drop_duplicates(subset=["Node index"])
    try:
        df_melt_crys = dat_to_dataframe(elmer_results + "/save_line_melt_crys_relaxed.dat")
    except FileNotFoundError:
        df_melt_crys = dat_to_dataframe(elmer_results + "/save_line_melt_crys.dat")
    df_melt_crys_finalstep = df_melt_crys.loc[
        df_melt_crys["Call count"] == df_melt_crys["Call count"][df_melt_crys.index[-1]]
    ].drop_duplicates(subset=["Node index"])

    coords_crc_melt = np.array(
        [
            df_crc_melt_finalstep["coordinate 1"].to_numpy(),
            df_crc_melt_finalstep["coordinate 2"].to_numpy(),
        ]
    )
    coords_melt_surf = np.array(
        [
            np.round(df_melt_surf_finalstep["coordinate 1"].to_numpy(), 100),
            np.round(df_melt_surf_finalstep["coordinate 2"].to_numpy(), 100),
        ]
    )
    coords_melt_crys = np.array(
        [
            np.round(df_melt_crys_finalstep["coordinate 1"].to_numpy(), 100),
            np.round(df_melt_crys_finalstep["coordinate 2"].to_numpy(), 100),
        ]
    )

    # BEGIN DIRTYFIX
    # this is required because Elmer sometimes slightly shifts coordinates

    # 1. neighboring boundaries must share one point
    i_surf = np.argmin(coords_melt_surf[0, :])
    i_crys = np.argmax(coords_melt_crys[0, :])
    # x_new = (coords_melt_crys[0, i_crys] + coords_melt_surf[0, i_surf]) * 0.5
    # y_new = (coords_melt_crys[1, i_crys] + coords_melt_surf[1, i_surf]) * 0.5
    x_new = coords_melt_crys[0, i_crys]
    y_new = coords_melt_crys[1, i_crys]
    coords_melt_surf[0, i_surf] = x_new
    coords_melt_surf[1, i_surf] = y_new
    coords_melt_crys[0, i_crys] = x_new
    coords_melt_crys[1, i_crys] = y_new

    i_surf = np.argmax(coords_melt_surf[0, :])
    i_cruc = np.argmax(coords_crc_melt[1, :])
    # x_new = (coords_crc_melt[0, i_cruc] + coords_melt_surf[0, i_surf]) * 0.5
    # y_new = (coords_crc_melt[1, i_cruc] + coords_melt_surf[1, i_surf]) * 0.5
    x_new = coords_crc_melt[0, i_cruc]
    y_new = coords_crc_melt[1, i_cruc]
    coords_melt_surf[0, i_surf] = x_new
    coords_melt_surf[1, i_surf] = y_new
    coords_crc_melt[0, i_cruc] = x_new
    coords_crc_melt[1, i_cruc] = y_new


    i = np.argmin(coords_crc_melt[0, :])
    if coords_crc_melt[0, i] < 1e-6:
        coords_crc_melt[0, i] = 0
    else:
        raise ValueError("Check the coordinates!")
    i = np.argmin(coords_melt_crys[0, :])
    if coords_melt_crys[0, i] < 1e-6:
        coords_melt_crys[0, i] = 0
    else:
        raise ValueError("Check the coordinates!")
    # END DIRTYFIX

    # fig, ax = plt.subplots()
    # ax.plot(coords_crc_melt[0, :], coords_crc_melt[1, :], "x-")
    # ax.plot(coords_melt_surf[0, :], coords_melt_surf[1, :], "x-")
    # ax.plot(coords_melt_crys[0, :], coords_melt_crys[1, :], "x-")
    # ax.set_aspect("equal")
    # plt.show()

    print("Mean node temperatures on boundary")
    print("Crucible:", df_crc_melt_finalstep["temperature"].mean())
    print("Surface:", df_melt_surf_finalstep["temperature"].mean())
    print("Crystal:", df_melt_crys_finalstep["temperature"].mean())

    # extract main dimensions
    r_crucible = coords_crc_melt[0, :].max()
    r_crystal = coords_melt_crys[0, :].max()
    h_melt = coords_crc_melt[1, :].max() - coords_crc_melt[1, :].min()
    z_min = coords_crc_melt[1, :].min()

    ########################################################################
    # input 2

    # Circumferential mesh size
    res_phi = 60

    # boundary layer
    smallest_element_crucible = 0.00009
    layer_thickness_crucible = 0.0005
    growth_rate_crucible = 1.4

    smallest_element_surface = 0.00009
    layer_thickness_surface = 0.0005
    growth_rate_surface = 1.4

    smallest_element_meniscus = 0.00009
    layer_thickness_meniscus = 0.0005
    growth_rate_meniscus = 1.5

    smallest_element_crystal = smallest_element_crucible
    layer_thickness_crystal = 0.0005
    growth_rate_crystal = 1.5

    if dimension == 3:
        # # medium resolution, approx. 120k cells
        # smallest_element_crucible *= 2
        # layer_thickness_crucible *= 4
        # smallest_element_surface *= 2
        # layer_thickness_surface *= 4
        # smallest_element_meniscus *= 2
        # layer_thickness_meniscus *= 4
        # smallest_element_crystal = 0.00012
        # growth_rate_crystal = 1.1
        # # layer_thickness_crystal *= 2

        # high resolution, approx. 870k cells
        res_phi = 120
        smallest_element_crucible *= 1
        layer_thickness_crucible *= 2.5
        smallest_element_surface *= 1
        layer_thickness_surface *= 2.5
        smallest_element_meniscus *= 1
        layer_thickness_meniscus *= 2.5
        smallest_element_crystal = 0.00008
        growth_rate_crystal = 1.1
        # layer_thickness_crystal *= 2
    ########################################################################
    # nemoblock meshing
    mesh = Mesh()

    ####################
    # Mesh structure  (left: symmetry axis)
    # fmt: off
    #
    # ^ z-direction
    # |
    # |
    # |------->
    #   r-direction
    #  __ __ __
    # |  |  |  |
    # .c1|r1|r2|
    # |__|__|__|
    # |  |  |  |
    # .c2|r3|r4|
    # |__|__|__|
    #
    # 2D coord. system:
    # ^ y-direction
    # |
    # |
    # |------->
    #   x-direction
    #
    # fmt: on
    ####################

    ####################
    # Surfaces defined by splines
    s_fs = spline(coords_melt_surf.T, "linear")  # free surface
    s_if = spline(coords_melt_crys.T, "linear")  # phase interface

    # # fig, ax = plot_spline(s_bt, [0, r_crucible])
    # fig, ax = plot_spline(s_fs, [r_crystal, r_crucible])
    # plot_spline(s_if, [0, r_crystal], fig, ax)
    # ax.axis('equal')
    # plt.show()
    ####################
    ####################

    # Geometry
    # cylinder c1
    # c1_r_top = r_crystal
    c1_z_top = z_min + h_melt / 2 + h_melt * 1 / 20
    # c1_r_bt = r_crystal
    # ring r1
    r1_r_top = r_crystal + (r_crucible - r_crystal) * 0.4  # adjust here!
    # r1_z_top = c1_z_top + h_melt / 8
    # r1_r_bt = r_crucible * 0.75
    # # ring r2
    # r2_r_top = r_crystal + (r_crucible - r_crystal) * 0.5

    ####################
    # Mesh sizes
    print("\nGrading bottom")
    res_z_bot, grading_bot = boundary_layer(
        h_melt / 2,
        "xmin",
        smallest_element_crucible,
        layer_thickness_crucible,
        growth_rate_crucible,
    )
    print("\nGrading top")
    res_z_top, grading_top = boundary_layer(
        h_melt / 2,
        "xmax",
        smallest_element_surface,
        layer_thickness_surface,
        growth_rate_surface,
    )
    print("\nGrading outside")
    res_r_out, grading_r_out = boundary_layer(
        r_crucible - r1_r_top,
        "xmax",
        smallest_element_crucible,
        layer_thickness_crucible,
        growth_rate_crucible,
    )
    print("\nGrading meniscus")
    res_r_men, grading_r_men = boundary_layer(
        r1_r_top - r_crystal,
        "xmin",
        smallest_element_meniscus,
        layer_thickness_meniscus,
        growth_rate_meniscus,
    )
    print("\nGrading crystal")
    if dimension == 3:
        res_r_crys, grading_r_crys = boundary_layer(
            r_crystal / 2,
            "xmax",
            smallest_element_crystal,
            layer_thickness_crystal,
            growth_rate_crystal,
        )
    else:
        res_r_crys, grading_r_crys = boundary_layer(
            r_crystal,
            "xmax",
            smallest_element_crystal,
            layer_thickness_crystal,
            growth_rate_crystal,
        )

    # res_z_bot = 10
    # res_z_top = 10
    # res_r_out = 10
    # res_r_men = 10
    # res_r_crys = 10

    ####################
    # Blocks (defined as cylinders & rings)

    def coords2D(x):  # to transform coordinate system
        return [x[0], x[2], -x[1]]

    def evaluate_spline(r_start, r_end, phi, spline, res=100):
        return [
            coords2D(cartesian(r, phi, spline(r)))
            for r in np.linspace(r_start, r_end, res, endpoint=False)
        ]

    phi1 = -2.5
    phi2 = 2.5
    if dimension == 3:
        c1 = create_cylinder(
            mesh,
            [r_crystal, s_fs(r_crystal)],
            [r_crystal, h_melt / 2],
            res_r_crys,
            res_phi,
            res_z_top,
        )
        c1.set_spline_surface(s_if, "top")
    else:
        c1 = Block(
            mesh,
            coords2D(cartesian(0, phi1, h_melt / 2)),
            coords2D(cartesian(r_crystal, phi1, h_melt / 2)),
            coords2D(cartesian(r_crystal, phi2, h_melt / 2)),
            coords2D(cartesian(0, phi2, h_melt / 2)),
            coords2D(cartesian(0, phi1, s_if(0))),
            coords2D(cartesian(r_crystal, phi1, s_if(r_crystal))),
            coords2D(cartesian(r_crystal, phi2, s_if(r_crystal))),
            coords2D(cartesian(0, phi2, s_if(0))),
        )
        c1.p3 = "p0"
        c1.p7 = "p4"
        c1.set_number_of_cells(res_r_crys, 1, res_z_top)
        c1.create()

        c1.e3.type = "spline"
        c1.e3.points = evaluate_spline(0, r_crystal, phi1, s_if)
        c1.e2.type = "spline"
        c1.e2.points = evaluate_spline(0, r_crystal, phi2, s_if)

    if dimension == 3:
        r1 = create_ring(
            mesh,
            [r1_r_top, s_fs(r1_r_top)],
            [r1_r_top, h_melt / 2],
            c1.surf_rad,
            res_r_men,
            res_phi,
            res_z_top,
        )
        r1.set_spline_surface(s_fs, "top", res=1000)
    else:
        r1 = Block(mesh)
        r1.set_connection(c1, "left")
        r1.p1 = coords2D(cartesian(r1_r_top, phi1, h_melt / 2))
        r1.p2 = coords2D(cartesian(r1_r_top, phi2, h_melt / 2))
        r1.p5 = coords2D(cartesian(r1_r_top, phi1, s_fs(r1_r_top)))
        r1.p6 = coords2D(cartesian(r1_r_top, phi2, s_fs(r1_r_top)))
        r1.cells_x1 = res_r_men
        r1.create()
        r1.e3.type = "spline"
        r1.e3.points = evaluate_spline(r_crystal, r1_r_top, phi1, s_fs)
        r1.e2.type = "spline"
        r1.e2.points = evaluate_spline(r_crystal, r1_r_top, phi2, s_fs)
    if dimension == 3:
        r2 = create_ring(
            mesh,
            [r_crucible, h_melt],
            [r_crucible, h_melt / 2],
            r1.surf_rad,
            res_r_out,
            res_phi,
            res_z_top,
        )
        r2.set_spline_surface(s_fs, "top", res=1000)
    else:
        r2 = Block(mesh)
        r2.set_connection(r1, "left")
        r2.p1 = coords2D(cartesian(r_crucible, phi1, h_melt / 2))
        r2.p2 = coords2D(cartesian(r_crucible, phi2, h_melt / 2))
        r2.p5 = coords2D(cartesian(r_crucible, phi1, s_fs(r_crucible)))
        r2.p6 = coords2D(cartesian(r_crucible, phi2, s_fs(r_crucible)))
        r2.cells_x1 = res_r_out
        r2.create()
        r2.e3.type = "spline"
        r2.e3.points = evaluate_spline(r1_r_top, r_crucible, phi1, s_fs)
        r2.e2.type = "spline"
        r2.e2.points = evaluate_spline(r1_r_top, r_crucible, phi2, s_fs)
    if dimension == 3:
        c2 = create_cylinder(
            mesh,
            [r_crystal, h_melt / 2],
            [r_crystal, 0],
            res_r_crys,
            res_phi,
            res_z_top,
            cylinder_on_top=c1,
        )
    else:
        c2 = Block(mesh)
        c2.set_connection(c1, "top")
        c2.p0 = coords2D(cartesian(0, phi1, 0))
        c2.p1 = coords2D(cartesian(r_crystal, phi1, 0))
        c2.p2 = coords2D(cartesian(r_crystal, phi2, 0))
        c2.p3 = "p0"
        c2.cells_x3 = res_z_bot
        c2.create()
    if dimension == 3:
        r3 = create_ring(
            mesh,
            [r1_r_top, h_melt / 2],
            [r1_r_top, 0],
            c2.surf_rad,
            res_r_men,
            res_phi,
            res_z_top,
            ring_on_top=r1,
        )
    else:
        r3 = Block(mesh)
        r3.set_connection(c2, "left")
        r3.set_connection(r1, "top")
        r3.p1 = coords2D(cartesian(r1_r_top, phi1, 0))
        r3.p2 = coords2D(cartesian(r1_r_top, phi2, 0))
        r3.create()
    if dimension == 3:
        r4 = create_ring(
            mesh,
            [r_crucible, h_melt / 2],
            [r_crucible, 0],
            r3.surf_rad,
            res_r_out,
            res_phi,
            res_z_top,
            ring_on_top=r2,
        )
    else:
        r4 = Block(mesh)
        r4.set_connection(r3, "left")
        r4.set_connection(r2, "top")
        r4.p1 = coords2D(cartesian(r_crucible, phi1, 0))
        r4.p2 = coords2D(cartesian(r_crucible, phi2, 0))
        r4.create()

    # ####################
    # # Grading
    if dimension == 3:
        c1.set_grading_axial(grading_top)
        r1.set_grading_axial(grading_top)
        r2.set_grading_axial(grading_top)
        c2.set_grading_axial(grading_bot)
        r3.set_grading_axial(grading_bot)
        r4.set_grading_axial(grading_bot)
        c1.set_grading_radial(grading_r_crys)
        c2.set_grading_radial(grading_r_crys)
        r1.set_grading_radial(grading_r_men)
        r3.set_grading_radial(grading_r_men)
        r2.set_grading_radial(grading_r_out)
        r4.set_grading_radial(grading_r_out)
    else:
        c1.grading = f"simpleGrading ({grading_r_crys} 1 {grading_top})"
        r1.grading = f"simpleGrading ({grading_r_men} 1 {grading_top})"
        r2.grading = f"simpleGrading ({grading_r_out} 1 {grading_top})"
        c2.grading = f"simpleGrading ({grading_r_crys} 1 {grading_bot})"
        r3.grading = f"simpleGrading ({grading_r_men} 1 {grading_bot})"
        r4.grading = f"simpleGrading ({grading_r_out} 1 {grading_bot})"

    # ####################
    # # Patches
    crucible_bot = Patch(mesh, "wall crucibleBot")
    if dimension == 3:
        crucible_bot.faces += c2.surf_bt
        crucible_bot.faces += r3.surf_bt
        crucible_bot.faces += r4.surf_bt
    else:
        crucible_bot.add_face(c2.face_bottom)
        crucible_bot.add_face(r3.face_bottom)
        crucible_bot.add_face(r4.face_bottom)

    crucible_side = Patch(mesh, "wall crucibleSide")
    if dimension == 3:
        crucible_side.faces += r2.surf_rad
        crucible_side.faces += r4.surf_rad
    else:
        crucible_side.add_face(r2.face_right)
        crucible_side.add_face(r4.face_right)

    free_surf = Patch(mesh, "wall freeSurf")
    if dimension == 3:
        free_surf.faces += r1.surf_top
        free_surf.faces += r2.surf_top
    else:
        free_surf.add_face(r1.face_top)
        free_surf.add_face(r2.face_top)

    top_surf = Patch(mesh, "wall crysInter")
    if dimension == 3:
        top_surf.faces += c1.surf_top
    else:
        top_surf.add_face(c1.face_top)

    if dimension == 2:
        mesh_front = Patch(mesh, "wedge front")
        mesh_front.add_face(c1.face_front)
        mesh_front.add_face(c2.face_front)
        mesh_front.add_face(r1.face_front)
        mesh_front.add_face(r2.face_front)
        mesh_front.add_face(r3.face_front)
        mesh_front.add_face(r4.face_front)

    if dimension == 2:
        mesh_back = Patch(mesh, "wedge back")
        mesh_back.add_face(c1.face_back)
        mesh_back.add_face(c2.face_back)
        mesh_back.add_face(r1.face_back)
        mesh_back.add_face(r2.face_back)
        mesh_back.add_face(r3.face_back)
        mesh_back.add_face(r4.face_back)

    try:
        shutil.copytree(of_template, sim_dir)
    except FileExistsError:
        print("could not copy OpenFOAM template, files already exists")


    try:
        shutil.copy(
            f"{elmer_results}/save_line_crc_melt_relaxed.dat",
            f"{sim_dir}/save_line_crc_melt_relaxed.dat",
        )
    except FileNotFoundError:
        shutil.copy(
            f"{elmer_results}/save_line_crc_melt.dat",
            f"{sim_dir}/save_line_crc_melt_relaxed.dat",
        )
    try:
        shutil.copy(
            f"{elmer_results}/save_line_melt_crys_relaxed.dat",
            f"{sim_dir}/save_line_melt_crys_relaxed.dat",
        )
    except FileNotFoundError:
        shutil.copy(
            f"{elmer_results}/save_line_melt_crys.dat",
            f"{sim_dir}/save_line_melt_crys_relaxed.dat",
        )
    try:
        shutil.copy(
            f"{elmer_results}/save_line_melt_surf_relaxed.dat",
            f"{sim_dir}/save_line_melt_surf_relaxed.dat",
        )
    except FileNotFoundError:
        shutil.copy(
            f"{elmer_results}/save_line_melt_surf.dat",
            f"{sim_dir}/save_line_melt_surf_relaxed.dat",
        )

    try:
        os.remove(f"{sim_dir}/system/blockMeshDict")
    except FileNotFoundError:
        print("No blockMeshDict to remove.")
    mesh.write(f"{sim_dir}/system/")

    # extract heater power scaling for scaling of Lorentz forces
    with open(f"{elmer_results}/../probes_result.yml") as f:
        power_scaling = yaml.safe_load(f)["res heater power scaling 1"]


    lorentz_force_to_txt(f"{elmer_results}/../case_melt_t0001.vtu", f"{sim_dir}/lorentz-force.csv", scaling=power_scaling)
    lorentz_force_to_txt(f"{elmer_results}/../case_melt_t0001.vtu", f"{sim_dir}/lorentz-force_with-header.csv", header=True, scaling=power_scaling)


if __name__ == "__main__":
    # prepare_simulation()
    prepare_simulation(
        3,
        "simdata_aluminum_of-3d_elmer-2d",
        "simdata/2023-09-28_11-47_eo_aluminum_transient_hf-bc/02_simulation/elmer_03/results",
        "openfoam_template_transient"
    )
