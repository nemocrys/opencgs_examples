import shutil
import matplotlib.pyplot as plt
import numpy as np
from nemoblock import *
from pyelmer.post import dat_to_dataframe
from extract_lorentz_force import lorentz_force_to_txt


def prepare_simulation(
    dimension=2,
    sim_dir="./simdata_openfoam",
    elmer_results="./simdata_elmer/results",
    of_template="./openfoam_template_transient",
    use_hardcoded_grading=False,
    highres=True,
):
        
    # read elmer output

    try:
        df_crc_melt = dat_to_dataframe(elmer_results + "/save_line_crc_melt_relaxed.dat")
    except FileNotFoundError:
        df_crc_melt = dat_to_dataframe(elmer_results + "/save_line_crc_melt.dat")
    df_crc_melt_finalstep = df_crc_melt.loc[
        df_crc_melt["Call count"] == df_crc_melt["Call count"][df_crc_melt.index[-1]]
    ]
    try:
        df_melt_surf = dat_to_dataframe(elmer_results + "/save_line_melt_surf_relaxed.dat")
    except FileNotFoundError:
        df_melt_surf = dat_to_dataframe(elmer_results + "/save_line_melt_surf.dat")
    df_melt_surf_finalstep = df_melt_surf.loc[
        df_melt_surf["Call count"] == df_melt_surf["Call count"][df_melt_surf.index[-1]]
    ]
    try:
        df_melt_crys = dat_to_dataframe(elmer_results + "/save_line_melt_crys_relaxed.dat")
    except FileNotFoundError:
        df_melt_crys = dat_to_dataframe(elmer_results + "/save_line_melt_crys.dat")
    df_melt_crys_finalstep = df_melt_crys.loc[
        df_melt_crys["Call count"] == df_melt_crys["Call count"][df_melt_crys.index[-1]]
    ]

    # DIRTYFIX: sort out points with duplicate x-coordinate to get a proper graph
    # coords_crc_melt = np.array(
    #     [
    #         df_crc_melt_finalstep["coordinate 1"].to_numpy(),
    #         df_crc_melt_finalstep["coordinate 2"].to_numpy(),
    #     ]
    # )
    x = np.round(df_crc_melt_finalstep["coordinate 1"].to_numpy(), 6)
    y = np.round(df_crc_melt_finalstep["coordinate 2"].to_numpy(), 6)
    _, indices = np.unique(x, return_index=True)
    x_cleaned = []
    y_cleaned = []
    for i in indices:
        x_cleaned.append(x[i])
        y_cleaned.append(y[i])
    y_cleaned[np.argmax(np.array(x_cleaned))] = max(y)
    coords_crc_melt = np.array([x_cleaned, y_cleaned])

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
    # this is required because mgdyn solver slightly shifts coordinates

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
    i_cruc = np.argmax(coords_crc_melt[0, :])
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
    # input

    # Circumferential mesh size
    res_phi = 180  # only used in 3D

    # boundary layer
    smallest_element_crucible = 0.00009
    layer_thickness_crucible = 0.005

    smallest_element_top = 0.00009
    layer_thickness_top = 0.005

    smallest_element_meniscus = 0.00025
    layer_thickness_meniscus = 0.005
    
    if highres:  # this is default
        growth_rate_crucible = 1.2
        growth_rate_top = 1.2
        growth_rate_meniscus = 1.2
    else:
        growth_rate_crucible = 1.4
        growth_rate_top = 1.4
        growth_rate_meniscus = 1.5

    smallest_element_crystal_side = smallest_element_meniscus
    layer_thickness_crystal_side = layer_thickness_meniscus
    growth_rate_crystal_side = growth_rate_meniscus

    ########################################################################
    # nemoblock meshing
    mesh = Mesh()

    ####################
    # Mes structure  (left: symmetry axis)
    # fmt: off
    #
    # coordinate system:
    #
    # ^ z-direction
    # |
    # |
    # |------->
    #   r-direction
    # .__
    # |  | \
    # .  |    \
    # |c2| r3 / \
    # .__|___/    \
    # |  |   \ r2 /
    # .c1| r1 \  /
    # |__|_____\/
    #
    # fmt: on

    ####################

    # create splines
    s_bt = spline(coords_crc_melt.T, "linear")  # crucible: bottom of c1, r1, r2  # linear required here together with dirtyfix from above
    s_fs = spline(coords_melt_surf.T)  # free surface: right of r2, top of r3
    s_ph = spline(coords_melt_crys.T)  # melting front: top of c2
    # fig, ax = plot_spline(s_bt, [0, r_crucible])
    # # plot_spline(s_bt, [0, r_crucible], fig, ax)
    # plot_spline(s_fs, [r_crystal, r_crucible], fig, ax)
    # plot_spline(s_ph, [0, r_crystal], fig, ax)
    # ax.axis('equal')
    # plt.show()

    # Geometry
    # cylinder c1
    c1_r_top = r_crystal * 0.75
    c1_z_top = z_min + h_melt / 2 + h_melt * 0.05
    c1_r_bt = r_crystal * 0.75
    # ring r1
    r1_r_top = r_crystal + (r_crucible - r_crystal) * 0.38  # adjust here!
    r1_z_top = c1_z_top + h_melt * 0.05  # adjust here!
    r1_r_bt = r_crucible * 0.79  # adjust here!
    # ring r2 = r3_r_out
    r2_r_top = r_crystal + (r_crucible - r_crystal) * 0.5

    ####################
    # Mesh sizes
    if use_hardcoded_grading:
        res_z_c1, grading_bottom  = 51, '( (0.1278222910381451 14 10.699320537907195) (0.8721777089618549 37 1) )'
        res_z_c2, grading_top = 51, '( (0.8721777089618549 37 1) (0.1278222910381451 14 0.09346387898717928) )'
        res_r_r1, grading_meniscus = 49, '( (0.0643953934740883 4 1.7279999999999998) (0.9356046065259117 45 1) )'
        res_r_c1, grading_crys_rad = 59, '( (0.9471653543307086 55 1) (0.05283464566929135 4 0.5787037037037038) )'
    else:
        res_z_c1, grading_bottom = boundary_layer(
            h_melt / 2,
            "xmin",
            smallest_element_crucible,
            layer_thickness_crucible,
            growth_rate_crucible,
        )
        res_z_c2, grading_top = boundary_layer(
            h_melt / 2, "xmax", smallest_element_top, layer_thickness_top, growth_rate_top
        )
        res_r_r1, grading_meniscus = boundary_layer(
            h_melt / 2,
            "xmin",
            smallest_element_meniscus,
            layer_thickness_meniscus,
            growth_rate_meniscus,
        )
        if dimension == 2:
            radius_c1 = r_crystal
        else:
            radius_c1 = r_crystal * 0.5
        res_r_c1, grading_crys_rad = boundary_layer(
            radius_c1,
            "xmax",
            smallest_element_crystal_side,
            layer_thickness_crystal_side,
            growth_rate_crystal_side,
        )
    # res_z_c1 = 10
    # res_z_c2 = 10
    # res_r_r1 = 10
    # res_r_c1 = 10

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
            [c1_r_top, c1_z_top],
            [c1_r_bt, s_bt(c1_r_bt)],
            res_r_c1,
            res_phi,
            res_z_c1
        )
        c1.set_spline_surface(s_bt, "bottom", 1000)
    else:
        c1 = Block(
            mesh,
            coords2D(cartesian(0, phi1, s_bt(0))),
            coords2D(cartesian(c1_r_bt, phi1, s_bt(c1_r_bt))),
            coords2D(cartesian(c1_r_bt, phi2, s_bt(c1_r_bt))),
            coords2D(cartesian(0, phi2, s_bt(0))),
            coords2D(cartesian(0, phi1, c1_z_top)),
            coords2D(cartesian(c1_r_top, phi1, c1_z_top)),
            coords2D(cartesian(c1_r_top, phi2, c1_z_top)),
            coords2D(cartesian(0, phi2, c1_z_top)),
        )
        c1.p3 = "p0"
        c1.p7 = "p4"
        c1.set_number_of_cells(res_r_c1, 1, res_z_c1)
        c1.create()

        c1.e0.type = "spline"
        c1.e0.points = evaluate_spline(0, c1_r_bt, phi1, s_bt)
        c1.e1.type = "spline"
        c1.e1.points = evaluate_spline(0, c1_r_bt, phi2, s_bt)
    if dimension == 3:
        c2 = create_cylinder(
            mesh,
            [r_crystal, s_ph(r_crystal)],
            [c1_r_top, c1_z_top],
            res_r_c1,
            res_phi,
            res_z_c2,
            cylinder_below=c1,
        )
        c2.set_spline_surface(s_ph, "top", 1000)
    else:
        c2 = Block(mesh)
        c2.set_connection(c1, "bottom")
        c2.p4 = coords2D(cartesian(0, phi1, s_ph(0)))
        c2.p5 = coords2D(cartesian(r_crystal, phi1, s_ph(r_crystal)))
        c2.p6 = coords2D(cartesian(r_crystal, phi2, s_ph(r_crystal)))
        c2.p7 = "p4"
        c2.cells_x3 = res_z_c2
        c2.create()
        c2.e3.type = "spline"
        c2.e3.points = evaluate_spline(0, r_crystal, phi1, s_ph)
        c2.e2.type = "spline"
        c2.e2.points = evaluate_spline(0, r_crystal, phi2, s_ph)

    if dimension == 3:
        r1 = create_ring(
            mesh,
            [r1_r_top, r1_z_top],
            [r1_r_bt, s_bt(r1_r_bt)],
            c1.surf_rad,
            res_r_r1,
            res_phi,
            res_z_c1,
        )
        r1.set_spline_surface(s_bt, "bottom", 1000)
    else:
        r1 = Block(mesh)
        r1.set_connection(c1, "left")
        r1.p1 = coords2D(cartesian(r1_r_bt, phi1, s_bt(r1_r_bt)))
        r1.p2 = coords2D(cartesian(r1_r_bt, phi2, s_bt(r1_r_bt)))
        r1.p5 = coords2D(cartesian(r1_r_top, phi1, r1_z_top))
        r1.p6 = coords2D(cartesian(r1_r_top, phi2, r1_z_top))
        r1.cells_x1 = res_r_r1
        r1.create()
        r1.e0.type = "spline"
        r1.e0.points = evaluate_spline(c1_r_bt, r1_r_bt, phi1, s_bt)
        r1.e1.type = "spline"
        r1.e1.points = evaluate_spline(c1_r_bt, r1_r_bt, phi2, s_bt)
    if dimension == 3:
        r2 = create_ring(
            mesh,
            [r2_r_top, s_fs(r2_r_top)],
            [r_crucible, s_bt(r_crucible)],
            r1.surf_rad,
            res_z_c2,
            res_phi,
            res_z_c1,
        )
        r2.set_spline_surface(s_bt, "bottom", 5000)  # high number required here together with dirtyfix from above
        r2.set_spline_surface(s_fs, "side", 1000)
    else: 
        r2 = Block(mesh)
        r2.set_connection(r1, "left")
        r2.p1 = coords2D(cartesian(r_crucible, phi1, s_bt(r_crucible)))
        r2.p2 = coords2D(cartesian(r_crucible, phi2, s_bt(r_crucible)))
        r2.p5 = coords2D(cartesian(r2_r_top, phi1, s_fs(r2_r_top)))
        r2.p6 = coords2D(cartesian(r2_r_top, phi2, s_fs(r2_r_top)))
        r2.cells_x1 = res_z_c2
        r2.create()
        r2.e0.type = "spline"
        r2.e0.points = evaluate_spline(r1_r_bt, r_crucible, phi1, s_bt, res=100000)
        r2.e1.type = "spline"
        r2.e1.points = evaluate_spline(r1_r_bt, r_crucible, phi2, s_bt, res=100000)
    if dimension == 3:
        r3 = create_ring(
            mesh,
            [r2_r_top, s_fs(r2_r_top)],
            [r1_r_top, r1_z_top],
            c2.surf_rad,
            res_r_r1,
            res_phi,
            res_z_c2,
            faces_outside=r2.surf_top,
        )
        r3.set_spline_surface(s_fs, "top", 1000)
    else:
        r3 = Block(mesh)
        r3.set_connection(c2, "left")
        r3.set_connection(r1, "bottom")
        r3.face_right = r2.face_top
        r3.create()
        r3.e3.type = "spline"
        r3.e3.points = evaluate_spline(r_crystal, r2_r_top, phi1, s_fs)
        r3.e2.type = "spline"
        r3.e2.points = evaluate_spline(r_crystal, r2_r_top, phi2, s_fs)

    # ####################
    # Grading
    if dimension == 3:
        c2.set_grading_axial(grading_top)
        r3.set_grading_axial(grading_top)
        r2.set_grading_radial(grading_top)

        c1.set_grading_axial(grading_bottom)
        r1.set_grading_axial(grading_bottom)
        r2.set_grading_axial(grading_bottom)

        r1.set_grading_radial(grading_meniscus)
        r3.set_grading_radial(grading_meniscus)

        c1.set_grading_radial(grading_crys_rad)
        c2.set_grading_radial(grading_crys_rad)
    else:
        c1.grading = f"simpleGrading ({grading_crys_rad} 1 {grading_bottom})"
        c2.grading = f"simpleGrading ({grading_crys_rad} 1 {grading_top})"
        r1.grading = f"simpleGrading ({grading_meniscus} 1 {grading_bottom})"
        r2.grading = f"simpleGrading ({grading_top} 1 {grading_bottom})"
        r3.grading = f"simpleGrading ({grading_meniscus} 1 {grading_top})"

    ####################
    # Patches
    crucible_bot = Patch(mesh, "wall crucibleBot")
    if dimension == 3:
        crucible_bot.faces += c1.surf_bt
        crucible_bot.faces += r1.surf_bt
    else:
        crucible_bot.add_face(c1.face_bottom)
        crucible_bot.add_face(r1.face_bottom)
    crucible_side = Patch(mesh, "wall crucibleSide")
    if dimension == 3:
        crucible_side.faces += r2.surf_bt
    else:
        crucible_side.add_face(r2.face_bottom)
    free_surf = Patch(mesh, "wall freeSurf")
    if dimension == 3:
        free_surf.faces += r2.surf_rad
        free_surf.faces += r3.surf_top
    else:
        free_surf.add_face(r2.face_right)
        free_surf.add_face(r3.face_top)
    top_surf = Patch(mesh, "wall crysInter")
    if dimension == 3:
        top_surf.faces += c2.surf_top
    else:
        top_surf.add_face(c2.face_top)

    if dimension == 2:
        mesh_front = Patch(mesh, "wedge front")
        mesh_front.add_face(c1.face_front)
        mesh_front.add_face(c2.face_front)
        mesh_front.add_face(r1.face_front)
        mesh_front.add_face(r2.face_front)
        mesh_front.add_face(r3.face_front)

        mesh_back = Patch(mesh, "wedge back")
        mesh_back.add_face(c1.face_back)
        mesh_back.add_face(c2.face_back)
        mesh_back.add_face(r1.face_back)
        mesh_back.add_face(r2.face_back)
        mesh_back.add_face(r3.face_back)

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

    lorentz_force_to_txt(f"{elmer_results}/../case_melt_t0001.vtu", f"{sim_dir}/lorentz-force.csv")
    lorentz_force_to_txt(f"{elmer_results}/../case_melt_t0001.vtu", f"{sim_dir}/lorentz-force_with-header.csv", header=True)


if __name__ == "__main__":
    prepare_simulation(
        dimension=2,
        sim_dir="simdata_debug_openfoam_linearUpwind_init-1685_bc-heatflux",
        elmer_results="simdata/2023-07-28_10-35_eo_exp_15kg_DC/02_simulation/elmer_01/results",
        # elmer_results="simdata_elmer/2023-06-02_13-01_ss_30-kg-melt_2heaters_50A_90deg/02_simulation/results",
        # use_hardcoded_grading=True
    )
    # prepare_simulation(
    #     dimension=3,
    #     sim_dir="simdata_openfoam_3D",
    #     elmer_results="simdata_elmer/2023-05-10_14-12_ss_2-heaters_with-joule-heat_90-deg-phaseshift/02_simulation/results",
    #     use_hardcoded_grading=False
    # )
