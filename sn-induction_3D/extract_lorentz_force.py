import meshio
import numpy as np
import pandas as pd

# NOTE: apply correction for direction of B in 2D only, otherwise own calculation of F_L has wrong sign in 3D
def lorentz_force_to_txt(vtu_file, output_file, plot=False, header=False, scaling=1, coordinate_permutation=False):
    meshvtu = meshio.read(vtu_file)

    # take element data
    elements = meshvtu.cells[0]
    joule_heat = meshvtu.cell_data["joule heating e"][0]
    jxb = meshvtu.cell_data["jxb re e"][0]
    j_re = meshvtu.cell_data["current density re e"][0]
    j_im = meshvtu.cell_data["current density im e"][0]
    b_re = meshvtu.cell_data["magnetic flux density re e"][0]
    b_im = meshvtu.cell_data["magnetic flux density im e"][0]

    data = []
    for i in range(len(elements.data)):  # for each element
        element = elements.data[i]
        # compute center of element
        x_mid = 0
        y_mid = 0
        z_mid = 0
        for node in element:
            x_mid += meshvtu.points[node][0]
            y_mid += meshvtu.points[node][1]
            z_mid += meshvtu.points[node][2]
        number_of_nodes = len(element)
        x_mid /= number_of_nodes
        y_mid /= number_of_nodes
        z_mid /= number_of_nodes
        # compute lorentz force
        f_l = 0.5 * (
            np.cross(
                [j_re[i][0], j_re[i][1], j_re[i][2]],
                [-b_re[i][0], -b_re[i][1], -b_re[i][2]],  # correct direction of B
            )
            + np.cross(
                [j_im[i][0], j_im[i][1], j_im[i][2]],
                [-b_im[i][0], -b_im[i][1], -b_im[i][2]],  # correct direction of B
            )
        )
        jxb_x = jxb[i][0]
        jxb_y = jxb[i][1]
        jxb_z = jxb[i][2]
        fl_x = f_l[0]
        fl_y = f_l[1]
        fl_z = f_l[2]
        if coordinate_permutation:  
            x_mid, y_mid, z_mid = x_mid, z_mid, y_mid
            jxb_x, jxb_y, jxb_z = jxb_x, jxb_z, jxb_y
            fl_x, fl_y, fl_z = fl_x, fl_z, fl_y
        data.append(
            [
                i,
                number_of_nodes,
                x_mid,
                y_mid,
                z_mid,
                joule_heat[i][0],
                jxb_x,
                jxb_y,
                jxb_z,
                fl_x,
                fl_y,
                fl_z,
            ]
        )
    df = pd.DataFrame(
        data,
        columns=[
            "Element",
            "Number_of_nodes",
            "x_mid",
            "y_mid",
            "z_mid",
            "joule_heat",
            "jxb_x",
            "jxb_y",
            "jxb_z",
            "F_l_x",
            "F_l_y",
            "F_l_z",
        ],
    )
    for column in ["joule_heat", "jxb_x", "jxb_y", "jxb_z", "F_l_x", "F_l_y", "F_l_z"]:
        df[column] *= scaling
    
    if plot:
        import matplotlib.pyplot as plt
        import matplotlib.tri as tri
        fig, ax = plt.subplots(3, 2, figsize=(10, 8))
        x = df["x_mid"].to_numpy()
        y = df["y_mid"].to_numpy()
        for i, component in zip([0, 1, 2], ["x", "y", "z"]):
            f_min = min([df[f"jxb_{component}"].min(), df[f"F_l_{component}"].min()])
            f_max = max([df[f"jxb_{component}"].max(), df[f"F_l_{component}"].max()])
            ax[i, 0].scatter(x, y, c=df[f"jxb_{component}"], cmap="turbo", vmin=f_min, vmax=f_max)
            ax[i, 1].scatter(x, y, c=df[f"F_l_{component}"], cmap="turbo", vmin=f_min, vmax=f_max)
            norm = plt.Normalize(vmin=f_min, vmax=f_max)
            scalarmap = plt.cm.ScalarMappable(norm=norm, cmap='turbo')
            cbar = fig.colorbar(scalarmap, ax=[ax[i, 0], ax[i, 1]], location='right', aspect=50)
            cbar.set_label(f"Lorentz force density {component}")
            ax[i, 0].set_aspect("equal")
            ax[i, 1].set_aspect("equal")
        ax[0, 0].set_title("Elmer jxb re")
        ax[0, 1].set_title("Own calculation")
        # fig.tight_layout()
        plt.show()
    df.to_csv(output_file, header=header)


if __name__ == "__main__":
    lorentz_force_to_txt(
        "./simdata_elmer/2023-05-10_14-12_ss_2-heaters_with-joule-heat_90-deg-phaseshift/02_simulation/case_melt_t0010.vtu",
        "./lorentz-force.csv",
        True
    )
