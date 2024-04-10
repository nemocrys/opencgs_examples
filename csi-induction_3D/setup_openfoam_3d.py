import os
import shutil
from setup_openfoam import prepare_simulation
from pyelmer.post import dat_to_dataframe


########################################################################
# input l = 10mm
sim_dir = "simdata_of-3D_l=10mm_steady"
elmer_2D_result = "C:/Users/enders-seidlitz/Documents/Results/2023-09-26_csi-exp-2_sim-test-cz-csi-2d_T-coupling/simdata/T-coupling_2023-09-28_06-02_qt_csi-exp-2_1mm-steps_relax-0.5_conv-5e-5_cluster/02_selected_simulations/eo_length=0.01_vpull=0.2/02_simulation/elmer_08/results"
elmer_3D_result = "C:/Users/enders-seidlitz/Documents/Results/2023-10-24_csi-3d/simdata/crystal_length_10mm/results"

# # input l = 87mm
# sim_dir = "simdata_of-3D_l=87mm_steady"
# elmer_2D_result = "C:/Users/enders-seidlitz/Documents/Results/2023-09-26_csi-exp-2_sim-test-cz-csi-2d_T-coupling/simdata/T-coupling_2023-09-28_06-02_qt_csi-exp-2_1mm-steps_relax-0.5_conv-5e-5_cluster/02_selected_simulations/eo_length=0.087_vpull=0.2/02_simulation/elmer_07/results"
# elmer_3D_result = "C:/Users/enders-seidlitz/Documents/Results/2023-10-24_csi-3d/simdata/crystal_length_87mm/results"

# Apply the scaling factor from effective heat conductivity to scale the temperature gradient
gradient_scaling = 2
# of_template = "openfoam_template_transient_3D"
of_template = "openfoam_template_steady_3D"
########################################################################

# create setup based on 2D data
prepare_simulation(3, sim_dir, elmer_2D_result, of_template)
# replace 2D with 3D input
for file in [
    "save_line_crc_melt_relaxed.dat",
    "save_line_melt_crys_relaxed.dat",
    "save_line_melt_surf_relaxed.dat",
    ]:
    os.remove(f"{sim_dir}/{file}")

for file in [
    "save_line_crc_melt.dat",
    "save_line_melt_crys.dat",
    "save_line_melt_surf.dat",
    "save_line_crc_melt.dat.names",
    "save_line_melt_crys.dat.names",
    "save_line_melt_surf.dat.names",
    ]:
    shutil.copy(f"{elmer_3D_result}/{file}", f"{sim_dir}/{file}")

for file in ["save_line_crc_melt.dat", "save_line_melt_crys.dat", "save_line_melt_surf.dat"]:
    df = dat_to_dataframe(f"{sim_dir}/{file}")
    for column in ["temperature grad 1", "temperature grad 2", "temperature grad 3"]:
        df[column] *= gradient_scaling
    df.to_csv(f"{sim_dir}/{file}", header=False, sep=" ", index=False)


