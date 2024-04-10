import os
import shutil
from setup_openfoam import prepare_simulation
from extract_lorentz_force import lorentz_force_to_txt
from pyelmer.post import dat_to_dataframe

########################################################################
# input graphite
# sim_dir = "simdata_of-3D_graphite_cluster_highres"
# elmer_2D_result = "simdata/2023-09-28_11-47_eo_graphite_transient_hf-bc/02_simulation/elmer_08/results"
# elmer_3D_result = "C:/Users/enders-seidlitz/Documents/Results/2023-10-05_sim-test-cz-induction-3D/simdata_graphite_crucible_restart-s2s/results"
# heater_power_scaling = 1.4719  # graphite crucible, radiation s2s

# input aluminum
sim_dir = "simdata_of-3D_aluminum_cluster_highres"
elmer_2D_result = "simdata/2023-09-28_11-47_eo_aluminum_transient_hf-bc/02_simulation/elmer_03/results"
elmer_3D_result = "C:/Users/enders-seidlitz/Documents/Results/2023-10-05_sim-test-cz-induction-3D/simdata_aluminum_crucible_restart-s2s/results"
heater_power_scaling = 4.7156  # graphite crucible, radiation s2s

# Apply the scaling factor from effective heat conductivity to scale the temperature gradient
gradient_scaling = 2.14
########################################################################

# create setup based on 2D data
prepare_simulation(3, sim_dir, elmer_2D_result, "openfoam_template_transient")

# replace 2D with 3D input
for file in [
    "lorentz-force.csv",
    "lorentz-force_with-header.csv",
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

lorentz_force_to_txt(f"{elmer_3D_result}/case_melt_t0001.vtu", f"{sim_dir}/lorentz-force.csv", scaling=heater_power_scaling, coordinate_permutation=False)
lorentz_force_to_txt(f"{elmer_3D_result}/case_melt_t0001.vtu", f"{sim_dir}/lorentz-force-with-header.csv", scaling=heater_power_scaling, coordinate_permutation=False, header=True)
