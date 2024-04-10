import opencgs.control as ctrl
from setup_elmer import geometry, simulation

    
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
    simulation(sim_config, model, sim_dir)

    # use this to execute without meshing (uses model.pkl in sim_dir)
    # simulation(sim_config, None, "simdata/crystal_length_87mm")

