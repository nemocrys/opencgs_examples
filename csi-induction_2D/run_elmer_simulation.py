"""
Main script to execute the simulation.
"""
import opencgs.control as ctrl
from opencgs.sim import ParameterStudy, SteadyStateSim, QuasiTransientSim
from setup_elmer import geometry, simulation


if __name__ == "__main__":
    try:
        git_metadata = ctrl.get_git_metadata()
    except:
        git_metadata = "not available"  # no git on cluster nodes

    geo_config = ctrl.load_config("./config_geo.yml")
    sim_config = ctrl.load_config("./config_sim.yml")
    mat_config = ctrl.load_config("./config_mat.yml")
    config = ctrl.load_config("./config.yml")
    config.update({"metadata": git_metadata})

    # This is used to run steady state simulations / parameter studies
    if "quasi_transient" in config:
        sim = QuasiTransientSim(
            geo=geometry,
            geo_config=geo_config,
            sim=simulation,
            sim_config=sim_config,
            mat_config=mat_config,
            **config)
    elif "study_params" in config:
        sim = ParameterStudy(SteadyStateSim, geometry, geo_config, simulation, sim_config, mat_config, **config)
    else:
        sim = SteadyStateSim(
            geometry, geo_config, simulation, sim_config, mat_config, **config
        )

    sim.execute()
