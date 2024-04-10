from opencgs.sim import CoupledSim, QuasiTransientSim
import opencgs.control as ctrl
from setup_elmer import geometry, simulation
from setup_openfoam import prepare_simulation

if __name__ == "__main__":
    try:
        git_metadata = ctrl.get_git_metadata()
    except:
        git_metadata = "not available"
    config = ctrl.load_config("./config.yml")
    config.update({"metadata": git_metadata})
    config_geo = ctrl.load_config("./config_geo.yml")
    config_sim = ctrl.load_config("./config_sim.yml")
    config_mat = ctrl.load_config("./config_mat.yml")
    config_coupling = ctrl.load_config("./config_coupling.yml")
    config_of = ctrl.load_config("./config_openfoam.yml")

    if "quasi_transient" in config:
        sim = QuasiTransientSim(
            geo=geometry,
            geo_config=config_geo,
            sim=simulation,
            setup_openfoam=prepare_simulation,
            setup_openfoam_config=config_of,
            coupling_config=config_coupling,
            sim_config=config_sim,
            mat_config=config_mat,
            simulation_class=CoupledSim,
            **config,
        )
    else:
        sim = CoupledSim(
            geometry,
            config_geo,
            simulation,
            config_sim,
            prepare_simulation,
            config_of,
            config_coupling,
            config_mat,
            **config,
        )
    sim.execute()
