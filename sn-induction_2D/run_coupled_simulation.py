from opencgs.sim import CoupledSim
import opencgs.control as ctrl

from setup_elmer import geometry, simulation
from setup_openfoam import prepare_simulation


if __name__ == "__main__":
    try:
        git_metadata = ctrl.get_git_metadata()
    except:
        git_metadata = "not available"
    config = ctrl.load_config(
        "./config.yml"
    )  # only sim_name, quasi_transient is used
    config.update({"metadata": git_metadata})
    config_geo = ctrl.load_config("./config_geo.yml")
    config_sim = ctrl.load_config("./config_sim.yml")
    config_mat = ctrl.load_config("./config_mat.yml")
    config_coupling = ctrl.load_config("./config_coupling.yml")
    config_of = ctrl.load_config("./config_openfoam.yml")

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
