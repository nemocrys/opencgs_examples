# sim-kristmag-si
2D-2D or 2D-3D simulation of silicon Czochralski growth using a TMF based on the KristMAG technology. Azimuthal averaging with [azimuthalAverage](azimuthalAverage) is based on the code provided in [nemoFoam-utils](https://github.com/nemocrys/nemoFoam-utils).

## Overview

An overview of the simulation setup can be found [here](figures/setup.png). The following result was obtained with an AC current:

![result-2D-simulation](figures/2D-3D_AC-global.png)

The corresponding flow veloctiy and temperature fields (left: snapshot, right: time-average) are:

![result-flow-velocity](figures/AC_flow-velocity.png)
![result-flow-temperature](figures/AC_temperature.png)

## Configuration, setup, and execution

- The configuration of the simulations is stored in the yaml-files:
  - Main configurations such as the simulation name are set in [config.yml](config.yml).
  - Geometry parameters are defined in [config_geo.yml](config_geo.yml). Note, that some parameters of the meshing are directly set in [setup_elmer.py](setup_elmer.py) and [setup_openfoam.py](setup_openfoam.py).
  - The global Elmer simulation is configured in [config_sim.yml](config_sim.yml).
  - The material properties used in the global Elmer simulation are configured in [config_mat.yml](config_mat.yml).
  - The coupling is configured in [config_coupling.yml](config_coupling.yml).
- The mesh of the global model is set up in [setup_elmer.py](setup_elmer.py), which contains also the setup of the global simulation.
- The mesh of the flow model is set up in [setup_openfoam.py](setup_openfoam.py).
- The flow model is configured in the corresponding templates in [openfoam_template_steady](openfoam_template_steady), [openfoam_template_transient](openfoam_template_transient), and [openfoam_template_transient_3D](openfoam_template_transient_3D). The most important settings can be found in constant/transportProperties (material parameters), system/controlDict (simulation time, output), system/changeDictionaryDict (coupling with heat fluxes or fixed temperatures, Marangoni coefficient, crystal and crucible rotation), and system/elmerToFoamDict (material parameters).
- The simulation is executed using the run scripts:
  - Simulations using the global Elmer model only are executed with the [run_elmer_simulation.py](run_elmer_simulation.py) script.
  - Simulations using the coupled Elmer-OpenFOAM model only are executed with the [run_coupled_simulation.py](run_coupled_simulation.py) script.

## Additional details

For a more detailed description including simulation results see:

> A. Wintzer, *Validation of multiphysical models for Czochralski crystal growth*. PhD thesis, Technische Universität Berlin, Berlin, 2024.
