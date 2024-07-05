# sim-test-cz-csi-3d
3D simulation of CsI Czochralski growth with induction heating using a coupling between a global Elmer model and a melt flow OpenFOAM model.

## Overview

An overview of the simulation setup can be found [here](../csi-induction_2D/figures/setup.png).  With the 3D model, the following result was obtained after a single coupling iteration:

![result-3D-simulation-T-EM](figures/l=87mm_T-control.png)

![result-3D-simulation-flow](figures/l=87mm_flow.png)


## Configuration, setup, and execution

This model is based on the 2D model provided [here](../csi-induction_2D).

- The configuration of the simulation is stored in the yaml-files:
  - Geometry parameters are defined in [config_geo.yml](config_geo.yml). Note, that some parameters of the meshing are directly set in [setup_elmer.py](setup_elmer.py) and [setup_openfoam.py](setup_openfoam.py).
  - The global Elmer simulation is configured in [config_sim.yml](config_sim.yml).
  - The material properties used in the global Elmer simulation are configured in [config_mat.yml](config_mat.yml).
- The mesh of the global model is set up in [setup_elmer.py](setup_elmer.py), which contains also the setup of the global simulation.
- The mesh of the flow model is set up in [setup_openfoam.py](setup_openfoam.py).
- The flow model is configured in the corresponding templates in [openfoam_template_steady_3D](openfoam_template_steady_3D) or [openfoam_template_transient_3D](openfoam_template_transient_3D). The most important settings can be found in constant/transportProperties (material parameters), system/controlDict (simulation time, output), system/changeDictionaryDict (coupling with heat fluxes or fixed temperatures) and 0.orig/U (crystal rotation).
- The flow model is set up using [setup_openfoam_3d.py](setup_openfoam_3d.py) with the path to the results of the global model defined in the script.
- The simulation is executed in the following way:
  - The global model is executed with Python using [run.py](run.py).
  - The flow model is executed using the `Allrun` script.
- Feedback from the flow model to the global model is possible, but not implemented in an automatized way. Manual coupling by modifying the Elmer sif-file was performed in the reference.

It should be noted that mesh generation for the global 3D model requires a significant amount of memory and was performed on a simulation cluster for the example shown above.

## Additional details

For a more detailed description including simulation results see:

> A. Wintzer, *Validation of multiphysical models for Czochralski crystal growth*. PhD thesis, Technische Universität Berlin, Berlin, 2024.
