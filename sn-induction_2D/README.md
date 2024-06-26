# sn-induction_2D
2D simulation of CsI Czochralski growth with induction heating using an iterative coupling between a global Elmer model and a melt flow OpenFOAM model. This model is based on https://github.com/nemocrys/test-cz-induction.

## Overview

An overview of the simulation setup can be found [here](https://camo.githubusercontent.com/c0b28a2c445645c2045e898809f38f6a29eb8656abb430bd06d42fff57a1543a/68747470733a2f2f6172732e656c732d63646e2e636f6d2f636f6e74656e742f696d6167652f312d73322e302d53303032323032343832323030323338582d6772325f6c72672e6a7067). The following result was obtained using a graphite crucible:

![result-2D-simulation](figures/evaluation-2D_graphite_x2.png)

## Configuration, setup, and execution

- The configuration of the simulations is stored in the yaml-files:
  - Main configurations such as the simulation name are set in [config.yml](config.yml).
  - Geometry parameters are defined in [config_geo.yml](config_geo.yml). Note, that some parameters of the meshing are directly set in [setup_elmer.py](setup_elmer.py) and [setup_openfoam.py](setup_openfoam.py).
  - The global Elmer simulation is configured in [config_sim.yml](config_sim.yml).
  - Iterative crystal diameter computation can be configured in [config_di.yml](config_sim.yml).
  - The material properties used in the global Elmer simulation are configured in [config_mat.yml](config_mat.yml).
  - The coupling is configured in [config_coupling.yml](config_coupling.yml).
- The mesh of the global model is set up in [setup_elmer.py](setup_elmer.py), which contains also the setup of the global simulation.
- The mesh of the flow model is set up in [setup_openfoam.py](setup_openfoam.py).
- The flow model is configured in the corresponding templates in [openfoam_template_steady](openfoam_template_steady) or [openfoam_template_transient](openfoam_template_transient). The most important settings can be found in constant/transportProperties (material parameters), system/controlDict (simulation time, output), system/changeDictionaryDict (coupling with heat fluxes or fixed temperatures), system/elmerToFoamDict (material parameters), and 0.orig/U (crystal rotation).
- The simulation is executed using the run scripts:
  - Simulations using the global Elmer model only are executed with the [run_elmer_simulation.py](run_elmer_simulation.py) script.
  - Simulations using the coupled Elmer-OpenFOAM model only are executed with the [run_coupled_simulation.py](run_coupled_simulation.py) script.

## Referencing
If you use this code in your research, please cite:

> A. Wintzer, *Validation of multiphysical models for Czochralski crystal growth*. PhD thesis, Technische Universität Berlin, Berlin, 2024.

**Note: this reference will be updated soon**
