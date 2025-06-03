# sn-resistance-test_2D
2D simulation of a heating test with resistance heating.


## Overview

The current resistance heating test is based on the previous setup, which can be found [here](https://github.com/nemocrys/opencgs_examples/tree/heating_tests_readme/sn-resistance_2D). The following results correspond to the heating test with insulation:



![result-2D-simulation](figures/resistive_T.png)


## Configuration, setup, and execution

- The configuration of the simulations is stored in the yaml-files:
    - Main configurations such as the simulation name are set in [config.yml](config.yml).
    - Geometry parameters are defined in [config_geo.yml](config_geo.yml). Note, that some parameters of the meshing are directly set in [setup.py](setup.py).
    - The global Elmer simulation is configured in [config_sim.yml](config_sim.yml).
    - The material properties used in the global Elmer simulation are configured in [config_mat.yml](config_mat.yml).

- The mesh of the global model is set up in [setup.py](setup.py), which contains also the setup of the global simulation.

- The simulation is executed using the run script:
  - Simulations using the global Elmer model only are executed with the [run.py](run.py) script.

## Heat flux estimation

- After the simulation execution the heat fluxes can be found in :
    - boundary-scalars.dat and  boundary-scalars.dat.names files within the similation folder (Elmer Output).
    - heat-fluxes.yml file at the result folder (OpenCGS post-processing).
    

The results are in SI units (W) and the sign points the direction according to the boundary definition in  [setup.py](setup.py).


Simulation heat flux results for T<sub>cr</sub> = 780 °C :

| Boundary Name                     | Value (W) | Description                                     |
|-----------------------------------|-----------|-------------------------------------------------|
| vessel_if_axbt_vessel             | 238       | Vessel-Bottom Axis interface                    |
| vessel_bnd_vessel_outside         | 1633      | Outer Vessel boundary (water-cooled boundary)    |
| heater_bnd_heater                 | 1637      | Heater-Atmosphere boundary                      |


The heating power is computed iteratively to match the target crucible temperature. From the Elmer output, the **heater power scaling** : 3.27. Multiplying the scaling factor by the initial power value declared in [config_sim.yml](config_sim.yml) (i.e 500 W)  gives a total power of **1637 W**.




## Additional details

- For a more detailed description including simulation results see:

> Sepehr Foroushani, Arved Wintzer, Frank-Michael Kiessling, and Kaspars Dadzis. “Heating efficiency and energy saving potential of Czochralski crystal growth furnaces.” *Journal of Crystal Growth*, vol. 662, 2025, Art. 128106. 
> DOI: https://doi.org/10.1016/j.jcrysgro.2025.128106  