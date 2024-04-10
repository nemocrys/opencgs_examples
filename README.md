# opencgs_examples

Collection of crystal growth simulations using [opencgs](https://github.com/nemocrys/opencgs).

The project is developed and maintained by the [**Model experiments group**](https://www.ikz-berlin.de/en/research/materials-science/section-fundamental-description#c486) at the Leibniz Institute for Crystal Growth (IKZ).

### Referencing
If you use this code in your research, please cite the reference provided in the respective directory or this open-access article:

> A. Enders-Seidlitz, J. Pal, and K. Dadzis, Development and validation of a thermal simulation for the Czochralski crystal growth process using model experiments *Journal of Crystal Growth*,  593 (2022) 126750. [https://doi.org/10.1016/j.jcrysgro.2022.126750](https://doi.org/10.1016/j.jcrysgro.2022.126750).

## Overview
This repository contains various 2D and 2D models of Czochralski crystal growth using [Elmer](http://www.elmerfem.org/blog/) and [OpenFOAM](https://www.openfoam.com/).

- [csi-induction_2D](csi-induction_2D): 2D simulation of CsI Czochralski growth
- [csi-induction_3D](csi-induction_3D): 3D simulation of CsI Czochralski growth
- [kristmag-si_2D-3D](kristmag-si_2D-3D): 2D-2D or 2D-3D simulation of Si Czochralski growth using a TMF based on the KristMAG technology
- [sn-induction_2D](sn-induction_2D): 2D simulation of Sn Czochralski growth with induction heating
- [sn-induction_3D](sn-induction_3D): 3D simulation of Sn Czochralski growth with induction heating
- [sn-resistance_2D](sn-resistance_2D): 2D simulation of Sn Czochralski growth with resistance

Additional examples can be found here:
- [2D simulation of Sn Czochralski growth with induction heating](https://github.com/nemocrys/test-cz-induction)
- [2D simulation of GaAs vertical gradient freeze](https://github.com/nemocrys/vertical-gradient-freeze)

## Computational setup (Docker)

The setup for the simulations is provided in form of a docker image, so just an installation of [Docker](https://docs.docker.com/get-docker/) is required on your system. The image nemocrys/opencgs:v1.0.0 is used (see [opencgs](https://github.com/nemocrys/opencgs) for more information).

On Windows, the container can be started with the following command:
```
docker run --rm -v ${PWD}:/home/workdir nemocrys/opencgs:v1.0.0 bas
```
On Linux, the container can be started with:
```
docker run -it --rm -v $PWD:/home/workdir -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) nemocrys/opencgs:v1.0.0 bash
```

This will map the current working directory (e.g., a copy of this repository) into the container and, on Linux, set the user's group and user id. The simulation can then be executed using the provided `python3` or `Allrun`scripts.

## Model description

*Add description or reference here*

## Acknowledgements

[This project](https://nemocrys.github.io/) has received funding from the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 851768).

<img src="https://github.com/nemocrys/test-cz-induction/blob/main/EU-ERC.png">