# COLLAPSO1D - Core Collapse Supernova in 1D <a href="https://github.com/pikarpov-LANL/COLLAPSO1D"><img src="https://github.com/pikarpov-LANL/COLLAPSO1D/blob/images/docs/images/collapso1d.png?raw=true"  alt="COLLAPSO1D logo" align="right" width="80"></a>

This is a 1D lagrangian code for modeling Core-Collapse Supernovae (CCSN). For progenitors, it takes KEPLER generated data. Currently, there is support for [Heger & Woosley, 2000:](https://2sn.org/stellarevolution/) and [Sukhbold et al, 2016](http://doi.org/10.17617/1.b).

Turbulence is treated through mixing length theory (MLT) and Machine Learning (ML) based models. The latter has been trained using the [Sapsan](https://github.com/pikarpov-LANL/Sapsan) ML pipeline.

PyTorch is implemented based on [pytorch-fortran](https://github.com/alexeedm/pytorch-fortran).


## [Documentation](https://pikarpov-LANL.github.io/COLLAPSO1D/)
For installation, quick start, and the details of the code please refer to the documentation.

-------

COLLAPSO1D has a BSD-style license, as found in the [LICENSE](https://github.com/pikarpov-LANL/COLLAPSO1D/blob/master/LICENSE) file.

© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.