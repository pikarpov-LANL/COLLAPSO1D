# 1D Core-Collapse Supernovae Code

1. [Dependencies](README.md#dependencies)
2. [Quick Start](README.md#quick-start)
3. [Setup Details](README.md#setup-details)
4. [Make Commands](README.md#make-commands)

This is a 1D lagrangian code to explore CCSN modeling. For progenitors, it takes KEPLER generated data. Currently, there is support for [Heger & Woosley, 2000:](https://2sn.org/stellarevolution/) and [Sukhbold et al, 2016](https://arxiv.org/abs/1510.04643).

Turbulence is treated through mixing length theory (MLT) and Machine Learning (ML) based models. The latter has been trained using the [Sapsan](https://github.com/pikarpov-LANL/Sapsan) ML pipeline.

PyTorch is implemented based on [pytorch-fortran](https://github.com/alexeedm/pytorch-fortran).

## Dependencies

### GFORTRAN
Both the main CCSN code and the PyTorch wrapper can be compiled with `gfortran >= 9.4.0`.

### PyTorch
To install the latest PyTorch for CPU, follow the official [instructions](https://pytorch.org/). Lastly, we need to make sure all compilers link correctly.

### CMake
Make sure you have cmake or install it by (tested on cmake==3.22.1)
```bash
sudo apt install cmake
```

### ~/.bashrc
Add the following to you `~/.bashrc` and then `source ~/.bashrc`:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
```

### Speed-up `make` commands
`CMAKE_PREFIX_PATH` is a variable set in `Makefile` pointing to the location of Torch cmake config files. By default, it is determined by importing torch which is fairly slow, but it can also be set explicitely. To get the torch `{Location}`:
```bash
pip show torch
```
Then set `CMAKE_PREFIX_PATH={Locaton}/torch/share/cmake` in `Makefile`.

---
## Quick Start
To process progenitor data and to compile the code with PyTorch included, all you need to do is run:
```shell
make
```
This will move the processed input data and executable to the project folder. To run the model:
```shell
cd project/1dmlmix
./1dmlmix
```
Congratulations! You just blew up a star!

## Setup details

### Progenitor Data Setup

You can setup parameters and choose which stellar progenitor data to prepare. In our case, we are using Stan Woosley's datasets.
```shell
vi prep_data/setup
```
to process and move the data to the project folder:
```shell
make data
```
Take note of the `Number of cells` (depend on the initial cell mass), since that determines the grid size (<Number of Cells (from Data)>) parameter for the main simulation setup.

### 1dccsn simulation setup
All of the simulation parameters can be adjust here (don't forget about `<Number of Cells (from Data)>`):
```shell
vi project/1dmlmix/setup
```
Next to compile:
```shell
make project
```
Lastly to run:
```shell
cd project/1dmlmix
./1dmlmix
```
to process and move the data to the project folder:
```shell
make data
```
Take note of the `Number of cells` (depend on the initial cell mass), since that determines the grid size (<Number of Cells (from Data)>) parameter for the main simulation setup.

### Setup Conversion from Binary to Readable Output
The setup can be done either by editing `setup_readout` or by cmd arguments. In the first case:
```shell
cd project/1dmlmix
vi setup_readout
./readout
```
or you can provide the same 3 arguemnts, (Input, Output, #Dumps), as arguments to the `readout` executable:
```shell
cd project/1dmlmix
./readout Input Output ndumps
```

---
### Make Commands
Data preparation:
```shell
make data
```
Model compilation
```shell
make project
```
To convert from binary output to readable tables, you need to run `readout`. Compile it by:
```shell
make readout
```
Combined data prep, model compilation, and to get a binary to readable output executable:
```shell
make
```
In addition, you can prepare examples:
```shell
make examples
```
Lastly, clean up everything:
```shell
make clean
```

## [Notes](https://www.overleaf.com/read/pgsnmxgdjkrq)

Notes on ML subgrid turbulence model implementation within this code can be found on [Overleaf](https://www.overleaf.com/read/pgsnmxgdjkrq).


-------
1dccsn code has a BSD-style license, as found in the [LICENSE](https://github.com/pikarpov-LANL/1dccsn/blob/master/LICENSE) file.

© 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.