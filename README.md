# 1D Core-Collapse Supernovae Code

This is a 1D lagrangian code to explore CCSN modeling. For progenitors, it takes KEPLER generated [data](https://2sn.org/stellarevolution/) of Alex Heger & Stan Woosley.

Turbulence is treated through mixing length theory (MLT) and Machine Learning (ML) based models.

PyTorch is implemented based on [pytorch-fortran](https://github.com/alexeedm/pytorch-fortran).

## Dependencies

### NVFORTRAN
While the base code can be compiled via `gfortran`, the PyTorch implimentation requires the [nvidia HPC toolkit 21.9](https://developer.nvidia.com/nvidia-hpc-sdk-219-downloads). You can register and download it for free via the link. Since it includes CUDA 11.4, while the maximum supported by PyTorch is 11.3, we also need to setup CUDA separately. 

### PyTorch
To install the latest PyTorch, follow the official [instructions](https://pytorch.org/). Lastly, we need to make sure all compilers link correctly.

### ~/.bashrc
Add the following to you `~/.bashrc` and then `source ~/.bashrc`:
```bash
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/compilers/bin/:$PATH
export CUDACXX=/usr/local/cuda/bin/nvcc
export CUDA_HOME=/usr/local/cuda
export PATH=/usr/local/cuda:$PATH
export PATH=/usr/local/cuda/bin:$PATH
export CUDA_HOME=/usr/local/cuda
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
```

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
vi read_data/setup
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
