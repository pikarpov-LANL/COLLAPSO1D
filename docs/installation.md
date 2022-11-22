# Installation

## Dependencies
### GFORTRAN
Both the main CCSN code and the PyTorch wrapper can be compiled with `gfortran >= 9.4.0`. If missing, install it via:
```
sudo apt install gfortran
```

### PyTorch
To install the latest PyTorch for <ins>CPU</ins>, follow the official [instructions](https://pytorch.org/). I would highly recommend installing it in a dedicated conda environemnt. As an example, here is how to create one and install PyTorch:
```
conda create -n py310 python=3.10
conda activate py310

# check the pytorch installation instructions!
pip install torch==1.11.0+cpu --extra-index-url https://download.pytorch.org/whl/cpu
```
!!! Warning
    `torch>=1.12.0` will work for inferencing, hence for the CCSN code will be fine, but loading and training the model will fail. Thus, `resnet_forward` will still work, but `polynomial` example will fail.

### CMake
Make sure you have cmake or install it by (tested on cmake==3.22.1)
```bash
sudo apt install cmake
```

### HDF5
EOS tables (SFHo by default) require an hdf5 installation. If missing, get it via:
```
sudo apt-get install libhdf5-dev
```
Then you need to provide paths to the hdf5 libraries in the `Makefile`. If yours differ from default, edit the following variables:
```
HDF5PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
HDF5INCS=-I/usr/include/hdf5/serial
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

## Troubleshoot
### `No CMAKE_CXX_COMPILER could be found`

CMake can't find your C compiler. Either check your GCC path or if you are on Ubuntu, run:
```bash
sudo apt-get update && sudo apt-get install build-essential
```