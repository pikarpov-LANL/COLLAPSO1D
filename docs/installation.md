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
conda create -n torch_cpu python=3.10
conda activate torch_cpu

# check the pytorch installation instructions!
pip3 install torch torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cpu
```

Lastly, we need to make sure all compilers link correctly.

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

## Troubleshoot
### `No CMAKE_CXX_COMPILER could be found`

CMake can't find your C compiler. Either check your GCC path or if you are on Ubuntu, run:
```bash
sudo apt-get update && sudo apt-get install build-essential
```