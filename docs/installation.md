# Installation

COLLAPSO1D supports `gfortran` and Intel's Classic Fortran `ifort`. The latter runs ~70% faster, but requires a more involved installation.

## Defendencies

### PyTorch
!!! Warning "PyTorch 2.X has <ins>not</ins> been tested."
    
To install the latest PyTorch 1.X for <ins>CPU</ins>, follow the official [instructions](https://pytorch.org/get-started/locally/). I would highly recommend installing it in a dedicated conda environemnt. As an example, here is how to create one and install PyTorch:
```
conda create -n py310 python=3.10
conda activate py310

# check the pytorch installation instructions!
pip install torch==1.11.0+cpu --extra-index-url https://download.pytorch.org/whl/cpu
```
!!! Warning
    `torch>=1.12.0` will work for inferencing, hence for the CCSN code will be fine, but loading and training the model will fail. Thus, `resnet_forward` will still work, but `polynomial` example will fail.

!!! Tip "Speed-up `make` commands"
    `CMAKE_PREFIX_PATH` is a variable set in `Makefile` pointing to the location of Torch cmake config files. By default, it is determined by importing torch which is fairly slow, but it can also be set explicitely. To get the torch `{Location}`:
    ```bash
    pip show torch
    ```
    Then set `CMAKE_PREFIX_PATH={Locaton}/torch/share/cmake` in `Makefile`.

### CMake
Make sure you have `cmake` or install it by (tested on `cmake==3.22.1`)
```bash
sudo apt install cmake
```

### EOS Tables

You will also need to download the SFHo EOS Table and put it in the executable directory.
```bash
cd project/1dmlmix
wget https://su.drive.sunet.se/index.php/s/FQkikyGcRnHTZNL/download/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.3_20190605.h5
```
Click to learn more about the [EOS Tables in COLLAPSO1D](eosdriver.md).

---

## GFORTRAN
Both the main CCSN code and the PyTorch wrapper can be compiled with `gfortran >= 9.4.0`. If missing, install it via:
```
sudo apt install gfortran
```

In the `Makefile`, set
```
COMPILER=gfortran
```

### HDF5
EOS tables (SFHo by default) require an `hdf5` installation. If missing, get it via:
```
sudo apt-get install libhdf5-dev
```
Then you need to provide paths to the `hdf5` libraries in the `Makefile`. If yours differ from default, edit:
```
HDF5PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
HDF5INCS=-I/usr/include/hdf5/serial
```
Add the following to you `~/.bashrc` and then `source ~/.bashrc`:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
```

### Test Installation
Run the following to test your installation with `gfortran`.
```bash
make test
```

---

## Intel `ifort`

I would highly recommend installing Intel's oneAPI through `APT`. Even though there are many options listed, a good number of them are bugged. In addition, we need to install a full basekit, instead of only the Fortran compiler, as there are issues with individual packages.

Please follow [APT Installation Guide](https://www.intel.com/content/www/us/en/develop/documentation/installation-guide-for-intel-oneapi-toolkits-linux/top/installation/install-using-package-managers/apt.html) on Intel's website to install `intel-hpckit`. After the pre-installation steps, you should run the following:
```bash
sudo apt install intel-hpckit
```
Next you will need to add executables to `PATH`. The easiest way is to run:
```bash
source /opt/intel/oneapi/setvars.sh
```
Add this line to your `~/.bashrc` to avoid re-running the above initialization on every start-up.

In the `Makefile`, set
```
COMPILER=ifort
```

!!! Warning
    The `readout.f90` to convert unformatted binary output to readable text will still be compiled with `gfortran` (hardcoded in the `Makefile`), since `ifort` is ~x20 slower at parsing unformatted binary files for some reason. These slow downs have no effect of the actual CCSN calculation, hence `ifort` remains vastly superior for the main code.

### HDF5 with `ifort`

To use COLLAPSO1D with EOS Tables, we will also need to re-compile the `hdf5` library with `ifort`. For that, we need to do a custom installation of `hdf5`. I tested on hdf5-1.12.2.

1. Download the latest [HDF5 Source Code](https://www.hdfgroup.org/downloads/hdf5/source-code/).
2. Configure
   ```
   FC=ifort FCFLAGS="-O3" CC=icc CFLAGS="-O3" ./configure --enable-fortran --enable-shared
   ```
3. Install
    ```
    make install
    ```
4. Check 
   ```
   make check-install
   ```
5. Add library path
   ```
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{path_to_hdf5}/hdf5/lib
   ```
    Add the above to your `~/.bashrc` to avoid re-running the above initialization on every start-up

Full installation instructions detailing all of the available commands and flag can be found inside your hdf5 directory in `release_docs/INSTALL`.

Lastly, don't forget to change in the COLLAPSO1D's `Makefile`.
```
HDF5PATH={path_to_hdf5}/hdf5/lib
HDF5INCS=-I{path_to_hdf5}/hdf5/include
```

---

## Troubleshoot
### `No CMAKE_CXX_COMPILER could be found`

CMake can't find your C compiler. Either check your GCC path or if you are on Ubuntu, run:
```bash
sudo apt-get update && sudo apt-get install build-essential
```

### Version Compatibility Issues

If you are having trouble getting the dependencies and COLLAPSO1D to work, please check the GitHub compilation test. The script sets up a clean Linux environment with explicitely defined package versions. It checks if COLLAPSO1D compiles correctly with `gfortran` on every push to the repo.

??? Quote "compilation_test.yml"
    ```yaml
    --8<-- ".github/workflows/compilation_test.yml"
    ```