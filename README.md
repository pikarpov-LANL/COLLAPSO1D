# 1D Core-Collapse Supernovae Code

This is a 1D lagrangian code to explore CCSN modeling. For progenitors, it takes KEPLER generated [data](https://2sn.org/stellarevolution/) of Alex Heger & Stan Woosley.

Turbulence is treated through mixing length theory (MLT), BHR, and Machine Learning (ML) based models.

## Quick Start
First we need to compile a script to read the stellar progenitor data. In our case, we are using Stan Woosley's datasets:
```shell
cd read_data
gfortran read_woosley.f -o a.out
```
when running `a.out`, it will ask for the innermost zone mass and then evolves the masses based on prescriptions that worked well a long time ago. Zone mass should be <1. Lastly, it will ask for the final binary file name that will be used as an input for the simulation. The name of the input and output files should not contain <ins>underscores</ins>.
```shell
./a.out
cp InputName ../1dcollapse/
```

Next we need to compile the code itself, including the EOS:
```shell
cd 1dcollapse
gfortran -O 1dmlmix.f90 ocean.f90 nse4c.f90 sleos.f90 -o model
```


## [Notes](https://www.overleaf.com/read/pgsnmxgdjkrq)

Notes on ML subgrid turbulence model implementation within this code can be found on [Overleaf](https://www.overleaf.com/read/pgsnmxgdjkrq).
