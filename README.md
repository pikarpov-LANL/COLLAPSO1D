# 1D Core-Collapse Supernovae Code

This is a 1D lagrangian code to explore CCSN modeling. For progenitors, it takes KEPLER generated [data](https://2sn.org/stellarevolution/) of Alex Heger & Stan Woosley.

Turbulence is treated through mixing length theory (MLT), BHR, and Machine Learning (ML) based models.

## Quick Start
First we need to compile a script to read the stellar progenitor data. In our case, we are using Stan Woosley's datasets:
```shell
cd 1dcollapse
gfortran read_woosley.f -o a.out
```
when running `a.out`, it will ask for the innermost zone mass and then evolves the masses based on prescriptions that worked well a long time ago. Zone mass should be <1. Lastly, it will ask for the final binary file name that will be used as an input for the simulation. The name of the input and output files should not contain <ins>underscores</ins>.
```shell
./a.out
```
Next we need to edit `inlahyc` to set the input and output file names, as well as some other initial conditions. See the details in section [Model Setup](#model-setup-inlahyc) below. Lastly, compile the code itself, including the EOS:
```shell
cd 1dcollapse
gfortran -O 1dburn.f ocean.f nse4c.f sleos.f -o RunName
```

## Model Setup (inlahyc)
here is an example of how `inlahyc` looks like:
```
Test
TestOut
1
.0001 1.
1.0 .25
1 4 2.d0
1206,1.d-5,1,0.01,1
1,1,1.0,1.0
```
Here is a description of each line:
| Line | Definition |
| -------- | ---------- |
| 1 | input file name |
| 2 | output file name |
| 3 | dump # to read the input file |
| 4 | initial timestep and total timestep |
| 5 | artificial viscosity values |
| 6 | external force (not default) <br/> equation of state option <br/> if > 1, include core mass |
| 7 | delp (not default) <br/> nups - number of steps per luminosity output <br/> damping term <br/> damping zones below this number |
| 8 | iflxlm - flux limiter option <br/> capture rate option <br/> changing the nuclear potential energy - shouldn't be altered <br/> yefact - not used|

## [Notes](https://www.overleaf.com/read/pgsnmxgdjkrq)

Notes on ML subgrid turbulence model implementation within this code can be found on [Overleaf](https://www.overleaf.com/read/pgsnmxgdjkrq).
