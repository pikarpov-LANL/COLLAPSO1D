# 1D Core-Collapse Supernovae Code (Legacy)

This is a 1D lagrangian code to explore CCSN modeling. For progenitors, it takes KEPLER generated [data](https://2sn.org/stellarevolution/) of Alex Heger & Stan Woosley.

Turbulence is treated through mixing length theory (MLT), BHR, and Machine Learning (ML) based models.

## Quick Start
First we need to compile a script to read the stellar progenitor data. In our case, we are using Stan Woosley's datasets:
```bash
cd setup/
gfortran read_woosley.f -o a.out
./a.out
```
when running `a.out`, it will ask for the *innermost zone mass* and then evolves the masses based on prescriptions that worked well a long time ago. Zone mass should be <1. Lastly, it will ask for the final *binary file name* that will be used as an input for the simulation. The name of the input and output files should not contain <ins>underscores</ins>.
```shell
./a.out
mv InputName ../legacy/1dcollapse/
```

Next we need to compile the code itself, including the EOS:
```bash
cd legacy/1dcollapse
gfortran -O 1dmlmix.f90 ocean.f90 nse4c.f90 sleos.f90 -o model
```

Lastly, update `inlahyc` input file to read `InputName` and run the model:
```bash
./model
```

## Data Defaults

Input file, i.e. progenitor mass, should be edited directly in the `read_woosley.f`. There are 3 available from Heger & Woosley, 2000:

1. s15presn 
2. s20presn (Default)
3. s25presn

For the prompts:

| Prompt              | Value  |
| ------------------- | ------ |
| Innermost zone mass | 0.0005 |
| Binary file name    | Data   |

This will result in `grid_szie=1706`

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
| 4 | dump time interval and total timestep |
| 5 | artificial viscosity values |
| 6 | external force (not default) <br/> equation of state option <br/> if > 1, include core mass |
| 7 | # of cells in the grid *(hint: copy it from `./a.out` output's last line)* <br/> delp (not default) <br/> nups - number of steps per luminosity output <br/> damping term <br/> damping zones below this number |
| 8 | iflxlm - flux limiter option <br/> capture rate option <br/> changing the nuclear potential energy - shouldn't be altered <br/> yefact - not used|

## [Notes](https://www.overleaf.com/read/pgsnmxgdjkrq)

Notes on ML subgrid turbulence model implementation within this code can be found on [Overleaf](https://www.overleaf.com/read/pgsnmxgdjkrq).
