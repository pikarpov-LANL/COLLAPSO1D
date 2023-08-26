
# EOS Tables

COLLAPSO1D supports a EOS tables through the integrated [EOSdriver](https://github.com/evanoconnor/EOSdriver) by Evan O'Connor. By default, the code is set up with [SFHo](http://adsabs.harvard.edu/abs/2012arXiv1207.2184S) tables in mind ([download link](https://su.drive.sunet.se/index.php/s/FQkikyGcRnHTZNL)). For a full list of supported tables, their format, and how they are interpolated upon read-in, please refer to Evan's [official website](https://ttt.astro.su.se/~eoco/eos.html) and the corresponding GitHub repos. You will be able to download the formatted tables from there as well.

## Installation

!!! Tip
    COLLAPSO1D already inlcudes EOSdriver, so you don't need to perform this step.

Here are the instructions to install standalone EOSdriver.

1. Install hdf5 on linux (make sure it is serial, which is default)
    ```
    sudo apt-get install libhdf5-dev
    ```
2. Get EOSdriver
   ```
   git clone https://github.com/evanoconnor/EOSdriver.git
   ```
3. make driver executables
    1. Edit `HDF5LIBS` and `HDF5INCS` in `make.inc` with correct hdf5 paths, e.g.
        ```
        HDF5LIBS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5 -lz
        HDF5INCS=-I/usr/include/hdf5/serial
        ```
    2. type `make`


## EOSdriver Variables

| nuc_eos_full | Units                    | Intent | Description                                        |
| :----------- | :----------------------- | :----: | :------------------------------------------------- |
| xrho         | g/cm^3                   | inout  | density                                            |
| xtemp        | MeV                      | inout  | temperature                                        |
| xye          | number fraction / baryon |   in   | electron fraction                                  |
| xenr         | erg/g                    | inout  | energy                                             |
| xprs         | dyn/cm^2                 |  out   | pressure                                           |
| xent         | k_B / baryon             | inout  | entropy                                            |
| xcs2         | cm^2/s^2                 |  out   | speed of sound squared (not relativistic)          |
| xdedt        | $erg/g/MeV$              |  out   | $C_{\nu}$                                          |
| xdpderho     | $dynes \; g/cm^2/erg$    |  out   | $dP/d\epsilon$ at constant $\rho$                  |
| xdpdrhoe     | $dynes \; cm^3/cm^2/g$   |  out   | $dP/d\rho$ at constant $\epsilon$                  |
| xxa          | mass fraction            |  out   | $\alpha$ particle mass fraction                    |
| xxh          | mass fraction            |  out   | average heavy nuclus mass fraction                 |
| xxn          | mass fraction            |  out   | neutron mass fraction                              |
| xxp          | mass fraction            |  out   | proton mass fraction                               |
| xabar        | A                        |  out   | average heavy nucleus mass number                  |
| xzbar        | Z                        |  out   | average heavy nucleus atomic number                |
| xmu_e        | $MeV$ (or / baryon?)     |  out   | electron chemical potential                        |
| xmu_n        | $MeV$                    |  out   | neutron chemical potential                         |
| xmu_p        | $MeV$                    |  out   | proton chemical potential                          |
| xmuhat       | $MeV$                    |  out   | mu_n - mu_p                                        |
| keytemp      | 0,1,2,3                  |   in   | primary value                                      |
| keyerr       |                          |  out   | error output; should be 0                          |
| rfeps        |                          |   in   | root finding relative accuracy, set around 1.0d-10 |


## I/O of EOSdriver

| keytemp | description                                     |
| ------- | ----------------------------------------------- |
| `0`     | coming in with rho,eps,ye (solve for temp)      |
| `1`     | coming in with rho,temperature,ye               |
| `2`     | coming in with rho,entropy,ye (solve for temp)  |
| `3`     | coming in with pressure,temp,ye (solve for rho) |


## EOSdriver vs. COLLAPSO1D

Comparison of the EOSdriver with COLLAPSO1D variables:

| nuc_eos_full | Units                    | COLLAPSO1D | Units            |
| :----------- | :----------------------- | :--------- | :--------------- |
| xrho         | $g/cm^3$                 | rho(k)     | 2.d6 $g/cm^3$    |
| xtemp        | $MeV$                    | temp(k)    | 1.d9 $K$         |
| xye          | number fraction / baryon | ye(k)      | same             |
| xenr         | $erg/g$                  | u(k)       | uergg            |
| xprs         | $dyn/cm^2$               | pr(k)      | 2.d22 $dyn/cm^2$ |
| xent         | $k_B / baryon$           | u2         | sfac             |
| xcs2         | $cm^2/s^2$               | vsound     | 1.d8 $cm/s$      |
| xdedt        | $erg/g/MeV$              | dusl(?)    | does it matter?  |
| xdpderho     | $dynes \; g/cm^2/erg$    |            | does it matter?  |
| xdpdrhoe     | $dynes \; cm^3/cm^2/g$   |            | does it matter?  |
| xxa          | mass fraction            | xalpha(k)  | same             |
| xxh          | mass fraction            | xheavy(k)  | same             |
| xxn          | mass fraction            | xn(k)      | same             |
| xxp          | mass fraction            | xp(k)      | same             |
| xabar        | A                        | abar(k)    | same             |
| xzbar        | Z                        | zbar(k)    | does it matter?  |
| xmu_e        | $MeV$                    | xmue(k)    | 1.d9 $K$         |
| xmu_n        | $MeV$                    |            | does it matter?  |
| xmu_p        | $MeV$                    |            | does it matter?  |
| xmuhat       | $MeV$                    | xmuhat(k)  | 1.d9 $K$         |
| keytemp      | 0,1,2,3                  |
| keyerr       |                          |
| rfeps        |                          |


## `slwrap` variables

Detailed description of COLLAPSO1D variables from `#!fortran subroutine slwrap`

| COLLAPSO1D | Description                                                            |
| :--------- | :--------------------------------------------------------------------- |
| inpvar(1)  | temperature (may not need the others)                                  |
| k          | iterator                                                               |
| yek        | electron franction (Ye)                                                |
| brydns     | density in sw units                                                    |
| pprev      | funky proton fraction (probs not needed)                               |
| psl        | pressure                                                               |
| usl        | energy or entropy (probs the latter, but depends on which eos is used) |
| dusl       | derivative of above?                                                   |
| gamsl      | effective gamma for sound speed (won't be needed)                      |
| etak       | degeneracy (can get it from electron chemical potential/T)             |
| xpk        | mass fractions                                                         |
| xnk        |
| xak        |
| xhk        |
| yehk       | Ye mass fraction                                                       |
| abark      | abar                                                                   |
| xmuh       | mu hat                                                                 |
| stot       | entropy                                                                |
