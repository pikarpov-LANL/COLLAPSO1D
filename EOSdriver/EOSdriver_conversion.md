
# COLLAPSO1D with SFHo EoS tables

### Installation

1. Install COLLAPSO1D (follow setup instructions in README)
   ```
   git clone https://github.com/pikarpov-LANL/COLLAPSO1D.git
   ```
2. Install hdf5 on linux (make sure it is serial, which is default)
    ```
    sudo apt-get install libhdf5-dev
    ```
3. Get EOSdriver
   ```
   git clone https://github.com/evanoconnor/EOSdriver.git
   ```
4. make driver executables
    1. Edit `HDF5LIBS` and `HDF5INCS` in `make.inc` with correct hdf5 paths, e.g.
        ```
        HDF5LIBS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5 -lz
        HDF5INCS=-I/usr/include/hdf5/serial
        ```
    2. type `make`

### Where to integrate?

Subroutines :
* energ
* eosflg
* eos3 (or create eos5?)
* printout - `s` and `uint` are swapped for eos4
  * probably want all entropy (s) as primary, but need to check
* step - a few `u` adjusting `if` statements based on `ieos`
* unit - units specific based on eos (check `!--3b) Conversion to the Swesty-Lattimer units from code units:  `)

New subroutines:
* rootemp5
* sfhowrap

### Variable compatability

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

keytemp: 
* 0 -> coming in with rho,eps,ye (solve for temp)
* 1 -> coming in with rho,temperature,ye
* 2 -> coming in with rho,entropy,ye (solve for temp)
* 3 -> coming in with pressure,temp,ye (solve for rho)


And here is a detailed description of COLLAPSO1D variables from `slwrap`

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

Comparison of the EOSdriver with COLLAPSO1D variables:

| nuc_eos_full | Units                    | COLLAPSO1D | Units             |
| :----------- | :----------------------- | :--------- | :---------------- |
| xrho         | $g/cm^3$                 | rho        | 2.d-6 $g/cm^3$    |
| xtemp        | $MeV$                    | inpvar(1)  | 1.d-9 $K$ (?)     |
| xye          | number fraction / baryon | yehk       | same              |
| xenr         | $erg/g$                  | usl(u)(?)  | ergmev*uergg/avo  |
| xprs         | $dyn/cm^2$               | psl        | 2.d-22 $dyn/cm^2$ |
| xent         | $k_B / baryon$           | stot       |
| xcs2         | $cm^2/s^2$               | vsound     | 1.d-8 $cm/s$      |
| xdedt        | $erg/g/MeV$              | dusl(?)    |
| xdpderho     | $dynes \; g/cm^2/erg$    |
| xdpdrhoe     | $dynes \; cm^3/cm^2/g$   |
| xxa          | mass fraction            | xak        | same              |
| xxh          | mass fraction            | xhk        | same              |
| xxn          | mass fraction            | xnk        | same              |
| xxp          | mass fraction            | xpk        | same              |
| xabar        | A                        | abark      |
| xzbar        | Z                        | zbark      |
| xmu_e        | $MeV$                    | xmue       |
| xmu_n        | $MeV$                    |
| xmu_p        | $MeV$                    |
| xmuhat       | $MeV$                    | xmuh       | $\eta T$, code    |
| keytemp      | 0,1,2,3                  |
| keyerr       |                          |
| rfeps        |                          |

| Name     | Value     | Comment       |
| -------- | --------- | ------------- |
| ggcgs    | 6.67e-8   |
| avo      | 6.02e23   |
| aradcgs  | 7.565e-15 |
| boltzk   | 1.381e-16 |
| hbar     | 1.055e-27 |
| ccgs     | 3e10      |
| emssmev  | 0.511e0   |
| boltzmev | 8.617e-11 |
| ergmev   | 6.2422e5  |
| sigma1   | 9d-44     |
| sigma2   | 5.6d-45   |
| c2cgs    | 6.15e-4   |
| c3cgs    | 5.04e-10  |
| fermi    | 1d-13     |
| uergg    | 1d16      | ergs per gram |

### Questions

* What is the actual difference between EoS 3 and 4? 
  * Former has two if statements in the code
  * !-- switch back to internal energy variable of state 
* Need help testing correctness of the table integration (units and etc.)
  * feed a range of densities [up to $10^{14}$] into EOS: check pressure, chemical potential, and other.
* What about proton and neutron chemical potential? Should I use them anywhere?
* Also, any use for `xdpderho` & `xdpdrhoe`?

call integrals


