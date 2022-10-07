
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
    1. Edit `HDF5LIBS` and `HDF5INCS` with correct hdf5 paths, e.g.
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
* unit - units specific based on eos

New subroutines:
* rootemp5
* sfhowrap

### Variable compatability

| EOSdriver    | EOSdriver Units             | COLLAPSO1D | COLLAPSO1D Units |
| :----------- | :-------------------------- | :--------- | :--------------- |
| pointsrho    | dimensionless               |
| pointstemp   | dimensionless               |
| pointsye     | dimensionless               |
| logrho       | $log_{10}(\rho[g/cm^3])$    |
| logrho       | $log_{10}(T[MeV])$          |
| logrho       | number fraction             |
| Abar         | A                           |
| Zbar         | Z                           |
| Xa           | mass fraction               | xak
| Xh           | mass fraction               | xhk
| Xn           | mass fraction               | xpk
| Xp           | mass fraction               | xpk
| cs2          | $cm^2/s^2$                  |
| dedt         | $erg/g/MeV$                 |
| dpderho      | $dynes \; g/cm^2/erg$       |
| dpdrhoe      | $dynes \; cm^3/cm^2/g$      |
| energy_shift | $erg/g$                     |
| entropy      | $k_B/baryon$                |
| gamma        | dimensionless               |
| logenergy    | $log_{10}(\epsilon[erg/g])$ |
| logpress     | $log_{10}(P[dynes/cm^2])$   |
| mu_e         | $MeV/baryon$                |
| mu_p         | $MeV/baryon$                |
| mu_n         | $MeV/baryon$                |
| muhat        | $MeV/baryon$                | xmuh
| munu         | $MeV/baryon$                |

And here is a detailed description of COLLAPSO1D variables

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

### Questions

* What is the actual difference between EoS 3 and 4? 
  * Former has two if statements in the code
  * !-- switch back to internal energy variable of state 
* Need help testing correctness of the table integration (units and etc.)
  * feed a range of densities [up to $10^{14}$] into EOS: check pressure, chemical potential, and other.
