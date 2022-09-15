
## General:

`netwinv4`, `*.atb`, `#.inc` are all eos/nse files

Note that there are some warnings in the code.  Most due to one of the old EOSs used that Nadyozhin and Blinnikov developed.

## Setup Code 

### setup
Compile data format conversion: 

```shell
gfortran readalexnew.f -o a.out
```

When it runs, it asks for the innermost zone mass and then evolves the masses based on prescriptions that worked well a long time ago.  

### run15f1

Core-Collapse Code:

| File     | Description             |
| -------- | ----------------------- |
| 1dburn.f | main core-collapse code |
| ocean.f  | the low-density eos     |
| nse4c.f  | the nse code            |
| sleos.f  | the Swesty-Lattimer EOS |

## to compile:

```shell
gfortran -O 1dburn.f ocean.f nse4c.f sleos.f -o goodname
```

It reads inlahyc that has the 

1. input file
2. output file
3. dump number to read the input file
4. initial timestep and total timestep
5. artificial viscosity values
6. options: 
   1. external force (not used in this code)
   2. equation of state option
   3. if >1, include core mass
7. options: number of zones
   1. delp (not used in current version)
   2. nups (number of steps per luminosity output)
   3. damping term
   4. damping zones below this number
8. options: 
   1. iflxlm (flux limiter option)
   2. capture rate option
   3. changing the nuclear potential energy (why?)
   4. yefact (not used)

readoutput codes read the binary output and put into ascii
