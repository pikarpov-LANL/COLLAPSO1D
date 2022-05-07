
## General:

netwinv4, *.atb, #.inc are all eos/nse files

Note that there are some warnings in the code.  Most due to one of the old EOSs used that Nadyozhin and Blinnikov developed.

## Setup Code
readalexnew.f 
compile gfortran readalexnew.f -o a.out
When it runs, it asks for the innermost zone mass and then evolves the masses based on prescriptions that worked well a long time ago.  

run15f1:
Core-Collapse Code:
1dburn.f is the main core-collapse code
ocean.f is the low-density eos
nse4c.f is the nse code
sleos.f is the Swesty-Lattimer EOS

## to compile:
gfortran -O 1dburn.f ocean.f nse4c.f sleos.f -o goodname

It reads inlahyc that has the 
input file
output file
dump number to read the input file
initial timestep and total timestep
artificial viscosity values
options: external force (not used in this code), equation of state option, if >1, include core mass
number of zones, delp (not used in current version), nups (number of steps per luminosity output), damping term, damping zones below this number
iflxlm (flux limiter option), capture rate option, changing the nuclear potential energy (why?), yefact (not used)

readoutput codes read the binary output and put into ascii
