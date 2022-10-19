This is the main core-collapse code. It contains the following subroutines:

| Subroutine | Description                                                                                                                                                                              |
| :--------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| hydro      | advances the system of hydro equations by one time step (print\_nuloss variable is a flag to avoid printing nu losses the 2nd time)                                                      |
| artvis     | updates the q-value for artificial viscosity                                                                                                                                             |
| coulomb    | computes Coulomb corrections as given in Shapiro and Teukolsky. p. 31 (2.4.9) and (2.4.11)                                                                                               |
| density    | calculates the density using the continuity equaiton                                                                                                                                     |
| energ      | computes the change in internal energy                                                                                                                                                   |
| eosflg     | determines what kind of eos to use depending on physical conditions:<br> eosflg = 1: freeze-out, just Ocean's eos + Coul corr. <br>eosflg = 2: NSE with Raph's routines, + Ocean eos + Coul <br>eosflg = 3: Swesty's eos and if e-neutrinos or x-neutrinos are trapped 
| eospg      | computes the pressure and sound speed for all particles on a list assuming a perfect gas equation of state                                                                               |
| eospgr     | computes the pressure, and sound speed according to an equation of state that includes gas and radiation pressure. This part of the code has not been debugged and most likely won't run |
| eos3       | compute pressure and temperatures with the Ocean eos assuming NSA                                                                                                                        |
| forces     | computes the force on the cells that need to have their forces evaluated. Neutrino diffusion is also performed here.                                                                     |
| gravity    | adds the gravitational force on the particles. This will have the option of a neutron star core or it can allow the particles to make up the core.                                       |
| mmw        | sets the mean molecular weight of the gas assuming complete ionization                                                                                                                   |
| nserho     | figures out the NSE eq. assuming that yp and yn were previously known at different density and ye, but \textbf{same temperature}                                                         |
| nsetemp    | figures out the NSE eq. assuming that yp and yn were previously known at the \textbf{same density and ye, but different temperatures}                                                    |
| nuabs      | computes the neutrino absorption by nucleons (all neutrino energies are in MeV)                                                                                                          |
| nuann      | computes the rate of neutrino anti- neutrino annihilation into e+/e- pairs (see Goodman, Dar, Nussinov, ApJ 314 L7)                                                                      |
| nubeta     | treats cases where beta eq. has occurred.In beta eq.: munue(beta)=mue-muhat, so we compute Ynue(munue(beta)) and unue(munue(beta)) assuming thermal distribution at matter temperature, compare with actual Ynue and unue, and move things in the right direction
| nucheck    |
| nuconv     |
| nudiff     |
| nuecap     |
| nuinit     |
| nulum      |
| nupp       |
| nupress    |
| nuscat     |
| nusphere   |
| nuwork     |
| pppb       |
| preset     | sets up all quantities needed before starting a simulation                                                                                                                               |
| rmp        |
| amp        |
| rootemp1   | computes temperature using a Newton-Raphson procedure found in the numerical recipes, p.254                                                                                              |
| rootemp2   |
| rootemp3   |
| rootemp4   |
| rooteta    |
| slwrap     | wrapper routine for the swesty-lattimer eos                                                                                                                                              |
| readini    | reads initial conditions                                                                                                                                                                 |
| printout   | prints out all the results                                                                                                                                                               |
| integrals  | numerical approximations to the fermi integrals (Takahashi et al, 1978)                                                                                                                  |
| epcapture  |
| step       | integrates the system of equations using a Runge-Kutta-Fehlberg integrator of second order. Particles are allowed to have individual time-steps. All particles are synchronized every dtime at which time the subroutine is exited |
| unit       | computes the transformation between the physical units (cgs) and the units used in the code, and the value of the physical constant in code units                                        |
| burn       | nuclear network subroutine - uses a 14 elements alpha network from fkt                                                                                                                   |
| derivn     | calculates time derivatives of abundances i.e. ff(i) is a system of differential equations                                                                                               |
| epsb       | change of total nuclear binding energy/mass units erg/g                                                                                                                                  |
| genpar     |
| lequb      |
| matr       |
| newab      |
| rates      | calculates rates <br> rrat(i) forward rates (captures 1-14) 15:CC 16:CO 17:OO <br> pf(7,i) coefficients <br> rlam(i) backward rates (1-14 photo-disintegrations) <br>pb(7,i) coefficients|
| rrate      | program reads coefficients of rates, pf forward rates, pb backward rates (=photodisintegrations)                                                                                         |