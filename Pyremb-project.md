src_dir: ./
output_dir: ./docs
docmark: >
predocmark: !
predocmark_alt: #
project: PYREMB
summary: A Subsystem DFT package able to combine arbitrary - periodic and non-periodic -  molecular systems
author: Christian Schwermann
author_description: PhD student in the [Doltsinis group](https://www.uni-muenster.de/Physik.FT/Forschung/agdoltsinis/index.html) at [WWU Muenster](https://www.uni-muenster.de)
email: c.schwermann@wwu.de
print_creation_date: true
year: 2020 Christian Schwermann <br/> ALL RIGHTS RESERVED
graph: true

The periodic embedding code PEREMB provides an interface between - in principle - arbitrary DFT codes to facilitate subsystem DFT calculations with isolated, periodic and mixed isolated-periodic systems.

The program uses provided densities and atomic positions to calculate the non-additive kinetic potential, the non-additive exchange-correlation potential, the Hartree potential and the nuclear potential. Those potentials together yield the embedding potential, which can be printed out and used in DFT calculations with external programs.

The source code is available at [https://github.com/cschwermann/pyremb](https://github.com/cschwermann/pyremb).

**Note** The program is incomplete and not in a usable state. 

## Features

 * The subroutines calculating the individual contributions to the embedding potential are clearly structured to be easily accessible in external programs by including specific modules. In particular, the subroutines could be called individually if not all contributions are needed (such as the electrostatic terms).
 * An interpolation of all data belonging to a "subsystem" type onto the arbitrary grid of another "subsystem" type was implemented to allow the combination of arbitrary subsystems on arbitrary grids. This first implementation uses a basic Shepard algorithm.
 * Integration subroutines for arbitrary grids were implemented, aiming at energy evaluations by integrating input potentials.
 * In order to allow the usage of GGA functionals, the density functional subroutines were generalized to in principle call the respective libXC subroutines for GGA functionals, if a density gradient is provided. Calls to the libXC energy subroutines were also implemented.    


## Non-Features

While the individual subroutines were tested and found to be working, the current implementation of PEREMB is incomplete. In particular, 

 * no python interface exists, as f2py is incompatible with the used Fortran types. Tests with F90Wrap were inconclusive.
 * The Shepard interpolation routine works well for small arrays, but scales quadratically with the number of grid points, which makes application to actual molecular systems unfeasible. 
 * The energy evaluation is incomplete. The most difficult part is the calculation of inter-subsystem nuclei-nuclei interaction, as this might require an immensely fine grid or some kind of Ewald summation. Attempts to implement an Ewald summation method failed due to convergence difficulties.
 *Also, the libXC interface is not entirely clear as to what exactly is the output of the energy subroutines. The output is an array containing a value for each grid point; however, integration as well as integration of the product of the array with the density did not yield any meaningful value -- not even scaled by factors such as multiples of pi or the cell volume.  


## Outlook 

Thus, further work is needed to develop a fully functional PEREMB version which enables FDE calculations with BAND:

 * Python wrappers might be created using "forpy". Otherwise, they have to be written manually (which should be possible).
 * More efficient implementations  of the Shepard interpolation exist  (e.g. within the NAG library), which reduce the quadratic scaling by only considering close-lying grid points. Otherwise, more sophisticated interpolation routines might be needed.
 * easily yield all the elctrostatic terms. 
 * libXC could be entirely avoided by using XCfun; however, the documentation also does not seem very clear.
