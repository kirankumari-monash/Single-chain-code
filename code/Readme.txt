==============  Compiling The Programme ====================
To compile the programme type
$ make sens
This will compile a single processor (multiple runs per directory)

=============== Status of the code =========================


Major Implementation: 
   1) Fraenkel spring 
   2) Modified LJ potential
   3) double precision using flag
   

Minor Modification:
   1) Corrected the array size of maxminev subroutine of utils.f90 (not properly defined earlier)
  
============== Running The Programme ========================

Input parameter description file: intputc.dat

This file takes the following parameters as input. Use the format
given in the sample file.


SpType (int) 
   is the spring type
   1 Hookean
   2 FENE
   3 ILC
   4 WLC

NBeads (int) 
   is the number of beads

h* (float) 
   is the HI parameter 

z* (float)
   is the dimensionless energy parameter in the Gaussian potential

d* (float)
    is the dimensionless radius parameter in the Gaussian potential

L0* (float)
    is the dimensionless maximum length of a spring
    for finitely extensible springs.  In terms of Birds
    notation. L0* = sqrt(b) 

gdot* (float)
    is the dimensionless shear rate

emax (float)
    this is the maximum strain units upto which the simulation is to be run
    Typical value = 10	  
    emax = gdot* x total_time

Nsamples (int)
    Number of samples to be taken in the production run, till end of emax
    Typical value = 100

ndelts (int)
     Number of time step sizes.
     Normally simulation is run over several time step sizes, and a 
     time step convergence to zero is obtained.  
     Parameters for each time step value is mentioned below
     Typical value 1 or 3

dtseq desne nblock ntot tol
     Parameters for each time step size, the number of lines corresponds
     to the number of ndelts mentioned above

dtseq (float)
     dimensionless timestep size for equilibrium
     typical value = 1

dtsne  (float)
     dimensionless time step for production (shear, elogation)

nblock  (int)
     number of trajectories per processor.
     Several processors can run in the same directory.  Each
     processor reads the number of completed trajectories from disk
     finds the remaining trajectories to be done and carries out
     the lesser of nblock or remaining trajectories.
     typical value = 1000

ntot (int)
     total number of trajectories.
     Each processor computes the completed trajectory (from gavgs.NN)
     and finds the remaining number of trajectories for each run.
     Typical value = 1000000

tol (float)
     tolerence for semi-implicit predictor corrector
     Typical value = 0.001



============== Interpreting results ========================

Results for each time step result is stored separately in a file called
output.NN,  where NN is a number from 01, 02, ... upto ndelts (number
of time steps)

The columns in the file are self explanatory, as given in the third
commented line.  # is treated by gnuplot as a comment line.

R^2 is the squared end to end distance
Rg^2 is the squared Radius of gyration
psi1 is the first normal stress coefficient
etap is the shear viscosity coefficient
