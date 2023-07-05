# LagrProb 
LagrRadBC is lagrangian model written on Fortran 90 to simulate bottom boundary layer formation  due to the scavenging process.

LagRadBC_main.f90 is the source file. It was tested under Windows 10 using Intel Fortran Compiler with VS2008 shell and MSMPI support. 
Module IFPORT is used for the DRAND () function that allows double precision random number generation.

For compilation use Default Real(8), Integer(8).
Before running create directory 'out\' for storing output files

File LagrRadBC.inp - text input file with parameters and decriptions
