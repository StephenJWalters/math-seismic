The four folders here contain code for analytic and numerical solutions to the seismic equations in continuous media.
The code here was used to produce the four figures in the paper "" which may be found at ... 
The Fortran code was used for all forward integration of seismic pulses, and the data stored in text files. 
This data was then read by the Matlab code to produce the figures. 
The Fortran code was compiled using the PGI Fortran compiler, but adaptation to GFortran or other is not difficult. 
Adaptation of the Matlab code to Python and matplotlib also should be straightforward, if desired.
In each case the files should be placed in a folder, and a subfolder "data" should be created.
The pgi compile commands are placed in comments at the start of each .f90 fortran file.
The support module provides subroutines which support the main file.
the mmul module provides a fast matrix multiplication routine which runs on the GPU. This was adapted from the PGI Compiler user guide. It requires that any matrices used have number of rows and columns in whole multiples of 16. This avoids "if" statments and make the routine simpler and faster, for a small cost in flexibility.
