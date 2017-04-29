RUSTICO   Rapid foUrier STatIstics COde

Author: Hector Gil Marin

Date: 1st May 2017

email: hector.gil.marin@gmail.com or hector.gilmarin@lpnhe.in2p3.fr

====Compilation====
icc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c cubature.c ps_write.c order_algorithm.c -O3   -lm -openmp -lfftw3_omp  -lfftw3 -lm -I/users/hectorgm/fftw3_threads/include/ -L/users/hectorgm/fftw3_threads/lib/ -o file.out

====Run====
./file.out param_file.txt

====Parameter file options====

#Type of Box (periodic/cutsky): 'periodic' for periodic boxes with boundary conditions. 'cutsky' for actual observations or mocks which include sky mask

#Type of file (ascii/gadget): 'ascii' is the option required for 'cutsky'. 'periodic' option allows 'ascii' files or 'gadget' files. See 'ascii file structure' below for the format of the file

#Number of gadget files(int): In case the gadget boxes are split in more than 1 gadget file

#RSD distorsion on gadget periodic box (yes/no): yes. For gadget boxes allow this option to distort particles along the z-axis for redshift space distortions

#Size of the Box (double/double): Low and Upper limits, respectively, of the cubic box where the galaxies are placed. 

#Type of Computation (DSE/DSY/FFT): Type of Computation for the Power Spectrum: Direct Sum Exact (DSE); Direct Sum Yamamoto (DSY); Fast Fourier Transform (FFT). For the bispectrum computation FFT is required.

#Binning for the Power Spectrum (linear/log10):

====Ascii File Structure=====

For 'periodic' option: 4-column entry
x-position (double), y-position (double), z-position (double), weight (double)

For example,
2.234479 8.066518 6.439348 1.0
1.406139 5.003024 3.664719 1.0
9.559097 7.672305 10.529753 1.0
16.126962 1.955036 1.807368 1.0
13.377709 2.506474 2.225270 1.0
11.555976 0.431673 3.149166 1.0
14.122240 16.198794 3.532218 1.0
14.214268 14.341477 4.141918 1.0

For 'skycut' option: data and random entry required
#data. 8-column entry:
Right Ascension (double), declination (double), redshift (double), weight_fkp (double), weight_colision (integer), weight_systematics (double), number_density (double), veto (integer)
For example
1.291761884477e+02 4.894649924515e+01 5.425299000000e-01 1.238355000000e-01 1 1.029931000000e+00 7.075228831797e-04 1
1.174169630420e+02 3.927675925314e+01 3.996815000000e-01 6.395705000000e-01 1 1.020382000000e+00 5.635492881551e-05 1
1.169127239443e+02 3.944331088030e+01 5.377024000000e-01 1.231396000000e-01 1 9.972784000000e-01 7.120864449779e-04 0
1.169501719710e+02 3.949076919700e+01 5.191724000000e-01 1.143980000000e-01 1 1.061827000000e+00 7.741411563139e-04 1
1.175284705313e+02 4.017649329509e+01 5.431913000000e-01 1.246438000000e-01 1 1.022436000000e+00 7.022861947405e-04 1
1.238161587640e+02 4.663678411199e+01 5.896081000000e-01 1.800036000000e-01 1 1.009149000000e+00 4.555444446667e-04 1
1.276012766986e+02 4.977593704841e+01 5.481966000000e-01 1.309460000000e-01 1 1.053429000000e+00 6.636735753669e-04 1

#randoms. 6-column entry:
Right Ascension (double), declination (double), redshift (double), weight_fkp (double), number_density (double), veto (int)
For example



====Citation====
If you use this code for your published or unpublished work, please refer it to Gil-Marin, Hector in prep. 
====Disclaimer====
blah blah blah
