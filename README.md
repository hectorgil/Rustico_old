RUSTICO: Rapid foUrier STatIstics COde
Author: Hector Gil Marin
Date: 4th Feb. 2019
Internal Version: 3.2
email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu
(c) All rights reserved

====Compilation Examples====

icc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c -O3 -lm -openmp -lfftw3_omp  -lfftw3 -lpthread -lfftw3_threads -I/users/hectorgm/fftw3_threads/include/ -L/users/hectorgm/fftw3_threads/lib/ -o file_icc.out

gcc main.c bispectrum.c functions.c mass_assignment.c fftw_compute.c read_positions.c read_line.c cubature.c ps_write.c order_algorithm.c -O3 -std=c99 -lm -fopenmp -lfftw3_omp  -lfftw3 -lpthread -lfftw3_threads -I/users/hectorgm/fftw3_threads_gcctest/include -L/users/hectorgm/fftw3_threads_gcctest/lib/ -o file_gcc.out

you will need to have the FFTW libraries installed with the following option in the configure file,

./configure  --enable-openmp --enable-threads  --prefix="/home/hector/fftw3_threads/" CC="gcc"

where CC="icc" in case you use the intel compiler. 

and provide the link to the lib and include paths in case they are not install authomatically in the root path. 

====Run====

./file_icc.out param_file.txt

====Parameter file options====

#Type of Box (periodic/cutsky): 'periodic' for periodic boxes with boundary conditions. 'cutsky' for actual observations or mocks which include sky mask

#Type of file (ascii/gadget): 'ascii' is the option required for 'cutsky'. 'periodic' option allows 'ascii' files or 'gadget' files. Gadget units assumed kpc/h. See 'ascii file structure' below for the format of the file

#Number of gadget files(int): In case the gadget boxes are split in more than 1 gadget file

#RSD distorsion on gadget periodic box (yes/no): yes. For gadget boxes allow this option to distort particles along the z-axis for redshift space distortions. The values of redhisft and Omega matter, will be taken from the header of the gadget file.

#Size of the Box (double/double): Low and Upper limits, respectively, of the cubic box where the galaxies are placed. 

#Type of Computation (DSE/DSY/FFT): Type of Computation for the Power Spectrum: Direct Sum Exact (DSE); Direct Sum Yamamoto (DSY); Fast Fourier Transform (FFT). For the bispectrum computation FFT is required.

#Binning for the Power Spectrum (linear/log10): Binning type for the power spectrum output. Linear or 10-base logarithmic. 

#Size of the bin for the power spectrum (double). Size of the bin for the power spectrum. In case log10 is choosen as binning, the interval is provided in log10-scale.

#k-range for computation (double/double): Low and Upper limits, respectively, of the k-values choosen for printing the power spectrum.

#Do Bispectrum (yes/no): Whether the bispectrum should be computed by the code

Do Multigrid (yes/no): Option for the bispectrum computation. If enable, the bispectrum triangles will be split according to their k-values and associated to different grid-sizes for a more optimal computation (large scale modes do not requires small grid cell ressolution). However, each grid-size computation will requires to re-associate the particles to the grid cells, which a potential lose of optimality. We recomend enable such option when many triangle shapes are required and when the datasets do not consists of many particles. In practice, each specific case will requires testing for determing the best performance option. When the multigrid option is enable, we require the interlacing option to be also enabled (see below), with at least 2 interlacing steps.

#Triangle Shapes (ALL/EQU/ISO/SQU). Triangle shapes to be computed. All (ALL), equilateral (EQU), Isosceles (ISO), squeezed (SQU). We define the squeezed triangles as those |k2-k3|<=k1 and K1<=0.1 K2; where by definition K1<=K2<=K3. Note that this condition is applied to the center of bin k-values and not to the effective k-values.

#Size of the bin for the bispectrum (double): Size of the bispectrum bin

#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): This determine the way how the Bispectrum is normalized: FFT is through performing Eq. xxx of ... using Fourier Transforms; APR is using the approximate solution of Eq. xxx 8pi^2k1k2k3Dk^3 and EXA is using the full and exact analytic solution (only available for Equilateral triangles at the moment). SUM option is computing APR and EXA for each triangle shape in the k-bin, whereas EFF is computing APR and EXA only for the effective k-values in the k-bin. 

#Write triangles in each bin(yes/no): This option enables writting each triangle shape inside each k-bin. This option is only available for squeezed triangles (SQU choice enabled).

#Path for triangles in each bin: Determines the path for writting the above triangles

#Path of data: path of data file

#Path of randoms: path for the random file. If periodic box enabled, write 'nothing'

#Path of output: path were the output files will be written

#Identifier of output: identification string for the output files

#Write header: option for writting the header for power spectrum and bispectrum output files. 

#Number of Grid Cells power (int): Number of grid-cells-per-side power input. If the input is n, the number of grid-cells per side will be 2^n. The input is an integer in the range 4<n<15

#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): mass interpolation assignment scheme: nearest-grid-point NGC, cloud-in-cell CIC, triangular-shaped-cloud TSC, piecewise-cubic-spline, piecewise-quartic-spline P4S, piecewise-quintic-spline P5S.

#Type of Yamamoto (GridCenter/GridAverage): Option only for skycut option. (details to be referred in a paper).

#Number of interlacing steps (int): Number of interlacing steps.

#Do Grid Correction? (yes/no): Grid correction option. (Jing et al 2005)

#Redshift Range (double/double):  Low and Upper limits, respectively, for the redshift cuts in the skycut option

#Omega matter value (double): Omega matter value used for converting redshifts to comoving distances in the skycut option.

#Area effective value in deg^2 (double): Value of the area used for the normalization and of the power spectrum and bispectrum in the skycut option

#Hexadecapole as (L4/L2L2): Hexadecapole type of computation for the skycut option. L4 refers to Eq. xxx of Bianchi et al. 2015 and L2L2 to Eq. xxx of Scoccimarro et al. 2015.

#Compute Normalization as (area/density): Compute the normalization of the power spectrum using either the area value of the number density column of the input. Only for skycut option

#Compute Normalization using (randoms/data):Compute the normalization of the power spectrum and bispectrum using either the data or random n(z) computed from the objects. Only for skycut option

#Compute Shot noise as (double): Shot noise factor parameter. See. Gil-Marin et al. 2014 Eq. ..... Only for skycut option


====Output Structure====

The power spectrum output have the following format:
k-eff, k-centerbin, Monopole-Pshotnoise, Quadrupole, Hexadecapole, number of modes, Pshotnoise

The bispectrum output have the following format:
k1-eff, k1-centerbin, k2-eff, k2-centerbin, k3-eff, k3-centerbin, B0-Bshotnoise, Bshotnoise, Reduced Bispectrum, Reduced Bispectrum shot noise, number of triangles

For the skycut option, the code automatically generates two extra file for the number density of objects as a function of redshift for the data and random files with the following format,

#data

#Interval: xxx Mpc/h

#z < nobs > < wc nobs> < wc wfkp nobs>

#random

#Interval: yyy  Mpc/h

#alpha< ns > alpha < wfkp ns >

where xxx and yyy are the interval choosen by the code to bin the data and randoms, respectively, < nobs > are the raw number density of observed objects, < wc nobs> is the number density of observed objects weighted by the collision weights, and < wc wfkp nobs> weighted by the collision weights and fkp weights. For the randoms, these number densities are scaled by alpha in order to match the data ones. 

====Citation====

If you use this code for your published or unpublished work, please refer it to Gil-Marin, Hector in prep. 2017

====Disclaimer====

....
