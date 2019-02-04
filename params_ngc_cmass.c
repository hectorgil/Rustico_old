#Main parameters
#Type of Box (periodic/cutsky): cutsky
#Type of file (ascii/gadget): ascii
#Number of gadget files(int): 1
#RSD distorsion on gadget periodic box (yes/no): no
#Size of the Box (double/double): -2600.0 +2400.0
#Type of Computation (DSE/DSY/FFT): FFT
#Binning for the Power Spectrum (linear/log10): linear
#Size of the bin for the power spectrum (double): 0.01
#k-range for computation (double/double): 0 0.32

#Bispectrum parameters
#Do Bispectrum (yes/no): no
#Do Multigrid (yes/no): no
#Triangle Shapes (ALL/EQU/ISO/SQU): ALL
#Size of the bin for the bispectrum (double): 0.01077117481
#Normalization of triangles(FFT/APR_SUM/APR_EFF,EXA_EFF): FFT
#Write triangles in each bin(yes/no): no
#Path for triangles in each bin: ./power_spectra/triangles

#Read inout parameters
#Path of data: /mnt/lustre/eboss/EZ_MOCKs/EZmock_LRG_v4.0/CMASS_NGC/EZmock_CMASS_LRG_NGC_DR12v5_0001.dat 
#Path of randoms: /mnt/lustre/eboss/EZ_MOCKs/EZmock_LRG_v4.0/RANDOM/random_20x_CMASS_LRG_NGC_DR12v5.dat
#Path of output: ./cmass
#Identifier of output: cmass_NGC_0001
#Write header: yes

#FFT parameters
#Number of Grid Cells power (int): 9
#Type of mass assingment (NGC/CIC/TSC/PCS/P4S/P5S): PCS
#Type of Yamamoto (GridCenter/GridAverage): GridCenter
#Number of interlacing steps (int): 2
#Do Grid Correction? (yes/no): yes


#Cutsky parameters
#Redshift Range (double/double): 0.6 1.0
#Omega matter value (double): 0.31
#Area effective value in deg^2 (double): 6969.68
#Hexadecapole as (L4/L2L2): L2L2
#Compute Normalization as (area/density): density
#Compute Normalization using (randoms/data): data
#Compute Shot noise as (double): 1.0
