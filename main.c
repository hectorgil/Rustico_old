/*
RUSTICO   Rapid foUrier STatIstics COde
Author: Hector Gil Marin
Date: 1st May 2017
email: hector.gil.marin@gmail.com or hector.gilmarin@lpnhe.in2p3.fr
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "mass_assignment.h"
#include "read_positions.h"
#include "functions.h"
#include "fftw_compute.h"
#include "ps_write.h"
//#include "order_algorithm.h"
#include "bispectrum.h"

typedef struct{
          double OMEGA_M;
}f_params;


int main(int argc, char *argv[])
{

FILE *f,*g;
int i,tid;

char name_file_ini[2000];//name of inizialization file
//Main Parameters read in
double L1,L2;//Box limit in Mpc/h
char type_of_survey[50];//type of survey: Periodic or Cutsky
char type_of_computation[10];//type of computation: DSY /  FFT / DSE
int power_grid;//power in number of grid cells: integer between 6 and 15 
int ngrid;//Number of grid cells: pow(2,power_grid)
char binning_type[10];//Type of binning for the power spectrum: linear or log
double bin_ps;//Size of the bin for the power spectrum in h/Mpc
char header[10];//write header yes/no
char type_of_file[10];//ascii or gadget(only for periodic boxes)
int gadget_files;//number of gadget files per realization
int snapshot_num;//individual gadget file
char RSD[10];//RSD distorsion on Gadget boxes?

//Bispectrum parameters
char do_bispectrum[10];//Do Bispectrum?
double Deltakbis;//Triangle Bin size in terms of k-fundamental
double kmin,kmax;//Minimum and maximum k
char triangles_num[2000];//FFT,APR_SUM,APR_EFF,EXA_EFF
char write_triangles[2000];
char path_for_triangles[2000];
char triangles_id[2000];
char do_multigrid[10];
char triangle_shapes[10];

//Read inout parameters
char name_data_in[2000];//Path and name of data source
char name_gadget_file[2000];//Full path for gadget like file
long int Ndata;//Number of lines of data source
long int Ndata2;//Number of lines of data used
char name_randoms_in[2000];//Path and name of random source
long int Nrand;//Number of lines of random source
long int Nrand2;//Number of lines of random used
char name_path_out[2000];//Path where to write output
char name_id[2000];//String to identify output
char name_ps_out[2000];//Final name of the output for the power spectrum
char name_bs_out[2000];//Final name of the output for the bispectrum
char name_den_out[2000];//Final name of the output for the density

//FFT parameters
char type_of_mass_assigment[10];//Type of mass assignment: NGC, CIC, TSC, PCS, P4S, P5S
int Ninterlacing;//Number of interlacing steps (1 for no-interlacing steps)
char grid_correction_string[10];// Do Grid Correction: yes/no
int mode_correction;//Correction factor power
char type_of_yamamoto[20];

//Cutsky parameters
double z_min,z_max;//Minimum and maximum (excluding) redshift cuts
double Omega_m;//Value of Omega matter;
double Area_survey;//Value of Area of the survey in deg^2
char Hexadecapole_type[20];//L4 or L2L2
char type_normalization_mode[20];//Normalize according to the area of the survey or the n(z) value: Area / Density
char type_normalization_mode2[20];//Normalize according to the n(z) of data or randoms file
double Shot_noise_factor;// Factor between 0 and 1. 0 Correspond ....

//Positions pointers
  double* pos_x;
  double* pos_y;
  double* pos_z;
  double* weight;
  double* pos_x_rand;
  double* pos_y_rand;
  double* pos_z_rand;
  double* weight_rand;
//  double DeltaR;  

 //Maximum and minimum values for positions of galaxies and randoms.
  double max,min;
//Parameters relative to shot noise, effective redshifts, number of particles and normalization
double Psn_1a, Psn_1b, Psn_2, z_effective_data,I_norm_data,z_effective_rand,I_norm_rand,alpha_data,alpha_rand,alpha,I22,I_norm_data2,I_norm_data3,I_norm_data4, I_norm_rand2,I_norm_rand3,I_norm_rand4;
double P_shot_noise1,P_shot_noise2,P_shot_noise;
double Bsn1,Bsn2,IN1,IN2,I3_norm_data2,I3_norm_data3,I3_norm_data4,I3_norm_rand2,I3_norm_rand3,I3_norm_rand4,I33,Bsn,IN;
int num_effective;
double num_effective2;

    int n_lines_parallel;//Number of parallel threads



//Read Inizialization parameters
sprintf(name_file_ini,argv[1]);
f=fopen(name_file_ini,"r");
if(f==NULL){printf("File %s not found...\t exiting now\n",name_file_ini);return 0;}
else{printf("Reading Inizialization file: %s\n\n",name_file_ini);}
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %s\n",&type_of_survey);
fscanf(f,"%*s %*s %*s %*s %s\n",&type_of_file);
fscanf(f,"%*s %*s %*s %*s %d\n",&gadget_files);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %s\n",RSD);
printf("== Inizialization Parameters ==\n");
printf("Type of Survey: %s\n",type_of_survey);
printf("Type of File: %s\n",type_of_file);
if(strcmp(type_of_file, "gadget") == 0){printf("Number of gadget files: %d\n",gadget_files);}
if(strcmp(type_of_file, "gadget") == 0){printf("RSD distortion: %s\n",RSD);}
fscanf(f,"%*s %*s %*s %*s %*s %lf %lf\n",&L1,&L2);
printf("Box edges at %lf Mpc/h and %lf Mpc/h; Size of the Box %lf Mpc/h\n",L1,L2,L2-L1);
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_computation);
printf("Type of Computation: %s\n",type_of_computation);
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n",binning_type);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %lf\n",&bin_ps);
fscanf(f,"%*s %*s %*s %*s %lf %lf\n\n",&kmin,&kmax);
printf("Binning for the Power Spectrum: %s; Size of bin: %lf; k-range %lf < k[h/Mpc] < %lf\n\n",binning_type,bin_ps,kmin,kmax);
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",do_bispectrum);
fscanf(f,"%*s %*s %*s %s\n",do_multigrid);
fscanf(f,"%*s %*s %*s %s\n",triangle_shapes);
fscanf(f,"%*s %*s %*s %*s %*s %*s %*s %*s %lf\n",&Deltakbis);
fscanf(f,"%*s %*s %*s %s\n",&triangles_num);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",&write_triangles);
fscanf(f,"%*s %*s %*s %*s %*s %*s %s\n\n",&path_for_triangles);

printf("== Bispectrum Parameters ==\n");
printf("Do Bispectrum? %s\n",do_bispectrum);
if(strcmp(do_bispectrum, "yes") == 0){printf("Bispectrum bins %lf\nMultigrid Computation:%s\nTriangle Shapes: %s\nTriangle normalization: %s\nWrite individual Triangles:%s \n\n",Deltakbis,do_multigrid,triangle_shapes,triangles_num,write_triangles);}

fscanf(f,"%*s %*s %*s\n");
fscanf(f,"%*s %*s %*s %s\n",name_data_in);
if(strcmp(type_of_file, "ascii") == 0)
{
g=fopen(name_data_in,"r");
if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_data_in);return 0;}
fclose(g);
Ndata=countlines(name_data_in);
}

printf("== Read in/out options ==\n");
printf("Reading data file %s; %d lines\n",name_data_in,Ndata);

fscanf(f,"%*s %*s %*s %s\n",name_randoms_in);

if(strcmp(type_of_survey, "cutsky") == 0)
{
g=fopen(name_randoms_in,"r");
if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_randoms_in);return 0;}
fclose(g);

Nrand=countlines(name_randoms_in);
printf("Reading randoms file %s; %d lines\n",name_randoms_in,Nrand);
}

if(strcmp(type_of_file, "gadget") == 0)
{ 
  Ndata=0;
  for(snapshot_num=0;snapshot_num<gadget_files;snapshot_num++)
  {
     sprintf(name_gadget_file, "%s.%d", name_data_in, snapshot_num);
     g=fopen(name_gadget_file,"r");
     if(g==NULL){printf("File %s does not exist. Exiting now...\n",name_data_in);return 0;}
     fclose(g);
     Ndata+=count_particles_gadget(name_gadget_file);
  }
P_shot_noise=pow(L2-L1,3)/Ndata*1.;
}

fscanf(f,"%*s %*s %*s %s\n",name_path_out);
printf("Output files at %s\n",name_path_out);
fscanf(f,"%*s %*s %*s %s\n",name_id);
printf("Output Id %s\n",name_id);
fscanf(f,"%*s %*s %s\n",&header);
printf("Write header? %s\n\n",header);
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %*s %*s %d\n",&power_grid);
ngrid=pow(2,power_grid);
fscanf(f,"%*s %*s %*s %*s %*s %s\n",type_of_mass_assigment);
printf("== FFT options ==\n");
printf("Number of k-modes and grid cells per side: %ld\n",ngrid);
printf("Type of mass assingment %s\n",type_of_mass_assigment);

if(strcmp(type_of_mass_assigment, "NGC") == 0){mode_correction=1;}
if(strcmp(type_of_mass_assigment, "CIC") == 0){mode_correction=2;}
if(strcmp(type_of_mass_assigment, "TSC") == 0){mode_correction=3;}
if(strcmp(type_of_mass_assigment, "PCS") == 0){mode_correction=4;}
if(strcmp(type_of_mass_assigment, "P4S") == 0){mode_correction=5;}
if(strcmp(type_of_mass_assigment, "P5S") == 0){mode_correction=6;}
fscanf(f,"%*s %*s %*s %*s %s\n",type_of_yamamoto);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Type of Yamamoto: %s\n",type_of_yamamoto);}
fscanf(f,"%*s %*s %*s %*s %*s %d\n",&Ninterlacing);
printf("Number of Interlacing steps %d\n",Ninterlacing);
fscanf(f,"%*s %*s %*s %*s %s\n",grid_correction_string);
printf("Do grid correction? %s\n\n",grid_correction_string);
if(strcmp(grid_correction_string, "no") == 0){mode_correction=0;}
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %lf %lf\n",&z_min,&z_max);
if(strcmp(type_of_survey, "cutsky") == 0){printf("== Cutsky options ==\n");}
if(strcmp(type_of_survey, "cutsky") == 0){printf("Redshift cuts: %lf < z < %lf\n",z_min,z_max);}
fscanf(f,"%*s %*s %*s %*s %lf\n",&Omega_m);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Omega_m: %lf\n",Omega_m);}
fscanf(f,"%*s %*s %*s %*s %*s %*s %lf\n",&Area_survey);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Area of the survey %lf deg2\n",Area_survey);}
fscanf(f,"%*s %*s %*s %s\n",Hexadecapole_type);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Do Hexadecapole as %s\n",Hexadecapole_type);}
fscanf(f,"%*s %*s %*s %*s %s\n",type_normalization_mode);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Compute normalization using %s\n",type_normalization_mode);}
fscanf(f,"%*s %*s %*s %*s %s\n",type_normalization_mode2);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Compute normalization using %s file n(z)\n",type_normalization_mode2);}
fscanf(f,"%*s %*s %*s %*s %*s %lf\n",&Shot_noise_factor);
if(strcmp(type_of_survey, "cutsky") == 0){printf("Shot noise factor set to %lf\n",Shot_noise_factor);}
fclose(f);


//Error conditions
if(strcmp(type_of_file, "gadget") != 0 && strcmp(type_of_file, "ascii") != 0){printf("File type must be either 'gadget' or 'ascii'. Entry read %s. Exiting now...\n",type_of_file);return 0;}
if( strcmp(type_of_survey, "cutsky") != 0 && strcmp(type_of_survey, "periodic") != 0){printf("Survey type must be either 'cutsky' or 'periodic'. Entry read %s. Exiting now...\n",type_of_survey);return 0;}
if( strcmp(type_of_survey, "cutsky") == 0 && strcmp(type_of_file, "gadget") == 0 ){printf("Warning. Cutsky+gadget option not available. Exiting now...\n");return 0;}
if(gadget_files<1){printf("Warning. gadget files entry must be >0. Entered value %d. Exiting now...\n",gadget_files);return 0;}
if( strcmp(RSD, "yes") != 0 && strcmp(RSD, "no") != 0 ){printf("Warning. RSD option only accepts either 'yes' or 'no' entries. Entry read %s. Exiting now...\n",RSD);return 0;}
if(L2<=L1){printf("Error: L2 parameter has to be large then L1. Exiting now...\n");return 0;}
if( strcmp(type_of_computation, "DSE") !=0 &&  strcmp(type_of_computation, "DSY") !=0 &&  strcmp(type_of_computation, "FFT") !=0){printf("Type of computation only accepts 'DSE', 'DSY' or 'FFT' options. Entry read %s. Exiting now...\n",type_of_computation);return 0;}
if( strcmp(type_of_computation, "FFT") !=0 && strcmp(do_bispectrum, "yes") ==0){printf("Warning. Bispectrum computation only accepts FFT computation option. Entry read %s. Exiting now...\n",do_bispectrum);return 0;}
if( strcmp(binning_type, "log10") !=0 && strcmp(binning_type, "linear") !=0){printf("Warning. Binning type must be either 'log10' or 'linear'. Entry read %s. Exiting now...\n",binning_type);return 0;}
if( strcmp(binning_type, "log10") ==0 && strcmp(do_bispectrum, "yes") ==0){printf("Warning. 'log10' type of binning not available for bispectrum computation on this version. Exiting now...\n"); return 0;}
if(bin_ps<=0){printf("Error: Size of the power spectrum bin has to be greater than 0: Entry read %lf. Exiting now...\n",bin_ps);return 0;}
if(kmin<0 || kmax<=0 || kmin>kmax){printf("Warning: Unusual values for maximum and/or minimum k-values for the bispectrum computation: kmin=%lf, kmax=%lf. Exiting now...\n",kmin, kmax);return 0;}
if(kmin==0 && strcmp(binning_type, "log10") == 0){printf("Cannot set kmin=0 and log-k binning. Exiting now...\n");return 0;}
if( strcmp(do_bispectrum, "yes") != 0 && strcmp(do_bispectrum, "no") != 0){printf("Error. Bispectrum entry must be either 'yes' or 'no'. Read entry %s. Exiting now...\n",do_bispectrum);return 0;}
if( strcmp(do_multigrid, "yes") != 0 && strcmp(do_multigrid, "no") != 0){printf("Error. Multigrid entry must be either 'yes' or 'no'. Read entry %s. Exiting now...\n",do_multigrid);return 0;}
if( strcmp(do_multigrid, "yes") == 0 && Ninterlacing<2 ){printf("Warning. Multigrid option requires a number of interlacing steps >1. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"NGC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"CIC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(do_multigrid, "yes") == 0 && strcmp(type_of_mass_assigment,"TSC") ==0 ){printf("Warning. Multigrid option requires a mass interpolation grid of at least PCS. Exiting now...\n");return 0;}
if( strcmp(triangle_shapes,"ALL") !=0 && strcmp(triangle_shapes,"EQU") !=0 && strcmp(triangle_shapes,"ISO") !=0 && strcmp(triangle_shapes,"SQU") !=0 ){printf("Error. Triangle shapes entry only accepts 'ALL', 'ISO', 'EQU' or 'SQU'. Read entry %s. Exiting now...\n",triangle_shapes);return 0;}
if(Deltakbis<=0){printf("Error: Size of the bispectrum bin has to be greater than 0: %lf kf. Exiting now...\n",Deltakbis);return 0;}
if( strcmp(triangles_num, "FFT") != 0 && strcmp(triangles_num, "APR_SUM") !=0 && strcmp(triangles_num, "EXA_SUM") !=0  && strcmp(triangles_num, "APR_EFF") !=0 && strcmp(triangles_num, "EXA_EFF") !=0){printf("Error. Number of Triangles normalization option only accepts: 'FFT', 'APR_SUM', 'EXA_SUM', 'APR_EFF', 'EXA_EFF'. Entry read %s. Exiting now...\n",triangles_num);return 0;}
if( strcmp(triangles_num, "EXA_SUM") == 0 &&  strcmp(triangle_shapes,"EQU") !=0){printf("Warning. 'EXA_SUM' triangle normalization option it is only available for equilateral triangles. Exiting now...\n");return 0;}//not available at the moment
if( strcmp(triangles_num, "EXA_EFF") == 0 &&  strcmp(triangle_shapes,"EQU") !=0){printf("Warning. 'EXA_EFF' triangle normalization option it is only available for equilateral triangles. Exiting now...\n");return 0;}//not available at the moment
if( strcmp(write_triangles, "yes") != 0 &&  strcmp(write_triangles, "no") !=0){printf("Error. write triangle entry must be either 'yes' or 'now'. Exiting now...\n");return 0;}
if( strcmp(write_triangles, "yes") == 0 && strcmp(triangle_shapes,"SQU") !=0){printf("Waring. Write triangle option 'yes' is only recomended for squeezed triangle shapes 'SQU'. Exiting now...\n");return 0;}
if( strcmp(header, "yes") !=0 && strcmp(header, "no") !=0){printf("Waring. Write header option must be either 'yes' or 'no'. Entry read %s. Exiting now...\n",header);return 0;}
if(power_grid<4 || power_grid>15){printf("Warning: Unusual value for number of grid cells per side: 2^%d=%ld. Exiting now...\n",power_grid,ngrid);return 0;}
if(strcmp(type_of_mass_assigment,"NGC") !=0 && strcmp(type_of_mass_assigment,"CIC") !=0 && strcmp(type_of_mass_assigment,"TSC") !=0 && strcmp(type_of_mass_assigment,"PCS") !=0 && strcmp(type_of_mass_assigment,"P4S") !=0 &&  strcmp(type_of_mass_assigment,"P5S") !=0){printf("Error. Type of mass assigment must be either 'NGC', 'CIC', 'TSC', 'PCS', 'P4S' or 'P5S'. Entry read %s. Exiting now...\n",type_of_mass_assigment);return 0;}
if( strcmp(type_of_yamamoto, "GridCenter") != 0 && strcmp(type_of_yamamoto, "GridAverage") != 0){printf("Error. Type of Yamamoto option must be either 'GridCenter' or 'GridAverage'. Entry read %s. Exiting now...\n",type_of_yamamoto);return 0;}
if(Ninterlacing<=0){printf("Error: Number of interglacing steps has to be equal or larger than 1. %d\n",Ninterlacing);return 0;}
if( strcmp(grid_correction_string, "yes") !=0 && strcmp(grid_correction_string, "no") !=0){printf("Grid correction input must be either 'yes' or 'no'. Entry read %s. Exiting now...\n",grid_correction_string);return 0;}

if(strcmp(type_of_survey, "cutsky") == 0){
if(z_min>=z_max){printf("Error. Minimum value for redshift is larger than the maximum: z_min=%lf; z_max=%lf. Exiting now...\n",z_min,z_max);return 0;}
if(Omega_m<=0 || Omega_m>1){printf("Warning. Unusual value for Omega_m, Omega_m=%lf. Exiting now...\n",Omega_m);return 0;}
if(Area_survey<=0){printf("Warning. Usual value for the Area of the survey: %lf. Exiting now...\n",Area_survey);return 0;}
if( strcmp(Hexadecapole_type, "L4") !=0 && strcmp(Hexadecapole_type, "L2L2") !=0){printf("Hexadecapole option must be either 'L2L2' or 'L4'. Entry read %s. Exiting now...\n",Hexadecapole_type);return 0;}
if( strcmp(type_normalization_mode, "area") !=0 && strcmp(type_normalization_mode, "density") !=0){printf("Error. Normalisation type must be either 'area' or 'density'. Entry read %s. Exiting now...\n",type_normalization_mode);return 0;}
if( strcmp(type_normalization_mode2, "data") !=0 && strcmp(type_normalization_mode2, "randoms") !=0){printf("Error. Normalisation type must be either 'data' or 'randoms'. Entry read %s. Exiting now...\n",type_normalization_mode2);return 0;}
if(Shot_noise_factor>1 || Shot_noise_factor<0){printf("Warning. Usual value for the Shot noise factor: %lf. Exiting now...\n",Shot_noise_factor);return 0;}
if( strcmp(do_bispectrum, "yes") == 0 && strcmp(type_normalization_mode, "density") == 0 ){printf("Warning. Bispectrum computation requires a normalization by 'area' and not by 'density'. Exiting now...\n");return 0;}
}
//etc strings.....

//Reading files.
double parameter_value[26];
parameter_value[0]=Omega_m;
parameter_value[1]=z_min;
parameter_value[2]=z_max;
parameter_value[3]=Ndata;
parameter_value[13]=Area_survey;

if(strcmp(type_of_survey, "cutsky") == 0)
{
                Ndata2=get_number_used_lines_data(name_data_in,parameter_value);
}
if(strcmp(type_of_survey, "periodic") == 0)
{
Ndata2=Ndata;
}

if(strcmp(type_of_file, "ascii") == 0)//these are only kept stored during all the process for ascii files. Gadget files keep the name of the file and read it each time they need
{
                pos_x = (double*) calloc(Ndata2, sizeof(double));
		pos_y = (double*) calloc(Ndata2, sizeof(double));
	        pos_z = (double*) calloc(Ndata2, sizeof(double));
		weight = (double*) calloc(Ndata2, sizeof(double));
}

if(strcmp(type_of_survey, "cutsky") == 0)
{
//pos_x,pos_y,pos_z,weight are loaded with the position of particles
//Ndata is uploaded to the number of particles used
//Psn_1a,Psn_1b,Psn_2 are uploaded with information on the shot noise
//z_efffective is uploaded with the effective redshift of the sample
//num_effective is uploaded with the effective number of particles
//I_norm is uploaded with information relative to the normalization based on density
//alpha is uploeaded with information relative to the effective ratio between data and randoms
printf("Reading %s...",name_data_in);
get_skycuts_data(name_data_in, pos_x, pos_y, pos_z, weight, parameter_value,type_normalization_mode);

Ndata2=parameter_value[3];
Psn_1a=parameter_value[4];
Psn_1b=parameter_value[5];
Psn_2=parameter_value[6];
z_effective_data=parameter_value[7];
num_effective=(int)(parameter_value[8]);
I_norm_data=parameter_value[9];
min=parameter_value[10];
max=parameter_value[11];
alpha_data=parameter_value[12];

I_norm_data2=parameter_value[14];
I_norm_data3=parameter_value[15];
I_norm_data4=parameter_value[16];
num_effective2=parameter_value[17];

I3_norm_data2=parameter_value[18];
I3_norm_data3=parameter_value[19];
I3_norm_data4=parameter_value[20];

Bsn1=parameter_value[21];
Bsn2=parameter_value[22];
IN1=parameter_value[23];
IN2=parameter_value[24];

parameter_value[3]=Ndata;
sprintf(name_den_out,"%s/Density_galaxies_%s.txt",name_path_out,name_id);
get_skycuts_write_density_data(name_data_in, parameter_value,name_den_out);
parameter_value[3]=Ndata2;

printf("Ok!\n");

parameter_value[0]=Omega_m;
parameter_value[1]=z_min;
parameter_value[2]=z_max;
parameter_value[3]=Nrand;
Nrand2=get_number_used_lines_randoms(name_randoms_in,parameter_value);

        pos_x_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_y_rand = (double*) calloc(Nrand2, sizeof(double));
        pos_z_rand = (double*) calloc(Nrand2, sizeof(double));
        weight_rand = (double*) calloc(Nrand2, sizeof(double));
								
printf("Reading %s...",name_randoms_in);
parameter_value[12]=alpha_data;
sprintf(name_den_out,"%s/Density_randoms_%s.txt",name_path_out,name_id);
get_skycuts_randoms(name_randoms_in, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, parameter_value,type_normalization_mode, type_normalization_mode2);
Nrand2=parameter_value[3];
z_effective_rand=parameter_value[7];
if(min>parameter_value[10]){min=parameter_value[10];}
if(max<parameter_value[11]){max=parameter_value[11];}
alpha_rand=parameter_value[12];

I_norm_rand=parameter_value[9];
I_norm_rand2=parameter_value[14];
I_norm_rand3=parameter_value[15];
I_norm_rand4=parameter_value[16];

I3_norm_rand2=parameter_value[17];
I3_norm_rand3=parameter_value[18];
I3_norm_rand4=parameter_value[19];

alpha=alpha_data/alpha_rand;
I_norm_rand=I_norm_rand*alpha;

if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "density") == 0 ){I22=I_norm_rand;}
if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "density") == 0){I22=I_norm_data;}

if(strcmp(type_normalization_mode2, "randoms") == 0 && strcmp(type_normalization_mode, "area") == 0 ){I22=I_norm_rand2;I33=I3_norm_rand2;}
if(strcmp(type_normalization_mode2, "data") == 0 && strcmp(type_normalization_mode, "area") == 0){I22=I_norm_data2;I33=I3_norm_data2;}


P_shot_noise1=(Psn_1a+alpha*Psn_2)/I22;
P_shot_noise2=(Psn_1b+alpha*Psn_2)/I22;
P_shot_noise=P_shot_noise1*Shot_noise_factor+P_shot_noise2*(1.-Shot_noise_factor);
Bsn=(Bsn1-alpha*alpha*Bsn2)/I33;
IN=(IN1*Shot_noise_factor+IN2*(1.-Shot_noise_factor))/I33;

parameter_value[3]=Nrand;
sprintf(name_den_out,"%s/Density_randoms_%s.txt",name_path_out,name_id);
get_skycuts_write_density_randoms(name_randoms_in, parameter_value,alpha,name_den_out);
parameter_value[3]=Nrand2;

printf("Ok!\n\n");


}


if(strcmp(type_of_survey, "periodic") == 0 && strcmp(type_of_file, "ascii") == 0)
{
printf("Reading %s...",name_data_in);
get_periodic_data(name_data_in, pos_x, pos_y, pos_z, weight, parameter_value);
P_shot_noise=pow(L2-L1,3)/Ndata*1.;
min=parameter_value[10];
max=parameter_value[11];
printf("Ok!\n\n");
}

if(max>L2 || min<L1){printf("Warning: Limits of the box are exceeded by the data or random galaxies: data particles found at the limits %lf and %lf. Exiting now...\n",min,max);return 0;}

sprintf(name_ps_out,"%s/Power_Spectrum_%s.txt",name_path_out,name_id);
sprintf(name_bs_out,"%s/Bispectrum_%s.txt",name_path_out,name_id);
sprintf(triangles_id,"%s/Triangles_%s",path_for_triangles,name_id);

f=fopen(name_ps_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_ps_out);return 0;}

//write header for the power spectrum file
if( strcmp(header, "yes") == 0)
{
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %d\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Number of random elements used: %d\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %d\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s\n", type_normalization_mode );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization %lf\n",I22);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Shot noise factor %lf\n",Shot_noise_factor);}
fprintf(f,"#Shot noise value %lf\n",P_shot_noise);
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %ld\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf\n",Omega_m);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Hexadecapole as %s\n",Hexadecapole_type);}

}
fclose(f);

//Write header for the Bispectrum file
if( strcmp(header, "yes") == 0 && strcmp(do_bispectrum, "yes") == 0)
{
f=fopen(name_bs_out,"w");
if(f==NULL){printf("Could not write %s\n. Exiting now...\n",name_bs_out);return 0;}
fprintf(f,"#Data file %s\n",name_data_in);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Random file %s\n",name_randoms_in);}
fprintf(f,"#Number of data elements used: %d\n",Ndata2);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Number of random elements used: %d\n",Nrand2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol: %d\n",num_effective);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective number of data elements weighted by wcol*wsys: %lf\n",num_effective2);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from data %lf\n",z_effective_data);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Effective redshift value from randoms %lf\n",z_effective_rand);}
fprintf(f,"#Size of the Box %lf Mpc/h\n",L2-L1);
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Area of the survey used %lf deg^2\n",Area_survey);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of alpha: %lf\n",alpha);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization using %s file\n", type_normalization_mode2 );}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Normalization %lf\n",I33);}
fprintf(f,"#Type of Computation: %s\n",type_of_computation);
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of grid cells: %ld\n",ngrid);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Type of Mass Assigment: %s\n",type_of_mass_assigment);}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Type of Yamamoto: %s\n",type_of_yamamoto);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Number of Interlacing steps: %d\n",Ninterlacing);}
if(strcmp(type_of_computation, "FFT") == 0){fprintf(f,"#Grid Correction: %s\n",grid_correction_string);}
if(strcmp(type_of_survey, "cutsky") == 0){fprintf(f,"#Value of Om=%lf\n",Omega_m);}
if(strcmp(type_of_file, "gadget") == 0){fprintf(f,"#RSD: %s\n",RSD); }
fprintf(f,"#Computation using multigrid: %s\n",do_multigrid);
if(strcmp(triangle_shapes, "EQU") == 0){fprintf(f,"#Shapes of the triangles: equilateral\n");}
if(strcmp(triangle_shapes, "ISO") == 0){fprintf(f,"#Shapes of the triangles: isosceles\n");}
if(strcmp(triangle_shapes, "SQU") == 0){fprintf(f,"#Shapes of the triangles: squeezed\n");}
if(strcmp(triangle_shapes, "ALL") == 0){fprintf(f,"#Shapes of the triangles: all\n");}
fclose(f);
}
//return 0;
//Start Power Spectrum Computation for Cutsky
printf("== Computing the Power Spectrum ==\n");

//Determine number of processors available for openmp  
        #pragma omp parallel for private(i,tid) shared(n_lines_parallel,ngrid)
	for(i=0;i<ngrid;i++)
	{
		tid=omp_get_thread_num();
		if(tid==0 && i==0){n_lines_parallel=omp_get_num_threads();}
	}
	i=fftw_init_threads();
	fftw_plan_with_nthreads(n_lines_parallel);
	printf("Number of processors used: %d\n",n_lines_parallel);

//Compute and write the Power Spectrum for FFT+Skycut type of survey.
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0)
{
parameter_value[0]=L1;
parameter_value[1]=L2;
check_box_for_yamamoto(parameter_value,ngrid);
L1=parameter_value[0];
L2=parameter_value[1];

if(strcmp(type_of_yamamoto, "GridCenter") == 0){loop_interlacing_skycut(kmin,kmax,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2, L1, L2, ngrid, P_shot_noise, bin_ps, I22, alpha, mode_correction, n_lines_parallel, binning_type, Hexadecapole_type, name_ps_out, type_of_mass_assigment,do_bispectrum);}

if(strcmp(type_of_yamamoto, "GridAverage") == 0){loop_interlacing_skycut2(kmin,kmax,Ninterlacing, pos_x, pos_y, pos_z, weight,Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand,Nrand2, L1, L2, ngrid, P_shot_noise, bin_ps, I22, alpha, mode_correction, n_lines_parallel, binning_type, Hexadecapole_type, name_ps_out, type_of_mass_assigment,do_bispectrum);}

}

if(strcmp(type_of_computation, "DSY") == 0 || strcmp(type_of_computation, "DSE") == 0)
{

     if( strcmp(type_of_survey, "cutsky") == 0)
     {
         loop_directsum_yamamoto_skycut_caller(kmin,kmax,pos_x, pos_y, pos_z, weight, Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2, L1, L2, P_shot_noise, bin_ps, I22, alpha, n_lines_parallel, binning_type, Hexadecapole_type, name_ps_out,type_of_computation);         
     }
     if(strcmp(type_of_survey, "periodic") == 0)
     {
          printf("No Direct Sum for periodic box at the moment. Exiting now...\n");return 0;
     }

}

//Compute and write the Power Spectrum for FFT+Box with constant line of sight along z.
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "periodic") == 0)
{
if(strcmp(type_of_file, "ascii") == 0)
{
loop_interlacing_periodic(kmin,kmax,Ninterlacing, pos_x, pos_y, pos_z, weight, Ndata, L1, L2, ngrid, P_shot_noise, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out, type_of_mass_assigment,do_bispectrum);
}
if(strcmp(type_of_file, "gadget") == 0)
{
loop_interlacing_periodic_gadget(kmin,kmax,Ninterlacing, name_data_in ,gadget_files, L1, L2, ngrid, bin_ps, mode_correction, n_lines_parallel, binning_type, name_ps_out, type_of_mass_assigment,Shot_noise_factor,grid_correction_string,RSD);
}
}

if(strcmp(do_bispectrum, "no") == 0){
printf("Computation of Power Spectrum finished sucessfully!\n\n");
return 0;
}

printf("== Computing the Bispectrum ==\n");

if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "cutsky") == 0)
{
//write bispectrum header for cutsky
loop_bispectrum_skycut_caller(kmin, kmax, Ninterlacing,  pos_x, pos_y, pos_z, weight, Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2, L1, L2, ngrid, P_shot_noise, Deltakbis, I33,I22, IN, Bsn, alpha, mode_correction, n_lines_parallel, binning_type, name_bs_out, type_of_mass_assigment,triangles_num,write_triangles,triangles_id, do_multigrid, triangle_shapes);
}
if(strcmp(type_of_computation, "FFT") == 0 && strcmp(type_of_survey, "periodic") == 0)
{
if(strcmp(type_of_file, "ascii") == 0)
{
//write bispectrum header for periodic
loop_bispectrum_skycut_caller(kmin, kmax, Ninterlacing,  pos_x, pos_y, pos_z, weight, Ndata2, pos_x_rand, pos_y_rand, pos_z_rand, weight_rand, Nrand2, L1, L2, ngrid, P_shot_noise, Deltakbis, 0,0, 0, 0, 0, mode_correction, n_lines_parallel, binning_type, name_bs_out, type_of_mass_assigment,triangles_num,write_triangles,triangles_id, do_multigrid, triangle_shapes);
}
if(strcmp(type_of_file, "gadget") == 0)
{
loop_bispectrum_periodic_for_gadget_caller(kmin,kmax,Ninterlacing, L1, L2, ngrid, Deltakbis, mode_correction, n_lines_parallel, binning_type, name_bs_out, type_of_mass_assigment, triangles_num, write_triangles, triangles_id, name_data_in,gadget_files, do_multigrid, triangle_shapes,RSD);

}

}
if(strcmp(type_of_computation, "FFT") != 0)
{
printf("No direct Sum for the Bispectrum yet!\n");
return 0;
}



//Improve: sorting algorithm. Avoid system calls

printf("Computation of the Bispectrum finished sucessfully!\n\n");
return 0;
}
