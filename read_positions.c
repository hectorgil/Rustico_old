#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cubature.h"

typedef struct {
	              double OMEGA_M;
}f_params;


void z_to_r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) 
{ 
     f_params params_function = *(f_params *) fdata;//cast from void to f_params 
     double omega= params_function.OMEGA_M; 
     double z=x[0]; 
     double c=299792.458;//speed of light 
     fval[0]=c/(100.*sqrt(omega*pow(1+z,3)+1.-omega)); 
}

long int get_number_used_lines_data(char *filename, double parameter_value[])
{
FILE *f;
long int npar=(int)(parameter_value[3]);
double z_min=parameter_value[1];
double z_max=parameter_value[2];
long int i;
int weight_col;
double redshift;
long int npar_used=0;
int veto;
f=fopen(filename,"r");
for(i=0;i<npar;i++)
{

fscanf(f,"%*f %*f %lf %*f %d %*f %*f\n",&redshift, &weight_col);veto=1;//data
//fscanf(f,"%*f %*f %lf %*f %*f %*f %d\n",&redshift,&veto);weight_col=1;//EZmocks
if(redshift>z_min && redshift<z_max && weight_col>0 && veto==1)
{
npar_used++;
}
}
fclose(f);
//free(f);
return npar_used;
}

long int get_number_used_lines_randoms(char *filename, double parameter_value[])
{
FILE *f;
double redshift;
long int npar=(int)(parameter_value[3]);
double z_min=parameter_value[1];
double z_max=parameter_value[2];
long int i;
int weight_col;
int veto;
long int npar_used=0;
f=fopen(filename,"r");
for(i=0;i<npar;i++)
{
        fscanf(f,"%*f %*f %lf %*f %*f\n",&redshift);veto=1;//data
//          fscanf(f,"%*f %*f %lf %*f %*f %*f %d\n",&redshift,&veto);//EZmocks

if(redshift>z_min && redshift<z_max && veto==1)
{
npar_used++;
}

}
fclose(f);

return npar_used;
}

void get_skycuts_write_density_data(char *filename, double parameter_value[], char *name_den_out)
{
double Pi=(4.*atan(1.));
double Omega_m=parameter_value[0];
double z_min=parameter_value[1];
double z_max=parameter_value[2];
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;

double MIN[1];
double MAX[1];
double r_min,r_max;
int npar,weight_col;
long int i;
FILE *f,*g;
double redshift,weight_fkp,weight_sys, radial,nuissance;
npar=(int)(parameter_value[3]);
double Area;
int n_bin_r;
double DeltaR=parameter_value[25];
int index_radial;
double *radial_cell_haloes, *radial_all_weight_cell_haloes,*radial_fkp_cell_haloes, *radial_weight_cell_haloes;
double *z_cell_haloes;
int veto;
Area=parameter_value[13]*pow(Pi/180.,2);

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);



  n_bin_r=(int)((r_max-r_min)/DeltaR);//printf("\n %lf %d\n",DeltaR,n_bin_r);
  radial_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  radial_all_weight_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  z_cell_haloes= (double*) calloc(n_bin_r, sizeof(double));

f=fopen(filename,"r");
for(i=0;i<npar;i++)
{
fscanf(f,"%*f %*f %lf %lf %d %lf %*f\n",&redshift, &weight_fkp, &weight_col, &weight_sys);veto=1;
//fscanf(f,"%*f %*f %lf %lf %*f %*f %d\n",&redshift,&weight_fkp,&veto);weight_col=1;weight_sys=1;//EZmocks

if(redshift>z_min && redshift<z_max && weight_col>0 && veto==1)
  {

  MAX[0]=redshift;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

      index_radial=(int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%d)\n",radial,redshift,r_min,r_max,i);break;}
      radial_cell_haloes[index_radial]+=1.;
      radial_fkp_cell_haloes[index_radial]+=weight_fkp;
      radial_weight_cell_haloes[index_radial]+=weight_col*weight_sys;
      radial_all_weight_cell_haloes[index_radial]+=weight_col*weight_sys*weight_fkp;
      z_cell_haloes[index_radial]+=redshift;

}

}
printf("\nWriting %s...",name_den_out);
        f=fopen(name_den_out,"w");
        fprintf(f,"#Interval: %lf Mpc/h\n",DeltaR);
        fprintf(f,"# z <nobs> <wc nobs> <wc wfkp nobs>\n");
        for(i=0;i<n_bin_r;i++)
        {
                if(radial_cell_haloes[i]!=0)
                {
                        z_cell_haloes[i]=z_cell_haloes[i]/radial_cell_haloes[i];
                        fprintf(f,"%lf %.16lf %.16lf %.16lf\n",z_cell_haloes[i],radial_cell_haloes[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.),  radial_weight_cell_haloes[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.), radial_all_weight_cell_haloes[i]/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.) );
                }
        }
        fclose(f);
printf("Ok!\n");
free(radial_cell_haloes);
free(radial_all_weight_cell_haloes);
free(radial_fkp_cell_haloes);
free(radial_weight_cell_haloes);
free(z_cell_haloes);

}


void get_skycuts_data(char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[], char *type_normalization_mode)
{
double Omega_m=parameter_value[0];
double z_min=parameter_value[1];
double z_max=parameter_value[2];

f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;

double MIN[1];
double MAX[1];
double Pi=(4.*atan(1.));
int npar,weight_col;
long int i;
FILE *f;
double RA,dec,redshift,weight_fkp,weight_sys,n_z,theta, radial,nuissance,max,min;
npar=(int)(parameter_value[3]);
double Area;
int n_bin_r;
double r_min,r_max;
Area=parameter_value[13]*pow(Pi/180.,2);
double normalization;
double zeff;
double num,num3;
int num2;
long int npar_used;
double Psn_1a;
double Psn_1b;
double Psn_2;
double Bsn_1,Bsn_2;
double alpha_data;
double I22,I22_w_data,I22_w_data_will;
double I33,I33_w_data,I33_w_data_will;
double IN1,IN2;
double DeltaR;
double DeltaR_min;
int index_radial;
double *radial_cell_haloes, *radial_all_weight_cell_haloes,*radial_fkp_cell_haloes, *radial_weight_cell_haloes;
double *radial_weight_cell_haloesN1, *radial_weight_cell_haloesN2, *z_cell;
int i_DeltaR;
double I22_min,I22_w_data_min,I22_w_data_will_min;
double I33_min,I33_w_data_min,I33_w_data_will_min;
double IN1_min,IN2_min;
int veto;
I22_min=0;
I22_w_data_min=0;
I22_w_data_will_min=0;
I33_min=0;
I33_w_data_min=0;
I33_w_data_will_min=0;
IN1_min=0;
IN2_min=0;

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);


i_DeltaR=0;
do
{
normalization=0;
zeff=0;
num=0;
num2=0;
num3=0;
npar_used=0;
Psn_1a=0;
Psn_1b=0;
Psn_2=0;
Bsn_1=0;
Bsn_2=0;
alpha_data=0;

i_DeltaR++;
DeltaR=i_DeltaR*0.5;

  n_bin_r=(int)((r_max-r_min)/DeltaR);//printf("\n %lf %d\n",DeltaR,n_bin_r);
  radial_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  radial_all_weight_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell_haloes = (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell_haloesN1= (double*) calloc(n_bin_r, sizeof(double));
  radial_weight_cell_haloesN2= (double*) calloc(n_bin_r, sizeof(double));
  z_cell= (double*) calloc(n_bin_r, sizeof(double));

  max=-9999999;
  min=9999999;
f=fopen(filename,"r");
for(i=0;i<npar;i++)
{
fscanf(f,"%lf %lf %lf %lf %d %lf %lf\n", &RA, &dec, &redshift, &weight_fkp, &weight_col, &weight_sys, &n_z);veto=1;
//fscanf(f,"%lf %lf %lf %lf %*f %lf %d\n",&RA,&dec,&redshift,&weight_fkp,&n_z,&veto);weight_col=1;weight_sys=1;//EZmocks

theta=90.-dec;
if(redshift>z_min && redshift<z_max && weight_col>0 && veto==1)
  {
  normalization+=n_z*weight_fkp*weight_fkp*weight_sys*weight_col;//effective normalization from data. Factor I22
  zeff+=redshift*weight_col*weight_fkp*weight_sys;//Effective redshift
  num+=weight_col*weight_sys;//Effective number of objects with wsys (real number)
  num2+=weight_col;//Effective number of objects without wsys (integer number)
  num3+=weight_col*weight_fkp*weight_sys;
//printf("%d %.10lf %.10lf %.10lf %d\n",i,normalization,zeff,num,num2);
  
//from deg to radiants
  RA=RA*Pi/180.;
  dec=dec*Pi/180.;
  theta=theta*Pi/180.;

  //Determination of comoving distance given the redshifts and the value of omegamatter
  MAX[0]=redshift;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

      index_radial=(int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%d)\n",radial,redshift,r_min,r_max,i);break;}   
      radial_cell_haloes[index_radial]+=1.;
      radial_fkp_cell_haloes[index_radial]+=weight_fkp;
      radial_weight_cell_haloes[index_radial]+=weight_col*weight_sys;
      radial_all_weight_cell_haloes[index_radial]+=weight_col*weight_sys*weight_fkp;
      radial_weight_cell_haloesN1[index_radial]+=pow(weight_col*weight_sys,2);
      radial_weight_cell_haloesN2[index_radial]+=weight_col*pow(weight_sys,2);
      z_cell[index_radial]+=redshift;


  //From polar to cartesian coordinates
  pos_x[npar_used]=radial*sin(theta)*cos(RA);
  pos_y[npar_used]=radial*sin(theta)*sin(RA);
  pos_z[npar_used]=radial*cos(theta);
 
  weight[npar_used]=weight_fkp*weight_col*weight_sys;//total weight
//  if(npar_used<6){printf("%lf %lf %lf %lf\n",pos_x[npar_used],pos_y[npar_used],pos_z[npar_used],weight[npar_used]);}

  Psn_1a+=pow(weight_fkp*weight_col*weight_sys,2);//false pairs
  Psn_1b+=weight_fkp*weight_col*weight_sys*weight_fkp*weight_sys;// true pairs
  Psn_2+=weight_fkp*weight_fkp*weight_col*weight_sys;//same quantity for both false and true pairs
  Bsn_1+=pow(weight_fkp*weight_col*weight_sys,3);
  Bsn_2+=pow(weight_fkp,3)*weight_col*weight_sys;


  if(pos_x[npar_used]>max || npar_used==0){max=pos_x[npar_used];}
  if(pos_y[npar_used]>max){max=pos_y[npar_used];}
  if(pos_z[npar_used]>max){max=pos_z[npar_used];}
  if(pos_x[npar_used]<min || npar_used==0){min=pos_x[npar_used];}
  if(pos_y[npar_used]<min){min=pos_y[npar_used];}
  if(pos_z[npar_used]<min){min=pos_z[npar_used];}
 
 alpha_data+=weight_col*weight_sys*weight_fkp; 
  npar_used++;

  }
  
  }
fclose(f);


        I22=0;I33=0;
        I22_w_data=0;I33_w_data=0;
        I22_w_data_will=0;I33_w_data_will=0;
        IN1=0;IN2=0;

        for(i=0;i<n_bin_r;i++)
        {
         I22+=Area*pow(radial_all_weight_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;
         I33+=Area*pow(radial_all_weight_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

          if(radial_cell_haloes[i]!=0)
          {
//Normalization of the Power Spectrum
             I22_w_data+=Area*pow(radial_fkp_cell_haloes[i]/radial_cell_haloes[i]*1.,2)*pow(radial_weight_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;
             I22_w_data_will+=Area*pow(radial_all_weight_cell_haloes[i]/radial_cell_haloes[i]*1.,2)*pow(radial_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

//Normalization of the Bispectrum
             I33_w_data+=Area*pow(radial_fkp_cell_haloes[i]/radial_cell_haloes[i]*1.,3)*pow(radial_weight_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

             I33_w_data_will+=Area*pow(radial_all_weight_cell_haloes[i]/radial_cell_haloes[i]*1.,3)*pow(radial_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

   
//Shot noise of the Bispectrum
                        IN2+=Area*pow(radial_fkp_cell_haloes[i]/radial_cell_haloes[i]*1.,3)*(radial_weight_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*(radial_weight_cell_haloesN2[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

                        IN1+=Area*pow(radial_fkp_cell_haloes[i]/radial_cell_haloes[i]*1.,3)*(radial_weight_cell_haloes[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*(radial_weight_cell_haloesN1[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ))*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;


                }
        }

free(radial_cell_haloes);
free(radial_all_weight_cell_haloes);
free(radial_fkp_cell_haloes);
free(radial_weight_cell_haloes);
free(radial_weight_cell_haloesN1);
free(radial_weight_cell_haloesN2);
free(z_cell);


if(i_DeltaR==1)//
{
I22_min=I22;
I22_w_data_min=I22_w_data;
I22_w_data_will_min=I22_w_data_will;
I33_min=I33;
I33_w_data_min=I33_w_data;
I33_w_data_will_min=I33_w_data_will;
IN1_min=IN1;
IN2_min=IN2;
}

if(I22<I22_min){I22_min=I22;DeltaR_min=DeltaR;}
if(I22_w_data<I22_w_data_min){I22_w_data_min=I22_w_data;}
if(I22_w_data_will<I22_w_data_will_min){I22_w_data_will_min=I22_w_data_will;}

if(I33<I33_min){I33_min=I33;}
if(I33_w_data<I33_w_data_min){I33_w_data_min=I33_w_data;}
if(I33_w_data_will<I33_w_data_will_min){I33_w_data_will_min=I33_w_data_will;}

if(IN1<IN1_min){IN1_min=IN1;}
if(IN2<IN2_min){IN2_min=IN2;}

//printf("%lf %lf %.10lf %lf %lf\n",DeltaR,I22,I33,IN1,IN2);
//printf("%lf %lf %.10lf %lf %lf\n",DeltaR,I22_min,I33_min,IN1_min,IN2_min);


}while(DeltaR<40 && strcmp(type_normalization_mode, "area") == 0);

if(strcmp(type_normalization_mode, "density") == 0){DeltaR_min=10.;}//Returns a fixed value if normalized by density

free(function_parameters);
//Copy needed information 
parameter_value[3]=npar_used*1.;
parameter_value[4]=Psn_1a;
parameter_value[5]=Psn_1b;
parameter_value[6]=Psn_2;
parameter_value[7]=zeff/num3*1.;
parameter_value[8]=num2;
parameter_value[9]=normalization;
parameter_value[10]=min;
parameter_value[11]=max;
parameter_value[12]=alpha_data;

parameter_value[14]=I22_min;
parameter_value[15]=I22_w_data_min;
parameter_value[16]=I22_w_data_will_min;
parameter_value[17]=num;
parameter_value[18]=I33_min;
parameter_value[19]=I33_w_data_min;
parameter_value[20]=I33_w_data_will_min;
parameter_value[21]=Bsn_1;
parameter_value[22]=Bsn_2;
parameter_value[23]=IN1_min;
parameter_value[24]=IN2_min;
parameter_value[25]=DeltaR_min;

//printf("tot %.10lf %.10lf %.10lf %d\n",normalization,zeff/num*1.,num,num2);


}

void get_skycuts_write_density_randoms(char *filename, double parameter_value[], double alpha, char *name_den_out)
{
double Omega_m=parameter_value[0];
double z_min=parameter_value[1];
double z_max=parameter_value[2];
f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;
double MIN[1];
double MAX[1];
double Pi=(4.*atan(1.));
double Area=parameter_value[13]*pow(Pi/180.,2);
double DeltaR=parameter_value[25];
long int i,npar;
FILE *f;
double redshift,weight_fkp, radial,nuissance;
double *radial_cell,*radial_fkp_cell,*z_cell;
double r_min,r_max;
int n_bin_r;
int index_radial;
int veto;
npar=(int)(parameter_value[3]);

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);

  n_bin_r=(int)((r_max-r_min)/DeltaR);//printf("\n %d\n",n_bin_r);
  radial_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell = (double*) calloc(n_bin_r, sizeof(double));
  z_cell = (double*) calloc(n_bin_r, sizeof(double));

f=fopen(filename,"r");
for(i=0;i<npar;i++)
{
        fscanf(f,"%*f %*f %lf %lf %*f\n",&redshift, &weight_fkp);veto=1;
//      fscanf(f,"%*f %*f %lf %lf %*f %*f %d\n",&redshift,&weight_fkp,&veto);//EZmocks

        if(redshift>z_min && redshift<z_max && veto==1)
        {
            MAX[0]=redshift;
            MIN[0]=0;
            adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

 index_radial=(int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%d)\n",radial,redshift,r_min,r_max,i);}
      radial_cell[index_radial]+=1.;
      radial_fkp_cell[index_radial]+=weight_fkp;
      z_cell[index_radial]+=redshift;

        }

}

         printf("\nWriting %s...",name_den_out);
         f=fopen(name_den_out,"w");
        fprintf(f,"#Interval: %lf Mpc/h\n",DeltaR);
        fprintf(f,"# alpha<ns> alpha<wfkp ns>\n");
        for(i=0;i<n_bin_r;i++)
        {
                if(radial_cell[i]!=0)
                {
                        z_cell[i]=z_cell[i]/radial_cell[i];
                        fprintf(f,"%lf %.16lf %.16lf \n",z_cell[i], radial_cell[i]*1.*alpha/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.), radial_fkp_cell[i]*1.*alpha/(Area*pow(r_min+i*(r_max-r_min)/n_bin_r*1.,2)*(r_max-r_min)/n_bin_r*1.));
                }
        }
        fclose(f);
printf("Ok!\n");


free(radial_cell);
free(radial_fkp_cell);
free(z_cell);

}

void get_skycuts_randoms(char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[],char *type_normalization_mode, char *type_normalization_mode2)
{
double Omega_m=parameter_value[0];
double z_min=parameter_value[1];
double z_max=parameter_value[2];

f_params *function_parameters;
function_parameters = (f_params *) malloc(sizeof(f_params));
(*function_parameters).OMEGA_M=Omega_m;

double MIN[1];
double MAX[1];
double Pi=(4.*atan(1.));
double Area=parameter_value[13]*pow(Pi/180.,2);

long int i,npar;
FILE *f;
double RA,dec,redshift,weight_fkp,n_z,theta, radial,nuissance,max,min;
npar=(int)(parameter_value[3]);

double normalization;
double zeff;
long int npar_used;
double alpha_data;
double alpha;//alpha count from data

double I22_w_randoms;
double I33_w_randoms;
double DeltaR;
double DeltaR_min;
int index_radial;
double *radial_cell,*radial_fkp_cell;
double r_min,r_max;
int i_DeltaR;
double I22_w_randoms_min;
double I33_w_randoms_min;
double numwfkp;
int n_bin_r;
int veto;
I22_w_randoms_min=0;
I33_w_randoms_min=0;

  MAX[0]=z_min;
  MIN[0]=0;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_min, &nuissance);
  MAX[0]=z_max;
  adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &r_max, &nuissance);



i_DeltaR=0;
do
{
normalization=0;
zeff=0;
numwfkp=0;
npar_used=0;
alpha_data=0;
alpha=parameter_value[12];

i_DeltaR++;
DeltaR=i_DeltaR*0.5;


  n_bin_r=(int)((r_max-r_min)/DeltaR);//printf("\n %d\n",n_bin_r);
  radial_cell = (double*) calloc(n_bin_r, sizeof(double));
  radial_fkp_cell = (double*) calloc(n_bin_r, sizeof(double));
  max=-9999999;
  min=9999999;
f=fopen(filename,"r");
for(i=0;i<npar;i++)
{
	fscanf(f,"%lf %lf %lf %lf %lf\n", &RA, &dec, &redshift, &weight_fkp, &n_z);veto=1;
//        fscanf(f,"%lf %lf %lf %lf %*f %lf %d\n",&RA,&dec,&redshift,&weight_fkp,&n_z,&veto);//EZmocks

	theta=90.-dec;
	if(redshift>z_min && redshift<z_max && veto==1)
		  {

		    normalization+=n_z*weight_fkp*weight_fkp;//effective normalization from data. Factor I22
          	zeff+=redshift*weight_fkp;//Effective redshift
                numwfkp+=weight_fkp;
			//From deg to rad
			RA=RA*Pi/180.;
			dec=dec*Pi/180.;
			theta=theta*Pi/180.;

            //Determination of comoving distance given the redshifts and the value of omegamatter
            MAX[0]=redshift;
            MIN[0]=0;
	        adapt_integrate(1, z_to_r , function_parameters, 1, MIN, MAX ,100000, 1e-6, 1e-6, &radial, &nuissance);

           //From polar to cartesian coordinates
		   pos_x[npar_used]=radial*sin(theta)*cos(RA);
		   pos_y[npar_used]=radial*sin(theta)*sin(RA);
		   pos_z[npar_used]=radial*cos(theta);

      index_radial=(int)(n_bin_r*(radial-r_min)/(r_max-r_min));
      if( index_radial<0 || index_radial>n_bin_r-1){printf("error bins (radial) radial=%lf (z=%lf)  (%lf,%lf) (line=%d)\n",radial,redshift,r_min,r_max,i);}
      radial_cell[index_radial]+=1.;
      radial_fkp_cell[index_radial]+=weight_fkp;

           weight[npar_used]=weight_fkp;//total weight
//           if(npar_used<6){printf("%lf %lf %lf %lf\n",pos_x[npar_used],pos_y[npar_used],pos_z[npar_used],weight[npar_used]);}


		    if(pos_x[npar_used]>max || npar_used==0){max=pos_x[npar_used];}
			if(pos_y[npar_used]>max){max=pos_y[npar_used];}
			if(pos_z[npar_used]>max){max=pos_z[npar_used];}
			if(pos_x[npar_used]<min || npar_used==0){min=pos_x[npar_used];}
			if(pos_y[npar_used]<min){min=pos_y[npar_used];}
			if(pos_z[npar_used]<min){min=pos_z[npar_used];}

			alpha_data+=weight_fkp;
		

					    npar_used++;
                               //printf("%ld %lf %lf %.20lf\n",npar_used,RA,dec, redshift);              
						  }

	  }
fclose(f);
alpha*=1./alpha_data;//printf("\n %lf \n",alpha);
I22_w_randoms=0;
I33_w_randoms=0;

        for(i=0;i<n_bin_r;i++)
        {         
                if(radial_cell[i]!=0)
                {
                   I22_w_randoms+=pow(radial_fkp_cell[i]*1./radial_cell[i]*1.,2)*Area*pow(alpha,2)*pow(radial_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),2)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;
                   I33_w_randoms+=pow(radial_fkp_cell[i]*1./radial_cell[i]*1.,3)*Area*pow(alpha,3)*pow(radial_cell[i]*n_bin_r/(r_max-r_min)*1./( (r_min+i*(r_max-r_min)/n_bin_r )*(r_min+i*(r_max-r_min)/n_bin_r)*Area  ),3)*pow(r_min+i*(r_max-r_min)/n_bin_r,2)*(r_max-r_min)/n_bin_r;

                }
         }

if(i_DeltaR==1)//
{
I22_w_randoms_min=I22_w_randoms;
I33_w_randoms_min=I33_w_randoms;
}

if(I22_w_randoms<I22_w_randoms_min){I22_w_randoms_min=I22_w_randoms;DeltaR_min=DeltaR;}
if(I33_w_randoms<I33_w_randoms_min){I33_w_randoms_min=I33_w_randoms;}

free(radial_cell);
free(radial_fkp_cell);

}while(DeltaR<40 && strcmp(type_normalization_mode, "area") == 0 && strcmp(type_normalization_mode2, "randoms") == 0);

if(strcmp(type_normalization_mode, "density") == 0){DeltaR_min=10.;}//If no area is used for normalization, keep deltaR by data
if(strcmp(type_normalization_mode, "area") == 0 && strcmp(type_normalization_mode2, "data") == 0){DeltaR_min=parameter_value[25];}//if no randoms are explored, keep deltaR by data


//Copy needed information 
parameter_value[3]=npar_used*1.;
parameter_value[7]=zeff/numwfkp*1.;
parameter_value[9]=normalization;
parameter_value[10]=min;
parameter_value[11]=max;
parameter_value[12]=alpha_data;

parameter_value[14]=I22_w_randoms_min;
parameter_value[15]=I22_w_randoms_min;
parameter_value[16]=I22_w_randoms_min;

parameter_value[17]=I33_w_randoms_min;
parameter_value[18]=I33_w_randoms_min;
parameter_value[19]=I33_w_randoms_min;
parameter_value[25]=DeltaR_min;

free(function_parameters);
//free(f);
}


void get_periodic_data(char *filename_data, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[])
{
long int i,npar;
double max,min;
FILE *f;
npar=(int)(parameter_value[3]);
f=fopen(filename_data,"r");
for(i=0;i<npar;i++)
{
fscanf(f,"%lf %lf %lf %lf\n", &pos_x[i], &pos_y[i], &pos_z[i], &weight[i]);

      if(pos_x[i]>max || i==0){max=pos_x[i];}
      if(pos_y[i]>max){max=pos_y[i];}
      if(pos_z[i]>max){max=pos_z[i];}
      if(pos_x[i]<min || i==0){min=pos_x[i];}
      if(pos_y[i]<min){min=pos_y[i];}
      if(pos_z[i]<min){min=pos_z[i];}

}
fclose(f);
parameter_value[10]=min;
parameter_value[11]=max;
//free(f);
}
