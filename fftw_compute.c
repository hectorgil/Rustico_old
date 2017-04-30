#include <stdio.h>
#include <stdlib.h>
#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>
#include <math.h>
#include "functions.h"
#define indexGG(n1,n2,Ngal) (n1)*(Ngal)+(n2)
#define indexGR(ngal,nrand,Nran) (ngal)*(Nran)+(nran)

struct io_header_1//Header for the reader
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int           flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;

struct particle_data
{
  float  Pos[3];
  float  Vel[3];
  float   Mass;
  int    Type;

} *P;

int load_snapshot(char *fname, int files, double *params)
{
        FILE *fd;
        char   buf[200];
        int    i,j,k,dummy,ntot_withmasses;
        int    t,n,off,pc,pc_new,pc_sph;
        int NumPart,Ngas;  
        int *Id;
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

        for(i=0, pc=0; i<files; i++, pc=pc_new)
    {
                if(files>1)
                        sprintf(buf,"%s.%ld",fname,i);
                else
                        sprintf(buf,"%s",fname);

                if(!(fd=fopen(buf,"r")))
                {
                        printf("can't open file `%s`\n",buf);
                        exit(EXIT_FAILURE);
                }


                printf("Reading `%s' ...\n",buf); fflush(stdout);

                fread(&dummy, sizeof(dummy), 1, fd);
                fread(&header1, sizeof(header1), 1, fd);
                fread(&dummy, sizeof(dummy), 1, fd);

                if(files==1)
                {
                        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                                NumPart+= header1.npart[k];
                        Ngas= header1.npart[0];
                }
                else
                {
                        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                                NumPart+= header1.npartTotal[k];//Suma del numero total de particulas
                        Ngas= header1.npartTotal[0];//Suma del numero total de particulas tipo 0 (gas)
                }

                for(k=0, ntot_withmasses=0; k<5; k++)
                {
                        if(header1.mass[k]==0)
                                ntot_withmasses+= header1.npart[k];
                }

                if(i==0)
{
        printf("allocating memory...\n");

        if( P!=NULL)
        {
        free(P);
        P=NULL;
        printf("P memory freed\n");
        }

        if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
                fprintf(stderr,"failed to allocate memory.\n");
                exit(EXIT_FAILURE);
    }


        Id=malloc(NumPart*sizeof(int));
        printf("done\n");

}
         SKIP;
                for(k=0,pc_new=pc;k<6;k++)
                {
                        for(n=0;n<header1.npart[k];n++)
                        {
                                fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);//Lee las posiciones
                                pc_new++;
                        }
                }
                SKIP;

                SKIP;
                for(k=0,pc_new=pc;k<6;k++)
                {
                        for(n=0;n<header1.npart[k];n++)
                        {
                                fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);//Lee las velocidades
                                pc_new++;
                        }

                }
                SKIP;

                SKIP;
                for(k=0,pc_new=pc;k<6;k++)
                {
                        for(n=0;n<header1.npart[k];n++)
                        {
                                fread(&Id[pc_new], sizeof(int), 1, fd);
                                pc_new++;
                        }
                }
                SKIP;

                fclose(fd);
    }

params[0]=NumPart;
params[1]=header1.redshift;
params[2]=header1.Omega0;
params[3]=header1.OmegaLambda;
free(Id);
}

long int count_particles_gadget(char *name_data_in)
{
        long int Ndata;
        FILE *fd;
        int dummy,k,ntot_withmasses;
         int NumPart,Ngas;

        fd=fopen(name_data_in,"r");
        if(fd==NULL){printf("Error reading %s file. Exiting now...\n",name_data_in);exit(EXIT_FAILURE);}
                fread(&dummy, sizeof(dummy), 1, fd);
                fread(&header1, sizeof(header1), 1, fd);
                fread(&dummy, sizeof(dummy), 1, fd);

                        for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
                                NumPart+= header1.npart[k];
                        Ngas= header1.npart[0];

                for(k=0, ntot_withmasses=0; k<5; k++)
                {
                        if(header1.mass[k]==0)
                                ntot_withmasses+= header1.npart[k];
                }


Ndata=NumPart;
return Ndata;
}

void loop_directsum_exact_skycut(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, int ngrid, double P_shot_noise, double Deltak, double I22, double alpha, int n_lines_parallel, char *binning_type, char *Hexadecapole_type, char *name_ps_out)
{
int tid;
long int j,i,i1,i2;
int xindex,yindex,zindex;
long int NGRID;//Number of k-modes between kmin and kmax
double KXX,KYY,KZZ;
double kdotx1,ckdotx1,skdotx1,musq,kampsq,kdotx2,ckdotx2,skdotx2;

double Pi=(4.*atan(1.));

double **monopole_real_rr;
double **monopole_real_gr;
double **monopole_real_gg;
double **quadrupole_real_rr;
double **quadrupole_real_gr;
double **quadrupole_real_gg;
double **hexadecapole_real_rr;
double **hexadecapole_real_gr;
double **hexadecapole_real_gg;
long int ngridtot=pow(ngrid,3);
double *P0, *P2, *P4;
double *kx;

double xam12;

double *KX;
double *KY;
double *KZ;
int hexadecapole_sw;
double SX,SY,SZ; 
if(strcmp(Hexadecapole_type, "L2L2") == 0){hexadecapole_sw=0;}
if(strcmp(Hexadecapole_type, "L4") == 0){hexadecapole_sw=1;}

kx=malloc(sizeof(double)*(ngrid));
for(i=0;i<ngrid;i++)
{
       if(i<ngrid/2+1)
       {
          kx[i]=i*1.0*(2.0*Pi/(L2-L1));
       }
       else
       {
          kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
       }
}
i=0;
for(j=0;j<ngridtot;j++)//Count the number of grid modes btween kmax and kmin
{
xindex=(int)(j/(ngrid*ngrid*1.));
yindex=(int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
SX=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(SX<=kmax*kmax && SX>kmin*kmin && SZ>=0)
{
i++;
}
}

NGRID=i;
printf("Computing aproximately %.0lf^3 (out of %d^3) number of k-modes in the range %lf < k < %lf (original range %lf < k < %lf) \n",cbrt(NGRID),ngrid,kmin,kmax, kx[0],kx[ngrid/2]);
KX=malloc(sizeof(double)*NGRID);
KY=malloc(sizeof(double)*NGRID);
KZ=malloc(sizeof(double)*NGRID);

i=0;
for(j=0;j<ngridtot;j++)
{
xindex=(int)(j/(ngrid*ngrid*1.));
yindex=(int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
SX=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(SX<=kmax*kmax && SX>kmin*kmin && SZ>=0)
{
KX[i]=kx[xindex];
KY[i]=kx[yindex];
KZ[i]=kx[zindex];
i++;
}
}
free(kx);

//monopole_real_rr= (double **)calloc(NGRID,sizeof(double*));
//monopole_real_gr= (double **)calloc(NGRID,sizeof(double*));
//monopole_real_gg= (double **)calloc(NGRID,sizeof(double*));

quadrupole_real_rr= (double **)calloc(NGRID,sizeof(double*));
quadrupole_real_gr= (double **)calloc(NGRID,sizeof(double*));
quadrupole_real_gg= (double **)calloc(NGRID,sizeof(double*));

if(hexadecapole_sw==1)
{
hexadecapole_real_rr= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_real_gr= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_real_gg= (double **)calloc(NGRID,sizeof(double*));

}

    for(i=0;i<NGRID;i++)
   {
  //  monopole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
  //  monopole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
  //  monopole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));

    quadrupole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));

if(hexadecapole_sw==1)
{   hexadecapole_real_rr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_real_gr[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_real_gg[i] = (double*)calloc(n_lines_parallel,sizeof(double));
}
    }


#pragma omp parallel for private(j,i1,i2,tid,KXX,KYY,KZZ,kampsq,kdotx1,kdotx2,ckdotx1,ckdotx2,skdotx1,skdotx2,xam12,musq) shared(NGRID,KX,KY,KZ,Ndata,Nrand,s_x,s_y,s_z,weight,s_x_ran,s_y_ran,s_z_ran,weight_ran,hexadecapole_sw,monopole_real_gg,quadrupole_real_gg,hexadecapole_real_gg,monopole_real_gr,quadrupole_real_gr,hexadecapole_real_gr,monopole_real_rr,quadrupole_real_rr,hexadecapole_real_rr)
        for(j=0;j<NGRID;j++)
        {
                tid=omp_get_thread_num();//thread number
                KXX=KX[j];
                KYY=KY[j];
                KZZ=KZ[j];
                kampsq = KXX*KXX+KYY*KYY+KZZ*KZZ;
                for(i1=0;i1<Ndata;i1++)
                {
                     kdotx1=KXX*s_x[i1]+KYY*s_y[i1]+KZZ*s_z[i1];
                     ckdotx1=cos(kdotx1)*weight[i1];
                     skdotx1=sin(kdotx1)*weight[i1];

                     for(i2=i1;i2<Ndata;i2++)
                     {
                       kdotx2=KXX*s_x[i2]+KYY*s_y[i2]+KZZ*s_z[i2];
                       ckdotx2=cos(kdotx2)*weight[i2];
                       skdotx2=sin(kdotx2)*weight[i2];

                       xam12=0.25*(pow(s_x[i1]+s_x[i2],2)+pow(s_y[i1]+s_y[i2],2)+pow(s_z[i1]+s_z[i2],2));
                       musq=0.25*(kdotx1*kdotx1+kdotx2*kdotx2+2.*kdotx2*kdotx1)/(xam12*kampsq);

    //                   monopole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2);
                       quadrupole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq;
                       if(hexadecapole_sw==1)
                       {
                          hexadecapole_real_gg[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq*musq;
                       }

                     }
                     for(i2=0;i2<Nrand;i2++)
                     {
                       kdotx2=KXX*s_x_ran[i2]+KYY*s_y_ran[i2]+KZZ*s_z_ran[i2];
                       ckdotx2=cos(kdotx2)*weight_ran[i2];
                       skdotx2=sin(kdotx2)*weight_ran[i2];
                   
                       xam12=0.25*(pow(s_x[i1]+s_x_ran[i2],2)+pow(s_y[i1]+s_y_ran[i2],2)+pow(s_z[i1]+s_z_ran[i2],2));
                       musq=0.25*(kdotx1*kdotx1+kdotx2*kdotx2+2.*kdotx2*kdotx1)/(xam12*kampsq);

      //                 monopole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2);
                       quadrupole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq;
                       if(hexadecapole_sw==1)
                       {
                          hexadecapole_real_gr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq*musq;
                       }
                }
             }
                   for(i1=0;i1<Nrand;i1++)
                   {
                     kdotx1=KXX*s_x_ran[i1]+KYY*s_y_ran[i1]+KZZ*s_z_ran[i1];
                     ckdotx1=cos(kdotx1)*weight_ran[i1];
                     skdotx1=sin(kdotx1)*weight_ran[i1];
                     for(i2=i1;i2<Nrand;i2++)
                     {
                       kdotx2=KXX*s_x_ran[i2]+KYY*s_y_ran[i2]+KZZ*s_z_ran[i2];
                       ckdotx2=cos(kdotx2)*weight_ran[i2];
                       skdotx2=sin(kdotx2)*weight_ran[i2];
     
                       xam12=0.25*(pow(s_x_ran[i1]+s_x_ran[i2],2)+pow(s_y_ran[i1]+s_y_ran[i2],2)+pow(s_z_ran[i1]+s_z_ran[i2],2));
                       musq=0.25*(kdotx1*kdotx1+kdotx2*kdotx2+2.*kdotx2*kdotx1)/(xam12*kampsq);

        //               monopole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2);
                       quadrupole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq;
                       if(hexadecapole_sw==1)
                       { 
                         hexadecapole_real_rr[j][tid]+=(ckdotx1*ckdotx2+skdotx1*skdotx2)*musq*musq;
                       }
                      }
                   }


            }


        for(i=0;i<NGRID;i++)
        {
               for(j=1;j<n_lines_parallel;j++)
                {
          //              monopole_real_gg[i][0]+=monopole_real_gg[i][j];
                        quadrupole_real_gg[i][0]+=quadrupole_real_gg[i][j];
                        if(hexadecapole_sw==1){
                        hexadecapole_real_gg[i][0]+=hexadecapole_real_gg[i][j];
                        }

            //            monopole_real_gr[i][0]+=monopole_real_gr[i][j];
                        quadrupole_real_gr[i][0]+=quadrupole_real_gr[i][j];
                        if(hexadecapole_sw==1){
                        hexadecapole_real_gr[i][0]+=hexadecapole_real_gr[i][j];
                        }

              //          monopole_real_rr[i][0]+=monopole_real_rr[i][j];
                        quadrupole_real_rr[i][0]+=quadrupole_real_rr[i][j];
                        if(hexadecapole_sw==1){
                        hexadecapole_real_rr[i][0]+=hexadecapole_real_rr[i][j];
                        }

                }

                //       monopole_real_gg[i][0]+=(alpha*alpha*monopole_real_rr[i][0]-alpha*monopole_real_gr[i][0]);
                       quadrupole_real_gg[i][0]+=(alpha*alpha*quadrupole_real_rr[i][0]-alpha*quadrupole_real_gr[i][0]);
                       if(hexadecapole_sw==1){
                       hexadecapole_real_gg[i][0]+=(alpha*alpha*hexadecapole_real_rr[i][0]-alpha*hexadecapole_real_gr[i][0]);
                       }


        }

//freeTokens(monopole_real_rr,NGRID);
//freeTokens(monopole_real_gr,NGRID);

freeTokens(quadrupole_real_rr,NGRID);
freeTokens(quadrupole_real_gr,NGRID);

if(hexadecapole_sw==1){
freeTokens(hexadecapole_real_rr,NGRID);
freeTokens(hexadecapole_real_gr,NGRID);
}

P0= (double *)calloc(NGRID,sizeof(double));
P2= (double *)calloc(NGRID,sizeof(double));
if(hexadecapole_sw==1){
P4= (double *)calloc(NGRID,sizeof(double));
}

#pragma omp parallel for private(i) shared(NGRID,hexadecapole_sw,P0,P2,P4,monopole_real_gg, quadrupole_real_gg,hexadecapole_real_gg)
        for(i=0;i<NGRID;i++)
        {
             P0[i]=0;//monopole_real_gg[i][0];
             P2[i]=quadrupole_real_gg[i][0];
             if(hexadecapole_sw==1){
             P4[i]=hexadecapole_real_gg[i][0];
             }

        }

//freeTokens(monopole_real_gg,NGRID);
freeTokens(quadrupole_real_gg,NGRID);
if(hexadecapole_sw==1){
freeTokens(hexadecapole_real_gg,NGRID);}


printf("Writing Power Spectrum output %s...",name_ps_out);
if(hexadecapole_sw==0){
write_power_spectrum_skyscuts_directsum_exactP2(kmin,kmax,KX, KY, KZ, P0,P2, Deltak,ngrid,NGRID, L1, L2, I22, name_ps_out, P_shot_noise,binning_type);
}
if(hexadecapole_sw==1){
write_power_spectrum_skyscuts_directsum_exactP4(kmin,kmax,KX, KY, KZ, P0,P2,P4, Deltak,ngrid,NGRID, L1, L2, I22, name_ps_out, P_shot_noise,binning_type);
}
printf("Ok!\n");

free(KX);
free(KY);
free(KZ);

}




void loop_directsum_yamamoto_skycut(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double  L1, double L2, int ngrid, double P_shot_noise, double  Deltak, double  I22, double  alpha, int  n_lines_parallel, char *binning_type, char  *Hexadecapole_type, char *name_ps_out)
{
int tid;
long int j,i;
int xindex,yindex,zindex;
long int NGRID;//Number of k-modes between kmin and kmax
double KXX,KYY,KZZ,XAM;
double kdotx,ckdotx,skdotx,musq,kampsq;
double Pi=(4.*atan(1.));

double **monopole_real;
double **monopole_imag;
double **quadrupole_real;
double **quadrupole_imag;
double **hexadecapole_real;
double **hexadecapole_imag;

double **monopole_real_rand;
double **monopole_imag_rand;
double **quadrupole_real_rand;
double **quadrupole_imag_rand;
double **hexadecapole_real_rand;
double **hexadecapole_imag_rand;
long int ngridtot=pow(ngrid,3);

double *deltak_re0, *deltak_im0, *deltak_re2, *deltak_im2, *deltak_re4, *deltak_im4;
double *kx;
double *xampsq;
double *xampsq_ran;
double *KX;
double *KY;
double *KZ;
int hexadecapole_sw;
double SX,SY,SZ;
if(strcmp(Hexadecapole_type, "L2L2") == 0){hexadecapole_sw=0;}
if(strcmp(Hexadecapole_type, "L4") == 0){hexadecapole_sw=1;}


        xampsq=(double*)malloc(sizeof(double)*Ndata);
        xampsq_ran=(double*)malloc(sizeof(double)*Nrand);
for(i=0;i<Ndata;i++){
SX=s_x[i]*s_x[i];
SY=s_y[i]*s_y[i];
SZ=s_z[i]*s_z[i];
xampsq[i]=SX+SY+SZ;
}
for(i=0;i<Nrand;i++){
SX=s_x_ran[i]*s_x_ran[i];
SY=s_y_ran[i]*s_y_ran[i];
SZ=s_z_ran[i]*s_z_ran[i];
xampsq_ran[i]=SX+SY+SZ;
}


kx=malloc(sizeof(double)*(ngrid));
for(i=0;i<ngrid;i++)
{
       if(i<ngrid/2+1)
       {
          kx[i]=i*1.0*(2.0*Pi/(L2-L1));
       }
       else
       {
          kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
       }
}
i=0;
for(j=0;j<ngridtot;j++)//Count the number of grid modes btween kmax and kmin
{
xindex=(int)(j/(ngrid*ngrid*1.));
yindex=(int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
SX=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(SX<=kmax*kmax && SX>kmin*kmin && SZ>=0)
{
i++;
}
}

NGRID=i;
printf("Computing aproximately %.0lf^3 (out of %d^3) number of k-modes in the range %lf < k < %lf (original range %lf < k < %lf) \n",cbrt(NGRID),ngrid,kmin,kmax, kx[0],kx[ngrid/2]);
KX=malloc(sizeof(double)*NGRID);
KY=malloc(sizeof(double)*NGRID);
KZ=malloc(sizeof(double)*NGRID);

i=0;
for(j=0;j<ngridtot;j++)
{
xindex=(int)(j/(ngrid*ngrid*1.));
yindex=(int)( (j-xindex*ngrid*ngrid)/(ngrid*1.));
zindex=j-xindex*ngrid*ngrid-yindex*ngrid;
SX=kx[xindex]*kx[xindex]+kx[yindex]*kx[yindex]+kx[zindex]*kx[zindex];
SZ=kx[zindex];
if(SX<=kmax*kmax && SX>kmin*kmin && SZ>=0)
{
KX[i]=kx[xindex];
KY[i]=kx[yindex];
KZ[i]=kx[zindex];
i++;
}
}
free(kx);

monopole_real= (double **)calloc(NGRID,sizeof(double*));
monopole_imag= (double **)calloc(NGRID,sizeof(double*));
quadrupole_real= (double **)calloc(NGRID,sizeof(double*));
quadrupole_imag= (double **)calloc(NGRID,sizeof(double*));
if(hexadecapole_sw==1)
{
hexadecapole_real= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_imag= (double **)calloc(NGRID,sizeof(double*));
}

monopole_real_rand= (double **)calloc(NGRID,sizeof(double*));
monopole_imag_rand= (double **)calloc(NGRID,sizeof(double*));
quadrupole_real_rand= (double **)calloc(NGRID,sizeof(double*));
quadrupole_imag_rand= (double **)calloc(NGRID,sizeof(double*));
if(hexadecapole_sw==1)
{
hexadecapole_real_rand= (double **)calloc(NGRID,sizeof(double*));
hexadecapole_imag_rand= (double **)calloc(NGRID,sizeof(double*));
}
    for(i=0;i<NGRID;i++)
   {
    monopole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    monopole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));
if(hexadecapole_sw==1)
{   hexadecapole_real[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_imag[i] = (double*)calloc(n_lines_parallel,sizeof(double));
}
    monopole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    monopole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    quadrupole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
if(hexadecapole_sw==1) 
{   hexadecapole_real_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
    hexadecapole_imag_rand[i] = (double*)calloc(n_lines_parallel,sizeof(double));
}
  }

        #pragma omp parallel for private(j,tid,KXX,KYY,KZZ,kampsq,i,XAM,kdotx,musq,ckdotx,skdotx) shared(NGRID,KX,KY,KZ,Ndata,Nrand,xampsq,s_x,s_y,s_z,monopole_real,monopole_imag,quadrupole_real,quadrupole_imag,hexadecapole_real,hexadecapole_imag,weight,xampsq_ran,s_x_ran,s_y_ran,s_z_ran,weight_ran, monopole_real_rand,monopole_imag_rand,quadrupole_real_rand,quadrupole_imag_rand,hexadecapole_real_rand,hexadecapole_imag_rand,hexadecapole_sw)
        for(j=0;j<NGRID;j++)
        {
                tid=omp_get_thread_num();//thread number
                KXX=KX[j];
                KYY=KY[j];
                KZZ=KZ[j];
                kampsq = KXX*KXX+KYY*KYY+KZZ*KZZ;
        for(i=0;i<Ndata;i++)
        {
                        XAM=xampsq[i];
                        kdotx=KXX*s_x[i]+KYY*s_y[i]+KZZ*s_z[i];
                        musq = kdotx*kdotx/(kampsq*XAM);
                        ckdotx=cos(kdotx)*weight[i];
                        skdotx=sin(kdotx)*weight[i];

                       monopole_real[j][tid]+=ckdotx;
                       monopole_imag[j][tid]+=skdotx;
                       quadrupole_real[j][tid]+=ckdotx*musq;
                       quadrupole_imag[j][tid]+=skdotx*musq;
                       if(hexadecapole_sw==1)
                       {hexadecapole_real[j][tid]+=ckdotx*musq*musq;
                       hexadecapole_imag[j][tid]+=skdotx*musq*musq;}

                        XAM=xampsq_ran[i];
                        kdotx=KXX*s_x_ran[i]+KYY*s_y_ran[i]+KZZ*s_z_ran[i];
                        musq = kdotx*kdotx/(kampsq*XAM);
                        ckdotx=cos(kdotx)*weight_ran[i];
                        skdotx=sin(kdotx)*weight_ran[i];
                      
                        monopole_real_rand[j][tid]+=ckdotx;
                        monopole_imag_rand[j][tid]+=skdotx;
                        quadrupole_real_rand[j][tid]+=ckdotx*musq;
                        quadrupole_imag_rand[j][tid]+=skdotx*musq;
                        if(hexadecapole_sw==1)
                        {hexadecapole_real_rand[j][tid]+=ckdotx*musq*musq;
                        hexadecapole_imag_rand[j][tid]+=skdotx*musq*musq;}
         }
         for(i=Ndata;i<Nrand;i++)
         {

                        XAM=xampsq_ran[i];
                        kdotx=KXX*s_x_ran[i]+KYY*s_y_ran[i]+KZZ*s_z_ran[i];
                        musq = kdotx*kdotx/(kampsq*XAM);
                        ckdotx=cos(kdotx)*weight_ran[i];
                        skdotx=sin(kdotx)*weight_ran[i];

                        monopole_real_rand[j][tid]+=ckdotx;
                        monopole_imag_rand[j][tid]+=skdotx;
                        quadrupole_real_rand[j][tid]+=ckdotx*musq;
                        quadrupole_imag_rand[j][tid]+=skdotx*musq;
                        if(hexadecapole_sw==1)
                        {hexadecapole_real_rand[j][tid]+=ckdotx*musq*musq;
                        hexadecapole_imag_rand[j][tid]+=skdotx*musq*musq;}
         }

        }



free(xampsq);
free(xampsq_ran);

alpha=-alpha;
 for(i=0;i<NGRID;i++)
        {
                for(j=1;j<n_lines_parallel;j++)
                {

                        monopole_real[i][0]+=monopole_real[i][j];
                        monopole_imag[i][0]+=monopole_imag[i][j];
                        quadrupole_real[i][0]+=quadrupole_real[i][j];
                        quadrupole_imag[i][0]+=quadrupole_imag[i][j];
                       if(hexadecapole_sw==1){
                        hexadecapole_real[i][0]+=hexadecapole_real[i][j];
                        hexadecapole_imag[i][0]+=hexadecapole_imag[i][j];}

                        monopole_real_rand[i][0]+=monopole_real_rand[i][j];
                        monopole_imag_rand[i][0]+=monopole_imag_rand[i][j];
                        quadrupole_real_rand[i][0]+=quadrupole_real_rand[i][j];
                        quadrupole_imag_rand[i][0]+=quadrupole_imag_rand[i][j];
                       if(hexadecapole_sw==1){
                        hexadecapole_real_rand[i][0]+=hexadecapole_real_rand[i][j];
                        hexadecapole_imag_rand[i][0]+=hexadecapole_imag_rand[i][j];}



                }
                       monopole_real[i][0]+=alpha*monopole_real_rand[i][0];
                       monopole_imag[i][0]+=alpha*monopole_imag_rand[i][0];
                       quadrupole_real[i][0]+=alpha*quadrupole_real_rand[i][0];
                       quadrupole_imag[i][0]+=alpha*quadrupole_imag_rand[i][0];
                       if(hexadecapole_sw==1){
                       hexadecapole_real[i][0]+=alpha*hexadecapole_real_rand[i][0];
                       hexadecapole_imag[i][0]+=alpha*hexadecapole_imag_rand[i][0];}


        }
freeTokens(monopole_real_rand,NGRID);
freeTokens(monopole_imag_rand,NGRID);
freeTokens(quadrupole_real_rand,NGRID);
freeTokens(quadrupole_imag_rand,NGRID);
if(hexadecapole_sw==1){
freeTokens(hexadecapole_real_rand,NGRID);
freeTokens(hexadecapole_imag_rand,NGRID);}

deltak_re0= (double *)calloc(NGRID,sizeof(double));
deltak_im0= (double *)calloc(NGRID,sizeof(double));
deltak_re2= (double *)calloc(NGRID,sizeof(double));
deltak_im2= (double *)calloc(NGRID,sizeof(double));
if(hexadecapole_sw==1){
deltak_re4= (double *)calloc(NGRID,sizeof(double));
deltak_im4= (double *)calloc(NGRID,sizeof(double));
}

#pragma omp parallel for private(i) shared(NGRID,hexadecapole_sw,deltak_re0,deltak_re2,deltak_re4, deltak_im0,deltak_im2,deltak_im4,monopole_real,monopole_imag,quadrupole_real, quadrupole_imag,hexadecapole_real,hexadecapole_imag)
 for(i=0;i<NGRID;i++)
        {
             deltak_re0[i]=monopole_real[i][0];
             deltak_im0[i]=monopole_imag[i][0];
             deltak_re2[i]=quadrupole_real[i][0];
             deltak_im2[i]=quadrupole_imag[i][0];
             if(hexadecapole_sw==1){
             deltak_re4[i]=hexadecapole_real[i][0];
             deltak_im4[i]=hexadecapole_imag[i][0];
             }

        }
freeTokens(monopole_real,NGRID);
freeTokens(monopole_imag,NGRID);
freeTokens(quadrupole_real,NGRID);
freeTokens(quadrupole_imag,NGRID);
if(hexadecapole_sw==1){
freeTokens(hexadecapole_real,NGRID);
freeTokens(hexadecapole_imag,NGRID);}


printf("Writing Power Spectrum output %s...",name_ps_out);
if(hexadecapole_sw==0){
write_power_spectrum_skyscuts_directsum_L2L2(kmin,kmax,KX, KY, KZ, deltak_re0, deltak_im0, deltak_re2, deltak_im2, Deltak,ngrid,NGRID, L1, L2, I22, name_ps_out, P_shot_noise,binning_type);
}
if(hexadecapole_sw==1){
write_power_spectrum_skyscuts_directsum_L4(kmin,kmax,KX, KY, KZ, deltak_re0, deltak_im0, deltak_re2, deltak_im2, deltak_re4, deltak_im4, Deltak, ngrid, NGRID, L1, L2, I22, name_ps_out, P_shot_noise,binning_type);
}
printf("Ok!\n");

free(KX);
free(KY);
free(KZ);

}


void loop_directsum_yamamoto_skycut_caller(double kmin,double kmax, double *s_x, double *s_y, double *s_z, double *weight, long int Ndata, double *s_x_ran, double *s_y_ran, double *s_z_ran, double *weight_ran, long int Nrand, double L1, double L2, double P_shot_noise, double Deltak, double I22, double alpha, int n_lines_parallel, char *binning_type, char *Hexadecapole_type, char *name_ps_out, char *type_of_computation)
{

double *kmin_i;
double *kmax_i;
int i;
int Nbins;
int ngrid;
double Pi=(4.*atan(1.));

if(strcmp(binning_type, "linear") == 0){Nbins=(int)((kmax-kmin)/Deltak)+1;}
if(strcmp(binning_type, "log10") == 0){Nbins=(int)((log10(kmax)-log10(kmin))/Deltak)+1;}

kmin_i= (double*)calloc(Nbins,sizeof(double));
kmax_i= (double*)calloc(Nbins,sizeof(double));

for(i=0;i<Nbins;i++)
{
//linear
if(strcmp(binning_type, "linear") == 0)
{
kmin_i[i]=kmin+Deltak*i;
kmax_i[i]=kmin+Deltak*(i+1);
}
//log
if(strcmp(binning_type, "log10") == 0)
{
kmin_i[i]=pow(10,log10(kmin)+Deltak*i);
kmax_i[i]=pow(10,log10(kmin)+Deltak*(i+1));
}

//Determine the appropiate value of ngrid
ngrid=1;
do
{

ngrid=ngrid*2;

}while(kmax_i[i]>Pi*ngrid/(L2-L1));

if(strcmp(type_of_computation, "DSY") == 0)
{
loop_directsum_yamamoto_skycut(kmin_i[i],kmax_i[i], s_x, s_y, s_z, weight, Ndata, s_x_ran, s_y_ran, s_z_ran, weight_ran, Nrand, L1, L2, ngrid, P_shot_noise, Deltak, I22, alpha, n_lines_parallel, binning_type, Hexadecapole_type, name_ps_out);
}
if(strcmp(type_of_computation, "DSE") == 0)
{
loop_directsum_exact_skycut(kmin_i[i],kmax_i[i], s_x, s_y, s_z, weight, Ndata, s_x_ran, s_y_ran, s_z_ran, weight_ran, Nrand, L1, L2, ngrid, P_shot_noise, Deltak, I22, alpha, n_lines_parallel, binning_type, Hexadecapole_type, name_ps_out);
}

}

free(s_x);
free(s_y);
free(s_z);
free(weight);
free(s_x_ran);
free(s_y_ran);
free(s_z_ran);
free(weight_ran);
free(kmin_i);
free(kmax_i);
}

void modify_input_for_yamamoto(int mode_yamamoto, double in[], int ngrid, double L1, double L2)
{
int i,j,k;
long int c;
double kmod2;
long int ngridtot=pow(ngrid,3);
        #pragma omp parallel for private(c,i,j,k,kmod2) shared(ngrid,ngridtot,L1,L2,in)
		for(c=0;c<ngridtot;c++)
		{
			i=(int)(c/(ngrid*ngrid*1.));
			j=(int)( (c-i*ngrid*ngrid)/(ngrid*1.));
			k=c-i*ngrid*ngrid-j*ngrid;
                        
            kmod2=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)+(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)+(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.);
        if(mode_yamamoto==1 && kmod2>0)
        {   in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/(kmod2);	    
        }
        if(mode_yamamoto==2 && kmod2>0)
        {   in[c]*=(L1+j*(L2-L1)/ngrid*1.)/(L1+i*(L2-L1)/ngrid*1.);
        }
        if(mode_yamamoto==3 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)/(L1+j*(L2-L1)/ngrid*1.);
        }
        if(mode_yamamoto==4 && kmod2>0)
        {   in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==5 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==6 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==7 && kmod2>0)
        {   in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/(kmod2*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        }
        if(mode_yamamoto==8 && kmod2>0)
        {   in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==9 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==10 && kmod2>0)
        {   in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==11 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==12 && kmod2>0)
        {   in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==13 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==14 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==15 && kmod2>0)
        {   in[c]*=(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==16 && kmod2>0) 
        {   in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==17 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==18 && kmod2>0)
        {   in[c]*=(L1+j*(L2-L1)/ngrid*1.)*(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==19 && kmod2>0)
        {   in[c]*=(L1+i*(L2-L1)/ngrid*1.)*(L1+i*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.)*(L1+k*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==20 && kmod2>0)
        {   in[c]*=(L1+j*(L2-L1)/ngrid*1.)/((L1+i*(L2-L1)/ngrid*1.));
        } 
        if(mode_yamamoto==21 && kmod2>0)
        {   in[c]*=(L1+k*(L2-L1)/ngrid*1.)/((L1+j*(L2-L1)/ngrid*1.));
        } 

              }

}


void fftw_yamamoto_skycut(int mode_yamamoto, double in[], double deltak_re[], double deltak_im[], int ngrid, double L1, double L2, int mode_mass_ass)
{
  fftw_complex *out;
  fftw_plan p;
  int i,j,k;
  long int c;
  double cx,cy,cz;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)

  long int index2;
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
			if(i<ngrid/2+1)
			{
				kx[i]=i*1.0*(2.0*Pi/(L2-L1));
			}
			else
			{
				kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
			}
	}

    out =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(ngridtotr2c));
    p =  fftw_plan_dft_r2c_3d(ngrid,ngrid,ngrid,in,out,FFTW_ESTIMATE);

    fftw_execute(p);//FFT
    fftw_destroy_plan(p);

   #pragma omp parallel for private(c,i,j,k,cx,cy,cz,index2) shared(kx,deltak_re,deltak_im,out,ngrid,ngridtotr2c,mode_mass_ass,mode_yamamoto,Pi)
  for(index2=0;index2<ngridtotr2c;index2++)
  {


i=(int)(index2/(ngrid*ngrid/2+ngrid));
j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

        cx=sin( kx[i]*Pi/(2.*kx[ngrid/2]) )/( kx[i]*Pi/(2.*kx[ngrid/2]) );
        cy=sin( kx[j]*Pi/(2.*kx[ngrid/2]) )/( kx[j]*Pi/(2.*kx[ngrid/2]) );
        cz=sin( kx[k]*Pi/(2.*kx[ngrid/2]) )/( kx[k]*Pi/(2.*kx[ngrid/2]) );
        if(kx[i]==0 || mode_mass_ass==0){cx=1.;}
        if(kx[j]==0 || mode_mass_ass==0){cy=1.;}
        if(kx[k]==0 || mode_mass_ass==0){cz=1.;}


if(mode_yamamoto==0){
          deltak_re[index2]=creal(out[index2])*pow(cx*cy*cz,-mode_mass_ass*1.);
          deltak_im[index2]=cimag(out[index2])*pow(cx*cy*cz,-mode_mass_ass*1.);           
}
if(mode_yamamoto==1){
           deltak_re[index2]=(creal(out[index2])*kx[i]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
           deltak_im[index2]=(cimag(out[index2])*kx[i]*kx[i])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==2){
deltak_re[index2]+=(creal(out[index2])*kx[i]*kx[j]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[i]*kx[j]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==3){
deltak_re[index2]+=(creal(out[index2])*kx[i]*kx[k]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[i]*kx[k]*2.)*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==4){
deltak_re[index2]+=(creal(out[index2])*kx[j]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[j]*kx[j])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==5){
deltak_re[index2]+=2.*(creal(out[index2])*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=2.*(cimag(out[index2])*kx[j]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==6){
deltak_re[index2]+=(creal(out[index2])*kx[k]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
deltak_im[index2]+=(cimag(out[index2])*kx[k]*kx[k])*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-1);
}
if(mode_yamamoto==7){
deltak_re[index2]=(creal(out[index2])*pow(kx[i],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]=(cimag(out[index2])*pow(kx[i],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
} 
if(mode_yamamoto==8){
deltak_re[index2]+=(creal(out[index2])*pow(kx[j],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=(cimag(out[index2])*pow(kx[j],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==9){
deltak_re[index2]+=(creal(out[index2])*pow(kx[k],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=(cimag(out[index2])*pow(kx[k],4))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==10){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[i],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[i],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==11){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[i],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[i],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==12){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[j],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[j],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==13){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[j],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[j],3))*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==14){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[k],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[k],3))*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==15){
deltak_re[index2]+=4.*(creal(out[index2])*pow(kx[k],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=4.*(cimag(out[index2])*pow(kx[k],3))*kx[j]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==16){
deltak_re[index2]+=6.*(creal(out[index2])*pow(kx[i]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=6.*(cimag(out[index2])*pow(kx[i]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==17){
deltak_re[index2]+=6.*(creal(out[index2])*pow(kx[i]*kx[k],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=6.*(cimag(out[index2])*pow(kx[i]*kx[k],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==18){
deltak_re[index2]+=6.*(creal(out[index2])*pow(kx[k]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=6.*(cimag(out[index2])*pow(kx[k]*kx[j],2))*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==19){
deltak_re[index2]+=12.*(creal(out[index2])*pow(kx[i],2))*kx[j]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=12.*(cimag(out[index2])*pow(kx[i],2))*kx[j]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==20){
deltak_re[index2]+=12.*(creal(out[index2])*pow(kx[j],2))*kx[i]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=12.*(cimag(out[index2])*pow(kx[j],2))*kx[i]*kx[k]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}
if(mode_yamamoto==21){
deltak_re[index2]+=12.*(creal(out[index2])*pow(kx[k],2))*kx[j]*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
deltak_im[index2]+=12.*(cimag(out[index2])*pow(kx[k],2))*kx[j]*kx[i]*pow(cx*cy*cz,-mode_mass_ass*1.)*pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],-2);
}


}

   
    fftw_free(out);
    free(kx);
fftw_cleanup();
//in not modified
//deltak_re filled
//deltak_im filled
//out is internally created and destroyed
} 




void loop_interlacing_skycut_for_bispectrum(int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand,long int Nrand, double L1, double L2, int ngrid, double alpha, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re0, double *deltak_im0)
{
long int i,j,k,i_inter;
long int c;
double phase_cos,phase_sin;
double new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b;
int i_yama;
double Sk;
long int index2;
double* delta_data;
double* delta_rand;
double* deltak_re0_b;
double* deltak_im0_b;
long int ngridtot=pow(ngrid,3);
long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)
double L2b,L1b;
double *kx;
double Pi=(4.*atan(1.));

        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }


  alpha=-alpha;

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

    delta_data = (double*) calloc(ngridtot, sizeof(double));
    printf("Assigning particles to the grid (Iteration %d) ...", i_inter);
    for(c=0; c<Ndata; c++)
    {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
    }

        delta_rand = (double*) calloc(ngridtot, sizeof(double));
       for(c=0; c<Nrand; c++)
        {
                if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
        }
    printf("Ok!\n");

    printf("Performing FFTs ...");
 #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]+=alpha*delta_rand[c];
        }
        free(delta_rand);

        if(i_inter==1){
       fftw_yamamoto_skycut(0, delta_data, deltak_re0, deltak_im0, ngrid,L1b,L2b, mode_correction);
  }
  else{
        deltak_re0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0, delta_data, deltak_re0_b, deltak_im0_b, ngrid,L1b,L2b, mode_correction);
#pragma omp parallel for private(i,j,k,c,index2,phase_cos,phase_sin,new_deltak_re0_b,new_deltak_im0_b) shared(ngrid,ngridtot,ngridtotr2c,deltak_re0,deltak_im0,deltak_re0_b,deltak_im0_b,kx,L2,L1,i_inter,Ninterlacing)
for(index2=0;index2<ngridtotr2c;index2++)
{
i=(int)(index2/(ngrid*ngrid/2+ngrid));
j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);
Sk=kx[i]+ kx[j]+ kx[k];

phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);
phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);

new_deltak_re0_b=deltak_re0_b[index2]*phase_cos-deltak_im0_b[index2]*phase_sin;
new_deltak_im0_b=deltak_re0_b[index2]*phase_sin+deltak_im0_b[index2]*phase_cos;

deltak_re0[index2]+=new_deltak_re0_b;
deltak_im0[index2]+=new_deltak_im0_b;


}

        free(deltak_re0_b);
        free(deltak_im0_b);
  }

        free(delta_data);
printf("Ok!\n");


}//loop interlacing
free(kx);
//deltak_re and deltak_im ready for bispectrum computation
}


void loop_interlacing_skycut(double kmin,double kmax,int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand,long int Nrand, double L1, double L2, int ngrid,  double P_shot_noise, double bin_ps, double I22, double alpha, int mode_correction, int n_lines_parallel, char *binning_type, char *Hexadecapole_type, char *name_ps_out, char *type_of_mass_assigment, char *do_bispectrum)
{
long int i,j,k,i_inter;
long int c;
double phase_cos,phase_sin;
double new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b;
int i_yama;
double Sk;
long int index2;
  //delta(k) pointers
  double* delta_data;
  double* delta_rand;
  double* deltak_re0;
  double* deltak_im0;
  double* deltak_re2;
  double* deltak_im2;
  double* deltak_re4;
  double* deltak_im4;
  double* deltak_re0_b;
  double* deltak_im0_b;
  double* deltak_re2_b;
  double* deltak_im2_b;
  double* deltak_re4_b;
  double* deltak_im4_b;
long int ngridtot=pow(ngrid,3);
long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));//N*N*(N/2+1)
int hexadecapole_sw;
if(strcmp(Hexadecapole_type, "L2L2") == 0){hexadecapole_sw=0;}
if(strcmp(Hexadecapole_type, "L4") == 0){hexadecapole_sw=1;}
  alpha=-alpha;

  double L2b,L1b;
  
  double *kx;
  double Pi=(4.*atan(1.));

        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }


for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

    delta_data = (double*) calloc(ngridtot, sizeof(double));
    printf("Assigning particles to the grid (Iteration %d) ...", i_inter);
    for(c=0; c<Ndata; c++)
    {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
    }

        delta_rand = (double*) calloc(ngridtot, sizeof(double));
       for(c=0; c<Nrand; c++)
        {
                if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c], L2b, L1b, ngrid);}
        }
    printf("Ok!\n");

    printf("Performing FFTs ...");
 #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand,hexadecapole_sw)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]+=alpha*delta_rand[c];
        }
        free(delta_rand);

  if(i_inter==1){
        deltak_re0= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0=(double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0, delta_data, deltak_re0, deltak_im0, ngrid,L1b,L2b, mode_correction);

         deltak_re2=(double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im2=(double*) calloc(ngridtotr2c,sizeof(double));
for(i_yama=1;i_yama<=6;i_yama++)
{
         modify_input_for_yamamoto(i_yama, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2, deltak_im2, ngrid,L1b,L2b, mode_correction);
}
if(hexadecapole_sw==1){

         deltak_re4= (double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im4= (double*) calloc(ngridtotr2c,sizeof(double));

for(i_yama=7;i_yama<=21;i_yama++)
{

         modify_input_for_yamamoto(i_yama, delta_data, ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4, deltak_im4, ngrid,L1b,L2b, mode_correction);
}

}


  }
  else{  
        deltak_re0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0, delta_data, deltak_re0_b, deltak_im0_b, ngrid,L1b,L2b, mode_correction);

         deltak_re2_b= (double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im2_b= (double*) calloc(ngridtotr2c,sizeof(double));

for(i_yama=1;i_yama<=6;i_yama++)
{
         modify_input_for_yamamoto(i_yama, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2_b, deltak_im2_b, ngrid,L1b,L2b, mode_correction);
}
if(hexadecapole_sw==1){

         deltak_re4_b= (double*) calloc(ngridtotr2c,sizeof(double));
         deltak_im4_b= (double*) calloc(ngridtotr2c,sizeof(double));
for(i_yama=7;i_yama<=21;i_yama++)
{
         modify_input_for_yamamoto(i_yama, delta_data,ngrid,L1,L2);
         fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4_b, deltak_im4_b, ngrid,L1b,L2b, mode_correction);
}

}

#pragma omp parallel for private(i,j,k,c,index2,phase_cos,phase_sin,new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b) shared(ngrid,ngridtot,ngridtotr2c,deltak_re0,deltak_im0,deltak_re2,deltak_im2,deltak_re4,deltak_im4,deltak_re0_b,deltak_im0_b,deltak_re2_b,deltak_im2_b,deltak_re4_b,deltak_im4_b,kx,L2,L1,i_inter,Ninterlacing,hexadecapole_sw)
for(index2=0;index2<ngridtotr2c;index2++)
{
i=(int)(index2/(ngrid*ngrid/2+ngrid));
j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

Sk=kx[i]+ kx[j]+ kx[k];

phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);
phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*Sk);

new_deltak_re0_b=deltak_re0_b[index2]*phase_cos-deltak_im0_b[index2]*phase_sin;
new_deltak_im0_b=deltak_re0_b[index2]*phase_sin+deltak_im0_b[index2]*phase_cos;

new_deltak_re2_b=deltak_re2_b[index2]*phase_cos-deltak_im2_b[index2]*phase_sin;
new_deltak_im2_b=deltak_re2_b[index2]*phase_sin+deltak_im2_b[index2]*phase_cos;

if(hexadecapole_sw==1){
new_deltak_re4_b=deltak_re4_b[index2]*phase_cos-deltak_im4_b[index2]*phase_sin;
new_deltak_im4_b=deltak_re4_b[index2]*phase_sin+deltak_im4_b[index2]*phase_cos;
}

deltak_re0[index2]+=new_deltak_re0_b;
deltak_im0[index2]+=new_deltak_im0_b;

deltak_re2[index2]+=new_deltak_re2_b;
deltak_im2[index2]+=new_deltak_im2_b;

if(hexadecapole_sw==1){
deltak_re4[index2]+=new_deltak_re4_b;
deltak_im4[index2]+=new_deltak_im4_b;
}

}

        free(deltak_re0_b);
        free(deltak_im0_b);
        free(deltak_re2_b);
        free(deltak_im2_b);
if(hexadecapole_sw==1){
         free(deltak_re4_b);
         free(deltak_im4_b);
}

  }


        free(delta_data);
printf("Ok!\n");

if(i_inter==Ninterlacing && strcmp(do_bispectrum, "no") == 0)
{
printf("positions freed\n");
free(pos_x);
free(pos_y);
free(pos_z);
free(weight);
free(pos_x_rand);
free(pos_y_rand);
free(pos_z_rand);
free(weight_rand);
}


}//loop interlacing

free(kx);
printf("Writing Power Spectrum output %s...",name_ps_out);
if(hexadecapole_sw==0){write_power_spectrum_skyscuts_L2L2(kmin,kmax,deltak_re0,deltak_im0,deltak_re2,deltak_im2,bin_ps, ngrid,L1,L2,I22,Ninterlacing,name_ps_out,P_shot_noise,binning_type);}
if(hexadecapole_sw==1){write_power_spectrum_skyscuts_L4(kmin,kmax,deltak_re0,deltak_im0,deltak_re2,deltak_im2,deltak_re4,deltak_im4,bin_ps, ngrid,L1,L2,I22,Ninterlacing,name_ps_out,P_shot_noise, binning_type);}
printf("Ok!\n");

}

void loop_interlacing_skycut2(double kmin, double kmax, int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double *pos_x_rand, double *pos_y_rand, double *pos_z_rand, double *weight_rand, long int Nrand, double L1, double L2, int ngrid,  double P_shot_noise, double bin_ps, double I22, double alpha, int mode_correction, int n_lines_parallel, char *binning_type, char *Hexadecapole_type, char *name_ps_out, char *type_of_mass_assigment, char *do_bispectrum)
{
long int i,j,k,i_inter;
long int c;
int i_yama_max;
i_yama_max;
double phase_cos,phase_sin;
double new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b;
int i_yama;

int hexadecapole_sw;
if(strcmp(Hexadecapole_type, "L2L2") == 0){hexadecapole_sw=0;i_yama_max=6;}
if(strcmp(Hexadecapole_type, "L4") == 0){hexadecapole_sw=1;i_yama_max=21;}
alpha=-alpha;

  //delta(k) pointers
  double* delta_data;
  double* delta_rand;
  double* deltak_re0;
  double* deltak_im0;
  double* deltak_re2;
  double* deltak_im2;
  double* deltak_re4;
  double* deltak_im4;
  double* deltak_re0_b;
  double* deltak_im0_b;
  double* deltak_re2_b;
  double* deltak_im2_b;
  double* deltak_re4_b;
  double* deltak_im4_b;
long int ngridtot=pow(ngrid,3);
long int ngridtotr2c=(pow(ngrid,3)/2+pow(ngrid,2));
long int index2;
double weight_pos;

  double L2b,L1b;
  
  double *kx;
  double Pi=(4.*atan(1.));

        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }


for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

if(i_inter==1)
{
        printf("Assigning particles to the grid and FFTing (Iteration %d) ...", i_inter);

          for(i_yama=0;i_yama<=i_yama_max;i_yama++)
          {
            
        delta_data = (double*) calloc(ngridtot, sizeof(double));
      
        for(c=0; c<Ndata; c++)
        {

if(i_yama==0){weight_pos=1;}
if(i_yama==1){weight_pos=pos_x[c]*pos_x[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==2){weight_pos=pos_x[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==3){weight_pos=pos_x[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==4){weight_pos=pos_y[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==5){weight_pos=pos_y[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==6){weight_pos=pos_z[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==7){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==8){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==9){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==10){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==11){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==12){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==13){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==14){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==15){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==16){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==17){weight_pos=pow(pos_x[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==18){weight_pos=pow(pos_y[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==19){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==20){weight_pos=pow(pos_y[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==21){weight_pos=pow(pos_z[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]==0){weight_pos=0;}
       weight_pos=weight_pos*weight[c];

       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight_pos, L2b, L1b, ngrid);}
    } 

        delta_rand = (double*) calloc(ngridtot, sizeof(double));

      for(c=0;c<Nrand;c++)
      {
if(i_yama==0){weight_pos=1.;}
if(i_yama==1){weight_pos=pos_x_rand[c]*pos_x_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==2){weight_pos=pos_x_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==3){weight_pos=pos_x_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==4){weight_pos=pos_y_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==5){weight_pos=pos_y_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==6){weight_pos=pos_z_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==7){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==8){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==9){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==10){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==11){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==12){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==13){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==14){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==15){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==16){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==17){weight_pos=pow(pos_x_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==18){weight_pos=pow(pos_y_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==19){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==20){weight_pos=pow(pos_y_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==21){weight_pos=pow(pos_z_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]==0){weight_pos=0;}
      weight_pos=weight_pos*weight_rand[c];

  
                if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}
                if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_pos, L2b, L1b, ngrid);}

        }

   

         
        #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand)
        for(c=0;c<ngridtot;c++) 
        { 
                                 delta_data[c]+=alpha*delta_rand[c];
        } 

        free(delta_rand);

if(i_yama==0){
        deltak_re0= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==1)
{
        deltak_re2= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im2= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==7)
{
        deltak_re4= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im4= (double*) calloc(ngridtotr2c,sizeof(double));
}


if(i_yama==0){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re0, deltak_im0, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=1 && i_yama<=6){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2, deltak_im2, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=7){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4, deltak_im4, ngrid,L1b,L2b, mode_correction);}

        free(delta_data);

}

    printf("Ok!\n");

  }
  else{  
        printf("Assigning particles to the grid and FFTing (Iteration %d) ...", i_inter);

        for(i_yama=0;i_yama<=i_yama_max;i_yama++)
          {

        delta_data = (double*) calloc(ngridtot, sizeof(double));

        for(c=0; c<Ndata; c++)
        {

if(i_yama==0){weight_pos=1;}
if(i_yama==1){weight_pos=pos_x[c]*pos_x[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==2){weight_pos=pos_x[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==3){weight_pos=pos_x[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==4){weight_pos=pos_y[c]*pos_y[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==5){weight_pos=pos_y[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==6){weight_pos=pos_z[c]*pos_z[c]/(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]);}
if(i_yama==7){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==8){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==9){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==10){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==11){weight_pos=pow(pos_x[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==12){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==13){weight_pos=pow(pos_y[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==14){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_x[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==15){weight_pos=pow(pos_z[c],2)*pos_z[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==16){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==17){weight_pos=pow(pos_x[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==18){weight_pos=pow(pos_y[c],2)*pos_z[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==19){weight_pos=pow(pos_x[c],2)*pos_y[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==20){weight_pos=pow(pos_y[c],2)*pos_x[c]*pos_z[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(i_yama==21){weight_pos=pow(pos_z[c],2)*pos_x[c]*pos_y[c]/pow(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c],2);}
if(pos_x[c]*pos_x[c]+pos_y[c]*pos_y[c]+pos_z[c]*pos_z[c]==0){weight_pos=0;}

       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c]*weight_pos, L2b, L1b, ngrid);}
    }

        delta_rand = (double*) calloc(ngridtot, sizeof(double));
       for(c=0; c<Nrand; c++)
        {

if(i_yama==0){weight_pos=1;}
if(i_yama==1){weight_pos=pos_x_rand[c]*pos_x_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==2){weight_pos=pos_x_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==3){weight_pos=pos_x_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==4){weight_pos=pos_y_rand[c]*pos_y_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==5){weight_pos=pos_y_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==6){weight_pos=pos_z_rand[c]*pos_z_rand[c]/(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]);}
if(i_yama==7){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==8){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==9){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==10){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==11){weight_pos=pow(pos_x_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==12){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==13){weight_pos=pow(pos_y_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==14){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_x_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==15){weight_pos=pow(pos_z_rand[c],2)*pos_z_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==16){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==17){weight_pos=pow(pos_x_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==18){weight_pos=pow(pos_y_rand[c],2)*pos_z_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==19){weight_pos=pow(pos_x_rand[c],2)*pos_y_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==20){weight_pos=pow(pos_y_rand[c],2)*pos_x_rand[c]*pos_z_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(i_yama==21){weight_pos=pow(pos_z_rand[c],2)*pos_x_rand[c]*pos_y_rand[c]/pow(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c],2);}
if(pos_x_rand[c]*pos_x_rand[c]+pos_y_rand[c]*pos_y_rand[c]+pos_z_rand[c]*pos_z_rand[c]==0){weight_pos=0;}

              if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
              if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_rand, pos_x_rand[c], pos_y_rand[c], pos_z_rand[c], weight_rand[c]*weight_pos, L2b, L1b, ngrid);}
        }


        #pragma omp parallel for private(c) shared(ngrid,ngridtot,alpha,delta_data,delta_rand)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]+=alpha*delta_rand[c];
        }
        free(delta_rand);

if(i_yama==0){
        deltak_re0_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im0_b= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==1)
{
        deltak_re2_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im2_b= (double*) calloc(ngridtotr2c,sizeof(double));
}
if(i_yama==7)
{
        deltak_re4_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im4_b= (double*) calloc(ngridtotr2c,sizeof(double));
}

if(i_yama==0){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re0_b, deltak_im0_b, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=1 && i_yama<=6){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re2_b, deltak_im2_b, ngrid,L1b,L2b, mode_correction);}
if(i_yama>=7){fftw_yamamoto_skycut(i_yama, delta_data, deltak_re4_b, deltak_im4_b, ngrid,L1b,L2b, mode_correction);}
free(delta_data);
}


#pragma omp parallel for private(i,j,k,c,index2,phase_cos,phase_sin,new_deltak_re0_b,new_deltak_im0_b,new_deltak_re2_b,new_deltak_im2_b,new_deltak_re4_b,new_deltak_im4_b) shared(ngrid,ngridtot,ngridtotr2c,deltak_re0,deltak_im0,deltak_re2,deltak_im2,deltak_re4,deltak_im4,deltak_re0_b,deltak_im0_b,deltak_re2_b,deltak_im2_b,deltak_re4_b,deltak_im4_b,kx,L2,L1,i_inter,Ninterlacing,hexadecapole_sw)
for(index2=0;index2<ngridtotr2c;index2++)
{
i=(int)(index2/(ngrid*ngrid/2+ngrid));
j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);


phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));

new_deltak_re0_b=deltak_re0_b[index2]*phase_cos-deltak_im0_b[index2]*phase_sin;
new_deltak_im0_b=deltak_re0_b[index2]*phase_sin+deltak_im0_b[index2]*phase_cos;

new_deltak_re2_b=deltak_re2_b[index2]*phase_cos-deltak_im2_b[index2]*phase_sin;
new_deltak_im2_b=deltak_re2_b[index2]*phase_sin+deltak_im2_b[index2]*phase_cos;

if(hexadecapole_sw==1){
new_deltak_re4_b=deltak_re4_b[index2]*phase_cos-deltak_im4_b[index2]*phase_sin;
new_deltak_im4_b=deltak_re4_b[index2]*phase_sin+deltak_im4_b[index2]*phase_cos;
}

deltak_re0[index2]+=new_deltak_re0_b;
deltak_im0[index2]+=new_deltak_im0_b;

deltak_re2[index2]+=new_deltak_re2_b;
deltak_im2[index2]+=new_deltak_im2_b;

if(hexadecapole_sw==1){
deltak_re4[index2]+=new_deltak_re4_b;
deltak_im4[index2]+=new_deltak_im4_b;
}


}

        free(deltak_re0_b);
        free(deltak_im0_b);
        free(deltak_re2_b);
        free(deltak_im2_b);
if(hexadecapole_sw==1){
         free(deltak_re4_b);
         free(deltak_im4_b);
}

printf("Ok!\n");

  }


if(i_inter==Ninterlacing && strcmp(do_bispectrum, "no") == 0)
{
//printf("positions freed\n");
free(pos_x);
free(pos_y);
free(pos_z);
free(weight);
free(pos_x_rand);
free(pos_y_rand);
free(pos_z_rand);
free(weight_rand);

}

}
free(kx);
printf("Writing Power Spectrum output %s...",name_ps_out);
if(hexadecapole_sw==0){write_power_spectrum_skyscuts_L2L2(kmin,kmax,deltak_re0,deltak_im0,deltak_re2,deltak_im2,bin_ps, ngrid,L1,L2,I22,Ninterlacing,name_ps_out,P_shot_noise,binning_type);}
if(hexadecapole_sw==1){write_power_spectrum_skyscuts_L4(kmin,kmax,deltak_re0,deltak_im0,deltak_re2,deltak_im2,deltak_re4,deltak_im4,bin_ps, ngrid,L1,L2,I22,Ninterlacing,name_ps_out,P_shot_noise,binning_type);}
printf("Ok!\n");

}



void loop_interlacing_periodic_for_bispectrum(int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double L1, double L2, int ngrid, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re, double *deltak_im)
{
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin;
  double* deltak_re_b;
  double* deltak_im_b;
  double new_deltak_re_b,new_deltak_im_b;
  double* delta_data;
  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
     L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
     L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));
     printf("Assigning particles to the grid (Iteration %d) ...", i_inter);
     for(c=0; c<Ndata; c++)
     {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}

     }

        #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,delta_data)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]=delta_data[c]/Ndata-pow(ngrid,-3);
        }



     printf("Ok!\n");
     printf("Performing FFTs ...");
     if(i_inter==1)
     {
        fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
        printf("Ok!\n");

     }
     else
     {
        deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
      #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
      for(index2=0;index2<ngridtotr2c;index2++)
      {

            i=(int)(index2/(ngrid*ngrid/2+ngrid));
            j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
            k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);
             phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
             new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
             deltak_re[index2]+=new_deltak_re_b;
             deltak_im[index2]+=new_deltak_im_b;
}

        free(deltak_re_b);
        free(deltak_im_b);

        printf("Ok!\n");


}

}//end of itteration loop
free(kx);

}

void loop_interlacing_periodic(double kmin, double kmax, int Ninterlacing, double *pos_x, double *pos_y, double *pos_z, double *weight, long int Ndata, double L1, double L2, int ngrid, double P_shot_noise, double bin_ps, int mode_correction, int n_lines_parallel, char *binning_type, char *name_ps_out, char *type_of_mass_assigment,char *do_bispectrum)
{
  
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin; 
  double* deltak_re;
  double* deltak_im;
  double* deltak_re_b;
  double* deltak_im_b;
  double new_deltak_re_b,new_deltak_im_b;
  double* delta_data;
  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{ 
     L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
     L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));
     printf("Assigning particles to the grid (Iteration %d) ...", i_inter);
     for(c=0; c<Ndata; c++)
     {
       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight[c], L2b, L1b, ngrid);}

     }


        #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,delta_data)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]=delta_data[c]/Ndata-pow(ngrid,-3);
        }
        


     printf("Ok!\n");
     printf("Performing FFTs ...");
     if(i_inter==1)
     {
        deltak_re= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
        printf("Ok!\n");

     }
     else
     {
        deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
      #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
      for(index2=0;index2<ngridtotr2c;index2++)
      {

            i=(int)(index2/(ngrid*ngrid/2+ngrid));
            j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
            k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

             phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
             new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
             deltak_re[index2]+=new_deltak_re_b;
             deltak_im[index2]+=new_deltak_im_b;

}

        free(deltak_re_b);
        free(deltak_im_b);
        
        printf("Ok!\n");

if(i_inter==Ninterlacing && strcmp(do_bispectrum, "no") == 0 )//only free them in the last interation and if no-bispectrum is computed
{
//printf("positions freed\n");
free(pos_x);
free(pos_y);
free(pos_z);
free(weight);
}


}

}//end of itteration loop
free(kx);
printf("Writing Power Spectrum output %s...",name_ps_out);
write_power_spectrum_periodic(kmin,kmax,deltak_re,deltak_im,bin_ps, ngrid,L1,L2,Ninterlacing,name_ps_out,P_shot_noise,binning_type);

printf("Ok!\n");
}


void loop_interlacing_periodic_gadget(double kmin,double kmax, int Ninterlacing, char *name_data_in ,int gadget_files, double L1, double L2, int ngrid, double bin_ps, int mode_correction, int n_lines_parallel, char *binning_type, char *name_ps_out, char *type_of_mass_assigment,double Shot_noise_factor,char *grid_correction_string, char *RSD)
{
  FILE *f;
  double *pos_x,*pos_y,*pos_z;
  double weight;
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin;
  double* deltak_re;
  double* deltak_im;
  double* deltak_re_b;
  double* deltak_im_b;
  double new_deltak_re_b,new_deltak_im_b;
  double* delta_data;
  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
  int snapshot_num;
  char name_gadget_file[500];
  long int NumPart_file;
  long int Ndata;
  double P_shot_noise;
  double scale_factor;
  double Omatter,Olambda;
  double params[4];
        kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

weight=1.0;
for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
     L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
     L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));
     printf("Assigning particles to the grid (Iteration %d) ...", i_inter);

Ndata=0;
for(snapshot_num=0;snapshot_num<gadget_files;snapshot_num++)
  {

  sprintf(name_gadget_file, "%s.%d", name_data_in, snapshot_num);

  load_snapshot(name_gadget_file, 1, params);
NumPart_file=(long int)(params[0]);
scale_factor=1./(1+params[1]);
Omatter=params[2];
Olambda=params[3];

Ndata=Ndata+NumPart_file;
 pos_x = (double*) calloc(NumPart_file, sizeof(double));
 pos_y = (double*) calloc(NumPart_file, sizeof(double));
 pos_z = (double*) calloc(NumPart_file, sizeof(double));

     for(c=0; c<NumPart_file; c++)
     {

     pos_x[c]=P[c].Pos[0]*0.001;
     pos_y[c]=P[c].Pos[1]*0.001;
if(strcmp(RSD, "no") == 0){pos_z[c]=P[c].Pos[2]*0.001;}
if(strcmp(RSD, "yes") == 0){pos_z[c]=(P[c].Pos[2]*0.001)+(P[c].Vel[2])*sqrt(scale_factor)/(100.*scale_factor*sqrt(Omatter*pow(scale_factor,-3)+Olambda));}


       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}

     }
free(pos_x);
free(pos_y);
free(pos_z);

printf("Ok!\n");
} 

        #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,delta_data)
     for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]=delta_data[c]/Ndata-pow(ngrid,-3);
        }



     printf("Ok!\n");
     printf("Performing FFTs ...");
   if(i_inter==1)
     {
        deltak_re= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
        printf("Ok!\n");

     }
     else
     {
        deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
      #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
    for(index2=0;index2<ngridtotr2c;index2++)
      {

            i=(int)(index2/(ngrid*ngrid/2+ngrid));
            j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
            k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);

             phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
             new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
             deltak_re[index2]+=new_deltak_re_b;
             deltak_im[index2]+=new_deltak_im_b;
} 

        free(deltak_re_b);
        free(deltak_im_b);

        printf("Ok!\n");


} 

}//end of itteration loop
free(kx);

printf("Writing Power Spectrum output %s...",name_ps_out);
P_shot_noise=pow(L2-L1,3)/Ndata;

write_power_spectrum_periodic(kmin,kmax,deltak_re,deltak_im,bin_ps, ngrid,L1,L2,Ninterlacing,name_ps_out,P_shot_noise,binning_type);

printf("Ok!\n");

}


void loop_interlacing_periodic_gadget_for_bispectrum(int Ninterlacing, double L1, double L2, int ngrid, int n_lines_parallel, char *type_of_mass_assigment, int mode_correction, double *deltak_re, double *deltak_im, char *name_data_in,int gadget_files, double *params_input, char *RSD)
{
  double *pos_x,*pos_y,*pos_z;
  long int Ndata;
  double weight=1.0;
  long int i,j,k,i_inter;
  long int c;
  double phase_cos,phase_sin;
  double* deltak_re_b;
  double* deltak_im_b;
  double new_deltak_re_b,new_deltak_im_b;
  double* delta_data;
  double L2b,L1b;
  double *kx;
  double Pi=(4.*atan(1.));
  long int ngridtot=pow(ngrid,3);
  long int ngridtotr2c=pow(ngrid,3)/2+pow(ngrid,2);
  long int index2;
  int snapshot_num;
  char name_gadget_file[500];
  int NumPart_file;
  double P_shot_noise;
  double scale_factor,Omatter,Olambda;
  double params[4];

       kx=malloc(sizeof(double)*(ngrid));
        for(i=0;i<ngrid;i++)
        {
                        if(i<ngrid/2+1)
                        {
                                kx[i]=i*1.0*(2.0*Pi/(L2-L1));
                        }
                        else
                        {
                                kx[i]=-(ngrid-i)*1.0*(2.0*Pi/(L2-L1));
                        }
        }

for(i_inter=1;i_inter<=Ninterlacing;i_inter++)
{
     L2b=L2-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);
     L1b=L1-(L2-L1)/ngrid*1.*1./Ninterlacing*1.*(i_inter-1);

     delta_data = (double*) calloc(ngridtot, sizeof(double));
     printf("Assigning particles to the grid (Iteration %d) ...", i_inter);

Ndata=0;
for(snapshot_num=0;snapshot_num<gadget_files;snapshot_num++)
  {

  sprintf(name_gadget_file, "%s.%d", name_data_in, snapshot_num);

  load_snapshot(name_gadget_file, 1, params);
NumPart_file=(long int)(params[0]);
scale_factor=1./(1+params[1]);
Omatter=params[2];
Olambda=params[3];

Ndata=Ndata+NumPart_file;
 pos_x = (double*) calloc(NumPart_file, sizeof(double));
 pos_y = (double*) calloc(NumPart_file, sizeof(double));
 pos_z = (double*) calloc(NumPart_file, sizeof(double));

printf("%ld\n",NumPart_file);
     for(c=0; c<NumPart_file; c++)
     {
     pos_x[c]=P[c].Pos[0]*0.001; //conversion to Mpc/h (originally in kpc/h)
     pos_y[c]=P[c].Pos[1]*0.001;
if(strcmp(RSD, "no") == 0){pos_z[c]=P[c].Pos[2]*0.001;}
if(strcmp(RSD, "yes") == 0){pos_z[c]=(P[c].Pos[2]*0.001)+(P[c].Vel[2])*sqrt(scale_factor)/(100.*scale_factor*sqrt(Omatter*pow(scale_factor,-3)+Olambda));}
        

       if(strcmp(type_of_mass_assigment, "NGC") == 0){ngc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "CIC") == 0){cic_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "TSC") == 0){tsc_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "PCS") == 0){pcs_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P4S") == 0){pq4s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}
       if(strcmp(type_of_mass_assigment, "P5S") == 0){pq5s_assingment(delta_data, pos_x[c], pos_y[c], pos_z[c], weight, L2b, L1b, ngrid);}

     }
free(pos_x);
free(pos_y);
free(pos_z);
printf("Ok!\n");
}
       #pragma omp parallel for private(c) shared(ngrid,ngridtot,Ndata,delta_data)
        for(c=0;c<ngridtot;c++)
        {
                                 delta_data[c]=delta_data[c]/Ndata-pow(ngrid,-3);
        }



     printf("Ok!\n");
     printf("Performing FFTs ...");
     if(i_inter==1)
     {
        fftw_yamamoto_skycut(0,delta_data, deltak_re, deltak_im, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
        printf("Ok!\n");

     }
     else
    {
        deltak_re_b= (double*) calloc(ngridtotr2c,sizeof(double));
        deltak_im_b= (double*) calloc(ngridtotr2c,sizeof(double));
        fftw_yamamoto_skycut(0,delta_data, deltak_re_b, deltak_im_b, ngrid,L1b,L2b, mode_correction);
        free(delta_data);
      #pragma omp parallel for private(index2,c,i,j,k,phase_cos,phase_sin,new_deltak_re_b,new_deltak_im_b) shared(ngrid,ngridtot,ngridtotr2c,kx,deltak_re,deltak_im,deltak_re_b,deltak_im_b,i_inter,Ninterlacing,L1,L2)
      for(index2=0;index2<ngridtotr2c;index2++)
      {

            i=(int)(index2/(ngrid*ngrid/2+ngrid));
            j=(int)( (index2-i*(ngrid*ngrid/2+ngrid))/(ngrid/2+1) );
            k=index2-i*(ngrid*ngrid/2+ngrid)-j*(ngrid/2+1);
             phase_cos=cos(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             phase_sin=sin(((L2-L1)/ngrid*1.)*(i_inter*1.-1.)/Ninterlacing*1.*(kx[i]+ kx[j]+ kx[k]));
             new_deltak_re_b=deltak_re_b[index2]*phase_cos-deltak_im_b[index2]*phase_sin;
             new_deltak_im_b=deltak_re_b[index2]*phase_sin+deltak_im_b[index2]*phase_cos;
             deltak_re[index2]+=new_deltak_re_b;
             deltak_im[index2]+=new_deltak_im_b;
}

        free(deltak_re_b);
        free(deltak_im_b);

        printf("Ok!\n");


}

}//end of itteration loop
free(kx);

params_input[0]=pow(L2-L1,3)/Ndata;

}

