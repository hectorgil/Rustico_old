#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "functions.h"

double Leg2(double x)
{
double f=0.5*(3.*x*x-1.);

return f;
}

double Leg4(double x)
{
double f=1./8.*(35.*x*x*x*x-30.*x*x+3.);

return f;
}

void write_power_spectrum_skyscuts_directsum_exactP2(double kmin,double kmax, double *kx, double *ky, double *kz, double *P0, double *P2, double Deltak,int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise,char *binning_type)
{
FILE *f;
long int i,j,k,l22;
int l,tid;
double **K;
double **Mono;
double **Quadru;
long int **nmodes;
double Pi=(4.*atan(1.));
int bintype_sw;
double kmineff;
double keff;
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}


       int nthreads;
  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)calloc(Nk,sizeof(double*));
        Quadru= (double **)calloc(Nk,sizeof(double*));
        K=(double **)calloc(Nk,sizeof(double*));
        nmodes=(long int **)calloc(Nk,sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc(nthreads,sizeof(double));
                Quadru[l] = (double*)calloc(nthreads,sizeof(double));
                nmodes[l] = (long int*)calloc(nthreads,sizeof(long int));
                K[l] = (double*)calloc(nthreads,sizeof(double));
        }
#pragma omp parallel for private(i,tid,l,kmineff,keff) shared(NGRID,K,kmin,kx,ky,kz,Mono,Quadru,P0,P2,nmodes,Deltak,Nk,bintype_sw)
        for(i=0;i<NGRID;i++)
        {
                tid=omp_get_thread_num();
if(bintype_sw==0){

l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){

l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}

                if(l<Nk)
                {
if(kz[i]==0)
{
               K[l][tid]+=pow(kx[i]*kx[i]+ky[i]*ky[i]+kz[i]*kz[i],0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=2.*P0[i];
               Quadru[l][tid]+=2.*P2[i];
               nmodes[l][tid]+=1;
}
else
{
               K[l][tid]+=2.*pow(kx[i]*kx[i]+ky[i]*ky[i]+kz[i]*kz[i],0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=2.*2.*P0[i];
               Quadru[l][tid]+=2.*2.*P2[i];
               nmodes[l][tid]+=2;
}
}
        }

free(P0);
free(P2);

for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
Quadru[l][0]+=Quadru[l][tid];
nmodes[l][0]+=nmodes[l][tid];
}
}

  f=fopen(name_ps_out,"a");
        for(l=0;l<Nk;l++)
        {
if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}

                if(nmodes[l][0]!=0 && keff>=kmin && keff<=kmax)
                {       K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=1./(nmodes[l][0]*1.*I22);
                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22);
                        

                        fprintf(f,"%lf %lf %lf %lf %i %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],nmodes[l][0], P_shot_noise);


                }
        }
        fclose(f);
freeTokens(K,Nk);
freeTokens(Mono,Nk);
freeTokens(Quadru,Nk);
freeTokensLInt(nmodes,Nk);
}

void write_power_spectrum_skyscuts_directsum_exactP4(double kmin,double kmax, double *kx, double *ky, double *kz, double *P0, double *P2, double *P4, double Deltak,int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise,char *binning_type)
{

}

void write_power_spectrum_skyscuts_directsum_L4(double kmin, double kmax, double kx[], double ky[], double kz[], double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[],double deltak_re4[], double deltak_im4[], double Deltak, int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise, char *binning_type)
{
FILE *f;
long int i,j,k,l22;
int l,tid;
double **K;
double **Mono;
double **Quadru;
double **Hexadeca;
long int **nmodes;
double Pi=(4.*atan(1.));
int bintype_sw;
double kmineff;
double keff;
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}


       
       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}


       int nthreads;
  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)calloc(Nk,sizeof(double*));
        Quadru= (double **)calloc(Nk,sizeof(double*));
        Hexadeca=(double **)calloc(Nk,sizeof(double*));
        K=(double **)calloc(Nk,sizeof(double*));
        nmodes=(long int **)calloc(Nk,sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc(nthreads,sizeof(double));
                Quadru[l] = (double*)calloc(nthreads,sizeof(double));
                Hexadeca[l] = (double*)calloc(nthreads,sizeof(double));
                nmodes[l] = (long int*)calloc(nthreads,sizeof(long int));
                K[l] = (double*)calloc(nthreads,sizeof(double));
        }



#pragma omp parallel for private(i,tid,l,kmineff,keff) shared(NGRID,K,kmin,kx,ky,kz,Mono,Quadru,Hexadeca,deltak_re0,deltak_re2,deltak_re4,deltak_im0,deltak_im2,deltak_im4,nmodes,Deltak,Nk,bintype_sw)
        for(i=0;i<NGRID;i++)
        {
                tid=omp_get_thread_num();

if(bintype_sw==0){

l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){

l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}

                if(l<Nk)
                {
if(kz[i]==0)
{
               K[l][tid]+=pow(kx[i]*kx[i]+ky[i]*ky[i]+kz[i]*kz[i],0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=(pow(deltak_re0[i],2)+pow(deltak_im0[i],2));
               Quadru[l][tid]+=deltak_re0[i]*0.5*(3.*deltak_re2[i]-deltak_re0[i])+deltak_im0[i]*0.5*(3.*deltak_im2[i]-deltak_im0[i]);
               Hexadeca[l][tid]+=deltak_re0[i]*(35.*deltak_re4[i]-30.*deltak_re2[i]+3.*deltak_re0[i])+deltak_im2[i]*(35.*deltak_im4[i]-30.*deltak_im2[i]+3.*deltak_im0[i]);
               nmodes[l][tid]+=1;
}
else
{
                K[l][tid]+=2.*pow(kx[i]*kx[i]+ky[i]*ky[i]+kz[i]*kz[i],0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=2.*(pow(deltak_re0[i],2)+pow(deltak_im0[i],2));
               Quadru[l][tid]+=2.*deltak_re0[i]*0.5*(3.*deltak_re2[i]-deltak_re0[i])+2.*deltak_im0[i]*0.5*(3.*deltak_im2[i]-deltak_im0[i]);
               Hexadeca[l][tid]+=2.*deltak_re0[i]*(35.*deltak_re4[i]-30.*deltak_re2[i]+3.*deltak_re0[i])+2.*deltak_im2[i]*(35.*deltak_im4[i]-30.*deltak_im2[i]+3.*deltak_im0[i]);
               nmodes[l][tid]+=2;
}


}
        }

free(deltak_re0);
free(deltak_re2);
free(deltak_re4);
free(deltak_im0);
free(deltak_im2);
free(deltak_im4);


for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
nmodes[l][0]+=nmodes[l][tid];
}
}


  f=fopen(name_ps_out,"a");

        for(l=0;l<Nk;l++)
        {
if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}

                if(nmodes[l][0]!=0 && keff>=kmin && keff<=kmax)
                {       K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=1./(nmodes[l][0]*1.*I22);
                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22);
                        Hexadeca[l][0]*=1./(nmodes[l][0]*1.*I22);

                        fprintf(f,"%lf %lf %lf %lf %lf %i %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);


                }
        }
        fclose(f);


freeTokens(K,Nk);
freeTokens(Mono,Nk);
freeTokens(Quadru,Nk);
freeTokens(Hexadeca,Nk);
freeTokensLInt(nmodes,Nk);

}

void write_power_spectrum_skyscuts_directsum_L2L2(double kmin, double kmax, double kx[], double ky[], double kz[], double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[], double Deltak, int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise, char *binning_type)
{

FILE *f;
long int i,j,k,l22;
int l,tid;
double **K;
double **Mono;
double **Quadru;
double **Hexadeca;
long int **nmodes;
double Pi=(4.*(atan(1.)));
int bintype_sw;
double kmineff;
double keff;
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}

      
       int nthreads;
  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)calloc(Nk,sizeof(double*));
        Quadru= (double **)calloc(Nk,sizeof(double*));
        Hexadeca=(double **)calloc(Nk,sizeof(double*));
        K=(double **)calloc(Nk,sizeof(double*));
        nmodes=(long int **)calloc(Nk,sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc(nthreads,sizeof(double));
                Quadru[l] = (double*)calloc(nthreads,sizeof(double));
                Hexadeca[l] = (double*)calloc(nthreads,sizeof(double));
                nmodes[l] = (long int*)calloc(nthreads,sizeof(long int));
                K[l] = (double*)calloc(nthreads,sizeof(double));
        }



#pragma omp parallel for private(i,tid,l,keff,kmineff) shared(NGRID,K,kmin,kx,ky,kz,Mono,Quadru,Hexadeca,deltak_re0,deltak_re2,deltak_im0,deltak_im2,nmodes,Deltak,Nk,bintype_sw)
        for(i=0;i<NGRID;i++)
        {
                tid=omp_get_thread_num();

if(bintype_sw==0){

l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){

l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}


                if(l<Nk)
                {
if(kz[i]==0)
{
               K[l][tid]+=pow(kx[i]*kx[i]+ky[i]*ky[i]+kz[i]*kz[i],0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=(pow(deltak_re0[i],2)+pow(deltak_im0[i],2));
               Quadru[l][tid]+=deltak_re0[i]*0.5*(3.*deltak_re2[i]-deltak_re0[i])+deltak_im0[i]*0.5*(3.*deltak_im2[i]-deltak_im0[i]);
               Hexadeca[l][tid]+=(315./8.*(pow(deltak_re2[i],2)+pow(deltak_im2[i],2))-135./4.*(deltak_re0[i]*deltak_re2[i]+deltak_im0[i]*deltak_im2[i])+27./8.*(pow(deltak_re0[i],2)+pow(deltak_im0[i],2)));
               nmodes[l][tid]+=1;

}
else
{
               K[l][tid]+=2.*pow(kx[i]*kx[i]+ky[i]*ky[i]+kz[i]*kz[i],0.5);//two modes coming from the condition kz>0
               Mono[l][tid]+=2.*(pow(deltak_re0[i],2)+pow(deltak_im0[i],2));
               Quadru[l][tid]+=2.*deltak_re0[i]*0.5*(3.*deltak_re2[i]-deltak_re0[i])+2.*deltak_im0[i]*0.5*(3.*deltak_im2[i]-deltak_im0[i]);
               Hexadeca[l][tid]+=2.*(315./8.*(pow(deltak_re2[i],2)+pow(deltak_im2[i],2))-135./4.*(deltak_re0[i]*deltak_re2[i]+deltak_im0[i]*deltak_im2[i])+27./8.*(pow(deltak_re0[i],2)+pow(deltak_im0[i],2)));
               nmodes[l][tid]+=2;
}

}

        }


free(deltak_re0);
free(deltak_re2);
free(deltak_im0);
free(deltak_im2);


for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
nmodes[l][0]+=nmodes[l][tid];
}
}


  f=fopen(name_ps_out,"a");
        for(l=0;l<Nk;l++)
        {
if(bintype_sw==0){keff=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){keff=pow(10,(l+0.5)*Deltak+log10(kmin));}
                if(nmodes[l][0]!=0 && keff>=kmin && keff<=kmax)
                {       K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=1./(nmodes[l][0]*1.*I22);
                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22);
                        Hexadeca[l][0]*=1.0/(nmodes[l][0]*1.*I22);
                        fprintf(f,"%lf %lf %lf %lf %lf %i %lf\n",keff, K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);


                }
        }
        fclose(f);


freeTokens(K,Nk);
freeTokens(Mono,Nk);
freeTokens(Quadru,Nk);
freeTokens(Hexadeca,Nk);
freeTokensLInt(nmodes,Nk);

}

void write_power_spectrum_skyscuts_L2L2(double kmin, double kmax, double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[], double Deltak, int ngrid, double L1, double L2, double I22, int N_interlacing, char *name_ps_out, double P_shot_noise, char *binning_type)
{//printf("hola ke ase\n");
double Pi=(4.*atan(1.));
double **K;
double **k_av;
double **Mono;
double **Quadru;
double **Hexadeca;
long int **nmodes;
int l,tid,i,j,k;
int i2,j2,k2;
long int l2;
FILE *f;
double *kx;
int bintype_sw;
double kmineff;
double keff;
long int index2;
long int ngridtot=pow(ngrid,3);
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

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


       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)+1;}



       int nthreads;

  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)calloc(Nk,sizeof(double*));
        Quadru= (double **)calloc(Nk,sizeof(double*));
        Hexadeca=(double **)calloc(Nk,sizeof(double*));
        K=(double **)calloc(Nk,sizeof(double*));
        k_av=(double **)calloc(Nk,sizeof(double*));
        nmodes=(long int **)calloc(Nk,sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc(nthreads,sizeof(double));
                Quadru[l] = (double*)calloc(nthreads,sizeof(double));
                Hexadeca[l] = (double*)calloc(nthreads,sizeof(double));
                nmodes[l] = (long int*)calloc(nthreads,sizeof(long int));
                K[l] = (double*)calloc(nthreads,sizeof(double));
                k_av[l] = (double*)calloc(nthreads,sizeof(double));
        }

#pragma  omp parallel for private(index2,l2,i,j,k,l,tid,i2,k2,j2,keff,kmineff) shared(ngrid,ngridtot,kx,Deltak,Nk,K,deltak_re0,deltak_im0,deltak_re2,deltak_im2,nmodes,Mono,Quadru,Hexadeca,N_interlacing,kmin,bintype_sw)
        for(l2=0;l2<ngridtot;l2++)
        {
                tid=omp_get_thread_num();
                i=(int)(l2/(ngrid*ngrid*1.));
                j=(int)( (l2-i*ngrid*ngrid)/(ngrid*1.));
                k=l2-i*ngrid*ngrid-j*ngrid;

if(bintype_sw==0){

l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){

l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}

                if(l<Nk)
                {
                K[l][tid]+=pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5);//averaged value

if(bintype_sw==0){k_av[l][tid]+=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){k_av[l][tid]+=pow(10,(l+0.5)*Deltak+log10(kmin));}

if(kx[k]>=0)
{
index2=((pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j);


if(index2>0){
          Mono[l][tid]+=pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2);
Quadru[l][tid]+=deltak_re0[index2]/N_interlacing*1.*0.5*(3.*deltak_re2[index2]/N_interlacing*1.-deltak_re0[index2]/N_interlacing*1.)+deltak_im0[index2]/N_interlacing*1.*0.5*(3.*deltak_im2[index2]/N_interlacing*1.-deltak_im0[index2]/N_interlacing*1.);
Hexadeca[l][tid]+=315./8.*(pow(deltak_re2[index2]/N_interlacing*1.,2)+pow(deltak_im2[index2]/N_interlacing*1.,2))-135./4.*(deltak_re0[index2]/N_interlacing*1.*deltak_re2[index2]/N_interlacing*1.+deltak_im0[index2]/N_interlacing*1.*deltak_im2[index2]/N_interlacing*1.)+27./8.*(pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2));
nmodes[l][tid]+=1;
}

}
else
{
k2=ngrid-k;
i2=ngrid-i;
j2=ngrid-j;
if(i2==ngrid){i2=0;}
if(j2==ngrid){j2=0;}
if(k2==ngrid){k2=0;}
index2=((pow(ngrid,2)*i2+ngrid*j2+2*k2)/2+i2*ngrid+j2);


if(index2>0){

Mono[l][tid]+=pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2);
Quadru[l][tid]+=deltak_re0[index2]/N_interlacing*1.*0.5*(3.*deltak_re2[index2]/N_interlacing*1.-deltak_re0[index2]/N_interlacing*1.)+deltak_im0[index2]/N_interlacing*1.*0.5*(3.*deltak_im2[index2]/N_interlacing*1.-deltak_im0[index2]/N_interlacing*1.);
Hexadeca[l][tid]+=315./8.*(pow(deltak_re2[index2]/N_interlacing*1.,2)+pow(deltak_im2[index2]/N_interlacing*1.,2))-135./4.*(deltak_re0[index2]/N_interlacing*1.*deltak_re2[index2]/N_interlacing*1.+deltak_im0[index2]/N_interlacing*1.*deltak_im2[index2]/N_interlacing*1.)+27./8.*(pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2));
                nmodes[l][tid]+=1;
}

}
               }
      }

for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
k_av[l][0]+=k_av[l][tid];
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
nmodes[l][0]+=nmodes[l][tid];
}
}


  f=fopen(name_ps_out,"a");
 
        for(l=0;l<Nk;l++)
        {
                if(nmodes[l][0]!=0)
                {       k_av[l][0]*=1.0/nmodes[l][0]*1.;
                        K[l][0]*=1.0/nmodes[l][0]*1.;
                        Mono[l][0]*=1./(nmodes[l][0]*1.*I22);
                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22);
                        Hexadeca[l][0]*=1.0/(nmodes[l][0]*1.*I22);
if(k_av[l][0]>=kmin && k_av[l][0]<=kmax ){

                        fprintf(f,"%lf %lf %lf %lf %lf %i %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);
}
                        
                }
        }
        fclose(f);



freeTokens(K,Nk);
freeTokens(k_av,Nk);
freeTokens(Mono,Nk);
freeTokens(Quadru,Nk);
freeTokens(Hexadeca,Nk);
freeTokensLInt(nmodes,Nk);
free(kx);

free(deltak_re0);
free(deltak_im0);
free(deltak_re2);
free(deltak_im2);


}


void write_power_spectrum_skyscuts_L4(double kmin, double kmax, double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[], double deltak_re4[], double deltak_im4[],  double Deltak, int ngrid, double L1, double L2, double I22, int N_interlacing, char *name_ps_out, double P_shot_noise, char *binning_type)
{
double Pi=(4.*atan(1.));
double **K;
double **k_av;
double **Mono;
double **Quadru;
double **Hexadeca;
long int **nmodes;
int l,tid,i,j,k;
int i2,j2,k2;
long int l2;
long int index2;
long int ngridtot=pow(ngrid,3);
FILE *f;
double *kx;
int bintype_sw;
double kmineff;
double keff;
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

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


       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}

       
       int nthreads;

  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)calloc(Nk,sizeof(double*));
        Quadru= (double **)calloc(Nk,sizeof(double*));
        Hexadeca=(double **)calloc(Nk,sizeof(double*));
        K=(double **)calloc(Nk,sizeof(double*));
        k_av=(double **)calloc(Nk,sizeof(double*));
        nmodes=(long int **)calloc(Nk,sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)calloc(nthreads,sizeof(double));
                Quadru[l] = (double*)calloc(nthreads,sizeof(double));
                Hexadeca[l] = (double*)calloc(nthreads,sizeof(double));
                nmodes[l] = (long int*)calloc(nthreads,sizeof(long int));
                K[l] = (double*)calloc(nthreads,sizeof(double));
                k_av[l] = (double*)calloc(nthreads,sizeof(double));
        }

#pragma  omp parallel for private(index2,l2,i,j,k,l,tid,i2,k2,j2,keff,kmineff) shared(ngrid,ngridtot,kx,Deltak,Nk,K,deltak_re0,deltak_im0,deltak_re2,deltak_im2,nmodes,Mono,Quadru,Hexadeca,N_interlacing,kmin,bintype_sw)
        for(l2=0;l2<ngridtot;l2++)
        {
tid=omp_get_thread_num();
                i=(int)(l2/(ngrid*ngrid*1.));
                j=(int)( (l2-i*ngrid*ngrid)/(ngrid*1.));
                k=l2-i*ngrid*ngrid-j*ngrid;

if(bintype_sw==0){

l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){

l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}
                if(l<Nk)
                {
                K[l][tid]+=pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5);//averaged value

if(bintype_sw==0){k_av[l][tid]+=(l+0.5)*Deltak+kmin;}
if(bintype_sw==1){k_av[l][tid]+=pow(10,(l+0.5)*Deltak+log10(kmin));}

if(kx[k]>=0)
{
index2=((pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j);
if(index2>0)
{
Mono[l][tid]+=pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2);
Quadru[l][tid]+=deltak_re0[index2]/N_interlacing*1.*0.5*(3.*deltak_re2[index2]/N_interlacing*1.-deltak_re0[index2]/N_interlacing*1.)+deltak_im0[index2]/N_interlacing*1.*0.5*(3.*deltak_im2[index2]/N_interlacing*1.-deltak_im0[index2]/N_interlacing*1.);
Hexadeca[l][tid]+=315./8.*(deltak_re0[index2]/N_interlacing*1.*deltak_re4[index2]/N_interlacing*1.+deltak_im0[index2]/N_interlacing*1.*deltak_im4[index2]/N_interlacing*1.)-135./4.*(deltak_re0[index2]/N_interlacing*1.*deltak_re2[index2]/N_interlacing*1.+deltak_im0[index2]/N_interlacing*1.*deltak_im2[index2]/N_interlacing*1.)+27./8.*(pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2));
                nmodes[l][tid]+=1;
}
}
else
{
k2=ngrid-k;
i2=ngrid-i;
j2=ngrid-j;
if(i2==ngrid){i2=0;}
if(j2==ngrid){j2=0;}
if(k2==ngrid){k2=0;}
index2=((pow(ngrid,2)*i2+ngrid*j2+2*k2)/2+i2*ngrid+j2);
if(index2>0)
{
Mono[l][tid]+=pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2);
Quadru[l][tid]+=deltak_re0[index2]/N_interlacing*1.*0.5*(3.*deltak_re2[index2]/N_interlacing*1.-deltak_re0[index2]/N_interlacing*1.)+deltak_im0[index2]/N_interlacing*1.*0.5*(3.*deltak_im2[index2]/N_interlacing*1.-deltak_im0[index2]/N_interlacing*1.);
Hexadeca[l][tid]+=315./8.*(deltak_re0[index2]/N_interlacing*1.*deltak_re4[index2]/N_interlacing*1.+deltak_im0[index2]/N_interlacing*1.*deltak_im4[index2]/N_interlacing*1.)-135./4.*(deltak_re0[index2]/N_interlacing*1.*deltak_re2[index2]/N_interlacing*1.+deltak_im0[index2]/N_interlacing*1.*deltak_im2[index2]/N_interlacing*1.)+27./8.*(pow(deltak_re0[index2]/N_interlacing*1.,2)+pow(deltak_im0[index2]/N_interlacing*1.,2));
}
                nmodes[l][tid]+=1;

}
               }
      }

for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
k_av[l][0]+=k_av[l][tid];
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
nmodes[l][0]+=nmodes[l][tid];
}
}


  f=fopen(name_ps_out,"a");
 
        for(l=0;l<Nk;l++)
        {
                if(nmodes[l][0]!=0)
                {       k_av[l][0]*=1./nmodes[l][0]*1.;
                        K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=1./(nmodes[l][0]*1.*I22);
                        Quadru[l][0]*=5./(nmodes[l][0]*1.*I22);
                        Hexadeca[l][0]*=1./(nmodes[l][0]*I22);
if(k_av[l][0]>=kmin && k_av[l][0]<=kmax){
                        fprintf(f,"%lf %lf %lf %lf %lf %i %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);
}
                        
                }
        }
        fclose(f);



freeTokens(K,Nk);
freeTokens(k_av,Nk);
freeTokens(Mono,Nk);
freeTokens(Quadru,Nk);
freeTokens(Hexadeca,Nk);
freeTokensLInt(nmodes,Nk);
free(deltak_re0);
free(deltak_im0);
free(deltak_re2);
free(deltak_im2);
free(deltak_re4);
free(deltak_im4);
free(kx);
}



void write_power_spectrum_periodic(double kmin, double kmax, double deltak_re[], double deltak_im[], double Deltak, int  ngrid, double L1, double L2, int Ninterlacing, char *name_ps_out, double P_shot_noise, char *binning_type)
{
double Pi=(4.*atan(1.));
double **K;
double **k_av;
double **Mono;
double **Quadru;
double **Hexadeca;
long int **nmodes;
int l,tid,i,j,k;
int i2,j2,k2;
long int l2;
FILE *f;
double *kx;
double argument;
int bintype_sw;
double kmineff;
double keff;
long int index2;
long int ngridtot=pow(ngrid,3);
if(strcmp(binning_type, "linear") == 0){bintype_sw=0;}
if(strcmp(binning_type, "log10") == 0){bintype_sw=1;}

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


       int Nk;
if(bintype_sw==0){Nk=(int)(sqrt(3.)*2.*Pi*ngrid/(Deltak*2.*(L2-L1)))+1;}
if(bintype_sw==1){Nk=(int)((log10(sqrt(3.)*2.*Pi*ngrid/(2.*(L2-L1)))-log10(2.*Pi/(L2-L1)))/Deltak)  +1;}

       int nthreads;

  #pragma omp parallel for private(i,tid) shared(nthreads,ngrid)
        for(i=0;i<ngrid;i++)
        {
                tid=omp_get_thread_num();
                if(tid==0 && i==0){nthreads=omp_get_num_threads();}
        }

        Mono= (double **)malloc((Nk)*sizeof(double*));
        Quadru= (double **)malloc((Nk)*sizeof(double*));
        Hexadeca=(double **)malloc((Nk)*sizeof(double*));
        K=(double **)malloc((Nk)*sizeof(double*));
        k_av=(double **)malloc((Nk)*sizeof(double*));
        nmodes=(long int **)malloc((Nk)*sizeof(long int*));

        for(l=0;l<Nk;l++)
        {
                Mono[l] = (double*)malloc((nthreads)*sizeof(double));
                Quadru[l] = (double*)malloc((nthreads)*sizeof(double));
                Hexadeca[l] = (double*)malloc((nthreads)*sizeof(double));
                nmodes[l] = (long int*)malloc((nthreads)*sizeof(long int));
                K[l] = (double*)malloc((nthreads)*sizeof(double));
                k_av[l] = (double*)malloc((nthreads)*sizeof(double));
        }
        for(l=0;l<Nk;l++)
        {
           for(tid=0;tid<nthreads;tid++)
           {
              k_av[l][tid]=0;
              Mono[l][tid]=0;
              Quadru[l][tid]=0;
              Hexadeca[l][tid]=0;
              K[l][tid]=0;
              nmodes[l][tid]=0;
           }
        }
#pragma  omp parallel for private(index2,l2,i,j,k,l,tid,i2,k2,j2,argument,kmineff,keff) shared(ngrid,ngridtot,kx,Deltak,Nk,K,deltak_re,deltak_im,nmodes,Mono,Quadru,Hexadeca,Ninterlacing,kmin,bintype_sw)
        for(l2=0;l2<ngridtot;l2++)
        {
                tid=omp_get_thread_num();
                i=(int)(l2/(ngrid*ngrid*1.));
                j=(int)( (l2-i*ngrid*ngrid)/(ngrid*1.));
                k=l2-i*ngrid*ngrid-j*ngrid;

if(bintype_sw==0){
l=(int)((pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5)-kmin)/Deltak);//-1.;
if(l<0){l=0;}
}

if(bintype_sw==1){
l=(int)(( log10(pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5))-log10(kmin))/Deltak);//-1.;
if(l<0){l=0;}
}

                if(l<Nk)
                {
                K[l][tid]+=pow(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k],0.5);
if(bintype_sw==0){k_av[l][tid]+=(l+0.5)*Deltak+(kmin);}
if(bintype_sw==1){k_av[l][tid]+=pow(10,(l+0.5)*Deltak+log10(kmin));}

argument=sqrt(kx[k]*kx[k]/(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]));

if(kx[i]*kx[i]+kx[j]*kx[j]+kx[k]*kx[k]==0){argument=0;}

if(kx[k]>=0)
{
index2=((pow(ngrid,2)*i+ngrid*j+2*k)/2+i*ngrid+j);
if(index2>0)
{
Mono[l][tid]+=pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2);
Quadru[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg2(argument);
Hexadeca[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg4(argument);
nmodes[l][tid]+=1;
}
}
else
{
k2=ngrid-k;
i2=ngrid-i;
j2=ngrid-j;
if(i2==ngrid){i2=0;}
if(j2==ngrid){j2=0;}
if(k2==ngrid){k2=0;}
index2=((pow(ngrid,2)*i2+ngrid*j2+2*k2)/2+i2*ngrid+j2);
if(index2>0)
{
Mono[l][tid]+=pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2);
Quadru[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg2(argument);
Hexadeca[l][tid]+=(pow(deltak_re[index2]/Ninterlacing*1.,2)+pow(deltak_im[index2]/Ninterlacing*1.,2))*Leg4(argument);
               nmodes[l][tid]+=1;
}
}
               }
      }
for(l=0;l<Nk;l++)
{
for(tid=1;tid<nthreads;tid++)
{
k_av[l][0]+=k_av[l][tid];
K[l][0]+=K[l][tid];
Mono[l][0]+=Mono[l][tid];
Quadru[l][0]+=Quadru[l][tid];
Hexadeca[l][0]+=Hexadeca[l][tid];
nmodes[l][0]+=nmodes[l][tid];
}
}


  f=fopen(name_ps_out,"a");

        for(l=0;l<Nk;l++)
        {
                if(nmodes[l][0]!=0)
                {       k_av[l][0]*=1./nmodes[l][0]*1.;
                        K[l][0]*=1./nmodes[l][0]*1.;
                        Mono[l][0]*=pow(L2-L1,3)/(nmodes[l][0]*1.);
                        Quadru[l][0]*=pow(L2-L1,3)*5./(nmodes[l][0]*1.);
                        Hexadeca[l][0]*=pow(L2-L1,3)*9/(nmodes[l][0]);
if(k_av[l][0]<=kmax && k_av[l][0]>=kmin){
                        fprintf(f,"%lf %lf %lf %lf %lf %i %lf\n",k_av[l][0],K[l][0],Mono[l][0]-P_shot_noise,Quadru[l][0],Hexadeca[l][0],nmodes[l][0], P_shot_noise);}


                }
        }
        fclose(f);



freeTokens(K,Nk);
freeTokens(k_av,Nk);
freeTokens(Mono,Nk);
freeTokens(Quadru,Nk);
freeTokens(Hexadeca,Nk);
freeTokensLInt(nmodes,Nk);
free(deltak_re);
free(deltak_im);
free(kx);
}
