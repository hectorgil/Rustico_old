#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double P_interpol(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N);
if(i==0)
{
//printf("Interpol %d %lf %lf %lf %lf %lf -> \t",i,k0,k[i],k[i+1],P[i],P[i+1]);
m=( log10(P[i]) - log10(P[i+1]) )/( log10(k[i]) - log10(k[i+1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
else
{
//printf("Interpol %d %lf %lf %lf %lf %lf- >\t",i,k0,k[i-1],k[i],P[i-1],P[i]);
m=( log10(P[i]) - log10(P[i-1]) )/( log10(k[i]) - log10(k[i-1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
//printf("%lf\n",P0);
return P0;
}


int countlines(char *filename)
{
FILE* myfile = fopen(filename, "r");
int ch, number_of_lines = -1;

do 
{
	    ch = fgetc(myfile);
		    if(ch == '\n')
				        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0) 
    number_of_lines++;

	fclose(myfile);
	return number_of_lines;
}

void check_box_for_yamamoto(double *parameters, int ngrid)
{
//printf("%lf %lf %d\n",parameters[0],parameters[1],ngrid);
int i,j;
int warning;
double epsilon,deltaepsilon,L1old,L2old;
double l;
double L1,L2;
L1old=parameters[0];
L2old=parameters[1];
L1=L1old;
L2=L2old;
j=0;
deltaepsilon=(L2-L1)/ngrid*1.*0.001;//0.1 percent of gridcell
//printf("delta=%lf\n",deltaepsilon);
do
{
warning=0;
epsilon+=j*deltaepsilon;
//printf("j=%d, epsilon=%lf\n",j,epsilon);
for(i=0;i<ngrid;i++)
{
l=L1+i*(L2-L1)/ngrid*1.;
//printf("i=%d, l=%lf\n",i,l);
if(l==0){warning=-1;break;}
}
L1+=epsilon;
L2+=epsilon;
j++;
}while(warning==-1);

if(L1!=L1old){printf("Box displaced by %lf (iterated %d times) from (L1,L2)=(%lf, %lf) to (%lf, %lf) for better performance\n",epsilon,j,L1old,L2old,L1,L2);}
parameters[0]=L1;
parameters[1]=L2;

}

void freeTokens(double **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);

}




void freeTokens2(double ***tokens, int N1, int *N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2[i];++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}



void freeTokensInt(int **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);
}


void freeTokensLInt(long int **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);
}

void freeTokensInt2(int ***tokens,int N1,int *N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2[i];++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}


