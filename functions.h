double P_interpol(double k0, double *k, double *P, int N);
int countlines(char *filename);

void freeTokens(double** tokens, int N);

void freeTokensInt(int** tokens, int N);

void freeTokensLInt(long int** tokens, int N);

void freeTokens2(double ***tokens, int N1, int *N2);

void freeTokensInt2(int ***tokens,int N1,int *N2);


void check_box_for_yamamoto(double *parametersL1, int ngrid);

//void free3D(double ***arr, int p, int c);

//double *** Create3D(int p, int c, int r);
