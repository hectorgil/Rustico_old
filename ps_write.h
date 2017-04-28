double Leg2(double x);

double Leg4(double x);

void write_power_spectrum_skyscuts_directsum_L4(double kmin, double kmax, double kx[], double ky[], double kz[], double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[],double deltak_re4[], double deltak_im4[], double Deltak, int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise, char *binning_type);

void write_power_spectrum_skyscuts_directsum_L2L2(double kmin, double kmax, double kx[], double ky[], double kz[], double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[], double Deltak, int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise, char *binning_type);

void write_power_spectrum_skyscuts_L2L2(double kmin, double kmax,double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[], double Deltak, int ngrid, double L1, double L2, double I22, int N_interlacing, char *name_ps_out, double P_shot_noise, char *binning_type);

void write_power_spectrum_skyscuts_L4(double kmin, double kmax, double deltak_re0[], double deltak_im0[], double deltak_re2[], double deltak_im2[], double deltak_re4[], double deltak_im4[],  double Deltak, long int ngrid, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise, char *binning_type);

void write_power_spectrum_periodic(double kmin, double kmax, double deltak_re[], double deltak_im[], double bin_ps, long int  ngrid, double L1, double L2, char *name_ps_out, double P_shot_noise, char *binning_type);

void write_power_spectrum_skyscuts_directsum_exactP2(double kmin,double kmax, double *KX, double *KY, double *KZ, double *P0, double *P2, double Deltak,int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise,char *binning_type);

void write_power_spectrum_skyscuts_directsum_exactP4(double kmin,double kmax, double *KX, double *KY, double *KZ, double *P0, double *P2, double *P4, double Deltak,int ngrid, long int NGRID, double L1, double L2, double I22, char *name_ps_out, double P_shot_noise,char *binning_type);

