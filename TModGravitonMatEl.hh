
#ifndef _TMODGRAVITONMATEL_HH_
#define _TMODGRAVITONMATEL_HH_

extern "C" {
  void __modgraviton_MOD_evalamp_gg_g_vv(double P[6][4], double *MReso,  double *GaReso, double *Gggcoupl, double *Gvvcoupl, int *MYIDUP, double *MatElSq);
  void __modgraviton_MOD_evalamp_qqb_g_vv(double P[6][4], double *MReso,  double *GaReso, double *Gqqcoupl, double *Gvvcoupl, int *MYIDUP, double *MatElSq);
}

#endif

