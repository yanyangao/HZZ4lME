
#ifndef _TMODZPRIMEMATEL_HH_
#define _TMODZPRIMEMATEL_HH_

extern "C" {
  void __modzprime_MOD_evalamp_qqb_zprime_vv(double P[6][4], double *MReso,  double *GaReso, double *Zqqcoupl, double *Zvvcoupl, int *MYIDUP, double *MatElSq);
  void __modzprime_MOD_evalamp_zprime_vv(double P[6][4], double *MReso,  double *GaReso, double *Zvvcoupl, int *MYIDUP, double *MatElSq);
}

#endif

