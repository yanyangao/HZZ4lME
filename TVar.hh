#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <TLorentzVector.h>
#include "TH2F.h"
#include "TH1F.h"

#define EBEAM 4000.00
#define fbGeV2 0.389379E12
#define smallnumber 1e-15
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define   eight_2Pi_to_5 7.83410393050320417e+04
#define    four_2Pi_to_2 39.478417604357432
class TVar{
public:
  enum VerbosityLevel {
    ERROR = 0,
    INFO = 1,
    DEBUG = 2
  };
  enum MatrixElement{
    MCFM,
    MadGraph,
    JHUGen
  };
  enum Process{
    ZZ_4l    =0,
    HZZ_4l   =1,
    PSHZZ_4l =2,
    Null
  };
  //---------------------------------
  // Function
  //---------------------------------
  static TString ProcessName(int temp){ 
    if(temp==TVar::ZZ_4l   ) 
      return TString("ZZ_4l");
    else if(temp==TVar::HZZ_4l   ) 
      return TString("HZZ_4l");
    else 
      return TString("UnKnown");
  };
  ClassDef(TVar,0)
};

struct branch_particle {
  int   PdgCode   ;
  int   Charge    ;
  double Px       ;
  double Py       ;
  double Pz       ;
  double E        ;
  double Eta      ;
  double Phi      ;

};
static const TString branch_format_particle =
 "PdgCode/I:"
 "Charge/I:"
 "Px/D:"
 "Py/D:"
 "Pz/D:"
 "E/D:"
 "Eta/D:"
 "Phi/D";

// in development
struct hzz4l_event_type{
  int PdgCode[4];
  TLorentzVector p[4];
  double Xsec   [10];
  double XsecErr[10];  
};
struct mcfm_event_type{
  int PdgCode[6];
  TLorentzVector p[6];
  double pswt;
};
struct event_type{
  TLorentzVector p1,p2,ep,em,nu,nb;
  double PSWeight;
};
struct anomcoup{
	   double delg1_z, delg1_g, lambda_g, lambda_z, delk_g, delk_z_,tevscale;
};

struct EffHist{
  TH2F* els_eff_mc;
  TH2F* mus_eff_mc;
};

#endif
