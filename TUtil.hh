// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)


#ifndef ZZ_COMMON
#define ZZ_COMMON
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <TFile.h>
#include <TF1.h>
#include "TVar.hh"
#include "TMCFM.hh"



using namespace std;
TString DbnEventLepSelName(int i);
void My_choose(TVar::Process process);
bool My_smalls(double s[][12], int npart);
double SumMatrixElementPDF(TVar::Process procees, mcfm_event_type* mcfm_event,double flavor_msq[][11],double* flux);

#endif
