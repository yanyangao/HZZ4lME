#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
// ME related
#include "TVar.hh"
#include "TEvtProb.hh"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double ERRORthreshold=1.0;
using namespace std;

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity);

//###################
//# main function
//###################
void runME_HZZ4l(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity=TVar::INFO){

  if (verbosity >= TVar::INFO) cout <<"=== Calculating differential cross-section ==========" <<endl;  
  xseccalc(inputDir, fileName, outputDir, maxevt, verbosity); 
 
}

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity){

  if (verbosity >= TVar::INFO) cout << "Input File: " << fileName << " \n";
  
  TFile* fin = new TFile(inputDir+fileName);
  TString outFileName = outputDir+fileName;
  outFileName.ReplaceAll(".root","_ME.root");
  cout << outFileName <<endl;
  TFile *newfile = new TFile(outFileName,"recreate");

  TTree* ch=(TTree*)fin->Get("passedEvents"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // Declare the matrix element related variables to be added to the existing ntuples
  double dXsecList;
  double dXsecErrList;
  
  evt_tree->Branch("dXsec"      ,&dXsecList           ,"dXsec/D");
  evt_tree->Branch("dXsecErr"   ,&dXsecErrList        ,"dXsecErr/D");
  
  // Initialize the branches to use to calculate the differential cross-sections
  Double_t EL1_ = 0.;
  Double_t pXL1_ = 0.;
  Double_t pYL1_ = 0.;
  Double_t pZL1_ = 0.;

  Double_t EL2_ = 0.;
  Double_t pXL2_ = 0.;
  Double_t pYL2_ = 0.;
  Double_t pZL2_ = 0.;
  
  Double_t EL3_ = 0.;
  Double_t pXL3_ = 0.;
  Double_t pYL3_ = 0.;
  Double_t pZL3_ = 0.;

  Double_t EL4_ = 0.;
  Double_t pXL4_ = 0.;
  Double_t pYL4_ = 0.;
  Double_t pZL4_ = 0.;
  
  ch->SetBranchAddress( "EL1"       , &EL1_      );   
  ch->SetBranchAddress( "pXL1"      , &pXL1_     );   
  ch->SetBranchAddress( "pYL1"      , &pYL1_     );   
  ch->SetBranchAddress( "pZL1"      , &pZL1_     );   

  ch->SetBranchAddress( "EL2"       , &EL2_      );   
  ch->SetBranchAddress( "pXL2"      , &pXL2_     );   
  ch->SetBranchAddress( "pYL2"      , &pYL2_     );   
  ch->SetBranchAddress( "pZL2"      , &pZL2_      );   

  ch->SetBranchAddress( "EL3"       , &EL3_      );   
  ch->SetBranchAddress( "pXL3"      , &pXL3_     );   
  ch->SetBranchAddress( "pYL3"      , &pYL3_     );   
  ch->SetBranchAddress( "pZL3"      , &pZL3_      );   

  ch->SetBranchAddress( "EL4"       , &EL4_      );   
  ch->SetBranchAddress( "pXL4"      , &pXL4_     );   
  ch->SetBranchAddress( "pYL4"      , &pYL4_     );   
  ch->SetBranchAddress( "pZL4"      , &pZL4_      );   
  
  

  // Create the instance of TEvtProb to calculate the differential cross-section
  TEvtProb Xcal2;  
  Xcal2.SetMatrixElement(TVar::MCFM);
  hzz4l_event_type hzz4l_event;
  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = maxevt > ch->GetEntries() ? ch->GetEntries() : maxevt; 
  
  if (verbosity >= TVar::INFO) printf("Total number of events = %d\n", Ntot);
  
  for(int ievt = 0; ievt < Ntot; ievt++){
    if (verbosity >= TVar::INFO && (ievt % 5 == 0)) 
      std::cout << "Doing Event: " << ievt << std::endl;
    
    dXsecList = 0.;
    dXsecErrList = 0.;
    
    ch->GetEntry(ievt);           
    
    hzz4l_event.p[0].SetXYZM(pXL1_, pYL1_, pZL1_, 0.);
    hzz4l_event.p[1].SetXYZM(pXL2_, pYL2_, pZL2_, 0.);
    hzz4l_event.p[2].SetXYZM(pXL3_, pYL3_, pZL3_, 0.);
    hzz4l_event.p[3].SetXYZM(pXL4_, pYL4_, pZL4_, 0.);
    
    double z1mass = (hzz4l_event.p[0]+hzz4l_event.p[1]).M();
    double z2mass = (hzz4l_event.p[1]+hzz4l_event.p[2]).M();

    // if ( TMath::Abs(z1mass - 91.1876) > 15. || TMath::Abs(z2mass - 91.1875) > 15. ) continue;
    
    if (verbosity >= TVar::DEBUG) {
      cout << "\n=========================================================\n";
      cout << "Entry: " << ievt << "\n";
      cout << "Input: ==================================================" <<endl;
      printf("lep1 p3 = (%4.4f, %4.4f, %4.4f)  lep2 p3 = (%4.4f, %4.4f, %4.4f)\n",
	     pXL1_, pYL1_, pZL1_, pXL2_, pYL2_, pZL2_); 
      printf("lep3 p3 = (%4.4f, %4.4f, %4.4f)  lep4 p3 = (%4.4f, %4.4f, %4.4f)\n",
	     pXL3_, pYL3_, pZL3_, pXL4_, pYL4_, pZL4_); 
      std::cout << "ZZ system (pX, pY, pZ, pZ) = ( " 
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Px() << ", "
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Py() << ", "
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Pz() << ", "
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Energy() << ")\n";
      std::cout << "Z1 mass = " << z1mass << "\tz2mass = " << z2mass << "\n";
      cout << "=========================================================\n";
    } 
    // finish loading event information
    
    // ==== Begin the differential cross-section calculation
    TVar::Process Proc = TVar::ZZ_4l;
    double Xsec = Xcal2.XsecCalc(Proc,hzz4l_event,verbosity);
    
    // fill in the differential cross-section and errors
    dXsecList = Xsec;  
    
    evt_tree->Fill();
    
  }//nevent
  
  if (verbosity >= TVar::INFO) cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
}  
