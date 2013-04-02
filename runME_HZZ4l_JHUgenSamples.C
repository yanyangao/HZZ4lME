// 
// This code is for private use only Yanyan Gao (ygao@fnal.gov)
// 

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TProfile.h"
#include <iostream>
#include "Math/LorentzVector.h"
#include "TLorentzRotation.h"
#include "Math/VectorUtil.h"
// ME related
#include "TVar.hh"
#include "TEvtProb.hh"

// analytic PDFs
#include "PDFs/RooXZsZs_5D.h"
#include "PDFs/RooSpinTwo_7D.h"
#include "PDFs/RooSpinOne_7D.h"
#include "PDFs/RooSpinOne_Decay.h"
#include "PDFs/RooSpinTwo_Decay.h"
#include "RooRealVar.h"
#include "math.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double ERRORthreshold=1.0;
using namespace std;

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity);

void checkZorder(float& z1mass, float& z2mass,
		 float& costhetastar, float& costheta1,
		 float& costheta2, float& phi, 
		 float& phistar1){
  
  float tempZ1mass=z1mass;
  float tempZ2mass=z2mass;
  float tempH1=costheta1;
  float tempH2=costheta2;
  float tempHs=costhetastar;
  float tempPhi1=phistar1;
  float tempPhi=phi;

  if(z2mass>z1mass){
    //cout<<"inverted"<<endl;
    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*TMath::Pi();
    if(phistar1<-3.1415)
      phistar1=phistar1+2*TMath::Pi();

  }else
    return;

}

void CalculateAnalyticalMELA(float mzz_,float m1_, float m2_,
			     float hs_, float h1_, float h2_, 
			     float phi_, float phi1_,
			     double &Psmh, 
			     double &Poneminus, double &Poneplus,
			     double &Poneminus_decay, double &Poneplus_decay,
			     double &Ptwomplus_gg, double &Ptwomplus_qq, double &Ptwomplus_decay,
			     double &Ptwohminus, double &Ptwohplus, double &Ptwobplus
			     ) {
  
  
    // - - - - - - - - - - - - - - -
    // calculating the analytical MELA
    // - - - - - - - - - - - - - - - 
    
    RooRealVar mzz("mzz","mzz",100.,1000.);
    RooRealVar m1("m1","m1",0.,1000.);
    RooRealVar m2("m2","m2",0.,1000.);
    
    RooRealVar hs("hs","hs",-1.,1.);
    RooRealVar h1("h1","h1",-1.,1.);
    RooRealVar h2("h2","h2",-1.,1.);
    
    RooRealVar phi("phi","phi",-3.2,3.2);
    RooRealVar phi1("phi1","phi1",-3.2,3.2);
    

    checkZorder(m1_,m2_,hs_,h1_,h2_,phi_,phi1_);
    mzz.setVal(mzz_);  m1.setVal(m1_);   m2.setVal(m2_);
    hs.setVal(hs_);  h1.setVal(h1_);   h2.setVal(h2_);  
    phi.setVal(phi_);  phi1.setVal(phi1_); 
    

    // - - - - - - - - - - - - - - - -
    // SM Higgs parameters
    // - - - - - - - - - - - - - - - - 
    
    RooRealVar a1("a1","a1",1.,-1000.,1000.);
    RooRealVar a2("a2","a2",0.,-1000.,1000.);
    RooRealVar a3("a3","a3",0.,-1000.,1000.);
    
    RooRealVar p1("p1","p1",0.,-1000.,1000.);
    RooRealVar p2("p2","p2",0.,-1000.,1000.);
    RooRealVar p3("p3","p3",0.,-1000.,1000.);
    
    // - - - - - - - - - - - - - - - -
    // common parameters
    // - - - - - - - - - - - - - - - - 
    
    RooRealVar mZ("mZ","mZ",91.188,0.,1000.);
    RooRealVar gamZ("gamZ","gamZ",2.5,0.,1000.);
    
    RooRealVar R1("R1","R1",.15,0.,1000.);  
    RooRealVar R2("R2","R2",.15,0.,1000.);
    
    // - - - - - - - - - - - - - - - -
    // SM Higgs PDF
    // - - - - - - - - - - - - - - - - 
    
    RooXZsZs_5D SMHiggs("SMHiggs","SMHiggs",m1,m2,h1,h2,phi,a1,p1,a2,p2,a3,p3,mZ,gamZ,mzz,R1,R2);
    Psmh = SMHiggs.getVal();
    
    //
    // Spin One parameters
    // 
    
    // 1-
    RooRealVar g1ValV("g1ValV","g1ValV",1.);
    RooRealVar g2ValV("g2ValV","g2ValV",0.);

    // this is the acceptance term associated with the production angles
    // the default setting is for setting no-acceptance
    RooRealVar aParam("aParam","aParam",0);

    RooSpinOne_7D oneminus("oneminus","oneminus", mzz,m1,m2,h1,h2,hs,phi,phi1,
			   g1ValV, g2ValV, 
			   R1, R2, aParam, mZ, gamZ);

    Poneminus = oneminus.getVal();
    
    RooSpinOne_Decay oneminus_decay("oneminus_decay","oneminus_decay", mzz,m1,m2,h1,h2,phi,
				  g1ValV, g2ValV, 
				  R1, R2, aParam, mZ, gamZ);
    Poneminus_decay = oneminus_decay.getVal();
    
    // 1+ 
    RooRealVar g1ValA("g1ValA","g1ValA",0.);
    RooRealVar g2ValA("g2ValA","g2ValA",1.);


    RooSpinOne_7D oneplus("oneplus","oneplus", mzz,m1,m2,h1,h2,hs,phi,phi1,
			   g1ValA, g2ValA, 
			   R1, R2, aParam, mZ, gamZ);
    
    Poneplus = oneplus.getVal();

    // 1+ decay
    
    RooSpinOne_Decay oneplus_decay("oneplus_decay","oneplus_decay", mzz,m1,m2,h1,h2,phi,
				  g1ValA, g2ValA, 
				  R1, R2, aParam, mZ, gamZ);

    Poneplus_decay = oneplus_decay.getVal();



    //
    // Spin 2 parameters
    // 
    
    
    RooRealVar fz1("fz1","fz1",0, 0.,1000.);  
    RooRealVar fz2("fz2","fz2",1, 0.,1000.);
    
    // - - - - - - - - - - - - - - - -
    // minimal coupling graviton parameters
    // - - - - - - - - - - - - - - - - 
    RooRealVar c1("c1","c1",1.,-1000.,1000.);
    RooRealVar c2("c2","c2",0.,-1000.,1000.);
    RooRealVar c3("c3","c3",0.,-1000.,1000.);
    RooRealVar c4("c4","c4",0.,-1000.,1000.);
    RooRealVar c5("c5","c5",0.,-1000.,1000.);
    RooRealVar c6("c6","c6",0.,-1000.,1000.);
    RooRealVar c7("c7","c7",0.,-1000.,1000.);

    RooRealVar useG("useG","useG",1.,-1000.,1000.);
    RooRealVar g1("g1","g1",1.,-1000.,1000.);
    RooRealVar g2("g2","g2",0.,-1000.,1000.);
    RooRealVar g3("g3","g3",0.,-1000.,1000.);
    RooRealVar g4("g4","g4",0.,-1000.,1000.);
    RooRealVar g5("g5","g5",1.,-1000.,1000.);
    RooRealVar g6("g6","g6",0.,-1000.,1000.);
    RooRealVar g7("g7","g7",0.,-1000.,1000.);
    RooRealVar g8("g8","g8",0.,-1000.,1000.);
    RooRealVar g9("g9","g9",0.,-1000.,1000.);
    RooRealVar g10("g10","g10",0.,-1000.,1000.);

    // 2m+ 
    
    RooSpinTwo_7D minGrav("minGrav","minGrav",mzz,m1,m2,hs,h1,h2,phi,phi1,
			  c1,c2,c3,c4,c5,c6,c7,
			  useG,
			  g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
			  fz1,fz2,
			  R1,R2,
			  mZ,gamZ);

    Ptwomplus_gg = minGrav.getVal();

    RooSpinTwo_Decay minGrav_decay("minGrav_decay","minGrav_decay",mzz,m1,m2,h1,h2,phi,
			  c1,c2,c3,c4,c5,c6,c7,
			  useG,
			  g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
			  fz1,fz2,
			  R1,R2,
			  mZ,gamZ);
    
    Ptwomplus_decay = minGrav_decay.getVal();

    // 2m+ (qqbar) 
    
    fz1.setVal(1.);
    fz2.setVal(0.);
    
    RooSpinTwo_7D minGravqq("minGravqq","minGravqq",mzz,m1,m2,hs,h1,h2,phi,phi1,
			  c1,c2,c3,c4,c5,c6,c7,
			  useG,
			  g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
			  fz1,fz2,
			  R1,R2,
			  mZ,gamZ);
    Ptwomplus_qq = minGravqq.getVal();
    
    // 2h-
    
    fz1.setVal(0.);
    fz2.setVal(0.);

    g1.setVal(0.);
    g2.setVal(0.);
    g3.setVal(0.);
    g4.setVal(0.);
    g5.setVal(0.);
    g6.setVal(0.);
    g7.setVal(0.);
    g8.setVal(1.);

    
    RooSpinTwo_7D twohminus("twohminus","twohminus",mzz,m1,m2,hs,h1,h2,phi,phi1,
			  c1,c2,c3,c4,c5,c6,c7,
			  useG,
			  g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
			  fz1,fz2,
			  R1,R2,
			  mZ,gamZ);
    Ptwohminus = twohminus.getVal();

    // 2h+
    
    fz1.setVal(0.);
    fz2.setVal(0.);
    
    g1.setVal(0.);
    g2.setVal(0.);
    g3.setVal(0.);
    g4.setVal(1.);
    g5.setVal(0.);
    g6.setVal(0.);
    g7.setVal(0.);
    g8.setVal(0.);

    
    RooSpinTwo_7D twohplus("twohplus","twohplus",mzz,m1,m2,hs,h1,h2,phi,phi1,
			  c1,c2,c3,c4,c5,c6,c7,
			  useG,
			  g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
			  fz1,fz2,
			  R1,R2,
			  mZ,gamZ);

    Ptwohplus  = twohplus.getVal();

    
    // 2b+
    
    fz1.setVal(0.);
    fz2.setVal(1.);

    g1.setVal(0.);
    g2.setVal(0.);
    g3.setVal(0.);
    g4.setVal(0.);
    g5.setVal(1.);
    g6.setVal(0.);
    g7.setVal(0.);
    g8.setVal(0.);
    
    
    RooSpinTwo_7D twobplus("twobplus","twobplus",mzz,m1,m2,hs,h1,h2,phi,phi1,
			  c1,c2,c3,c4,c5,c6,c7,
			  useG,
			  g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
			  fz1,fz2,
			  R1,R2,
			  mZ,gamZ);
    
    Ptwobplus  = twobplus.getVal();
    
    /*
    RooAbsReal* minGrav_decay_val = minGrav.createIntegral(RooArgSet(hs,phi1));
    Ptwomplus_decay_andrew = minGrav_decay_val->getVal();
    */

}

vector<TLorentzVector> Calculate4Momentum(double Mx,double M1,double M2,double theta,double theta1,double theta2,double Phi1,double Phi)
{
    double phi1,phi2;
    phi1=TMath::Pi()-Phi1;
    phi2=Phi1+Phi;
    
    
    double gamma1,gamma2,beta1,beta2;
    
    gamma1=(Mx*Mx+M1*M1-M2*M2)/(2*Mx*M1);
    gamma2=(Mx*Mx-M1*M1+M2*M2)/(2*Mx*M2);
    beta1=sqrt(1-1/(gamma1*gamma1));
    beta2=sqrt(1-1/(gamma2*gamma2));
    
    
    //gluon 4 vectors
    TLorentzVector p1CM(0,0,Mx/2,Mx/2);
    TLorentzVector p2CM(0,0,-Mx/2,Mx/2);
    
    //vector boson 4 vectors
    TLorentzVector kZ1(gamma1*M1*sin(theta)*beta1,0, gamma1*M1*cos(theta)*beta1,gamma1*M1*1);   
    TLorentzVector kZ2(-gamma2*M2*sin(theta)*beta2,0, -gamma2*M2*cos(theta)*beta2,gamma2*M2*1);
    
    //Rotation and Boost matrices. Note gamma1*beta1*M1=gamma2*beta2*M2.
    
    TLorentzRotation Z1ToZ,Z2ToZ;
    
    Z1ToZ.Boost(0,0,beta1);
    Z2ToZ.Boost(0,0,beta2);
    Z1ToZ.RotateY(theta);
    Z2ToZ.RotateY(TMath::Pi()+theta);
    
    
    //fermons 4 vectors in vector boson rest frame
    
    TLorentzVector p3Z1((M1/2)*sin(theta1)*cos(phi1),(M1/2)*sin(theta1)*sin(phi1),(M1/2)*cos(theta1),(M1/2)*1);
       
    TLorentzVector p4Z1(-(M1/2)*sin(theta1)*cos(phi1),-(M1/2)*sin(theta1)*sin(phi1),-(M1/2)*cos(theta1),(M1/2)*1);
      
    TLorentzVector p5Z2((M2/2)*sin(theta2)*cos(phi2),(M2/2)*sin(theta2)*sin(phi2),(M2/2)*cos(theta2),(M2/2)*1);
    
    TLorentzVector p6Z2(-(M2/2)*sin(theta2)*cos(phi2),-(M2/2)*sin(theta2)*sin(phi2),-(M2/2)*cos(theta2),(M2/2)*1);
      
    
    // fermions 4 vectors in CM frame
    
    TLorentzVector p3CM,p4CM,p5CM,p6CM;
    
    p3CM=Z1ToZ*p3Z1;
    p4CM=Z1ToZ*p4Z1;
    p5CM=Z2ToZ*p5Z2;
    p6CM=Z2ToZ*p6Z2;
    
    vector<TLorentzVector> p;
    
    p.push_back(p3CM);
    p.push_back(p4CM);
    p.push_back(p5CM);
    p.push_back(p6CM);

    return p;
}



//###################
//# main function
//###################
void runME_HZZ4l_JHUgenSamples(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity=TVar::INFO){

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

  TTree* ch=(TTree*)fin->Get("SelectedTree"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  evt_tree->SetName("newTree");
  
  // Declare the matrix element related variables to be added to the existing ntuples
  double dXsec_ZZ_DECAY_MCFM = 0.;
  double dXsec_ZZ_MCFM = 0.;
  double dXsec_HZZ_MCFM = 0.;
  double dXsec_HZZ_JHU = 0.;
  double dXsec_PSHZZ_JHU = 0.;
  double dXsec_HDHZZ_JHU = 0.;
  double dXsec_TZZ_JHU = 0.;
  double dXsec_VZZ_JHU = 0.;
  double dXsec_AVZZ_JHU = 0.;
  double dXsec_QQB_TZZ_JHU = 0.;
  double dXsec_TZZ_DECAY_JHU = 0.;
  double dXsec_VZZ_DECAY_JHU = 0.;
  double dXsec_AVZZ_DECAY_JHU = 0.;
  double dXsec_PTZZ_2hminus_JHU = 0.;
  double dXsec_TZZ_2hplus_JHU = 0.;
  double dXsec_TZZ_2bplus_JHU = 0.;

  double pseudoME = 0.;
  double graviME = 0.;
  
  evt_tree->Branch("dXsec_ZZ_DECAY_MCFM" , &dXsec_ZZ_DECAY_MCFM   ,"dXsec_ZZ_DECAY_MCFM/D");
  evt_tree->Branch("dXsec_ZZ_MCFM"   , &dXsec_ZZ_MCFM   ,"dXsec_ZZ_MCFM/D");
  evt_tree->Branch("dXsec_HZZ_MCFM"  , &dXsec_HZZ_MCFM   ,"dXsec_HZZ_MCFM/D");
  evt_tree->Branch("dXsec_HZZ_JHU"   , &dXsec_HZZ_JHU   ,"dXsec_HZZ_JHU/D");
  evt_tree->Branch("dXsec_PSHZZ_JHU" , &dXsec_PSHZZ_JHU ,"dXsec_PSHZZ_JHU/D");
  evt_tree->Branch("dXsec_HDHZZ_JHU" , &dXsec_HDHZZ_JHU ,"dXsec_HDHZZ_JHU/D");
  evt_tree->Branch("dXsec_TZZ_JHU"   , &dXsec_TZZ_JHU   ,"dXsec_TZZ_JHU/D");
  evt_tree->Branch("dXsec_VZZ_JHU"   , &dXsec_VZZ_JHU   ,"dXsec_VZZ_JHU/D");
  evt_tree->Branch("dXsec_AVZZ_JHU"  , &dXsec_AVZZ_JHU  ,"dXsec_AVZZ_JHU/D");
  evt_tree->Branch("dXsec_QQB_TZZ_JHU", &dXsec_QQB_TZZ_JHU   ,"dXsec_QQB_TZZ_JHU/D");
  evt_tree->Branch("dXsec_TZZ_DECAY_JHU", &dXsec_TZZ_DECAY_JHU   ,"dXsec_TZZ_DECAY_JHU/D");
  evt_tree->Branch("dXsec_VZZ_DECAY_JHU", &dXsec_VZZ_DECAY_JHU   ,"dXsec_VZZ_DECAY_JHU/D");
  evt_tree->Branch("dXsec_AVZZ_DECAY_JHU", &dXsec_AVZZ_DECAY_JHU   ,"dXsec_AVZZ_DECAY_JHU/D");
  evt_tree->Branch("dXsec_PTZZ_2hminus_JHU" , &dXsec_PTZZ_2hminus_JHU   ,"dXsec_PTZZ_2hminus_JHU/D");
  evt_tree->Branch("dXsec_TZZ_2hplus_JHU"   , &dXsec_TZZ_2hplus_JHU     ,"dXsec_TZZ_2hplus_JHU/D");
  evt_tree->Branch("dXsec_TZZ_2bplus_JHU"   , &dXsec_TZZ_2bplus_JHU     ,"dXsec_TZZ_2bplus_JHU/D");
  evt_tree->Branch("pseudoME"        , &pseudoME        ,   "pseudoME/D"  );
  evt_tree->Branch("graviME"         , &graviME         ,   "graviME/D"  );

  
  // 
  // analytical variables
  // 
  double Psmh_, Poneminus_, Poneplus_;
  double Poneminus_decay_, Poneplus_decay_;
  double Ptwomplus_gg_, Ptwomplus_qq_, Ptwomplus_decay_;
  double Ptwohminus_, Ptwohplus_, Ptwobplus_;

  evt_tree->Branch("Psmh"                , &Psmh_               ,"Psmh/D");
  evt_tree->Branch("Poneminus"           , &Poneminus_          ,"Poneminus/D");
  evt_tree->Branch("Poneplus"            , &Poneplus_           ,"Poneplus/D");
  evt_tree->Branch("Poneminus_decay"     , &Poneminus_decay_    ,"Poneminus_decay/D");
  evt_tree->Branch("Poneplus_decay"      , &Poneplus_decay_     ,"Poneplus_decay/D");
  evt_tree->Branch("Ptwomplus_gg"        , &Ptwomplus_gg_       ,"Ptwomplus_gg/D");
  evt_tree->Branch("Ptwomplus_qq"        , &Ptwomplus_qq_       ,"Ptwomplus_qq/D");
  evt_tree->Branch("Ptwomplus_decay"     , &Ptwomplus_decay_    ,"Ptwomplus_decay/D");
  evt_tree->Branch("Ptwohminus"          , &Ptwohminus_         ,"Ptwohminus/D");
  evt_tree->Branch("Ptwohplus"           , &Ptwohplus_          ,"Ptwohplus/D");
  evt_tree->Branch("Ptwobplus"           , &Ptwobplus_          ,"Ptwobplus/D");

 
  float m1,m2,h1,h2,hs,phi,phi1,mzz;  
  int mflavor = 3; // by default it is ee/mm 
  /*
  ch->SetBranchAddress( "z1mass"        , &m1      );   
  ch->SetBranchAddress( "z2mass"        , &m2      );   
  ch->SetBranchAddress( "costheta1"     , &h1      );   
  ch->SetBranchAddress( "costheta2"     , &h2      );   
  ch->SetBranchAddress( "costhetastar"  , &hs      );   
  ch->SetBranchAddress( "phi"           , &phi     );   
  ch->SetBranchAddress( "phistar1"      , &phi1    );   
  ch->SetBranchAddress( "zzmass"        , &mzz     );   
  */
  ch->SetBranchAddress( "Z1Mass"        , &m1      );   
  ch->SetBranchAddress( "Z2Mass"        , &m2      );   
  ch->SetBranchAddress( "helcosthetaZ1" , &h1      );   
  ch->SetBranchAddress( "helcosthetaZ2" , &h2      );   
  ch->SetBranchAddress( "costhetastar"  , &hs      );   
  ch->SetBranchAddress( "helphi"        , &phi     );   
  ch->SetBranchAddress( "phistarZ1"     , &phi1    );   
  ch->SetBranchAddress( "ZZMass"        , &mzz     );   
  if ( ch->GetBranchStatus("flavortype") ) 
    ch->SetBranchAddress( "flavortype"   , &mflavor);


 

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

  if ( ch->GetBranchStatus("EL1") ) 
    ch->SetBranchAddress( "EL1"       , &EL1_      );   
  if ( ch->GetBranchStatus("pXL1") ) 
    ch->SetBranchAddress( "pXL1"      , &pXL1_     );   
  if ( ch->GetBranchStatus("pYL1") ) 
    ch->SetBranchAddress( "pYL1"      , &pYL1_     );   
  if ( ch->GetBranchStatus("pZL1") ) 
    ch->SetBranchAddress( "pZL1"      , &pZL1_     );   
  if ( ch->GetBranchStatus("EL2") ) 
    ch->SetBranchAddress( "EL2"       , &EL2_      );   
  if ( ch->GetBranchStatus("pXL2") ) 
    ch->SetBranchAddress( "pXL2"      , &pXL2_     );   
  if ( ch->GetBranchStatus("pYL2") ) 
    ch->SetBranchAddress( "pYL2"      , &pYL2_     );   
  if ( ch->GetBranchStatus("pZL2") ) 
    ch->SetBranchAddress( "pZL2"      , &pZL2_      );   

  if ( ch->GetBranchStatus("EL3") ) 
    ch->SetBranchAddress( "EL3"       , &EL3_      );   
  if ( ch->GetBranchStatus("pXL3") ) 
    ch->SetBranchAddress( "pXL3"      , &pXL3_     );   
  if ( ch->GetBranchStatus("pYL3") ) 
    ch->SetBranchAddress( "pYL3"      , &pYL3_     );   
  if ( ch->GetBranchStatus("pZL3") ) 
    ch->SetBranchAddress( "pZL3"      , &pZL3_      );   

  if ( ch->GetBranchStatus("EL4") ) 
    ch->SetBranchAddress( "EL4"       , &EL4_      );   
  if ( ch->GetBranchStatus("pXL4") ) 
    ch->SetBranchAddress( "pXL4"      , &pXL4_     );   
  if ( ch->GetBranchStatus("pYL4") ) 
    ch->SetBranchAddress( "pYL4"      , &pYL4_     );   
  if ( ch->GetBranchStatus("pZL4") ) 
    ch->SetBranchAddress( "pZL4"      , &pZL4_      );   


  // Create the instance of TEvtProb to calculate the differential cross-section
  TEvtProb Xcal2;  
  hzz4l_event_type hzz4l_event;

  //==========================================
  // Loop All Events
  //==========================================
  

  int Ntot = maxevt > ch->GetEntries() ? ch->GetEntries() : maxevt; 
  if ( maxevt < 0. ) Ntot =  ch->GetEntries();

  pair<double,double> prob;

  if (verbosity >= TVar::INFO) printf("Total number of events = %d\n", Ntot);
  
  for(int ievt = 0; ievt < Ntot; ievt++){
    if (verbosity >= TVar::INFO && (ievt % 100 == 0)) 
      std::cout << "Doing Event: " << ievt << std::endl;
    
    // if ( ievt == 43167 ) continue;
    
    // 
    // initialise the differential cross-sections
    // 
    dXsec_ZZ_MCFM = 0.;
    dXsec_HZZ_MCFM = 0.;
    dXsec_HZZ_JHU = 0.;
    dXsec_PSHZZ_JHU = 0.;
    dXsec_HDHZZ_JHU = 0.;
    dXsec_TZZ_JHU = 0.;
    dXsec_VZZ_JHU = 0.;
    dXsec_AVZZ_JHU = 0.;
    dXsec_QQB_TZZ_JHU = 0.;
    dXsec_TZZ_DECAY_JHU = 0.;
    dXsec_VZZ_DECAY_JHU = 0.;
    dXsec_AVZZ_DECAY_JHU = 0.;
    dXsec_PTZZ_2hminus_JHU = 0.;
    dXsec_TZZ_2hplus_JHU = 0.;
    dXsec_TZZ_2bplus_JHU = 0.;
    dXsec_ZZ_DECAY_MCFM = 0.;
    pseudoME = 0.;
    graviME = 0.;
    
    ch->GetEntry(ievt);           
    
    /*
    cout << "mzz: " << mzz << " m1: " << m1 << " m2: " << m2 << endl;
    cout << "h1: " << h1 << " h2: " << h2 << " hs: " << hs << endl;
    cout << "phi: " << phi << " phi1: " << phi1 << endl;
    */
    
    // std::cout << "isnan(h2) = " << isnan(h2) << "\n";
    CalculateAnalyticalMELA(mzz, m1, m2, hs, h1, h2, phi, phi1,
			    Psmh_, Poneminus_, Poneplus_,
			    Poneminus_decay_, Poneplus_decay_,
			    Ptwomplus_gg_, Ptwomplus_qq_, Ptwomplus_decay_,
			    Ptwohminus_, Ptwohplus_, Ptwobplus_);

    // set four momenta
    vector<TLorentzVector> p;
    p=Calculate4Momentum(mzz,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi);
    // p=Calculate4Momentum(mzz,m1,m2,acos(0),acos(h1),acos(h2),0,phi);
   
    TLorentzVector Z1_minus = p[0];
    TLorentzVector Z1_plus  = p[1];
    TLorentzVector Z2_minus = p[2];
    TLorentzVector Z2_plus  = p[3];

    hzz4l_event.p[0].SetXYZM(Z1_minus.Px(), Z1_minus.Py(), Z1_minus.Pz(), 0.);
    hzz4l_event.p[1].SetXYZM(Z1_plus.Px(), Z1_plus.Py(), Z1_plus.Pz(), 0.);
    hzz4l_event.p[2].SetXYZM(Z2_minus.Px(), Z2_minus.Py(), Z2_minus.Pz(), 0.);
    hzz4l_event.p[3].SetXYZM(Z2_plus.Px(), Z2_plus.Py(), Z2_plus.Pz(), 0.);


    // flavor 1 for 4e, 2 for 4m, 3 for 2e2mu  
    if ( mflavor == 1 ) {
      hzz4l_event.PdgCode[0] = 11;
      hzz4l_event.PdgCode[1] = -11;
      hzz4l_event.PdgCode[2] = 11;
      hzz4l_event.PdgCode[3] = -11;
    }
    if ( mflavor == 2 ) {
      hzz4l_event.PdgCode[0] = 13;
      hzz4l_event.PdgCode[1] = -13;
      hzz4l_event.PdgCode[2] = 13;
      hzz4l_event.PdgCode[3] = -13;
    }
    if ( mflavor == 3 ) {
      hzz4l_event.PdgCode[0] = 11;
      hzz4l_event.PdgCode[1] = -11;
      hzz4l_event.PdgCode[2] = 13;
      hzz4l_event.PdgCode[3] = -13;
    }

    /*
    hzz4l_event.p[0].SetXYZM(pXL1_, pYL1_, pZL1_, 0.);
    hzz4l_event.p[1].SetXYZM(pXL2_, pYL2_, pZL2_, 0.);
    hzz4l_event.p[2].SetXYZM(pXL3_, pYL3_, pZL3_, 0.);
    hzz4l_event.p[3].SetXYZM(pXL4_, pYL4_, pZL4_, 0.);
    */

    double z1mass = (hzz4l_event.p[0]+hzz4l_event.p[1]).M();
    double z2mass = (hzz4l_event.p[2]+hzz4l_event.p[3]).M();
    double zzmass = (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).M();

    if (verbosity >= TVar::DEBUG) {
      cout << "\n=========================================================\n";
      cout << "Entry: " << ievt << "\n";
      cout << "Input: ==================================================" <<endl;
      printf("lep1 p3 = (%4.4f, %4.4f, %4.4f)  lep2 p3 = (%4.4f, %4.4f, %4.4f)\n",
	     p[0].Px(), p[0].Py(), p[0].Pz(), p[1].Px(), p[1].Py(), p[1].Pz()); 
      printf("lep3 p3 = (%4.4f, %4.4f, %4.4f)  lep4 p3 = (%4.4f, %4.4f, %4.4f)\n",
	     p[2].Px(), p[2].Py(), p[2].Pz(), p[3].Px(), p[3].Py(), p[3].Pz()); 
      std::cout << "ZZ system (pX, pY, pZ, E, mass) = ( " 
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Px() << ", "
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Py() << ", "
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Pz() << ", "
		<< (hzz4l_event.p[0]+hzz4l_event.p[1]+hzz4l_event.p[2]+hzz4l_event.p[3]).Energy()  << ", "
		<< zzmass << ")\n";
      std::cout << "Z1 mass = " << z1mass << "\tz2mass = " << z2mass << "\n";
      cout << "=========================================================\n";
    } 
    // finish loading event information
    
    // ==== Begin the differential cross-section calculation
    Xcal2.SetHiggsMass(zzmass);
    // calculate the ZZ using MCFM
    Xcal2.SetMatrixElement(TVar::MCFM);
    if ( mflavor < 3  )
      dXsec_ZZ_MCFM = Xcal2.XsecCalc(TVar::ZZ_4e,hzz4l_event,verbosity);
    else 
      dXsec_ZZ_MCFM = Xcal2.XsecCalc(TVar::ZZ_2e2m,hzz4l_event,verbosity);

    // trying to implement the integral
    
    dXsec_HZZ_MCFM = Xcal2.XsecCalc(TVar::HZZ_4l,hzz4l_event,verbosity);

    // calculate X->ZZ using JHUGen
    Xcal2.SetMatrixElement(TVar::JHUGen);
    
    // 0+ 
    dXsec_HZZ_JHU =  Xcal2.XsecCalc(TVar::HZZ_4l,hzz4l_event,verbosity);

    // 0-
    dXsec_PSHZZ_JHU = Xcal2.XsecCalc(TVar::PSHZZ_4l,hzz4l_event,verbosity);
    pseudoME =  dXsec_HZZ_JHU / ( dXsec_HZZ_JHU + 6*dXsec_PSHZZ_JHU );
    
    // 0h+
    dXsec_HDHZZ_JHU = Xcal2.XsecCalc(TVar::HDHZZ_4l,hzz4l_event,verbosity);

    // 1-
    dXsec_VZZ_JHU = Xcal2.XsecCalc(TVar::VZZ_4l,hzz4l_event,verbosity);    
    
    // 1+
    dXsec_AVZZ_JHU = Xcal2.XsecCalc(TVar::AVZZ_4l,hzz4l_event,verbosity);   

    // 1- decay only
    dXsec_VZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::VZZ_DECAY_4l,hzz4l_event,verbosity); 

    // 1+ decay only 
    dXsec_AVZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::AVZZ_DECAY_4l,hzz4l_event,verbosity); 
    
    // 2m+ 
    dXsec_TZZ_JHU = Xcal2.XsecCalc(TVar::TZZ_4l,hzz4l_event,verbosity);
    dXsec_QQB_TZZ_JHU = Xcal2.XsecCalc(TVar::QQB_TZZ_4l,hzz4l_event,verbosity);
    graviME =  dXsec_HZZ_JHU / ( dXsec_HZZ_JHU + 5*dXsec_QQB_TZZ_JHU );

    // 2m+ decay
    dXsec_TZZ_DECAY_JHU = Xcal2.XsecCalc(TVar::TZZ_DECAY_4l,hzz4l_event,verbosity);

    // 2h-
    dXsec_PTZZ_2hminus_JHU = Xcal2.XsecCalc(TVar::PTZZ_2hminus_4l,hzz4l_event,verbosity);

    // 2h+
    dXsec_TZZ_2hplus_JHU = Xcal2.XsecCalc(TVar::TZZ_2hplus_4l,hzz4l_event,verbosity);

    // 2b+
    dXsec_TZZ_2bplus_JHU = Xcal2.XsecCalc(TVar::TZZ_2bplus_4l,hzz4l_event,verbosity);

    // use the same constants defined in 
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CJLST/ZZMatrixElement/MELA/src/Mela.cc?revision=1.40&view=markup
    
    // add constants for discriminants
    // ---------------------------------
    if(mflavor==3){
      dXsec_PSHZZ_JHU*=6.0;
      dXsec_HDHZZ_JHU*=2.1;
      dXsec_VZZ_JHU*=16;
      dXsec_AVZZ_JHU*=13;
      dXsec_TZZ_JHU*=.6;
      dXsec_QQB_TZZ_JHU*=13;
    }else{
      dXsec_PSHZZ_JHU*=7.0;
      dXsec_HDHZZ_JHU*=2.3;
      dXsec_VZZ_JHU*=38;
      dXsec_AVZZ_JHU*=28;
      dXsec_TZZ_JHU*=1.4;
      dXsec_QQB_TZZ_JHU*=30;
    }

    // 
    // WARNING: tuned on the emu only at this point 
    // 

    // constants for JHUGen
    dXsec_TZZ_DECAY_JHU    *= 1.6e+9;
    dXsec_VZZ_DECAY_JHU    *= 1e+10;
    dXsec_AVZZ_DECAY_JHU   *= 1e+10;
    dXsec_PTZZ_2hminus_JHU *= 1e+10;
    dXsec_TZZ_2hplus_JHU   *= 1e+10;

    // constants for the analytical MELA
    Ptwomplus_gg_        *= 2e-07; 
    Ptwomplus_qq_        *= 2e-07; 
    Poneplus_            *= 5e-03;
    Poneminus_           *= 5e-03;
    Ptwomplus_decay_     *= 1e-08;
    Poneplus_decay_      *= 3.7e-04;
    Poneminus_decay_     *= 3.7e-04;
    Ptwohminus_           *= 0.4;
    Ptwohplus_           *= 0.4;
    Ptwobplus_           *= 2.3e-07;
    // ---------------------------------

    // 
    // calculate the last ZZ MCFM decay only MCFM
    // 
    
    int gridsize_hs = 10; 
    double hs_min = -1.;
    double hs_max = 1.;
    double hs_step =( hs_max - hs_min ) / double (gridsize_hs); 
    
    int gridsize_phi1 = 10; 
    double phi1_min = -TMath::Pi();
    double phi1_max = TMath::Pi();
    double phi1_step =( phi1_max - phi1_min ) / double (gridsize_phi1); 
    
    for ( int i_hs = 0; i_hs < gridsize_hs + 1; i_hs ++ ) {
      
      double hs_val = hs_min + i_hs * hs_step; 
      
      for ( int i_phi1 = 0; i_phi1 < gridsize_phi1 +1 ; i_phi1 ++ ) {
	double phi1_val = phi1_min + i_phi1 * phi1_step; 
	vector<TLorentzVector> p_decay;
	p_decay=Calculate4Momentum(mzz,m1,m2,acos(hs_val),acos(h1),acos(h2),phi1_val,phi);
	
	Z1_minus = p_decay[0];
	Z1_plus  = p_decay[1];
	Z2_minus = p_decay[2];
	Z2_plus  = p_decay[3];
	
	hzz4l_event.p[0].SetXYZM(Z1_minus.Px(), Z1_minus.Py(), Z1_minus.Pz(), 0.);
	hzz4l_event.p[1].SetXYZM(Z1_plus.Px(), Z1_plus.Py(), Z1_plus.Pz(), 0.);
	hzz4l_event.p[2].SetXYZM(Z2_minus.Px(), Z2_minus.Py(), Z2_minus.Pz(), 0.);
	hzz4l_event.p[3].SetXYZM(Z2_plus.Px(), Z2_plus.Py(), Z2_plus.Pz(), 0.);

	// calculate the ZZ using MCFM
	Xcal2.SetMatrixElement(TVar::MCFM);
	if ( mflavor < 3  )
	  dXsec_ZZ_DECAY_MCFM += Xcal2.XsecCalc(TVar::ZZ_4e,hzz4l_event,verbosity);
	else 
	  dXsec_ZZ_DECAY_MCFM += Xcal2.XsecCalc(TVar::ZZ_2e2m,hzz4l_event,verbosity);
      }
    }
    dXsec_ZZ_DECAY_MCFM =  dXsec_ZZ_DECAY_MCFM /  double ( (gridsize_hs + 1) * (gridsize_phi1 +1 )); 
    
    evt_tree->Fill();
    
  }//nevent
  
  if (verbosity >= TVar::INFO) cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
}  
