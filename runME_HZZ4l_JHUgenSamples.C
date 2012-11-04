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
#include "RooRealVar.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

double ERRORthreshold=1.0;
using namespace std;

void xseccalc(TString inputDir, TString fileName, TString outputDir, int maxevt, TVar::VerbosityLevel verbosity);

void checkZorder(double& z1mass, double& z2mass,
		 double& costhetastar, double& costheta1,
		 double& costheta2, double& phi, 
		 double& phistar1){
  
  double tempZ1mass=z1mass;
  double tempZ2mass=z2mass;
  double tempH1=costheta1;
  double tempH2=costheta2;
  double tempHs=costhetastar;
  double tempPhi1=phistar1;
  double tempPhi=phi;

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

pair<double,double> calculateGraviMela(double mzz_,double m1_, double m2_,
				       double hs_, double h1_, double h2_, 
				       double phi_, double phi1_){

  // - - - - - - - - - - - - - - -
  // measurables
  // - - - - - - - - - - - - - - - 

  RooRealVar mzz("mzz","mzz",100.,1000.);
  RooRealVar m1("m1","m1",0.,1000.);
  RooRealVar m2("m2","m2",0.,1000.);

  RooRealVar hs("hs","hs",-1.,1.);
  RooRealVar h1("h1","h1",-1.,1.);
  RooRealVar h2("h2","h2",-1.,1.);

  RooRealVar phi("phi","phi",-3.2,3.2);
  RooRealVar phi1("phi1","phi1",-3.2,3.2);

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

  RooRealVar fz1("fz1","fz1",0.,0.,1000.);  
  RooRealVar fz2("fz2","fz2",1.,0.,1000.);

  // - - - - - - - - - - - - - - - -
  // SM Higgs PDF
  // - - - - - - - - - - - - - - - - 
  
  RooXZsZs_5D SMHiggs("SMHiggs","SMHiggs",m1,m2,h1,h2,phi,a1,p1,a2,p2,a3,p3,mZ,gamZ,mzz,R1,R2);

  // - - - - - - - - - - - - - - - -
  // PS Higgs PDF
  // - - - - - - - - - - - - - - - -

  RooRealVar a1PS("a1PS","a1PS",0.,-1000.,1000.);
  RooRealVar a2PS("a2PS","a2PS",0.,-1000.,1000.);
  RooRealVar a3PS("a3PS","a3PS",1.,-1000.,1000.);

  RooXZsZs_5D PSHiggs("PSHiggs","PSHiggs",m1,m2,h1,h2,phi,a1PS,p1,a2PS,p2,a3PS,p3,mZ,gamZ,mzz,R1,R2);
  


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

  RooSpinTwo_7D minGrav("minGrav","minGrav",mzz,m1,m2,hs,h1,h2,phi,phi1,
			c1,c2,c3,c4,c5,c6,c7,
			useG,
			g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,
			fz1,fz2,
			R1,R2,
			mZ,gamZ);

  //
  // Spin One parameters
  // 
  
  // 1-
  RooRealVar g1ValV("g1ValV","g1ValV",1.);
  RooRealVar g2ValV("g2ValV","g2ValV",0.);
  // this is the acceptance term associated with the production angles
  // the default setting is for setting no-acceptance
  RooRealVar aParam("aParam","aParam",0);
  RooSpinOne_7D zprime("zprime","zprime", mzz,m1,m2,h1,h2,hs,phi,phi1,
		                          g1ValV, g2ValV, 
		                          R1, R2, aParam, mZ, gamZ);
  

  checkZorder(m1_,m2_,hs_,h1_,h2_,phi_,phi1_);

  mzz.setVal(mzz_);  m1.setVal(m1_);   m2.setVal(m2_);
  hs.setVal(hs_);  h1.setVal(h1_);   h2.setVal(h2_);  
  phi.setVal(phi_);  phi1.setVal(phi1_);  
  
  pair<double,double> result;
  result.first=SMHiggs.getVal();
  //result.second=minGrav.getVal();
  //result.second=zprime.getVal();
  result.second=PSHiggs.getVal();

  return result;
  
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

  TTree* ch=(TTree*)fin->Get("angles"); 
  if (ch==0x0) return; 
  TTree* evt_tree=(TTree*) ch->CloneTree(0, "fast");
  
  // Declare the matrix element related variables to be added to the existing ntuples
  double dXsec_ZZ = 0.;
  double dXsecErr_ZZ = 0.;

  double dXsec_HZZ = 0.;
  double dXsecErr_HZZ = 0.;

  double dXsec_HZZ_JHU = 0.;
  double dXsecErr_HZZ_JHU = 0.;

  double dXsec_PSHZZ_JHU = 0.;
  double dXsecErr_PSHZZ_JHU = 0.;

  double dXsec_TZZ_JHU = 0.;
  double dXsecErr_TZZ_JHU = 0.;

  double dXsec_VZZ_JHU = 0.;
  double dXsecErr_VZZ_JHU = 0.;

  
  evt_tree->Branch("dXsec_HZZ"   , &dXsec_HZZ      ,"dXsec_HZZ/D");
  evt_tree->Branch("dXsecErr_HZZ", &dXsecErr_HZZ  ,"dXsecErr_HZZ/D");

  evt_tree->Branch("dXsec_HZZ_JHU"   , &dXsec_HZZ_JHU      ,"dXsec_HZZ_JHU/D");
  evt_tree->Branch("dXsecErr_HZZ_JHU", &dXsecErr_HZZ_JHU   ,"dXsecErr_HZZ_JHU/D");

  evt_tree->Branch("dXsec_PSHZZ_JHU"   , &dXsec_PSHZZ_JHU      ,"dXsec_PSHZZ_JHU/D");
  evt_tree->Branch("dXsecErr_PSHZZ_JHU", &dXsecErr_PSHZZ_JHU   ,"dXsecErr_PSHZZ_JHU/D");

  evt_tree->Branch("dXsec_TZZ_JHU"   , &dXsec_TZZ_JHU      ,"dXsec_TZZ_JHU/D");
  evt_tree->Branch("dXsecErr_TZZ_JHU", &dXsecErr_TZZ_JHU   ,"dXsecErr_TZZ_JHU/D");

  evt_tree->Branch("dXsec_VZZ_JHU"   , &dXsec_VZZ_JHU      ,"dXsec_VZZ_JHU/D");
  evt_tree->Branch("dXsecErr_VZZ_JHU", &dXsecErr_VZZ_JHU   ,"dXsecErr_VZZ_JHU/D");
  
  double psig_new, pbkg_new, graviMela;

  evt_tree->Branch("Psmh_new"   , &psig_new      ,"Psmh_new/D");
  evt_tree->Branch("Pgrav_new", &pbkg_new   ,"Pgrav_new/D");
  evt_tree->Branch("graviMela"   , &graviMela      ,"graviMela/D");
 
  double m1,m2,h1,h2,hs,phi,phi1,mzz;  
  int mflavor = 3; // by default it is ee/mm 
  ch->SetBranchAddress( "z1mass"        , &m1      );   
  ch->SetBranchAddress( "z2mass"        , &m2      );   
  ch->SetBranchAddress( "costheta1"     , &h1      );   
  ch->SetBranchAddress( "costheta2"     , &h2      );   
  ch->SetBranchAddress( "costhetastar"  , &hs      );   
  ch->SetBranchAddress( "phi"           , &phi     );   
  ch->SetBranchAddress( "phistar1"      , &phi1    );   
  ch->SetBranchAddress( "zzmass"        , &mzz     );   
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
  hzz4l_event_type hzz4l_event;
  //==========================================
  // Loop All Events
  //==========================================
  
  int Ntot = maxevt > ch->GetEntries() ? ch->GetEntries() : maxevt; 
  
  pair<double,double> prob;

  if (verbosity >= TVar::INFO) printf("Total number of events = %d\n", Ntot);
  
  for(int ievt = 0; ievt < Ntot; ievt++){
    if (verbosity >= TVar::INFO && (ievt % 1000 == 0)) 
      std::cout << "Doing Event: " << ievt << std::endl;
    
    // 
    // initialise the differential cross-sections
    // 
    dXsec_ZZ = 0.;
    dXsecErr_ZZ = 0.;

    dXsec_HZZ = 0.;
    dXsecErr_HZZ = 0.;
    
    dXsec_HZZ_JHU = 0.;
    dXsecErr_HZZ_JHU = 0.;

    dXsec_PSHZZ_JHU = 0.;
    dXsecErr_PSHZZ_JHU = 0.;

    dXsec_TZZ_JHU = 0.;
    dXsecErr_TZZ_JHU = 0.;

    dXsec_VZZ_JHU = 0.;
    dXsecErr_VZZ_JHU = 0.;

    ch->GetEntry(ievt);           

    /*
    cout << "mzz: " << mzz << " m1: " << m1 << " m2: " << m2 << endl;
    cout << "h1: " << h1 << " h2: " << h2 << " hs: " << hs << endl;
    cout << "phi: " << phi << " phi1: " << phi1 << endl;
    */

    prob = calculateGraviMela(mzz,m1,m2,hs,h1,h2,phi,phi1);

    psig_new = prob.first;
    pbkg_new = prob.second;
    graviMela = 1./(1.+pbkg_new/psig_new);

    vector<TLorentzVector> p;
    p=Calculate4Momentum(mzz,m1,m2,acos(hs),acos(h1),acos(h2),phi1,phi);
   
    TLorentzVector Z1_minus = p[0];
    TLorentzVector Z1_plus  = p[1];
    TLorentzVector Z2_minus = p[2];
    TLorentzVector Z2_plus  = p[3];

    //if((Z1_minus+Z1_plus+ Z2_minus+Z2_plus).M()>130 || (Z1_minus+Z1_plus+ Z2_minus+Z2_plus).M()<120) continue; // || (Z1_minus+Z1_plus).M()<40 || (Z2_minus+Z2_plus).M()<12) continue;

    hzz4l_event.p[0].SetXYZM(p[0].Px(), p[0].Py(), p[0].Pz(), 0.);
    hzz4l_event.p[1].SetXYZM(p[1].Px(), p[1].Py(), p[1].Pz(), 0.);
    hzz4l_event.p[2].SetXYZM(p[2].Px(), p[2].Py(), p[2].Pz(), 0.);
    hzz4l_event.p[3].SetXYZM(p[3].Px(), p[3].Py(), p[3].Pz(), 0.);

    // flavor 1 for 4e, 2 for 4m, 3 for 2e2mu  
    if ( mflavor == 1 ) {
      hzz4l_event.PdgCode[0] = 13;
      hzz4l_event.PdgCode[1] = -13;
      hzz4l_event.PdgCode[2] = 13;
      hzz4l_event.PdgCode[3] = -13;
    }
    if ( mflavor == 2 ) {
      hzz4l_event.PdgCode[0] = 11;
      hzz4l_event.PdgCode[1] = -11;
      hzz4l_event.PdgCode[2] = 11;
      hzz4l_event.PdgCode[3] = -11;
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

    // if ( TMath::Abs(z1mass - 91.1876) > 15. || TMath::Abs(z2mass - 91.1875) > 15. ) continue;
    
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
    dXsec_ZZ = Xcal2.XsecCalc(TVar::ZZ_4l,hzz4l_event,verbosity);
    dXsec_HZZ = Xcal2.XsecCalc(TVar::HZZ_4l,hzz4l_event,verbosity);
    // calculate X->ZZ using JHUGen
    // 0+ 
    Xcal2.SetMatrixElement(TVar::JHUGen);
    dXsec_HZZ_JHU = Xcal2.XsecCalc(TVar::HZZ_4l,hzz4l_event,verbosity);
    // 0-
    Xcal2.SetMatrixElement(TVar::JHUGen);
    dXsec_PSHZZ_JHU = Xcal2.XsecCalc(TVar::PSHZZ_4l,hzz4l_event,verbosity);
    // spin 2
    Xcal2.SetMatrixElement(TVar::JHUGen);
    dXsec_TZZ_JHU = Xcal2.XsecCalc(TVar::TZZ_4l,hzz4l_event,verbosity);
    // spin 1
    Xcal2.SetMatrixElement(TVar::JHUGen);
    dXsec_VZZ_JHU = Xcal2.XsecCalc(TVar::VZZ_4l,hzz4l_event,verbosity);    
    
    evt_tree->Fill();
    
  }//nevent
  
  if (verbosity >= TVar::INFO) cout << "TotCalculatedEvts: " << evt_tree->GetEntries() <<endl; 
  
  newfile->cd(); 
  evt_tree->Write(); 
  newfile->Close();
  
}  
