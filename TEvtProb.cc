//-----------------------------------------------------------------------------
//
// Class EventProb Module
//
//   EventProb Module
//
// March 21 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "TEvtProb.hh"
#include "TVar.hh"


ClassImp(TEvtProb)

    using namespace std;

    //-----------------------------------------------------------------------------
    // Constructors and Destructor
    //-----------------------------------------------------------------------------
TEvtProb::TEvtProb() {
  mcfm_init_();
  SetEwkCoupligParameters();
  coupling_();
  myCSW_ = new HiggsCSandWidth("../txtFiles");
}

TEvtProb::~TEvtProb() {
    delete myCSW_;
}

//
// Directly calculate the ZZ->4l differential cross-section 
// WARNING: in development
// 
double TEvtProb::XsecCalc(TVar::Process proc, const hzz4l_event_type &hzz4l_event,
			  TVar::VerbosityLevel verbosity){

    //Initialize Process
    SetProcess(proc);
    
    if ( _matrixElement == TVar::MCFM) 
      My_choose(proc);
    
    //constants
    double sqrts = 2.*EBEAM;
    double W=sqrts*sqrts;
    
    //Weight calculation
    double flux=1.;
    double dXsec=0.;
    
    mcfm_event_type mcfm_event; 
    // assign the right initial momentum
    // assumes the events are boosted to have 0 transverse momenta
    double sysPz= ( hzz4l_event.p[0] + hzz4l_event.p[1] + hzz4l_event.p[2] + hzz4l_event.p[3] ).Pz();
    double sysE = ( hzz4l_event.p[0] + hzz4l_event.p[1] + hzz4l_event.p[2] + hzz4l_event.p[3] ).Energy();
    double pz0 = (sysE+sysPz)/2.; 
    double pz1 = -(sysE-sysPz)/2.;
    mcfm_event.p[0].SetPxPyPzE   (0., 0., pz0, TMath::Abs(pz0));
    mcfm_event.p[1].SetPxPyPzE   (0., 0., pz1, TMath::Abs(pz1));
    mcfm_event.p[2].SetPxPyPzE   (hzz4l_event.p[0].Px(), hzz4l_event.p[0].Py(), hzz4l_event.p[0].Pz(), hzz4l_event.p[0].Energy());
    mcfm_event.p[3].SetPxPyPzE   (hzz4l_event.p[1].Px(), hzz4l_event.p[1].Py(), hzz4l_event.p[1].Pz(), hzz4l_event.p[1].Energy());
    mcfm_event.p[4].SetPxPyPzE   (hzz4l_event.p[2].Px(), hzz4l_event.p[2].Py(), hzz4l_event.p[2].Pz(), hzz4l_event.p[2].Energy());
    mcfm_event.p[5].SetPxPyPzE   (hzz4l_event.p[3].Px(), hzz4l_event.p[3].Py(), hzz4l_event.p[3].Pz(), hzz4l_event.p[3].Energy());
    
    /*
    for ( int i = 0; i < 6; i++ ) {
      std::cout << "Particle " << i << " (Px, Py, Pz, E): " <<  mcfm_event.p[i].Px() << ", " << mcfm_event.p[i].Py() 
		<< ", " << mcfm_event.p[i].Pz() << ", " << mcfm_event.p[i].Energy() <<  "\n";
    }
    */
    //Matrix Element evaluation in qX=qY=0 frame
    //Evaluate f(x1)f(x2)|M(q)|/x1/x2 
    // 
    double qX=mcfm_event.p[0].Px()+mcfm_event.p[1].Px();
    double qY=mcfm_event.p[0].Py()+mcfm_event.p[1].Py();
    
    if((qX*qX+qY*qY)>0){
      double qE = mcfm_event.p[0].Energy()+mcfm_event.p[1].Energy();
      TVector3 boostV(qX/qE,qY/qE,0);
      for(int ipt=0;ipt<6;ipt++) mcfm_event.p[ipt].Boost(-boostV);
    }
    //event selections in Lab Frame
    double flavor_msq[nmsq][nmsq];
    double msqjk; 
    if ( _matrixElement == TVar::MCFM ) 
      msqjk = SumMatrixElementPDF(proc, &mcfm_event, flavor_msq, &flux);
    if ( _matrixElement == TVar::JHUGen ) {
      //
      // spin 0 
      //
      if ( proc == TVar::HZZ_4l || proc == TVar::PSHZZ_4l ) {
	// By default set the Spin 0 couplings for SM case
	// H->Glu Glu coupling constants 
	double hggcoupl[3] = {1.0, 0.0, 0.0};  
	// H->ZZ coupling constants 
	double hzzcoupl[4] = {1.0, 0.0, 0.0, 0.0}; 
	if ( proc == TVar::PSHZZ_4l ) {
	  hzzcoupl[0] = 0.0;
	  hzzcoupl[1] = 0.0;
	  hzzcoupl[2] = 0.0;
	  hzzcoupl[3] = 1.0;
	}
	msqjk = JHUGenMatEl(proc, &mcfm_event, _hmass, _hwidth, hggcoupl, hzzcoupl);
      }

      //
      // spin 2 
      // 
      if ( proc == TVar::TZZ_4l ) {
	// Graviton->Glu Glu coupling constants 
	double Gggcoupl[5] = {1.0, 0.0, 0.0, 0.0, 0.0}; // 2m+
	// double Gggcoupl[5] = {0.0, 0.0, 0.0, 1.0, 0.0}; // 2h+
	// double Gggcoupl[5] = {0.0, 1.0, 1.0, 0.0, 0.0}; // 2L+
	// double Gggcoupl[5] = {0.0, 0.0, 0.0, 0.0, 1.0}; // 2h-
	// Graviton->ZZ coupling constants 
	double Gvvcoupl[10]; 
	Gvvcoupl[0]=1.0; // 2m+
	Gvvcoupl[1]=0.0; 
	Gvvcoupl[2]=0.0;  
	Gvvcoupl[3]=0.0; // 2h+
	Gvvcoupl[4]=1.0; // 2m+ 
	Gvvcoupl[5]=0.0;
	Gvvcoupl[6]=0.0;
	Gvvcoupl[7]=0.0; // 2h-
	Gvvcoupl[8]=0.0;
	Gvvcoupl[9]=0.0;
	msqjk = JHUGenMatEl(proc, &mcfm_event, _hmass, _hwidth, Gggcoupl, Gvvcoupl);
      }
    }


    if(msqjk<=0){ mcfm_event.pswt=0; }
    
    flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
    //    dXsec=msqjk*flux;
    dXsec=msqjk;

    
    if (verbosity >= TVar::DEBUG)
      {
	cout <<" TEvtProb::XsecCalc(): dXsec=" << dXsec
	     <<" Msq="<<msqjk 
	     <<" flux="<<flux 
	     <<endl;
      }

    return dXsec;

}

// this appears to be some kind of 
// way of setting MCFM parameters through
// an interface defined in TMCFM.hh
void TEvtProb::SetHiggsMass(double mass){
    masses_mcfm_.hmass=mass;
    masses_mcfm_.hwidth=myCSW_->HiggsWidth(0, mass);
    _hmass = mass;
    _hwidth = myCSW_->HiggsWidth(0, mass);
    /*
    //
    // get higgs width for 125 and 250 GeV
    // 
    std::cout << "H125 width " << myCSW_->HiggsWidth(0, 125);
    std::cout << "H250 width " << myCSW_->HiggsWidth(0, 250);
    */
}
