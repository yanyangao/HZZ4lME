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
  myCSW_ = new HiggsCSandWidth();
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

    // This is bad.  It seems to be setting some kind of 
    // global BW function somewhere.
    // This function, My_choose, is defined in TUtil
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
    double msqjk = SumMatrixElementPDF(proc, &mcfm_event, flavor_msq, &flux);
    if(msqjk<=0){ mcfm_event.pswt=0; }
    
    flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
    dXsec=msqjk*flux;
    
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
}


