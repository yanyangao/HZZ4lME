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
}

void TEvtProb::TestModHiggsMatEl()
{
// input unit = GeV/100 such that 125GeV is 1.25 in the code
  double MReso = 125.0/100.0;
  double GaReso= 0.1/100.0;
  double P[6][4];
  double MatElSq;
  int MYIDUP[4];

// particle ID: +7=e+,  -7=e-,  +8=mu+,  -8=mu-
  MYIDUP[0]=+7;
  MYIDUP[1]=-7;
  MYIDUP[2]=+7;
  MYIDUP[3]=-7;

// p(i,0:3) = (E(i),px(i),py(i),pz(i))
// i=0,1: glu1,glu2 (outgoing convention)
// i=2,3: correspond to MY_IDUP(1),MY_IDUP(0)
// i=4,5: correspond to MY_IDUP(3),MY_IDUP(2)
  P[0][0]=-0.1210830252829298;
  P[0][1]=0.0;
  P[0][2]=0.0;
  P[0][3]=-0.1210830252829298;

  P[1][0]=-3.2260921717742228;
  P[1][1]=0.0;
  P[1][2]=0.0;
  P[1][3]=3.2260921717742228;

  P[2][0]=2.3786499379469044;
  P[2][1]=0.4825158868534503;
  P[2][2]=0.1249853079677575;
  P[2][3]=-2.3258401963636808;

  P[3][0]=0.5705642134010264;
  P[3][1]=-0.3042379113840483;
  P[3][2]=-0.0252235463240039;
  P[3][3]=-0.4820234305523382;

  P[4][0]=0.1627239095197414;
  P[4][1]=-0.0193830498689781;
  P[4][2]=0.0081768728141093;
  P[4][3]=-0.1613583182180193;

  P[5][0]=0.2352371361894799;
  P[5][1]=-0.1588949256004239;
  P[5][2]=-0.1079386344578630;
  P[5][3]=-0.1357872013572543;


  __modhiggs_MOD_evalamp_gg_h_vv(P, &MReso,  &GaReso, MYIDUP, &MatElSq);

 printf("\n ");
 printf("Matr.el. squared: %20.17e \n ",MatElSq);
 printf("result should be: %20.17e \n ",0.0112810555699413);
 printf("ratio: %20.17e \n ",MatElSq/0.0112810555699413);

}


