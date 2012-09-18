#ifndef _TEVTPROB_HH_
#define _TEVTPROB_HH_
//-----------------------------------------------------------------------------
// Description: Class TEvtProb: EvtProb base class
// ------------
//
//      Event Probability Density Calculation
//
// Feb 21 2011
// Sergo Jindariani
// Yanyan Gao
//-----------------------------------------------------------------------------
#include <sstream>
#include <string>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "assert.h"
#include "TROOT.h"
// ME related
#include "TMCFM.hh"
#include "TVar.hh"
#include "TUtil.hh"


//----------------------------------------
// Class TEvtProb
//----------------------------------------
class TEvtProb : public TObject {
  
public:
  //--------------------
  // Variables
  //--------------------
  TVar::Process _process;
  TVar::MatrixElement _matrixElement;
  
  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb();
  ~TEvtProb();
  
  //----------------------
  // Function
  //----------------------
  void SetProcess(TVar::Process tmp) { _process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }
  void XsecCalc(TVar::Process proc,
		         const hzz4l_event_type &hzz4l_event,
			 double Xsec, double XsecErr, TVar::VerbosityLevel verbosity);
  ClassDef(TEvtProb,0);
};

#endif

