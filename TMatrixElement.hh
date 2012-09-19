#ifndef _TMATRIXELEMENT_HH_
#define _TMATRIXELEMENT_HH_
//-----------------------------------------------------------------------------
// Description: Class TMatrixElement: MatrixElement base class
//---------------------------------------------------------------- ------------
//
//      Interface to Different Matrix Element
//
// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "TObject.h"

#include "TVar.hh"
#include "TMCFM.hh"

//----------------------------------------
// Class TMatrixElement
//----------------------------------------

class TMatrixElement : public TObject {

  private:
     static TMatrixElement* globalObj_;

  public:

  //--------------------
  // Static Object
  //--------------------

   static TMatrixElement* inst(){
     if(! globalObj_)
          globalObj_ = new TMatrixElement;
     return globalObj_;
   }

  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TMatrixElement();
 ~TMatrixElement();

  //--------------------
  // Variables
  //--------------------
  TVar::MatrixElement _matrixElement;

  //--------------------
  // Function
  //--------------------
  double SumMatrixElementPDF(TVar::Process process, double p[][mxpart],double flavor_msq[][11],double* flux);

  int     MCFM_MatrixElement(TVar::Process process,double p[][mxpart],double flavor_msq[][11]);
  int MadGraph_MatrixElement(TVar::Process process,double p[][mxpart],double flavor_msq[][11]);


  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }


  ClassDef(TMatrixElement,0)

};


#endif
