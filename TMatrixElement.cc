//-----------------------------------------------------------------------------
//
// Class MatrixElement Module
//
//   MatrixElement Module
//
// March 28 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------
#include "iostream"
#include "TMatrixElement.hh"
#include "TUtil.hh"

ClassImp(TMatrixElement)

//#########################
//## Static Object
//#########################

TMatrixElement* TMatrixElement::globalObj_;

//--------------------------------
//   Constructor
//--------------------------------

TMatrixElement::TMatrixElement() {
  _matrixElement = TVar::MCFM;
}

TMatrixElement::~TMatrixElement() {}

//
double TMatrixElement::SumMatrixElementPDF(TVar::Process process, double p4[][mxpart],double flavor_msq[][nmsq],double* flux){


           int NPart=npart_.npart+2;
 	   double fx1[nmsq];
           double fx2[nmsq];
           double msq[nmsq][nmsq];
 

           double xx[2]={-p4[3][0]/EBEAM,-p4[3][1]/EBEAM};
           if(xx[0] > 1.0 || xx[0]<=xmin_.xmin) return 0.0;
	   if(xx[1] > 1.0 || xx[1]<=xmin_.xmin) return 0.0;

    
	   //calculate invariant masses between partons/final state particles
	   double s[mxpart][mxpart];
	   for(int jdx=0;jdx< NPart ;jdx++){
	     s[jdx][jdx]=0;
	     for(int kdx=jdx+1;kdx<NPart;kdx++){
	       s[jdx][kdx]=2*(p4[3][jdx]*p4[3][kdx]-p4[2][jdx]*p4[2][kdx]-p4[1][jdx]*p4[1][kdx]-p4[0][jdx]*p4[0][kdx]);
	       s[kdx][jdx]=s[jdx][kdx];
	     }
	   }

    //remove events has small invariant mass
     if(My_smalls(s,npart_.npart)) return 0.0;

	   
           //Calculate Pdf
           //Always pass address through fortran function
           fdist_ (&density_.ih1, &xx[0], &scale_.scale, fx1); //P+=> W+->e+nu
           fdist_ (&density_.ih2, &xx[1], &scale_.scale, fx2); //P-=> W-->e-nu

           //Calculate Mateix Elements
           // Matrix Elements calculation based on qt boosted four vectors
                if (_matrixElement==TVar::MCFM    )         MCFM_MatrixElement (process,p4,msq);
           else if (_matrixElement==TVar::MadGraph)     MadGraph_MatrixElement (process,p4,msq);
           
          
	   double msqjk=0;
	   for (int ii=0;ii<nmsq;ii++){
	     for (int jj=0;jj<nmsq;jj++){
             //2-D matrix is reversed in fortran
             // msq[ parton2 ] [ parton1 ]
	     flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];
             msqjk+=flavor_msq[jj][ii];
	     }
	   }
	   
	   (*flux)=fbGeV2/(8*xx[0]*xx[1]*EBEAM*EBEAM);

	   return msqjk;

}

//---------------------------------------------
//
// MCFM MatrixElement routine
//
//---------------------------------------------
/* int TMatrixElement::MCFM_MatrixElement(TVar::Process process,double p4[][mxpart],double msq[][nmsq]){
 
   cout<<"Current process is: "<<process<<"\n";
                 if(process==TVar::WW)      qqb_ww_  (p4[0],msq[0]);
            else if(process==TVar::WZ)      qqb_wz_  (p4[0],msq[0]);
            else if(process==TVar::ZZ)      qqb_zz_  (p4[0],msq[0]);
            else if(process==TVar::ZZ_4l)   qqb_zz_  (p4[0],msq[0]);
            else if(process==TVar::Wp_gamma)qqb_wgam_(p4[0],msq[0]);
            else if(process==TVar::Wm_gamma)qqb_wgam_(p4[0],msq[0]);
            else if(process==TVar::Wp_1jet )qqb_w_g_ (p4[0],msq[0]);
            else if(process==TVar::Wm_1jet )qqb_w_g_ (p4[0],msq[0]);
            else if(process==TVar::HWW)     qqb_hww_ (p4[0],msq[0]);
            else if(process==TVar::HZZ)     qqb_hzz_ (p4[0],msq[0]);
            else if(process==TVar::Z_2l)    qqb_z_   (p4[0],msq[0]);

   return 0;
}
*/

