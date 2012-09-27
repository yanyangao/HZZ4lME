#include "fstream"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TUtil.hh"
#include "TCanvas.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "TProfile.h"

using namespace std;

void SetEwkCoupligParameters(){
  
  ewinput_.Gf_inp=1.16639E-05;
  ewinput_.aemmz_inp=7.81751E-03;
  ewinput_.wmass_inp=79.956049884402844;
  ewinput_.zmass_inp=91.1876;
  ewinput_.xw_inp=0.23116864;

}


void My_choose(TVar::Process process){
 
//ZZ_4l
if(process==TVar::ZZ_4l ){ 
 
    //81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
    //86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'

    npart_.npart=4;
    nqcdjets_.nqcdjets=0;

    vsymfact_.vsymfact=1.0;                                                                                                               
    interference_.interference=false;

    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;

       
    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    zcouple_.q2=-1.;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;

 }  else if ( process == TVar::HZZ_4l) {

  //  114 '  f(p1)+f(p2) --> H(--> Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))' 'N'
     npart_.npart=4;
     nqcdjets_.nqcdjets=0;

     bveg1_mcfm_.ndim=10;
     masses_mcfm_.mb=0;

     breit_.n2=1;
     breit_.n3=1;

     breit_.mass2 =masses_mcfm_.zmass;
     breit_.width2=masses_mcfm_.zwidth;
     breit_.mass3 =masses_mcfm_.zmass;
     breit_.width3=masses_mcfm_.zwidth;
     
     zcouple_.l1=zcouple_.le;
     zcouple_.r1=zcouple_.re;
     
     zcouple_.l2=zcouple_.le;
     zcouple_.r2=zcouple_.re;
 } 
 else{
     std::cerr <<"[My_choose]: Can't identify Process: " << process <<endl;
 } 
}

bool My_masscuts(double s[][12],TVar::Process process){

 double minZmassSqr=10*10;

 if(process==TVar::ZZ_4l){
   if(s[2][3]< minZmassSqr) return true;
   if(s[4][5]< minZmassSqr) return true;
 }
 return false;	 

}


bool My_smalls(double s[][12],int npart){

// Reject event if any s(i,j) is too small
// cutoff is defined in technical.Dat
	
      if ( 
       npart == 3 &&
       (
        (-s[5-1][1-1]< cutoff_.cutoff)  //gamma p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //gamma p2
     || (-s[4-1][1-1]< cutoff_.cutoff)  //e+    p1
     || (-s[4-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[3-1][1-1]< cutoff_.cutoff)  //nu    p1
     || (-s[3-1][2-1]< cutoff_.cutoff)  //nu    p2
     || (+s[5-1][4-1]< cutoff_.cutoff)  //gamma e+
     || (+s[5-1][3-1]< cutoff_.cutoff)  //gamma nu
     || (+s[4-1][3-1]< cutoff_.cutoff)  //e+    nu
	)	 
      ) 
        return true;
     
     else if (
       npart == 4 &&     
      (
        (-s[5-1][1-1]< cutoff_.cutoff)  //e-    p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[6-1][1-1]< cutoff_.cutoff)  //nb    p1
     || (-s[6-1][2-1]< cutoff_.cutoff)  //nb    p2
     || (+s[6-1][5-1]< cutoff_.cutoff)  //e-    nb
       )

     )
       
      return true;
     
     return false;
}




//Make sure
// 1. tot Energy Sum < 2EBEAM
// 2. PartonEnergy Fraction minimum<x0,x1<1
// 3. number of final state particle is defined
//
double SumMatrixElementPDF(TVar::Process process, mcfm_event_type* mcfm_event,double flavor_msq[][nmsq],double* flux){

  int NPart=npart_.npart+2;
  double p4[4][12];
  double fx1[nmsq];
  double fx2[nmsq];
  double msq[nmsq][nmsq];
  
  
  //Parton Density Function is always evalualted at pT=0 frame
  //Make sure parton Level Energy fraction is [0,1]
  //phase space function already makes sure the parton energy fraction between [min,1]
  //  x0 EBeam =>   <= -x1 EBeam
  
  double sysPz=mcfm_event->p[0].Pz()    +mcfm_event->p[1].Pz();
  double sysE =mcfm_event->p[0].Energy()+mcfm_event->p[1].Energy();
  
  //Ignore the Pt doesn't make significant effect
  //double sysPt_sqr=sysPx*sysPx+sysPy*sysPy;
  //if(sysPt_sqr>=1.0E-10)  sysE=TMath::Sqrt(sysE*sysE-sysPt_sqr);
  
  double xx[2]={(sysE+sysPz)/EBEAM/2,(sysE-sysPz)/EBEAM/2};
  if(xx[0] > 1.0 || xx[0]<=xmin_.xmin) return 0.0;
  if(xx[1] > 1.0 || xx[1]<=xmin_.xmin) return 0.0;
  
  //Convert TLorentzVector into 4x12 Matrix
  //reverse sign of incident partons back
  for(int ipar=0;ipar<2;ipar++){    
    if(mcfm_event->p[ipar].Energy()>0){
      p4[0][ipar] = -mcfm_event->p[ipar].Px();
      p4[1][ipar] = -mcfm_event->p[ipar].Py();
      p4[2][ipar] = -mcfm_event->p[ipar].Pz();
      p4[3][ipar] = -mcfm_event->p[ipar].Energy();
    }
  }
  //initialize decayed particles
  for(int ipar=2;ipar<NPart;ipar++){
    
    p4[0][ipar] = mcfm_event->p[ipar].Px();
    p4[1][ipar] = mcfm_event->p[ipar].Py();
    p4[2][ipar] = mcfm_event->p[ipar].Pz();
    p4[3][ipar] = mcfm_event->p[ipar].Energy();
    
  }
  
  //calculate invariant masses between partons/final state particles
  double s[12][12];
  for(int jdx=0;jdx< NPart ;jdx++){
    s[jdx][jdx]=0;
    for(int kdx=jdx+1;kdx<NPart;kdx++){
      s[jdx][kdx]=2*(p4[3][jdx]*p4[3][kdx]-p4[2][jdx]*p4[2][kdx]-p4[1][jdx]*p4[1][kdx]-p4[0][jdx]*p4[0][kdx]);
      s[kdx][jdx]=s[jdx][kdx];
    }
  }
  
  
  //remove events has small invariant mass
  if(My_masscuts(s,process)) return 0.0;
  if(My_smalls(s,npart_.npart)) return 0.0;
  
  
  //Calculate Pdf
  //Always pass address through fortran function
  fdist_ (&density_.ih1, &xx[0], &scale_.scale, fx1); 
  fdist_ (&density_.ih2, &xx[1], &scale_.scale, fx2); 
  
  if( process==TVar::ZZ_4l)      qqb_zz_  (p4[0],msq[0]);
  if( process==TVar::HZZ_4l)     qqb_hzz_ (p4[0],msq[0]);
  
  double msqjk=0;
  for(int ii=0;ii<nmsq;ii++){
    for(int jj=0;jj<nmsq;jj++){
      
      //2-D matrix is reversed in fortran
      // msq[ parton2 ] [ parton1 ]
      //      flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];

      flavor_msq[jj][ii] = msq[jj][ii];

      msqjk+=flavor_msq[jj][ii];
    }//ii
  }//jj
  
  (*flux)=fbGeV2/(8*xx[0]*xx[1]*EBEAM*EBEAM);
  
  if(msqjk != msqjk || flux!=flux ){
    cout << "SumMatrixPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << " flux="<< flux <<endl;
    msqjk=0;
    flux=0;
  }
  return msqjk;
  
}

double HiggsWidth(double mass){

     double width=0.004;
     if (mass==125.) width=0.41650E-02; // obtained by running MCFM standalone
     
 return width;


}
