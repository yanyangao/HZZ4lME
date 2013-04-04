//
// Run by root -l -b MELArocCurve.C
// 

//#include "../interface/TVar.hh"
#include "TVar.hh"


// this flag switches between the samples from 
// set to true if the files are from CJLST ZZMatrixElement/MELA 
// set to false if they are from YGao/HZZ4lME/ 

bool fromMELAPackage = false;

TGraph* makeROCcurve(char* drawVar, const int bins, float start, float end,
		     int lineColor, int lineStyle, int lineWidth, int flavor,
		     TString sigFile, TString bkgFile, TString testName
		     ){

  TString treeName("newTree");
  
  if ( fromMELAPackage ) treeName = "SelectedTree";

  TChain* SMHtree = new TChain(treeName);
  SMHtree->Add(sigFile.Data());

  TChain* PStree = new TChain(treeName);
  PStree->Add(bkgFile.Data());
  
  TH1F *SMHhisto, *PShisto;
  
  char drawString[150];
  
  TString cutstring;

  if ( flavor < 4 ) 
    cutstring = Form("flavortype==%i", flavor);
  else 
    cutstring = "";

  sprintf(drawString,"%s>>SMHhisto(%i,%f,%f)",drawVar,bins,start,end);
  // std::cout << drawString << "\n";
  SMHtree->Draw(drawString,cutstring);
  sprintf(drawString,"%s>>PShisto(%i,%f,%f)",drawVar,bins,start,end);
  // std::cout << drawString << "\n";
  PStree->Draw(drawString,cutstring);
  
  SMHhisto = (TH1F*) gDirectory->Get("SMHhisto");
  SMHhisto->Scale(1/SMHhisto->Integral(0,bins+1));
  PShisto = (TH1F*) gDirectory->Get("PShisto");
  PShisto->Scale(1/PShisto->Integral(0,bins+1));

  double effSMH[bins],effPS[bins];

  for(int i=0; i<bins; i++){

    effSMH[i] = SMHhisto->Integral(i+1, bins);
    effPS[i] = PShisto->Integral(i+1, bins);

  }

  TGraph* ROC = new TGraph(bins,effSMH,effPS);
  ROC->SetLineColor(lineColor);
  ROC->SetLineStyle(lineStyle);
  ROC->SetLineWidth(lineWidth);
  ROC->GetXaxis()->SetTitle("#epsilon_{0^{+}}");
  ROC->GetYaxis()->SetTitle(Form("#epsilon_{%s}", testName.Data()));
  
  delete SMHtree;
  delete PStree;
  
  return ROC;

}


void MELArocCurve(){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  tdrstyle();

  TString fileAppendix = "ME";
  
  if ( fromMELAPackage ) 
    fileAppendix = "withProbabilities_ABtest";
  
  TString sigFile(Form("SMHiggs_comb_withDiscriminants_%s.root", fileAppendix.Data()));
  
  std::vector<TString> bkgFiles;
  std::vector<int> tests;

  // 2m+ gg
  bkgFiles.push_back(Form("minGrav_comb_withDiscriminants_%s.root", fileAppendix.Data()));
  tests.push_back(TVar::TZZ_4l);

  // 2m+ qq
  bkgFiles.push_back(Form("minGravqq_comb_withDiscriminants_%s.root", fileAppendix.Data()));
  tests.push_back(TVar::QQB_TZZ_4l);

  // 1-
  bkgFiles.push_back(Form("1minus_8-126_comb_%s.root", fileAppendix.Data()));
  tests.push_back(TVar::VZZ_4l);

  // 1+
  bkgFiles.push_back(Form("1plus_8-126_comb_%s.root", fileAppendix.Data()));
  tests.push_back(TVar::AVZZ_4l);

  // 2h- gg
  bkgFiles.push_back(Form("2hminus_8_126_comb_%s.root", fileAppendix.Data()));
  tests.push_back(TVar::PTZZ_2hminus_4l);

  // 2h+ gg
  bkgFiles.push_back(Form("2hplus_8_126_comb_%s.root", fileAppendix.Data()));
  tests.push_back(TVar::TZZ_2hplus_4l);

  // 2b+ gg
  bkgFiles.push_back(Form("2bplus_8_126_comb_%s.root", fileAppendix.Data()));
  tests.push_back(TVar::TZZ_2bplus_4l);


  for (int test = 0 ; test < bkgFiles.size(); test ++ ) {

    bool drawdecay = false;
    
    if (    tests.at(test) == TVar::TZZ_4l 
	 || tests.at(test) == TVar::QQB_TZZ_4l 
	 || tests.at(test) == TVar::VZZ_4l 
	 || tests.at(test) == TVar::AVZZ_4l 	 ) 
      drawdecay = true;

    TString bkgFile = bkgFiles.at(test);
    TString testName = TVar::ProcessName(tests.at(test));
    
    // default for ggX-> 2m+
    
    TString melaVar = "Psmh/(Psmh+Ptwomplus_gg*2e-07)";
    TString jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_TZZ_JHU)";
    TString melaDecayVar = "Psmh/(Psmh+Ptwomplus_decay*1e-08)";
    TString jhugenDecayVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_TZZ_DECAY_JHU)";
    
    if ( fromMELAPackage )  {
      melaVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p2_mela_NEW*2e-07)";
      jhugenVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p2_VAJHU_NEW)";
      melaDecayVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p2_decay_mela_NEW*1e-08)";
      jhugenDecayVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p2_decay_VAJHU_NEW)";
    }
    
    // qqbar->2m+
    if ( tests.at(test) == TVar::QQB_TZZ_4l ) {
      melaVar = "Psmh/(Psmh+Ptwomplus_qq*2e-07)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_QQB_TZZ_JHU)";
      if ( fromMELAPackage ) {
	melaVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p2qqb_mela_NEW*2e-07)";
	jhugenVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p2qqb_VAJHU_NEW)";
      }
    }

    // qqbar-> 1-
    if (tests.at(test) == TVar::VZZ_4l ) {
      melaVar = "Psmh/(Psmh+Poneminus*5e-03)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_VZZ_JHU)";
      melaDecayVar = "Psmh/(Psmh+Poneminus_decay*3.7e-04)";
      jhugenDecayVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_VZZ_DECAY_JHU)";
      if ( fromMELAPackage ) {
	melaVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p1_mela_NEW*5e-03)";
	jhugenVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p1_VAJHU_NEW)";
	melaDecayVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p1_decay_mela_NEW*3.7e-04)";
	jhugenDecayVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p1_decay_VAJHU_NEW)";
      }
    }

    // qqbar-> 1+
    if (tests.at(test) == TVar::AVZZ_4l ) {
      melaVar = "Psmh/(Psmh+Poneplus*5e-03)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_AVZZ_JHU)";
      melaDecayVar = "Psmh/(Psmh+Poneplus_decay*3.7e-04)";
      jhugenDecayVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_AVZZ_DECAY_JHU)";
      
      if ( fromMELAPackage ) {
	melaVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p1plus_mela_NEW*5e-03)";
	jhugenVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p1plus_VAJHU_NEW)";
	melaDecayVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p1plus_decay_mela_NEW*3.7e-04)";
	jhugenDecayVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p1plus_decay_VAJHU_NEW)";
      }
    }

    // gg->2h-
    if (tests.at(test) == TVar::PTZZ_2hminus_4l ) {
      melaVar = "Psmh/(Psmh+Ptwohminus*0.4)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_PTZZ_2hminus_JHU)";
      
      if ( fromMELAPackage ) {
	melaVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p2hminus_mela_NEW*0.4)";
	jhugenVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p2hminus_VAJHU_NEW)";
      }
    }
    
    // gg->2h+
    if (tests.at(test) == TVar::TZZ_2hplus_4l ) {
      melaVar = "Psmh/(Psmh+Ptwohplus*0.4)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_TZZ_2hplus_JHU)";
      if ( fromMELAPackage ) {
	melaVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p2hplus_mela_NEW*0.4)";
	jhugenVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p2hplus_VAJHU_NEW)";
      }
    }

    // gg->2b+
    if (tests.at(test) == TVar::TZZ_2bplus_4l ) {
      melaVar = "Psmh/(Psmh+Ptwobplus*2.3e-07)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_TZZ_2bplus_JHU)";
      if ( fromMELAPackage ) {
	melaVar = "p0plus_mela_NEW/(p0plus_mela_NEW+p2bplus_mela_NEW*2.3e-07)";
	jhugenVar = "p0plus_VAJHU_NEW/(p0plus_VAJHU_NEW+p2bplus_VAJHU_NEW)";
      }
    }
    
    for ( int flavor = 1; flavor < 5; flavor ++ ) {

      if ( flavor != 3 ) continue;
      
      int nbins = 200;
      
      TGraph* MELA = makeROCcurve(melaVar,nbins,0,1,4,1,3,flavor,sigFile,bkgFile,testName);
      TGraph* JHUgen = makeROCcurve(jhugenVar,nbins,0,1,1,2,3,flavor,sigFile,bkgFile,testName);
      TGraph* MELA_decay; 
      TGraph* JHUgen_decay;
      
      if ( drawdecay) {
	MELA_decay = makeROCcurve(melaDecayVar,nbins,0,1,2,1,3,flavor,sigFile,bkgFile,testName);
	JHUgen_decay = makeROCcurve(jhugenDecayVar,nbins,0,1,6,2,3,flavor,sigFile,bkgFile,testName);
      }
      MELA->Draw("AC");
      if ( drawdecay ) 
	MELA_decay->Draw("sameC");
      JHUgen->Draw("sameC");
      if ( drawdecay ) 
	JHUgen_decay->Draw("sameC");
      
      TLegend* leg = new TLegend(.2,.6,.5,.9);
      leg->SetFillColor(0);
      leg->AddEntry(MELA,"analytical(prod+decay)","l");
      if ( drawdecay ) 
	leg->AddEntry(MELA_decay,"analytical(decay)","l");
      leg->AddEntry(JHUgen,"JHUgen(prod+decay)","l");
      if ( drawdecay ) 
	leg->AddEntry(JHUgen_decay,"JHUgen(decay)","l");
      
      leg->Draw();

      TString flavortype = "all";
      if (flavor == 1 ) flavortype = "4e";
      if (flavor == 2 ) flavortype = "4mu";
      if (flavor == 3 ) flavortype = "2e2mu";
      
      c1->SaveAs(Form("plots/ROCcurve_analytical_vs_JHUgen_%s_%s.eps", testName.Data(), flavortype.Data()));
      c1->SaveAs(Form("plots/ROCcurve_analytical_vs_JHUgen_%s_%s.png", testName.Data(), flavortype.Data()));
      
      
      TFile *file = new TFile(Form("plots/ROCcurve_JHUgen_%s_%s.root", testName.Data(), flavortype.Data()), "RECREATE");
      file->cd();
      
      JHUgen->SetName(Form("JHUGen_%s_%s", testName.Data(), flavortype.Data()));
      JHUgen->Write();
      MELA->SetName(Form("MELA_%s_%s", testName.Data(), flavortype.Data()));
      MELA->Write();
      if ( drawdecay ) {
	JHUgen_decay->SetName(Form("JHUGen_decay_%s_%s", testName.Data(), flavortype.Data()));
	JHUgen_decay->Write();
	MELA_decay->SetName(Form("MELA_decay_%s_%s", testName.Data(), flavortype.Data()));
	MELA_decay->Write();
      }
      file->Close();
    }
  }
}
