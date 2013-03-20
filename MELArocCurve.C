//
// Run by root -l -b MELArocCurve.C
// 

#include "TVar.hh"

TGraph* makeROCcurve(char* drawVar, const int bins, float start, float end,
		     int lineColor, int lineStyle, int lineWidth, int flavor,
		     TString sigFile, TString bkgFile, TString testName
		     ){

  TChain* SMHtree = new TChain("newTree");
  SMHtree->Add(sigFile.Data());

  TChain* PStree = new TChain("newTree");
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
  
  TString sigFile = "SMHiggs_comb_withDiscriminants_ME.root";
  
  std::vector<TString> bkgFiles;
  std::vector<int> tests;
  
  // 2m+ gg
  bkgFiles.push_back("minGrav_comb_withDiscriminants_ME.root");
  tests.push_back(TVar::TZZ_4l);

  // 2m+ qq
  bkgFiles.push_back("minGravqq_comb_withDiscriminants_ME.root");
  tests.push_back(TVar::QQB_TZZ_4l);

  // 1-
  bkgFiles.push_back("1minus_8-126_ME.root");
  tests.push_back(TVar::VZZ_4l);

  // 1+
  bkgFiles.push_back("1plus_8-126_ME.root");
  tests.push_back(TVar::AVZZ_4l);
  
  for (int test = 0 ; test < bkgFiles.size(); test ++ ) {

    TString bkgFile = bkgFiles.at(test);
    TString testName = TVar::ProcessName(tests.at(test));
    
    TString melaVar = "Psmh/(Psmh+2e-7*Ptwomplus_gg)";
    TString jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_TZZ_JHU)";
    TString melaDecayVar = "Psmh/(Psmh+Ptwomplus_decay)";
    TString jhugenDecayVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_TZZ_DECAY_JHU)";
    
    
    if ( tests.at(test) == TVar::QQB_TZZ_4l ) {
      melaVar = "Psmh/(Psmh+2e-07*Ptwomplus_qq)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_QQB_TZZ_JHU)";
    }


    if (tests.at(test) == TVar::VZZ_4l ) {
      melaVar = "Psmh/(Psmh+5e-03*Poneminus)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_VZZ_JHU)";
      melaDecayVar = "Psmh/(Psmh+Poneminus_decay)";
      jhugenDecayVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_VZZ_DECAY_JHU)";
    }

    if (tests.at(test) == TVar::AVZZ_4l ) {
      melaVar = "Psmh/(Psmh+5e-03*Poneplus)";
      jhugenVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_AVZZ_JHU)";
      melaDecayVar = "Psmh/(Psmh+Poneplus_decay)";
      jhugenDecayVar = "dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_AVZZ_DECAY_JHU)";
    }
    
    for ( int flavor = 1; flavor < 5; flavor ++ ) {

      TGraph* MELA = makeROCcurve(melaVar,500,0,1,4,1,3,flavor,sigFile,bkgFile,testName);
      TGraph* JHUgen = makeROCcurve(jhugenVar,500,0,1,1,2,3,flavor,sigFile,bkgFile,testName);
      TGraph* MELA_decay = makeROCcurve(melaDecayVar,500,0,1,2,1,3,flavor,sigFile,bkgFile,testName);
      TGraph* JHUgen_decay = makeROCcurve(jhugenDecayVar,500,0,1,6,2,3,flavor,sigFile,bkgFile,testName);
      
      MELA->Draw("AC");
      MELA_decay->Draw("sameC");
      JHUgen->Draw("sameC");
      JHUgen_decay->Draw("sameC");
      
      TLegend* leg = new TLegend(.2,.6,.5,.9);
      leg->SetFillColor(0);
      leg->AddEntry(MELA,"analytical(prod+decay)","l");
      leg->AddEntry(MELA_decay,"analytical(decay)","l");
      leg->AddEntry(JHUgen,"JHUgen(prod+decay)","l");
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
      
      JHUgen_decay->SetName(Form("JHUGen_decay_%s_%s", testName.Data(), flavortype.Data()));
      JHUgen_decay->Write();
      MELA_decay->SetName(Form("MELA_decay_%s_%s", testName.Data(), flavortype.Data()));
      MELA_decay->Write();

      file->Close();
    }
  }
}
