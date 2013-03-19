//
// Run by root -l -b MELArocCurve.C
// very very messy...
// YY will clean it if she is in a good cleaning mood
// 

//char* fileName = "1plus";
//char* fileName = "1minus";
char* fileName = "minGrav";
//char* fileName = "minGravqq";

TGraph* makeROCcurve(char* drawVar="melaLD", 
		     const int bins=100, float start=0, float end=1,
		     int lineColor=1, int lineStyle=1, int lineWidth=2,	int flavor = 4, 
		     char* sigFile="SMHiggs_comb_withDiscriminants_ME.root",
		     char* bkgFile="minGrav_comb_withDiscriminants_ME.root",
		     // char* bkgFile="minGravqq_comb_withDiscriminants_ME.root",
		     // char* bkgFile="minGrav_comb_withDiscriminants_ME.root",
		     // char* sigFile="fromAndrew/SMHiggs_comb_withDiscriminants_withDiscriminants_ME.root",
		     // char* bkgFile="fromAndrew/minGravqq_comb_withDiscriminants_withDiscriminants_ME.root",
		     //char* sigFile="wresolution/SMHiggs_126GeV_comb_w4muResolution_withDiscriminants_ME.root",
		     //char* bkgFile="wresolution/qqMinGrav_126GeV_comb_w4muResolution_withDiscriminants_ME.root",
		     //char* bkgFile="fromAndrew/minGrav_comb_withDiscriminants_withDiscriminants_ME.root",
		     // char* bkgFile="1minus_8-126_ME.root",
		     //char* bkgFile="1plus_8-126_ME.root",
		     // char* sigFile="fromYY/SMHiggs_8-126_withDiscriminants_ME.root",
		     // char* bkgFile="fromYY/minGrav_8-126_withDiscriminants_ME.root",
		     // char* bkgFile="fromYY/minGravqq_8-126_withDiscriminants_ME.root",
		     ){

  TChain* SMHtree = new TChain("newTree");
  SMHtree->Add(sigFile);

  TChain* PStree = new TChain("newTree");
  PStree->Add(bkgFile);

  // std::cout << bkgFile<< "\n";
  
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
  //   ROC->GetYaxis()->SetTitle("qq#rightarrow#epsilon_{2^{+}_{min}}");
  ROC->GetYaxis()->SetTitle(Form("#epsilon_{%s}", fileName));
  delete SMHtree;
  delete PStree;

  return ROC;

}


void MELArocCurve(){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  tdrstyle();

  for ( int flavor = 1; flavor < 5; flavor ++ ) {
    
    // if ( flavor != 3) continue;
    
    //
    // for spin 2 
    // 

    // gg->X
    // TGraph* MELA_check = makeROCcurve("spinTwoMinimalMelaLD",500,0,1,4,1,3,flavor);
    TGraph* MELA_check = makeROCcurve("Psmh/(2e-7*Ptwomplus_gg+Psmh)",500,0,1,4,1,3,flavor);
    TGraph* JHUgen = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+0.6*dXsec_TZZ_JHU)",500,0,1,1,2,3,flavor);
    // TGraph* MELA_check_decay = makeROCcurve("spinTwoDecayMinimalMelaLD",500,0,1,2,1,3,flavor);
    TGraph* MELA_check_decay = makeROCcurve("Psmh/(Psmh+Ptwomplus_decay*1e-8)",500,0,1,2,1,3,flavor);
    //TGraph* MELA_check_decay = makeROCcurve("spinTwoDecay",500,0,1,2,1,3,flavor);
    TGraph* JHUgen_decay = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+1e11*dXsec_TZZ_DECAY_JHU)",500,0,1,6,2,3,flavor);


    /*
    // qqbar->X 
    // TGraph* MELA_check = makeROCcurve("spinTwoqqMinimalMelaLD",500,0,1,4,1,3,flavor);
    TGraph* MELA_check = makeROCcurve("Psmh/(2e-7*Ptwomplus_qq+Psmh)",500,0,1,4,1,3,flavor);
    // TGraph* MELA_check = makeROCcurve("spinTwoQQ",500,0,1,4,1,3,flavor);
    TGraph* JHUgen = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+20*dXsec_QQB_TZZ_JHU)",500,0,1,1,2,3,flavor);
    // TGraph* MELA_check_decay = makeROCcurve("spinTwoDecayMinimalMelaLD",1000,0,0.15,2,1,3,flavor);
    // TGraph* MELA_check_decay = makeROCcurve("spinTwoDecay",500,0,0.2,2,1,3,flavor);
    TGraph* MELA_check_decay = makeROCcurve("Psmh/(Psmh+Ptwomplus_decay*1e-8)",500,0,1,2,1,3,flavor);
    TGraph* JHUgen_decay = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+1e10*dXsec_TZZ_DECAY_JHU)",1000,0,0.6,6,2,3,flavor);
    */

    MELA_check->Draw("AC");
    JHUgen->Draw("sameC");
    MELA_check_decay->Draw("sameC");
    JHUgen_decay->Draw("sameC");

    
    /*
    //
    // for spin 1
    //

    // 1-
    TGraph* MELA_check = makeROCcurve("Psmh/(Psmh+0.005*Poneminus)",100,0,1,4,1,3,flavor);
    TGraph* MELA_check_decay = makeROCcurve("Psmh/(Psmh+0.0005*Poneminus_decay)",100,0,1,1,1,3,flavor);
    TGraph* JHUgen = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+10*dXsec_VZZ_JHU)",100,0,1,2,2,3,flavor);
    TGraph* JHUgen_decay = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+1e10*dXsec_VZZ_DECAY_JHU)",100,0,1,6,2,3,flavor);

    // 1+
    TGraph* MELA_check = makeROCcurve("Psmh/(Psmh+0.005*Poneplus)",100,0,1,4,1,3,flavor);
    TGraph* MELA_check_decay = makeROCcurve("Psmh/(Psmh+0.0005*Poneplus_decay)",100,0,1,1,1,3,flavor);
    TGraph* JHUgen = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+10*dXsec_AVZZ_JHU)",100,0,1,2,2,3,flavor);
    TGraph* JHUgen_decay = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+1e10*dXsec_AVZZ_DECAY_JHU)",100,0,1,6,2,3,flavor);

    MELA_check->Draw("AC");
    MELA_check_decay->Draw("sameC");
    JHUgen->Draw("sameC");
    JHUgen_decay->Draw("sameC");
    */    
    

    TLegend* leg = new TLegend(.2,.6,.5,.9);
    leg->SetFillColor(0);
    leg->AddEntry(MELA_check,"analytical(prod+decay)","l");
    leg->AddEntry(MELA_check_decay,"analytical(decay)","l");
    leg->AddEntry(JHUgen,"JHUgen(prod+decay)","l");
    leg->AddEntry(JHUgen_decay,"JHUgen(decay)","l");
    
    leg->Draw();

    TString flavortype = "all";
    if (flavor == 1 ) flavortype = "4e";
    if (flavor == 2 ) flavortype = "4mu";
    if (flavor == 3 ) flavortype = "2e2mu";
    
    c1->SaveAs(Form("ROCcurve_analytical_vs_JHUgen_%s_%s.eps", fileName, flavortype.Data()));
    c1->SaveAs(Form("ROCcurve_analytical_vs_JHUgen_%s_%s.png", fileName, flavortype.Data()));

    
    TFile *file = new TFile(Form("ROCcurve_JHUgen_%s_%s.root", fileName, flavortype.Data()), "RECREATE");
    file->cd();
    
    JHUgen_decay->SetName(Form("JHUGen_%s_%s", fileName, flavortype.Data()));
    JHUgen_decay->Write();
    MELA_check_decay->SetName(Form("MELA_%s_%s", fileName, flavortype.Data()));
    MELA_check_decay->Write();
    
    // JHUgen->SetName(Form("JHUGen_%s_%s", fileName, flavortype.Data()));
    // JHUgen->Write();
    // MELA_check->SetName(Form("MELA_%s_%s", fileName, flavortype.Data()));
    // MELA_check->Write();
    file->Close();

  }
}
