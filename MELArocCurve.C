




TGraph* makeROCcurve(char* drawVar="melaLD", 
		     const int bins=100, float start=0, float end=1,
		     int lineColor=1, int lineStyle=1, int lineWidth=2,
		     char* sigFile="SMHiggsZZ_250_JHU_ME.root",
		     char* bkgFile="TZZ_2mplus_250_JHU_ME.root"
		     ){

  TChain* SMHtree = new TChain("angles");
  SMHtree->Add(sigFile);

  TChain* PStree = new TChain("angles");
  PStree->Add(bkgFile);
  
  TH1F *SMHhisto, *PShisto;
  
  char drawString[150];

  sprintf(drawString,"%s>>SMHhisto(%i,%f,%f)",drawVar,bins,start,end);
  SMHtree->Draw(drawString,"");
  sprintf(drawString,"%s>>PShisto(%i,%f,%f)",drawVar,bins,start,end);
  PStree->Draw(drawString,"");
  
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
  ROC->GetYaxis()->SetTitle("#epsilon_{2^{+}_{min}}");
  delete SMHtree;
  delete PStree;

  return ROC;

}


void MELArocCurve(){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  setTDRStyle();

  TGraph* MELA_check = makeROCcurve("Psmh_new/(Psmh_new+1e-7*Pgrav_new)",1000,0,1,  4,2);
  TGraph* JHUgen = makeROCcurve("dXsec_HZZ_JHU/(dXsec_HZZ_JHU+dXsec_TZZ_JHU)",1000,-10,10,                    6,1);

  MELA_check->Draw("AC");
  JHUgen->Draw("sameC");

  TLegend* leg = new TLegend(.2,.6,.5,.9);
  leg->SetFillColor(0);
  
  leg->AddEntry(MELA_check,"analytic","l");
  leg->AddEntry(JHUgen,"JHUgen","l");

  leg->Draw();

  c1->SaveAs("ROCcurve_MELA_vs_JHUgen.eps");
  c1->SaveAs("ROCcurve_MELA_vs_JHUgen.png");

}
