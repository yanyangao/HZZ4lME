

void overlayroc(){

  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("tdrstyle();");
  gROOT->ForceStyle();

  TString flavor="2e2mu";

  TFile *file1 = new TFile(Form("ROCcurve_JHUgen_minGrav_%s.root", flavor.Data()));
  gROOT->cd();
  TGraph *g1 = (TGraph*)file1->Get(Form("MELA_minGrav_%s", flavor.Data()));
  g1->SetLineColor(kBlue);
  g1->SetLineWidth(2);
  TGraph *g2 = (TGraph*)file1->Get(Form("JHUGen_minGrav_%s", flavor.Data()));
  g2->SetLineColor(kRed);
  g2->SetLineWidth(2);

  TFile *file2 = new TFile(Form("ROCcurve_JHUgen_minGravqq_%s.root", flavor.Data()));
  gROOT->cd();

  TGraph *g3 = (TGraph*)file2->Get(Form("MELA_minGravqq_%s", flavor.Data()));
  g3->SetLineColor(kMagenta);
  g3->SetLineWidth(2);

  TGraph *g4 = (TGraph*)file2->Get(Form("JHUGen_minGravqq_%s", flavor.Data()));
  g4->SetLineColor(kBlack);
  g4->SetLineWidth(2);


  TLegend* leg = new TLegend(.2,.6,.5,.9);
  leg->SetFillColor(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(g1,"Analytical(ggX)","l");
  leg->AddEntry(g2,"JHUGen(ggX)","l");
  leg->AddEntry(g3,"Analytical(qqbarX)","l");
  leg->AddEntry(g4,"JHUGen(qqbarX)","l");

  TCanvas *c1 = new TCanvas();
  g1->Draw("AC");
  g2->Draw("same");
  g3->Draw("same");
  g4->Draw("same");
  leg->Draw("same");
  c1->SaveAs(Form("ROCcurve_analytical_vs_JHUgen_decayonly_%s.png", flavor.Data()));
  c1->SaveAs(Form("ROCcurve_analytical_vs_JHUgen_decayonly_%s.eps", flavor.Data()));
    		

  

  
}
