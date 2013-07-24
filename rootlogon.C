{
  gSystem->SetIncludePath("-I$ROOTSYS/include -I../Higgs/Higgs_CS_and_Width/include");
  gSystem->Load("libgfortran.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("./libmcfm_6p6.so");
  gSystem->Load("./libME.so");
  gSystem->AddIncludePath("-I$ROOFITSYS/include/");
  gROOT->ProcessLine(".L PDFs/RooSpinTwo_7D.cc+");
  gROOT->ProcessLine(".L PDFs/RooXZsZs_5D.cc+");
  gROOT->ProcessLine(".L PDFs/RooSpinOne_7D.cc+");
  gROOT->ProcessLine(".L PDFs/RooSpinOne_Decay.cc+");
  gROOT->ProcessLine(".L PDFs/RooSpinTwo_Decay.cc+");
  gROOT->ProcessLine(".L PDFs/RooSpinZero_7DComplex.cc+");
}
