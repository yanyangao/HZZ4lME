{
  gSystem->SetIncludePath("-I$ROOTSYS/include -I../Higgs/Higgs_CS_and_Width/include");
  gSystem->Load("libgfortran.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("./libmcfm.so");
  gSystem->Load("./libME.so");
}
