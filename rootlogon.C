void rootlogon() {
  gSystem->SetBuildDir("$HOME/.root_build_dir");
  gPrintViaErrorHandler = kTRUE;
  gErrorIgnoreLevel = kWarning;
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gSystem->Load("libTree");
  gSystem->Load("libTreePlayer");
  gSystem->Load("libHist");
}
