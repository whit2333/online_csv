#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include "TCanvas.h"
int readsimroot(){
  std::string rootfile = std::string("/u/home/shuojia/simc_gfortran/worksim/csvelas7822.root");
  ROOT::RDataFrame d("h666",rootfile);
  auto h_Emiss = d.Histo1D({"hmsgtrdp","emiss",100,-1,1},"Em");

  auto *c = new TCanvas;
  h_Emiss->DrawCopy();

  return 0;
}
