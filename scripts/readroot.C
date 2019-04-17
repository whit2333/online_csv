#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include "TCanvas.h"
int readroot(){
 std::string rootfile = std::string("ROOTfiles/coin_replay_production_7593_-1.root");
 ROOT::RDataFrame d("T",rootfile);
  auto h_P_etottracknorm = d.Histo1D({"hmsgtrdp","dp",100,0.1,2},"H.gtr.dp");
  //auto all = h_P_etottracknorm->Integral();
  auto *c = new TCanvas;
  std::cout<<"check"<<"\n";
 h_P_etottracknorm->DrawCopy("E");

  return 0;
}
