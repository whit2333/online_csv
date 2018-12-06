#include <iostream>
#include <cmath>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "TCanvas.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TStyle.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libGenVector.so)

#include "THStack.h"
#include "TBufferJSON.h"

#include "nlohmann/json.hpp"

#include "THcParmList.h"
R__LOAD_LIBRARY(libPodd.so)
R__LOAD_LIBRARY(libHallA.so)
R__LOAD_LIBRARY(libdc.so)
R__LOAD_LIBRARY(libHallC.so)

// fmt - string formatting library
#include "fmt/core.h"
#include "fmt/ostream.h"
R__LOAD_LIBRARY(libfmt.so)

using Pvec3D   = ROOT::Math::XYZVector;
using Pvec4D   = ROOT::Math::PxPyPzMVector;

// VecOps::RVec is like std::vector with some extra bells and whistles 
using inters   = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;

void hms_s1_times(int RunNumber = 5893, const char* codatype = "COIN", int nevents = 100000) {

  std::string coda_type = codatype;
  //std::string coda_type = "SHMS";
  //std::string coda_type = "COIN";

  //std::string rootfile  = "/net/cdaq/cdaql3data/cdaq/sidis-fall-18/ROOTfiles/";
  std::string rootfile    = std::string("ROOTfiles/");
  rootfile    += std::string("coin_replay_production_");

  if (coda_type == "SHMS") {
    rootfile = std::string(
      "ROOTfiles/shms_replay_production_all_");
    // + std::to_string(RunNumber) + "_50000.root";
  }
  rootfile += std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";

  std::string run_list_json  = "DBASE/run_list.json";

  auto file = new TFile(rootfile.c_str());
  //new TBrowser;

  //auto ptheta_lab = hc_parms->Find("ptheta_lab");
  //ptheta_lab->Print();
  //auto hpcentral = hc_parms->Find("hpcentral");
  //hpcentral->Print();
  //auto ppcentral = hc_parms->Find("ppcentral");
  //ppcentral->Print();
  //std::cout << j.dump()  << "\n";
  //double P0_SHMS = hcana::json::FindVarValueOr(hc_parms,"ppcentral",0.0) ;
  //std::cout <<  "SHMS P0 = " << P0_SHMS << "\n";
  //double THETA_SHMS = hcana::json::FindVarValueOr(hc_parms,"ptheta_lab",0.0) ;

  //nlohmann::json j2;
  //{
  //  ifstream  in_run_file(run_list_json);
  //  in_run_file >> j2;
  //}
  //if( j2.find(std::to_string(RunNumber)) == j2.end() ) {
  //}
  //j2[std::to_string(RunNumber)] = j[std::to_string(RunNumber)];

  //{
  //  //std::cout << j.dump(2) << "\n";
  //  // write prettified JSON to another file
  //  std::ofstream o(run_list_json);
  //  o << std::setw(4) << j2 << std::endl;
  //}
  //std::cout << j2.dump() << "\n";

  ROOT::EnableImplicitMT(4);

  Pvec4D  Pbeam(0,0,10.598,0.000511);

  ROOT::RDataFrame d("T",rootfile);
  ROOT::RDataFrame d_sh("TSH",rootfile);

  std::string hpdelta = "P.gtr.dp > -10 && P.gtr.dp < 20 && H.gtr.dp > -10 && H.gtr.dp < 10";
  std::string epiCut  = "P.aero.npeSum > 1.0 && P.cal.eprtracknorm < 0.2 "
                        "&& H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 "
                        "&& H.cal.etottracknorm < 2.0 && H.cal.eprtracknorm  > 0.2"; 


  //tt->SetBranchAddress("H.gtr.dp", &hdelta);   
  //tt->SetBranchAddress("P.gtr.dp", &pdelta);  
  //tt->SetBranchAddress("H.gtr.th", &HgtrTh);   
  //tt->SetBranchAddress("P.gtr.th", &PgtrTh);  
  //tt->SetBranchAddress("H.gtr.ph", &HgtrPh);   
  //tt->SetBranchAddress("P.gtr.ph", &PgtrPh);  
  const double SHMSpartMass       = 0.1395704; // pion mass in GeV/c^2
  const double HMSpartMass        = 0.000510998; // electron mass in GeV/c^2
  const double speedOfLight       = 29.9792; // in cm/ns
  const double pOffset            = 1.5 + 17. + 6.; //9.5 + 10;  // in ns
  const double hOffset            = 335;
  const double SHMScentralPathLen = 18.1*100;  // SHMS Target to focal plane path length converted to cm
  const double SHMSpathLength     = 18.1*100;      // For now assume that it's same as SHMScentralPathLen
  const double HMScentralPathLen  = 22*100;     // HMS Target to focal plane path length converted to cm
  const double SHMScoinCorr       = 0.0;
  const double HMScoinCorr        = 0.0;
  const double DeltaHMSpathLength = HMScentralPathLen;// For now assume that it's same as HMScentralPathLen 

  //auto d_timing = d
  //.Define("DeltaHMSpathLength",
  //        [](double th, double ph, double x ) {
  //          double DeltaHMSpathLength = 12.462*th + 0.1138*th*x - 0.0154*x - 72.292*th*th - 0.0000544*x*x - 116.52*ph*ph;
  //          return DeltaHMSpathLength;
  //        } ,{"P.gtr.th","H.gtr.ph","H.gtr.x"})
  //.Define("PgtrBetaCalc",[=](double p){return p/sqrt(p*p + SHMSpartMass*SHMSpartMass);} ,{"P.gtr.p"})
  //.Define("HgtrBetaCalc",[=](double p){return p/sqrt(p*p + HMSpartMass*HMSpartMass);} ,{"H.gtr.p"});
  //.Define("SHMScoinCorr",[=](double beta){
  //  double SHMScoinCorr = SHMScentralPathLen / (speedOfLight*beta) + (SHMSpathLength - SHMScentralPathLen) / (speedOfLight*beta) + (PhodoStartTimeMean - PhodfpHitsTime);                                                                   
  //  return p/sqrt(p*p + SHMSpartMass*SHMSpartMass);} ,{"PgtrBetaCalc"})
  //SHMScoinCorr = SHMScentralPathLen / (speedOfLight*PgtrBetaCalc) + (SHMSpathLength - SHMScentralPathLen) / (speedOfLight*PgtrBetaCalc) + (PhodoStartTimeMean - PhodfpHitsTime);                                                                   
  //HMScoinCorr = HMScentralPathLen / (speedOfLight*HgtrBetaCalc) + DeltaHMSpathLength / (speedOfLight*HgtrBetaCalc) + (HhodoStartTimeMean - HhodfpHitsTime);      

  //SHMScoinCorr = SHMScentralPathLen / (speedOfLight*PgtrBetaCalc) + (SHMSpathLength - SHMScentralPathLen) / (speedOfLight*PgtrBetaCalc) + (PhodoStartTimeMean - PhodfpHitsTime);                                                                   
  //HMScoinCorr = HMScentralPathLen / (speedOfLight*HgtrBetaCalc) + DeltaHMSpathLength / (speedOfLight*HgtrBetaCalc) + (HhodoStartTimeMean - HhodfpHitsTime);      

  //SHMScorrCoinTimeROC1 = (TcoinpTRIG1_ROC1_tdcTimeRaw*0.1 - SHMScoinCorr) - (TcoinpTRIG4_ROC1_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset; // 0.1 to convert to ns 

  //SHMScorrCoinTimeROC2 = (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - SHMScoinCorr) - (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset; 

  //SHMScorrCoinTimeROC2 = 
  //(TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - 
  // SHMScoinCorr) - 
  //(TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - 
  // HMScoinCorr) - pOffset; 

  //ctimepi = SHMScorrCoinTimeROC2 ;


  //h_reactz_diff_cointtime->Fill(Preactz+Hreactz,SHMScorrCoinTimeROC2);
  //h1_epi_PcointimeROC2->Fill(SHMScorrCoinTimeROC2);    
  //h1_epi_PcointimeROC2_R->Fill(SHMScorrCoinTimeROC2);
  //h1_epi_PcointimeROC2_L->Fill(SHMScorrCoinTimeROC2);
  //h1_epi_PcointimeROC2_C->Fill(SHMScorrCoinTimeROC2);


  auto d0 = d
  .Define("hms_e_EoverP",
          [](doublers& EOverP) {
            return EOverP[EOverP > 0.5 && EOverP < 1.8];
          }, {"H.cal.etottracknorm"})
  .Define("hms_e_EoverP_nGood",
          [](doublers& EOverP) {
            return (int)EOverP.size();
          }, {"hms_e_EoverP"})
  .Define("shms_pi_preshower",
          [](doublers& EOverP) {
            return EOverP[EOverP < 0.2];
          }, {"P.cal.eprtracknorm"})
  .Define("shms_pi_preshower_nGood",
          [](doublers& EOverP) {
            return (int)EOverP.size();
          }, {"shms_pi_preshower"})
  ;

  // Find the coin peak 
  auto   h_coin_time = d0.Histo1D({"coin_time","coin_time",400,0,100},"CTime.ePiCoinTime_ROC2");
  int    coin_peak_bin = h_coin_time->GetMaximumBin();
  double coin_peak_center = h_coin_time->GetBinCenter(coin_peak_bin);
  std::cout << " coin peak at " << coin_peak_center << " ns\n";
    
  auto d2 = d0
  .Filter([=](double coin_time) {return std::abs(coin_time - coin_peak_center)< 1.5;},{"CTime.ePiCoinTime_ROC2"})
  .Filter([](double npe){return npe>1.0;},{"H.cer.npeSum"})
  .Filter(
    [](int n_good) {
        if( n_good >0 ) {
          return true;
        }
      return  false;
    },{"hms_e_EoverP_nGood"});
  auto d3 = d2
  .Filter(
    [](int n_good) {
        if( n_good >0 ) {
          return true;
        }
      return  false;
    },{"shms_pi_preshower_nGood"});


  auto d_beans  = d2
  .Filter(hpdelta)
  .Filter(epiCut);

  auto d2_2 = d2.Filter("fEvtHdr.fEvtType == 2");
  auto d2_4 = d2.Filter("fEvtHdr.fEvtType == 4");

  // Histograms lazy evaluated -> Nothing happens until they are used below (eg Draw() is called) 
  auto h_EOverP_0       = d0.Histo1D<doublers>({"hms_e_EoverP_0","SHMS total shower; SHMS E/P",    100,0.05,1.8},"P.cal.etottracknorm");
  auto h_EOverP_2       = d2.Histo1D<doublers>({"hms_e_EoverP_2","SHMS total shower; SHMS E/P",    100,0.05,1.8},"P.cal.etottracknorm");
  auto h_EOverP_3       = d3.Histo1D<doublers>({"hms_e_EoverP_3","SHMS total shower; SHMS E/P",    100,0.05,1.8},"P.cal.etottracknorm");

  auto h_EprOverP_0       = d0.Histo1D<doublers>({"hms_e_EoverP_0","SHMS pre-shower; SHMS E/P",    100,0.05,1.8},"P.cal.eprtracknorm");
  auto h_EprOverP_2       = d2.Histo1D<doublers>({"hms_e_EoverP_2","SHMS pre-shower; SHMS E/P",    100,0.05,1.8},"P.cal.eprtracknorm");
  auto h_EprOverP_3       = d3.Histo1D<doublers>({"hms_e_EoverP_3","SHMS pre-shower; SHMS E/P",    100,0.05,1.8},"P.cal.eprtracknorm");

  //auto h_EOverP_1       = d0.Histo1D<doublers>({"hms_e_EoverP_1","E/P cut; HMS E/P",100,0.05,1.8},"hms_e_EoverP");
  //auto h_EOverP_nGood_0 = d0.Histo1D({"h_EOverP_nGood_0","h_EOverP_nGood_0",10,0,10},"P.cal.ntracks");
  //auto h_EOverP_nGood_1 = d0.Histo1D({"h_EOverP_nGood_1","h_EOverP_nGood_0",10,0,10},"hms_e_EoverP_nGood");

  auto h_event_type = d0.Histo1D({"event_type","event_type",10,0,10},"fEvtHdr.fEvtType");

  
  auto total_charge = d_sh.Sum("H.BCM4A.scalerCurrent");

  auto hms_electron_counts = d2.Count();
  auto hms_electron_counts2 = d2_2.Count();
  auto hms_electron_counts4 = d2_4.Count();

  auto bean_count = d_beans.Count();
  //auto disp = d2.Display("fEvtHdr.fEvtType");

  //disp->Print();
  double hms_e_yield = (*hms_electron_counts)/(*total_charge);
  double hms_e_yield2 = (*hms_electron_counts2)/(*total_charge);
  double hms_e_yield4 = (*hms_electron_counts4)/(*total_charge);


  std::cout << " beans  : " << *bean_count << "\n";
  std::cout << " counts : " << *hms_electron_counts << "\n";
  std::cout << " charge : " << *total_charge << " uC\n";
  std::cout << " yield  : " << hms_e_yield << "\n";
  std::cout << " yield2 : " << hms_e_yield2 << "\n";
  std::cout << " yield4 : " << hms_e_yield4 << "\n";

  //auto s_dc_x_fp           = d2.Histo1D({"s_dc_x_fp ","xy fp; x",100,-50,50}, "P.dc.x_fp");
  //auto s_dc_y_fp           = d2.Histo1D({"s_dc_y_fp ","xy fp; y",100,-50,50}, "P.dc.y_fp");
  //auto s_dc_xp_fp          = d2.Histo1D({"s_dc_xp_fp","xy fp; xp",100,-50,50},"P.dc.xp_fp");
  //auto s_dc_yp_fp          = d2.Histo1D({"s_dc_yp_fp","xy fp; xp",100,-50,50},"P.dc.yp_fp");
  
  // -----------------------------------------------------------
  //
  TCanvas* c = nullptr;
  int b1 = 0;
  int b2 = 0;
  double hmax = 0.0;
  THStack*  hs = nullptr;
  //TLatex latex;

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);

  // This call starts the loop over the data.
  // A DrawCopy is used so that the histogram is not deleted at the end of 
  // scope, and thus stays visible on the canvas.
  h_EOverP_0->DrawCopy();

  c->cd(2);

  h_EOverP_2->DrawCopy();
  c->cd(3);
  h_event_type->DrawCopy();
  c->cd(4);
  h_coin_time->DrawCopy();
  gPad->BuildLegend();

  gSystem->mkdir("results/df_example", true);
  c->SaveAs((std::string("results/df_example/c1_")+std::to_string(RunNumber)+".pdf").c_str());
  c->SaveAs((std::string("results/df_example/c1_")+std::to_string(RunNumber)+".png").c_str());
  std::cout << "derp\n";
  
  // ---------------------------------------------------------
  //
  c = new TCanvas();
  hs = new THStack("SHMS_cal","SHMS calorimeter; E/P");

  h_EOverP_0->SetLineColor(1);
  h_EOverP_2->SetLineColor(4);
  h_EOverP_3->SetLineColor(2);
  h_EprOverP_0->SetLineColor(8);
  h_EprOverP_2->SetLineColor(6);
  h_EprOverP_3->SetLineColor(9);

  hs->Add((TH1*)h_EOverP_0->Clone());
  hs->Add((TH1*)h_EOverP_2->Clone());
  hs->Add((TH1*)h_EOverP_3->Clone());

  hs->Add((TH1*)h_EprOverP_0->Clone());
  hs->Add((TH1*)h_EprOverP_2->Clone());
  hs->Add((TH1*)h_EprOverP_3->Clone());

  hs->Draw("nostack");
  gPad->SetLogy(true);
  gPad->BuildLegend();

  c->SaveAs((std::string("results/df_example/c2_")+std::to_string(RunNumber)+".pdf").c_str());
  c->SaveAs((std::string("results/df_example/c2_")+std::to_string(RunNumber)+".png").c_str());

}
