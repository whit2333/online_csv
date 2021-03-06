#include "nlohmann/json.hpp"
#include <cmath>
#include <iostream>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libGenVector.so)

#include "THStack.h"

//#include "THcParmList.h"
// R__LOAD_LIBRARY(libPodd.so)
// R__LOAD_LIBRARY(libHallA.so)
// R__LOAD_LIBRARY(libdc.so)
// R__LOAD_LIBRARY(libHallC.so)

// fmt - string formatting library
//#include "fmt/core.h"
//#include "fmt/ostream.h"
// R__LOAD_LIBRARY(libfmt.so)

using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// VecOps::RVec is like std::vector with some extra bells and whistles
using inters = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;
using nlohmann::json;

void good_hms_counter(int RunNumber = 6018, int nevents = -1) {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list.json");
    json_input_file >> j;
  }

  auto runnum_str = std::to_string(RunNumber);
  if (j.find(runnum_str) == j.end()) {
    std::cout << "Run " << RunNumber << " not found in ddb2/run_list.json\n";
    std::cout << "Check that run number and replay exists. \n";
    std::cout << "If problem persists please contact Whit (717-341-1080)\n";
    std::cout << "In the meantime use: good_coin_counter_old.cxx \n";
  }
  double P0_shms_setting =
      j[runnum_str]["spectrometers"]["shms_momentum"].get<double>();
  double P0_shms = std::abs(P0_shms_setting);

  int ps7 = 7;
  if (j[runnum_str].find("daq") != j[runnum_str].end()) {
    ps7 = j[runnum_str]["daq"]["ps4"].get<int>();
  } else {
    std::cout << " using default ps4 = 7 \n";
  }
  //    The way the input rates are prescaled follows:
  //         input-rate/(2^{val - 1} + 1)
  double singles_ps_value = std::pow(2.0,ps7-1.0);
  std::cout << "prescale value " << singles_ps_value << "\n";

  std::string coda_type = "COIN";
  std::string rootfile = "ROOTfiles_online/";
  rootfile += std::string("coin_replay_production_");
  rootfile +=
      std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";

  // auto file = new TFile(rootfile.c_str());
  // new TBrowser;

  ROOT::EnableImplicitMT(24);

  Pvec4D Pbeam(0, 0, 10.598, 0.000511);

  // Detector tree 
  ROOT::RDataFrame d("T", rootfile);

  // HMS Scaler tree
  ROOT::RDataFrame d_sh("TSH", rootfile);
  //int N_scaler_events = *(d_sh.Count());
  
  std::string hpdelta = "H.gtr.dp > -10 && H.gtr.dp < 10";
  std::string epiCut = "H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 && "
                       "H.cal.etottracknorm < 2.0 && P.cal.etottracknorm<1.0";
  // && H.cal.eprtracknorm  > 0.2

  auto c_n_events_total = d.Count();

  auto c_n_events_hms = d.Filter("(fEvtHdr.fEvtType == 6) || (fEvtHdr.fEvtType == 2)").Count();


  auto d0 = d.Define("hms_e_EoverP",
                     [](doublers& EOverP) { return EOverP[EOverP > 0.5 && EOverP < 1.8]; },
                     {"H.cal.etottracknorm"})
                .Define("hms_e_EoverP_nGood", [](doublers& EOverP) { return (int)EOverP.size(); },
                        {"hms_e_EoverP"});

  auto d0_coin_only = d.Filter("fEvtHdr.fEvtType == 6");
  auto d0_hms_only  = d.Filter("fEvtHdr.fEvtType == 2");

  // Apply the electron cuts
  auto d_spec_cuts = d0.Filter(hpdelta).Filter(epiCut);
  auto d_spec_cuts_coin_only = d_spec_cuts.Filter("fEvtHdr.fEvtType == 6");
  auto d_spec_cuts_hms_only  = d_spec_cuts.Filter("fEvtHdr.fEvtType == 2");

  auto h_EOverP_0 = d0.Histo1D<doublers>(
      {"hms_e_EoverP_0", "HMS total shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.etottracknorm");
  auto h_EOverP_2 = d_spec_cuts.Histo1D<doublers>(
      {"hms_e_EoverP_2", "HMS total shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.etottracknorm");

  auto h_EprOverP_0 = d0.Histo1D<doublers>(
      {"hms_e_EoverP_0", "HMS pre-shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.eprtracknorm");
  auto h_EprOverP_2 = d_spec_cuts.Histo1D<doublers>(
      {"hms_e_EoverP_2", "HMS pre-shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.eprtracknorm");

  auto h_event_type   = d0.Histo1D({"EvtType", "EvtType", 10, 0, 10}, "fEvtHdr.fEvtType");
  auto h_event_type_2 = d_spec_cuts.Histo1D({"EvtType2", "EvtType", 10, 0, 10}, "fEvtHdr.fEvtType");

  auto h_beta_0 = d_spec_cuts.Histo1D({"h_beta_0", "h_beta_0", 150, 0, 1.5}, "H.gtr.beta");

  auto bcm4b_charge        = d_sh.Max("H.BCM4B.scalerChargeCut");
  auto el_real_scaler      = d_sh.Max("H.hEL_REAL.scaler");
  auto time_1MHz           = d_sh.Max("H.1MHz.scalerTime");
  //auto hTRIG1_ROC1_npassed = d_sh.Max("H.hTRIG1_ROC1.npassed");
  //auto H_hTRIG1_scaler     = d_sh.Max("H.hTRIG1.scaler");

  //{(hTRIG1_ROC1.npassed / H.hTRIG1.scaler)*100.0:%3.4f}
  //H.hEL_REAL.scaler/H.1MHz.scalerTime)/1000
  auto total_charge = bcm4b_charge;

  auto c_e_yield_raw  = d.Count();
  auto c_e_yield      = d_spec_cuts.Count();
  auto c_T2_yield_raw = d0_hms_only.Count();
  auto c_T6_yield_raw = d0_coin_only.Count();
  auto c_T2_yield     = d_spec_cuts_hms_only.Count();
  auto c_T6_yield     = d_spec_cuts_coin_only.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------

  double good_total_charge = *bcm4b_charge / 1000.0; // mC
  //double hms_live_time = double(*hTRIG1_ROC1_npassed) / double(*H_hTRIG1_scaler);

  double hms_scaler_yield     = ((*el_real_scaler) / good_total_charge);
  double hms_scaler_yield_unc = (std::sqrt(*el_real_scaler) / good_total_charge);

  double hms_e_yield      = (*c_e_yield_raw) / (*total_charge);
  double hms_e_yield2     = (*c_e_yield) / (*total_charge);

  double ps_cor_hms_e_yield  = hms_e_yield *singles_ps_value;
  double ps_cor_hms_e_yield2 = hms_e_yield2*singles_ps_value;

  std::cout << " good_total_charge " << good_total_charge << " \n";
  std::cout << " hms_scaler_yield " << hms_scaler_yield << " \n";


  // Update counts list
  json jruns;
  {
    std::ifstream input_file("db2/run_count_list.json");
    input_file >> jruns;
  }
  std::string run_str                = std::to_string(RunNumber);
  jruns[run_str]["hms e raw counts"] = int(*c_T2_yield_raw + *c_T6_yield_raw);
  jruns[run_str]["hms e counts"]     = int(*c_T2_yield + *c_T6_yield);
  jruns[run_str]["ps cor. hms e raw counts"] = int((*c_T2_yield_raw)*singles_ps_value + (*c_T6_yield_raw));
  jruns[run_str]["ps cor. hms e counts"]     = int((*c_T2_yield_raw)*singles_ps_value + (*c_T6_yield_raw));
  jruns[run_str]["charge bcm4b 2u cut"]     = good_total_charge;
  jruns[run_str]["hms e raw yield"] = double((*c_T2_yield_raw)*singles_ps_value + (*c_T6_yield_raw))/good_total_charge;
  jruns[run_str]["hms e yield"]     = double((*c_T2_yield)*singles_ps_value + (*c_T6_yield))/good_total_charge;
  jruns[run_str]["hms e yield unc."]     = double(std::sqrt(*c_T2_yield)*singles_ps_value + std::sqrt(*c_T6_yield))/good_total_charge;
  jruns[run_str]["hms ps4 factor"]        = singles_ps_value;
  jruns[run_str]["hms scaler yield"]      = hms_scaler_yield;
  jruns[run_str]["hms scaler yield unc."] = hms_scaler_yield_unc;
  //jruns[run_str]["hms live time"]         = hms_live_time;

  std::ofstream json_output_file("db2/run_count_list.json");
  json_output_file << std::setw(4) << jruns << "\n";

  std::cout << " ----------------------------------------------    \n";
  std::cout << " # of good coin  = "
            << ps_cor_hms_e_yield2
            << "    \n";
  std::cout << " ----------------------------------------------    \n";
  std::cout << " of  " << *c_n_events_total << " total triggers\n";
  std::cout << " and " << *c_n_events_hms << " hms triggers\n";
  // std::cout << " pions+kaons   : " << *coin_counts << "\n";
  // std::cout << " pions         : " << *pion_count << "\n";
  // std::cout << " random        : " << random_bg << "\n";
  // std::cout << " ratio         : "  << pi_K_ratio << " \n";

  // std::cout << " counts  : " << *hms_electron_counts << "\n";
  // std::cout << " charge  : " << *total_charge << " uC\n";
  // std::cout << " yield   : " << (*hms_electron_counts)  << " cnts, " <<
  // hms_e_yield  << " cnts/uC\n"; std::cout << " singles : " <<
  // (*hms_electron_counts2) << " cnts, " << hms_e_yield2 << " cnts/uC\n";
  // std::cout << " coin    : " << (*coin_counts) << " cnts, " << coin_yield <<
  // " cnts/uC\n";

  // auto s_dc_x_fp           = d2.Histo1D({"s_dc_x_fp ","xy fp; x",100,-50,50},
  // "P.dc.x_fp"); auto s_dc_y_fp           = d2.Histo1D({"s_dc_y_fp ","xy fp;
  // y",100,-50,50}, "P.dc.y_fp"); auto s_dc_xp_fp          =
  // d2.Histo1D({"s_dc_xp_fp","xy fp; xp",100,-50,50},"P.dc.xp_fp"); auto
  // s_dc_yp_fp          = d2.Histo1D({"s_dc_yp_fp","xy fp;
  // xp",100,-50,50},"P.dc.yp_fp");

  // -----------------------------------------------------------
  //
  TCanvas *c = nullptr;
  int b1 = 0;
  int b2 = 0;
  double hmax = 0.0;
  THStack *hs = nullptr;
  // TLatex latex;

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  c->Divide(2, 2);
  c->cd(1);

  // This call starts the loop over the data.
  // A DrawCopy is used so that the histogram is not deleted at the end of
  // scope, and thus stays visible on the canvas.
  h_EOverP_0->DrawCopy();

  c->cd(2);
  h_beta_0->DrawCopy();

  c->cd(3);
  gPad->SetLogy(true);
  h_event_type->DrawCopy();
  h_event_type_2->SetLineColor(2);
  h_event_type_2->DrawCopy("same");

  c->cd(4);
  //h_coin_time->DrawCopy();
  //h_coin_time2->SetLineColor(4);
  //h_coin_time2->DrawCopy("same");
  //h_coin_time3->SetLineColor(2);
  //h_coin_time3->DrawCopy("same");

  gPad->BuildLegend();

  gSystem->mkdir("results/df_example", true);
  c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) +
             ".pdf")
                .c_str());
  c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) +
             ".png")
                .c_str());

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  hs = new THStack("SHMS_cal", "SHMS calorimeter; E/P");

  h_EOverP_0->SetLineColor(1);
  h_EOverP_2->SetLineColor(4);
  h_EprOverP_0->SetLineColor(8);
  h_EprOverP_2->SetLineColor(6);

  hs->Add((TH1 *)h_EOverP_0->Clone());
  hs->Add((TH1 *)h_EOverP_2->Clone());

  hs->Add((TH1 *)h_EprOverP_0->Clone());
  hs->Add((TH1 *)h_EprOverP_2->Clone());

  hs->Draw("nostack");
  gPad->SetLogy(true);
  gPad->BuildLegend();

  c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
             ".pdf")
                .c_str());
  c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
             ".png")
                .c_str());
}
