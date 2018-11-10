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

void good_coin_counter2(int RunNumber = 6018, int nevents = -1) {

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
  auto bcm4b = d_sh.Max("H.BCM4B.scalerChargeCut");
  std::string hpdelta = "P.gtr.dp > -10 && P.gtr.dp < 20 && "
                        "H.gtr.dp > -10 && H.gtr.dp < 10";
  std::string epiCut = "P.aero.npeSum > 1.0 && P.cal.eprtracknorm < 0.2 && "
                       "H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 && "
                       "H.cal.etottracknorm < 2.0 && P.cal.etottracknorm<1.0";
  // && H.cal.eprtracknorm  > 0.2

  auto c_n_events_total = d.Count();

  auto c_n_events_coin =
      d.Filter("(fEvtHdr.fEvtType == 6) || (fEvtHdr.fEvtType == 4)").Count();

  auto d0 = d.Define("hms_e_EoverP",
                     [](doublers &EOverP) {
                       return EOverP[EOverP > 0.5 && EOverP < 1.8];
                     },
                     {"H.cal.etottracknorm"})
                .Define("hms_e_EoverP_nGood",
                        [](doublers &EOverP) { return (int)EOverP.size(); },
                        {"hms_e_EoverP"})
                .Define("shms_pi_preshower",
                        [](doublers &EOverP) { return EOverP[EOverP < 0.3]; },
                        {"P.cal.eprtracknorm"})
                .Define("shms_pi_preshower_nGood",
                        [](doublers &EOverP) { return (int)EOverP.size(); },
                        {"shms_pi_preshower"});

  auto d2 = d0.Filter([](double npe) { return npe > 1.0; }, {"H.cer.npeSum"});
  auto d3 = d2.Filter([](int n_good) { return n_good > 0; },
                      {"shms_pi_preshower_nGood"});

  auto d_spec_cuts = d2.Filter(hpdelta).Filter(epiCut);

  // Find the coin peak
  auto h_coin_time = d_spec_cuts.Histo1D(
      {"coin_time", "coin_time", 800, 0, 100}, "CTime.ePiCoinTime_ROC2");
  int coin_peak_bin = h_coin_time->GetMaximumBin();
  double coin_peak_center = h_coin_time->GetBinCenter(coin_peak_bin);
  // std::cout << " coin peak at " << coin_peak_center << " ns\n";

  // Cut on the coin time
  auto d_coin_time_cut = d_spec_cuts.Filter(
      [=](double coin_time) {
        return std::abs(coin_time - coin_peak_center) < 2.0;
      },
      {"CTime.ePiCoinTime_ROC2"});
  auto h_coin_time2 = d_coin_time_cut.Histo1D(
      {"coin_time2", "coin_time2", 800, 0, 100}, "CTime.ePiCoinTime_ROC2");

  auto d_random = d_spec_cuts.Filter(
      [=](double coin_time) {
        return std::abs(coin_time - coin_peak_center - 28) < 10.0;
      },
      {"CTime.ePiCoinTime_ROC2"});
  auto h_coin_time3 = d_random.Histo1D(
      {"coin_time2", "coin_time2", 800, 0, 100}, "CTime.ePiCoinTime_ROC2");

  auto d2_2 = d_coin_time_cut.Filter("fEvtHdr.fEvtType == 2");
  auto d2_coin = d_coin_time_cut.Filter(
      "(fEvtHdr.fEvtType == 6) || (fEvtHdr.fEvtType == 4)");

  // Histograms lazy evaluated -> Nothing happens until they are used below (eg
  // Draw() is called)
  auto h_EOverP_0 = d0.Histo1D<doublers>(
      {"hms_e_EoverP_0", "SHMS total shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.etottracknorm");
  auto h_EOverP_2 = d2.Histo1D<doublers>(
      {"hms_e_EoverP_2", "SHMS total shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.etottracknorm");
  auto h_EOverP_3 = d3.Histo1D<doublers>(
      {"hms_e_EoverP_3", "SHMS total shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.etottracknorm");

  auto h_EprOverP_0 = d0.Histo1D<doublers>(
      {"hms_e_EoverP_0", "SHMS pre-shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.eprtracknorm");
  auto h_EprOverP_2 = d2.Histo1D<doublers>(
      {"hms_e_EoverP_2", "SHMS pre-shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.eprtracknorm");
  auto h_EprOverP_3 = d3.Histo1D<doublers>(
      {"hms_e_EoverP_3", "SHMS pre-shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.eprtracknorm");

  // auto h_EOverP_1       = d0.Histo1D<doublers>({"hms_e_EoverP_1","E/P cut;
  // HMS E/P",100,0.05,1.8},"hms_e_EoverP"); auto h_EOverP_nGood_0 =
  // d0.Histo1D({"h_EOverP_nGood_0","h_EOverP_nGood_0",10,0,10},"P.cal.ntracks");
  // auto h_EOverP_nGood_1 =
  // d0.Histo1D({"h_EOverP_nGood_1","h_EOverP_nGood_0",10,0,10},"hms_e_EoverP_nGood");

  auto h_event_type =
      d0.Histo1D({"event_type", "event_type", 10, 0, 10}, "fEvtHdr.fEvtType");

  auto total_charge = d_sh.Sum("H.BCM4A.scalerCurrent");

  auto d_hgc_cut = d2_coin.Filter(
      [=](double npe, double dp) {
        double p_track = P0_shms * (100.0 + dp) / 100.0;
        // no cerenkov cut needed when momentum is below 2.8 GeV/c
        if (p_track < 2.8) {
          return true;
        }
        return npe > 1.0;
      },
      {"P.hgcer.npeSum", "P.gtr.dp"});

  auto hms_electron_counts = d2.Count();
  auto hms_electron_counts2 = d2_2.Count();
  auto coin_counts = d2_coin.Count();

  auto rand_bg_count = d_random.Count();

  auto pion_count = d_hgc_cut.Count();

  double random_bg    = double(*rand_bg_count) / 5.0;
  double hms_e_yield  = (*hms_electron_counts) / (*total_charge);
  double hms_e_yield2 = (*hms_electron_counts2) / (*total_charge);
  double coin_yield   = (*coin_counts) / (*total_charge);

  double pi_K_ratio =
      double(*pion_count) / (double(*coin_counts) - double(*pion_count));
  double kaon_counts = (double(*coin_counts) - double(*pion_count));

  json jruns;
  {
    std::ifstream input_file("db2/run_count_list.json");
    input_file >> jruns;
  }

  // std::ifstream input_file("db2/run_count_list.json");
  // input_file >> jruns;
  // input_file.close();
  // std::cout << jruns << std::endl;;
  // Open json db file to get current values (if run already exists, replace)
  // fill the run values
  json j_current_run;
  std::string run_str = std::to_string(RunNumber);
  j_current_run["total trigger events"] = int(*c_n_events_total);
  j_current_run["coin trigger events"] = int(*c_n_events_coin);
  j_current_run["pi+K counts"] = int(*coin_counts);
  j_current_run["random background"] = double(random_bg);
  j_current_run["pion counts"] = int(*pion_count);
  j_current_run["pi/K ratio"] = pi_K_ratio;
  j_current_run["kaon counts"] = kaon_counts;
  j_current_run["pion bg sub. counts"] =
      double(double(*pion_count) -
             random_bg * double(*pion_count) / (double(*coin_counts)));
  j_current_run["kaon bg sub. counts"] = double(
      double(kaon_counts) - random_bg * kaon_counts / (double(*coin_counts)));
  // std::cout << std::setw(4) << j_current_run << "\n";
  //
  jruns[run_str] = j_current_run;

  // jruns[run_str] = nlohmann::json::parse(j_current_run.dump());
  // jruns[run_str] = nlohmann::json::parse(j_current_run.dump());
  std::ofstream json_output_file("db2/run_count_list.json");
  json_output_file << std::setw(4) << jruns << "\n";

  std::cout << " ----------------------------------------------    \n";
  std::cout << " # of good coin  = "
            << int(double(*pion_count) -
                   random_bg * double(*pion_count) / (double(*coin_counts)))
            << "    \n";
  std::cout << " ----------------------------------------------    \n";
  std::cout << " of  " << *c_n_events_total << " total triggers\n";
  std::cout << " and " << *c_n_events_coin << " coin triggers\n";
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
  h_EOverP_2->DrawCopy();

  c->cd(3);
  h_event_type->DrawCopy();

  c->cd(4);
  h_coin_time->DrawCopy();
  h_coin_time2->SetLineColor(4);
  h_coin_time2->DrawCopy("same");
  h_coin_time3->SetLineColor(2);
  h_coin_time3->DrawCopy("same");

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
  h_EOverP_3->SetLineColor(2);
  h_EprOverP_0->SetLineColor(8);
  h_EprOverP_2->SetLineColor(6);
  h_EprOverP_3->SetLineColor(9);

  hs->Add((TH1 *)h_EOverP_0->Clone());
  hs->Add((TH1 *)h_EOverP_2->Clone());
  hs->Add((TH1 *)h_EOverP_3->Clone());

  hs->Add((TH1 *)h_EprOverP_0->Clone());
  hs->Add((TH1 *)h_EprOverP_2->Clone());
  hs->Add((TH1 *)h_EprOverP_3->Clone());

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
