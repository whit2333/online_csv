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

#ifdef __cpp_lib_filesystem
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

constexpr const double M_P  = .938272;
constexpr const double M_P2 = M_P * M_P;
constexpr const double M_J  = 3.096916;
constexpr const double M_J2 = M_J * M_J;

using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// VecOps::RVec is like std::vector with some extra bells and whistles
using inters   = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;
using nlohmann::json;

void good_jpsi_counter(int RunNumber = 7146, int nevents = -1, int prompt = 0, int update = 1) {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list_coin.json");
    try {
      json_input_file >> j;
    } catch (json::parse_error) {
      std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }

  auto runnum_str = std::to_string(RunNumber);
  if (j.find(runnum_str) == j.end()) {
    std::cout << "Run " << RunNumber << " not found in db2/run_list_coin.json\n";
    std::cout << "Check that run number and replay exists. \n";
    std::cout << "If problem persists please contact Sylvester (217-848-0565)\n";
  }
  double P0_shms_setting = j[runnum_str]["spectrometers"]["shms_momentum"].get<double>();
  double P0_shms         = std::abs(P0_shms_setting);

  int ps1 = 0;
  if (j[runnum_str].find("daq") != j[runnum_str].end()) {
    ps1 = j[runnum_str]["daq"]["ps1"].get<int>();
    std::cout << "ps1 = " << ps1 << "\n";
  } else {
    std::cout << " using default ps1 = 0 \n";
  }
  //  The way the input rates are prescaled follows:
  //       input-rate/(2^{val - 1} + 1)
  double shms_singles_ps_value = (ps1 >= 0) ? std::pow(2.0, ps1) : 0.;
  std::cout << "prescale value " << shms_singles_ps_value << "\n";
  int ps4 = 0;
  if (j[runnum_str].find("daq") != j[runnum_str].end()) {
    ps4 = j[runnum_str]["daq"]["ps4"].get<int>();
    std::cout << "ps4 = " << ps4 << "\n";
  } else {
    std::cout << " using default ps4 = 0 \n";
  }
  double hms_singles_ps_value = (ps4 >= 0) ? std::pow(2.0, ps4) : 0.;
  std::cout << "prescale value " << hms_singles_ps_value << "\n";

  std::string coda_type = "COIN";
  std::string rootfile  = "ROOTfiles/";
  rootfile += std::string("coin_replay_production_");
  rootfile += std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";

  bool found_good_file = false;
  if (!gSystem->AccessPathName(rootfile.c_str())) {
    TFile file(rootfile.c_str());
    if (file.IsZombie()) {
      std::cout << rootfile << " is a zombie.\n";
      std::cout
          << " Did your replay finish?  Check that the it is done before running this script.\n";
      // return;
    } else {
      found_good_file = true;
    }
  }
  if (!found_good_file) {
    rootfile = "ROOTfiles_online/";
    rootfile += std::string("shms_replay_production_");
    rootfile += std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";

    if (!gSystem->AccessPathName(rootfile.c_str())) {
      TFile file(rootfile.c_str());
      if (file.IsZombie()) {
        std::cout << rootfile << " is a zombie.\n";
        std::cout
            << " Did your replay finish?  Check that the it is done before running this script.\n";
      } else {
        found_good_file = true;
      }
    }
  }
  if (!found_good_file) {
    std::cout << " Error: suitable root file not found\n";
    return;
  }

  ROOT::EnableImplicitMT(24);

  //---------------------------------------------------------------------------
  // Detector tree
  ROOT::RDataFrame d("T", rootfile);

  // SHMS Scaler tree
  ROOT::RDataFrame d_sh("TSP", rootfile);
  // int N_scaler_events = *(d_sh.Count());

  std::string goodTrackSHMS = "P.gtr.dp > -10 && P.gtr.dp < 22 && P.tr.n == 1 && "
                              "TMath::Abs(P.gtr.th) < 0.04 && TMath::Abs(P.gtr.ph) < 0.025"
                              "&& TMath::Abs(P.gtr.y) < 5.0";
  std::string goodTrackHMS = "H.gtr.dp > -8 && H.gtr.dp < 8 && H.tr.n == 1&&"
                             "H.gtr.th < 0.04 && -0.01 < H.gtr.ph && H.gtr.ph < 0.02"
                             "&& TMath::Abs(H.gtr.y) < 2.5";

  std::string eCutSHMS = "P.cal.etottracknorm > 0.8 && "
                         "P.ngcer.npeSum > 2";
  std::string eCutHMS = "H.cal.etottracknorm > 0.6 && H.cal.etottracknorm < 2.&&"
                        "H.cer.npeSum > 0.5";

#if 0
  auto d_coin = d.Filter("fEvtHdr.fEvtType == 6");
  auto d_shms = d.Filter("fEvtHdr.fEvtType == 1");
  auto d_hms  = d.Filter("fEvtHdr.fEvtType == 4");
#endif
  auto& d_coin = d;
  auto& d_shms = d;
  auto& d_hms  = d;

  // Good track cuts
  auto dHMSGoodTrack_hms   = d_hms.Filter(goodTrackHMS);
  auto dHMSGoodTrack_coin  = d_coin.Filter(goodTrackHMS);
  auto dSHMSGoodTrack_shms = d_shms.Filter(goodTrackSHMS);
  auto dSHMSGoodTrack_coin = d_coin.Filter(goodTrackSHMS);
  auto dCOINGoodTrack      = dHMSGoodTrack_coin.Filter(goodTrackSHMS);
  // PID cuts
  // auto dHMSEl_hms   = dHMSGoodTrack_hms.Filter(eCutHMS);
  // auto dHMSEl_coin  = dHMSGoodTrack_coin.Filter(eCutHMS);
  auto dHMSEl_hms   = dHMSGoodTrack_hms;
  auto dHMSEl_coin  = dHMSGoodTrack_coin;
  auto dSHMSEl_shms = dSHMSGoodTrack_shms.Filter(eCutSHMS);
  auto dSHMSEl_coin = dSHMSGoodTrack_coin.Filter(eCutSHMS);
  auto dCOINEl      = dCOINGoodTrack.Filter(eCutSHMS + " && " + eCutSHMS);
  // Timing cuts
  // Find the timing peak
  // Find the coin peak
  auto h_coin_time =
      dCOINEl.Histo1D({"coin_time", "coin_time", 8000, 0, 1000}, "CTime.ePositronCoinTime_ROC2");
  h_coin_time->DrawClone();
  int    coin_peak_bin    = h_coin_time->GetMaximumBin();
  double coin_peak_center = h_coin_time->GetBinCenter(coin_peak_bin);
  // timing cut lambdas
  auto timing_cut = [=](double coin_time) { return std::abs(coin_time - coin_peak_center) < 2.; };
  auto anti_timing_cut = [=](double coin_time) {
    return std::abs(coin_time - coin_peak_center - 28.) < 2.;
  };
  // J/psi reconstruction
  auto p_electron = [](double px, double py, double pz) { return Pvec4D{px, py, pz, .000511}; };
  auto p_jpsi     = [](const Pvec4D& e1, const Pvec4D& e2) { return e1 + e2; };
  auto E_gamma    = [](const Pvec4D& jpsi) {
    return (M_J2 - 2. * jpsi.E() * M_P) / (2. * (jpsi.E() - M_P * jpsi.P() * cos(jpsi.Theta())));
  };
  auto E_gamma_free = [](const Pvec4D& jpsi) {
    return (jpsi.M2() - 2. * jpsi.E() * M_P) /
           (2. * (jpsi.E() - M_P * jpsi.P() * cos(jpsi.Theta())));
  };
  auto t = [](const double Egamma, Pvec4D& jpsi) {
    Pvec4D beam{0, 0, Egamma, 0};
    return (beam - jpsi).M2();
  };
  auto dCOINEl2 = dCOINEl.Define("p_electron", p_electron, {"P.gtr.px", "P.gtr.py", "P.gtr.pz"})
                      .Define("p_positron", p_electron, {"H.gtr.px", "H.gtr.py", "H.gtr.pz"})
                      .Define("p_jpsi", p_jpsi, {"p_electron", "p_positron"})
                      .Define("E_gamma", E_gamma, {"p_jpsi"})
                      .Define("E_gamma_free", E_gamma_free, {"p_jpsi"})
                      .Define("t", t, {"E_gamma", "p_jpsi"});
  // timing counts
  auto dHMSElInTime  = dHMSEl_coin.Filter(timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dHMSElRandom  = dHMSEl_coin.Filter(anti_timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dSHMSElInTime = dSHMSEl_coin.Filter(timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dSHMSElRandom = dSHMSEl_coin.Filter(anti_timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dCOINElInTime = dCOINEl2.Filter(timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dCOINElRandom = dCOINEl2.Filter(anti_timing_cut, {"CTime.ePositronCoinTime_ROC2"});

  std::string Jpsi_cut    = "p_jpsi.M() < 3.096916 + .05 && p_jpsi.M() < 3.096916 + .05";
  auto        dJpsi       = dCOINElInTime.Filter(Jpsi_cut);
  auto        dJpsiRandom = dCOINElRandom.Filter(Jpsi_cut);

  // scalers
  auto total_charge        = d_sh.Max("P.BCM4B.scalerChargeCut");
  auto shms_el_real_scaler = d_sh.Max("P.pEL_REAL.scaler");
  auto hms_el_real_scaler  = d_sh.Max("P.hEL_REAL.scaler");
  auto time_1MHz           = d_sh.Max("P.1MHz.scalerTime");
  auto time_1MHz_cut       = d_sh.Max("P.1MHz.scalerTimeCut");

  auto yield_all = d.Count();
  // 5 timing cut widths worth of random backgrounds
  auto yield_shms               = d_shms.Count();
  auto yield_hms                = d_hms.Count();
  auto yield_coin               = d_coin.Count();
  auto yield_HMSGoodTrack_hms   = dHMSGoodTrack_hms.Count();
  auto yield_HMSGoodTrack_coin  = dHMSGoodTrack_coin.Count();
  auto yield_SHMSGoodTrack_shms = dSHMSGoodTrack_shms.Count();
  auto yield_SHMSGoodTrack_coin = dSHMSGoodTrack_coin.Count();
  auto yield_COINGoodTrack      = dCOINGoodTrack.Count();
  auto yield_HMSEl_hms          = dHMSEl_hms.Count();
  auto yield_HMSEl_coin         = dHMSEl_coin.Count();
  auto yield_SHMSEl_shms        = dSHMSEl_shms.Count();
  auto yield_SHMSEl_coin        = dSHMSEl_coin.Count();
  auto yield_COINEl             = dCOINEl.Count();
  auto yield_HMSElInTime        = dHMSElInTime.Count();
  auto yield_HMSElRandom        = dHMSElRandom.Count();
  auto yield_SHMSElInTime       = dSHMSElInTime.Count();
  auto yield_SHMSElRandom       = dSHMSElRandom.Count();
  auto yield_COINElInTime       = dCOINElInTime.Count();
  auto yield_COINElRandom       = dCOINElRandom.Count();
  auto yield_jpsi_raw           = dJpsi.Count();
  auto yield_jpsi_random        = dJpsiRandom.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------
  auto n_jpsi_good  = *yield_jpsi_raw - *yield_jpsi_random / 5.;
  auto n_HMSElGood  = *yield_HMSElInTime - *yield_HMSElRandom / 5;
  auto n_SHMSElGood = *yield_SHMSElInTime - *yield_SHMSElRandom / 5;
  auto n_COINElGood = *yield_COINElInTime - *yield_COINElRandom / 5;

  double good_total_charge = *total_charge / 1000.0; // mC
  double good_time         = *time_1MHz_cut;         // s

#if 0
  // Not sure what we want from this TODO
  double shms_scaler_yield     = ((*shms_el_real_scaler) / good_total_charge);
  double hms_scaler_yield      = ((*hms_el_real_scaler) / good_total_charge);
  double shms_scaler_yield_unc = (std::sqrt(*shms_el_real_scaler) / good_total_charge);
  double shms_scaler_yield_unc = (std::sqrt(*hms_el_real_scaler) / good_total_charge);
#endif

  std::map<std::string, double> counts = {
      {"shms_raw_yield_coin", (*yield_SHMSGoodTrack_coin) / (good_total_charge)},
      {"shms_raw_yield_shms ", (*yield_SHMSGoodTrack_shms) / (good_total_charge)},
      {"hms_raw_yield_coin", (*yield_HMSGoodTrack_coin) / (good_total_charge)},
      {"hms_raw_yield_hms", (*yield_HMSGoodTrack_hms) / (good_total_charge)},
      {"coin_raw_yield", (*yield_COINGoodTrack) / (good_total_charge)},
      {"shms_e_yield_coin", (*yield_SHMSEl_coin) / (good_total_charge)},
      {"shms_e_yield_shms", (*yield_SHMSEl_shms) / (good_total_charge)},
      {"hms_e_yield_coin", (*yield_HMSEl_coin) / (good_total_charge)},
      {"hms_e_yield_hms", (*yield_HMSEl_hms) / (good_total_charge)},
      {"coin_ee_yield", (*yield_COINEl) / (good_total_charge)},
      {"hms_e_intime", (*yield_HMSElInTime) / (good_total_charge)},
      {"shms_e_intime", (*yield_SHMSElInTime) / (good_total_charge)},
      {"coin_e_intime", (*yield_COINElInTime) / (good_total_charge)},
      {"hms_e_random", (*yield_HMSElRandom) / (good_total_charge) / 5.},
      {"shms_e_random", (*yield_SHMSElRandom) / (good_total_charge) / 5.},
      {"coin_e_random", (*yield_COINElRandom) / (good_total_charge) / 5.},
      {"hms_e_good", (n_HMSElGood) / (good_total_charge)},
      {"shms_e_good", (n_SHMSElGood) / (good_total_charge)},
      {"coin_e_good", (n_COINElGood) / (good_total_charge)},
      {"ps_cor_shms_e_good", (n_SHMSElGood / good_total_charge) * shms_singles_ps_value},
      {"ps_cor_hms_e_good", (n_HMSElGood / good_total_charge) * hms_singles_ps_value},
      {"J/psi count", *yield_jpsi_raw},
      {"J/psi yield", *yield_jpsi_raw / (good_total_charge)},
      {"J/psi random background count", *yield_jpsi_random / 5.},
      {"J/psi random background", *yield_jpsi_random / (good_total_charge) / 5.},
      {"J/psi Good event count", n_jpsi_good},
      {"J/psi Good event yield", n_jpsi_good / (good_total_charge)},
      {"good_total_charge", good_total_charge}};

  // Update counts list
  json jruns;
  {
    std::ifstream input_file("db2/jpsi_run_count_list.json");
    try {
      input_file >> jruns;
    } catch (json::parse_error) {
      std::cerr << "error: json file is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }
  std::string run_str = std::to_string(RunNumber);
  std::cout << "----------------------------------------------------------" << std::endl;
  for (const auto& kv : counts) {
    std::cout << " " << kv.first;
    if (kv.first.find("yield") != std::string::npos) {
      std::cout << " (counts / mC)";
    }
    std::cout << ": " << kv.second;
    std::cout << "\n";
    jruns[run_str][kv.first] = kv.second;
  }

  jruns[run_str]["charge bcm4b 2u cut"] = good_total_charge;
  jruns[run_str]["time 1MHz 2u cut"]    = good_time;
  jruns[run_str]["shms ps1 factor"]     = shms_singles_ps_value;
  jruns[run_str]["hms ps4 factor"]      = hms_singles_ps_value;

  if (update) {
    std::cout << "Updating db2/jpsi_run_count_list.json with shms counts\n";
    std::ofstream json_output_file("db2/jpsi_run_count_list.json");
    json_output_file << std::setw(4) << jruns << "\n";
  }

#if 0
  // -----------------------------------------------------------
  //
  TCanvas* c    = nullptr;
  int      b1   = 0;
  int      b2   = 0;
  double   hmax = 0.0;
  THStack* hs   = nullptr;
  // TLatex latex;

  gSystem->mkdir("results/good_shms_counter", true);

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  c->Divide(2, 2);
  c->cd(1);

  // This call starts the loop over the data.
  // A DrawCopy is used so that the histogram is not deleted at the end of
  // scope, and thus stays visible on the canvas.
  // h_EOverP_0->DrawCopy();

  c->cd(2);
  h_beta_0->DrawCopy();

  c->cd(3);
  gPad->SetLogy(true);
  h_event_type->DrawCopy();
  h_event_type_2->SetLineColor(2);
  h_event_type_2->DrawCopy("same");

  c->cd(4);
  // h_coin_time->DrawCopy();
  // h_coin_time2->SetLineColor(4);
  // h_coin_time2->DrawCopy("same");
  // h_coin_time3->SetLineColor(2);
  // h_coin_time3->DrawCopy("same");

  gPad->BuildLegend();

  gSystem->mkdir("results/df_example", true);
  c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) + ".pdf").c_str());
  c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) + ".png").c_str());

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  c->Divide(2, 2);
  c->cd(1);
  hs = new THStack("ngcer", "ngcer; npe");
  hs->Add((TH1*)h_ngc_0->Clone());
  hs->Add((TH1*)h_ngc_1->Clone());

  hs->Draw("nostack");
  gPad->SetLogy(true);
  gPad->BuildLegend();

  c->cd(2);
  hs = new THStack("hgcer", "hgcer; npe");
  hs->Add((TH1*)h_hgc_0->Clone());
  hs->Add((TH1*)h_hgc_1->Clone());

  hs->Draw("nostack");
  gPad->SetLogy(true);
  gPad->BuildLegend();

  //// ---------------------------------------------------------
  ////
  // c = new TCanvas();
  // hs = new THStack("SHMS_cal", "SHMS calorimeter; E/P");

  // h_EOverP_0->SetLineColor(1);
  // h_EOverP_2->SetLineColor(4);
  // h_EprOverP_0->SetLineColor(8);
  // h_EprOverP_2->SetLineColor(6);

  // hs->Add((TH1 *)h_EOverP_0->Clone());
  // hs->Add((TH1 *)h_EOverP_2->Clone());

  // hs->Add((TH1 *)h_EprOverP_0->Clone());
  // hs->Add((TH1 *)h_EprOverP_2->Clone());

  // hs->Draw("nostack");
  // gPad->SetLogy(true);
  // gPad->BuildLegend();

  // c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
  //           ".pdf")
  //              .c_str());
  // c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
  //           ".png")
  //              .c_str());
#endif
}
