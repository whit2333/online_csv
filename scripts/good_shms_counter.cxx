#include "nlohmann/json.hpp"
#include <cmath>
#include <iostream>

#include "ROOT/RDataFrame.hxx"

#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

#ifdef __cpp_lib_filesystem
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

#include "monitor/DetectorDisplay.h"
#include "monitor/DisplayPlots.h"
#include "monitor/MonitoringDisplay.h"
R__LOAD_LIBRARY(libScandalizer.so)

using RDFNode = decltype(ROOT::RDataFrame{0}.Filter(""));
using Histo1DProxy =
    decltype(ROOT::RDataFrame{0}.Histo1D(ROOT::RDF::TH1DModel{"", "", 128u, 0., 0.}, ""));

struct RDFInfo {
  RDFNode&          df;
  const std::string title;
  RDFInfo(RDFNode& df, std::string_view title) : df{df}, title{title} {}
};

// =================================================================================
// Cuts
// =================================================================================
std::string goodTrack = "P.gtr.dp > -10 && P.gtr.dp < 22 && P.tr.n == 1&&"
                        "-0.05 < P.gtr.th && P.gtr.th < 0.05 && "
                        "-0.035 < P.gtr.ph && P.gtr.ph < 0.025"
                        "&& P.gtr.y > -2.0 && P.gtr.y < 3.0";
std::string eCut = "P.cal.etottracknorm > 0.80 && P.cal.etottracknorm < 2.&&"
                   "P.ngcer.npeSum > 2.";
void good_shms_counter(int RunNumber = 7146, int nevents = -1, const std::string& mode = "coin",
                       int update = 1) {

  if (mode != "coin" && mode != "shms") {
    std::cerr << "error: please specify a valid mode, either 'coin' or 'shms'" << std::endl;
    std::quick_exit(-127);
  }
  const std::string run_list_fname = "db2/run_list_" + mode + ".json";

  // ===============================================================================================
  // Initialization
  // ===============================================================================================
  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file(run_list_fname);
    try {
      json_input_file >> j;
    } catch (json::parse_error) {
      std::cerr << "error: json file, " << run_list_fname
                << ", is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }

  auto runnum_str = std::to_string(RunNumber);
  if (j.find(runnum_str) == j.end()) {
    std::cout << "Run " << RunNumber << " not found in " << run_list_fname << "\n";
    std::cout << "Check that run number and replay exists. \n";
    std::cout << "If problem persists please contact Sylvester (217-848-0565)\n";
  }

  // was this data taken with ps1 or ps2 or coin?
  int  ps1             = -1;
  int  ps2             = -1;
  int  ps5             = -1;
  int  ps6             = -1;
  bool singles_trigger = true;
  if (j[runnum_str].find("daq") != j[runnum_str].end()) {
    ps1 = j[runnum_str]["daq"]["ps1"].get<int>();
    ps2 = j[runnum_str]["daq"]["ps2"].get<int>();
    ps5 = j[runnum_str]["daq"]["ps5"].get<int>();
    ps6 = j[runnum_str]["daq"]["ps6"].get<int>();
    std::cout << "ps1 = " << ps1 << " and ps2 = " << ps2 << "\n";
    std::cout << "ps5 = " << ps1 << " and ps6 = " << ps2 << "\n";
  } else {
    std::cout << "Error: pre-scaler unspecified in " << run_list_fname << std::endl;
    std::quick_exit(-127);
  }

  int ps = std::max(ps1, ps2);
  if (ps1 > 0 && ps2 > 0) {
    std::cout << "Incorrect values for ps1 and ps2, only one should be set to be > 0\n";
    std::cout << "(ps1 = " << ps1 << " and ps2 = " << ps2 << ")" << std::endl;
    std::quick_exit(-127);
  } else if (ps == -1) {
    std::cout
        << "Warning: no data with singles pre-scaler taking, using coincidence trigger instead\n";
    singles_trigger = false;
    ps              = std::max(ps5, ps6);
    if (ps5 > 0 && ps6 > 0) {
      std::cout << "Incorrect values for ps5 and ps6, only one should be set to be > 0\n";
      std::cout << "(ps5 = " << ps5 << " and ps6 = " << ps6 << ")" << std::endl;
      std::quick_exit(-127);
    }
    std::cout << "Selected " << ((ps5 > ps6) ? "ps5" : "ps6") << std::endl;
  } else {
    std::cout << "Selected " << ((ps1 > ps2) ? "ps1" : "ps2") << std::endl;
  }

  if (ps == -1) {
    std::cout << "ERROR: no pre-scaler was set for the HMS, unable to proceed." << std::endl;
    std::quick_exit(-127);
  }
  const double ps_factor = (ps == 0) ? 1. : (std::pow(2, ps - 1) + 1);
  std::cout << "Using prescale factor " << ps_factor << std::endl;

  std::string rootfile = "ROOTfiles/" + mode + "_replay_production_" + std::to_string(RunNumber) +
                         "_" + std::to_string(nevents) + ".root ";
  bool found_good_file = false;
  std::cout << "Loading " << rootfile << "\n";
  TFile file(rootfile.c_str());
  if (file.IsZombie()) {
    std::cout << rootfile << " is a zombie.\n";
    std::cout
        << " Did your replay finish?  Check that the it is done before running this script.\n";
  } else {
    found_good_file = true;
  }
  if (!found_good_file) {
    std::cout << " Error: suitable root file not found\n";
    return;
  }
  // ===============================================================================================
  // Dataframe
  // ===============================================================================================

  ROOT::EnableImplicitMT(24);

  //---------------------------------------------------------------------------
  // Detector tree
  ROOT::RDataFrame d("T", rootfile);
  // SSHMS Scaler tree
  ROOT::RDataFrame d_sh("TSP", rootfile);

  // Select SHMS singles only
  auto dSHMS = d.Filter("fEvtHdr.fEvtType == 1");

  // Good track cuts
  auto dGoodTrack = dSHMS.Filter(goodTrack);

  // PID cuts
  auto dEl = dGoodTrack.Filter(eCut);

  // Data frame index
  std::vector<std::pair<std::string, RDFInfo>> dfs = {{"raw", {dSHMS, "SHMS"}},
                                                      {"tracked", {dGoodTrack, "Cuts: tracking"}},
                                                      {"identified", {dEl, "Cuts: tracking+PID"}}};

  // =========================================================================================
  // Histograms
  // =========================================================================================
  // P.cal.etottracknorm
  using HistoMap = std::map<std::string, Histo1DProxy>;
  HistoMap hPcalEP;
  HistoMap hPcerNphe;
  HistoMap hPdp;
  HistoMap hPth;
  HistoMap hPph;
  HistoMap hPy;
  for (auto& kval : dfs) {
    std::string name{kval.first};
    RDFInfo&    df_info{kval.second};
    // Calorimeter
    hPcalEP[name] = df_info.df.Histo1D({("P.cal.etottracknorm_" + name).c_str(),
                                        (df_info.title + ";SHMS E/P;counts").c_str(), 200, -.5, 2.},
                                       "P.cal.etottracknorm");
    // P.cer.npeSum
    hPcerNphe[name] =
        df_info.df.Histo1D({("P.ngcer.npeSum_" + name).c_str(),
                            (df_info.title + "; SHMS Cer #phe; counts ").c_str(), 200, -1, 15},
                           "P.ngcer.npeSum");
    // P.gtr.dp
    hPdp[name] = df_info.df.Histo1D({("P.gtr.dp_" + name).c_str(),
                                     (df_info.title + ";#deltap [%];counts").c_str(), 200, -30, 40},
                                    "P.gtr.dp");
    // P.gtr.th
    hPth[name] =
        df_info.df.Histo1D({("P.gtr.th_" + name).c_str(),
                            (df_info.title + ";#theta_{SHMS};counts ").c_str(), 200, -0.1, 0.1},
                           "P.gtr.th");
    // P.gtr.ph
    hPph[name] =
        df_info.df.Histo1D({("P.gtr.ph_" + name).c_str(),
                            (df_info.title + ";#phi_{SHMS};counts ").c_str(), 200, -0.1, 0.1},
                           "P.gtr.ph");
    // P.gtr.y
    hPy[name] = df_info.df.Histo1D({("P.gtr.y_" + name).c_str(),
                                    (df_info.title + ";ytar_{SHMS};counts ").c_str(), 200, -10, 10},
                                   "P.gtr.y");
  }

  // scalers
  auto total_charge  = d_sh.Max("P.BCM4B.scalerChargeCut");
  auto time_1MHz_cut = d_sh.Max("P.1MHz.scalerTimeCut");

  auto count_raw        = dSHMS.Count();
  auto count_tracked    = dGoodTrack.Count();
  auto count_identified = dEl.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------
  double good_total_charge = *total_charge / 1000.0; // mC
  double good_time         = *time_1MHz_cut;         // s

  std::map<std::string, double> counts = {
      {"count_e", (*count_identified)},
      {"count_tracked", (*count_tracked)},
      {"count_raw", (*count_raw)},
      {"count_identified_pscorr", (*count_identified) * ps_factor},
      {"count_tracked_pscorr", (*count_tracked) * ps_factor},
      {"count_raw_pscorr", (*count_raw) * ps_factor},
      {"good_total_charge", good_total_charge},
      {"good_time", good_time}};

  // Update counts list
  json jruns;
  {
    std::ifstream input_file("db2/shms_run_count_list.json");
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
    if (kv.first.find("pscorr") != std::string::npos && ps_factor == 1) {
      continue;
    }
    std::cout << " " << kv.first << ": " << kv.second;
    if (kv.first.find("count") != std::string::npos) {
      if (kv.first.find("pscorr") == std::string::npos) {
        std::cout << " --> yield (counts / mC): " << kv.second / good_total_charge << " +- "
                  << sqrt(kv.second) / good_total_charge;
      } else {
        std::cout << " --> yield (counts / mC): " << kv.second / good_total_charge << " +- "
                  << (1 / sqrt(kv.second / ps_factor)) * kv.second / good_total_charge;
      }
    }
    std::cout << "\n";
    jruns[run_str][kv.first] = kv.second;
  }

  jruns[run_str]["charge bcm4b 2u cut"] = good_total_charge;
  jruns[run_str]["time 1MHz 2u cut"]    = good_time;
  jruns[run_str]["ps_factor"]           = ps_factor;

  if (update) {
    std::cout << "Updating db2/shms_run_count_list.json with shms counts\n";
    std::ofstream json_output_file("db2/shms_run_count_list.json");
    json_output_file << std::setw(4) << jruns << "\n";
  }
  // =====================================================================================
  // Display
  // =====================================================================================
  // This is a naked pointer that we 'leak' on purpose so the connection stays alive
  // as long as the root session is running
  auto ddisplay    = new hallc::MonitoringDisplay(RunNumber);
  auto TrackingHdp = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.dp",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPdp["raw"]->SetLineColor(kGreen + 2);
        hPdp["raw"]->SetLineColor(kGreen + 2);
        hPdp["raw"]->SetLineWidth(2);
        hPdp["raw"]->DrawClone();
        hPdp["tracked"]->SetLineColor(kMagenta + 2);
        hPdp["tracked"]->SetLineWidth(2);
        hPdp["tracked"]->DrawClone("same");
        hPdp["identified"]->SetLineColor(kBlue + 1);
        hPdp["identified"]->SetLineWidth(2);
        hPdp["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingHth = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.th",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPth["raw"]->SetLineColor(kGreen + 2);
        hPth["raw"]->SetLineWidth(2);
        hPth["raw"]->DrawClone();
        hPth["tracked"]->SetLineColor(kMagenta + 2);
        hPth["tracked"]->SetLineWidth(2);
        hPth["tracked"]->DrawClone("same");
        hPth["identified"]->SetLineColor(kBlue + 1);
        hPth["identified"]->SetLineWidth(2);
        hPth["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingHph = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.ph",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPph["raw"]->SetLineColor(kGreen + 2);
        hPph["raw"]->SetLineWidth(2);
        hPph["raw"]->DrawClone();
        hPph["tracked"]->SetLineColor(kMagenta + 2);
        hPph["tracked"]->SetLineWidth(2);
        hPph["tracked"]->DrawClone("same");
        hPph["identified"]->SetLineColor(kBlue + 1);
        hPph["identified"]->SetLineWidth(2);
        hPph["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingHy = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.y",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPy["raw"]->SetLineColor(kGreen + 2);
        hPy["raw"]->SetLineWidth(2);
        hPy["raw"]->DrawClone();
        hPy["tracked"]->SetLineColor(kMagenta + 2);
        hPy["tracked"]->SetLineWidth(2);
        hPy["tracked"]->DrawClone("same");
        hPy["identified"]->SetLineColor(kBlue + 1);
        hPy["identified"]->SetLineWidth(2);
        hPy["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDcalEP = ddisplay->CreateDisplayPlot(
      "PID", "P.cal.etottracknorm",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPcalEP["raw"]->SetLineColor(kGreen + 2);
        hPcalEP["raw"]->SetLineWidth(2);
        hPcalEP["raw"]->DrawClone();
        hPcalEP["tracked"]->SetLineColor(kMagenta + 2);
        hPcalEP["tracked"]->SetLineWidth(2);
        hPcalEP["tracked"]->DrawClone("same");
        hPcalEP["identified"]->SetLineColor(kBlue + 1);
        hPcalEP["identified"]->SetLineWidth(2);
        hPcalEP["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDcerNphe = ddisplay->CreateDisplayPlot(
      "PID", "P.ngcer.nphe",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPcerNphe["raw"]->SetLineColor(kGreen + 2);
        hPcerNphe["raw"]->SetLineWidth(2);
        hPcerNphe["raw"]->DrawClone();
        hPcerNphe["tracked"]->SetLineColor(kMagenta + 2);
        hPcerNphe["tracked"]->SetLineWidth(2);
        hPcerNphe["tracked"]->DrawClone("same");
        hPcerNphe["identified"]->SetLineColor(kBlue + 1);
        hPcerNphe["identified"]->SetLineWidth(2);
        hPcerNphe["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/good_shms_counter/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
