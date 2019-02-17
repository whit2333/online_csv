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
std::string goodTrack = "H.gtr.dp > -8 && H.gtr.dp < 8 && H.tr.n == 1&&"
                        "-0.08 < H.gtr.th && H.gtr.th < 0.06 && "
                        "-0.03 < H.gtr.ph && H.gtr.ph < 0.04"
                        "&& TMath::Abs(H.gtr.y) < 2.0";
std::string eCut = "H.cal.etottracknorm > 0.80 && H.cal.etottracknorm < 2.&&"
                   "H.cer.npeSum > 1.";
void good_hms_counter(int RunNumber = 7146, int nevents = -1, const std::string& mode = "coin",
                      int update = 1) {

  if (mode != "coin" && mode != "hms") {
    std::cerr << "error: please specify a valid mode, either 'coin' or 'hms'" << std::endl;
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
    std::cout << "Run " << RunNumber << " not found in " << run_list_fname << ".json\n";
    std::cout << "Check that run number and replay exists. \n";
    std::cout << "If problem persists please contact Sylvester (217-848-0565)\n";
  }

  // was this data taken with ps3 or ps4?
  int ps3 = -1;
  int ps4 = -1;
  if (j[runnum_str].find("daq") != j[runnum_str].end()) {
    ps3 = j[runnum_str]["daq"]["ps3"].get<int>();
    ps4 = j[runnum_str]["daq"]["ps4"].get<int>();
    std::cout << "ps3 = " << ps3 << "\n";
    std::cout << "ps4 = " << ps4 << "\n";
  } else {
    std::cout << "Error: pre-scaler unspecified in " << run_list_fname << std::endl;
    std::quick_exit(-127);
  }

  if (ps3 == ps4 || (ps3 > 0 && ps4 > 0) || (ps3 < 0 && ps4 < 0)) {
    std::cout << "Incorrect values for ps3 and ps4, only one should be set to be > 0\n";
    std::cout << "(ps3 = " << ps3 << " and ps4 = " << ps4 << ")" << std::endl;
    std::quick_exit(-127);
  }
  const int    ps        = std::max(ps3, ps4);
  const double ps_factor = (ps == 0) ? 1. : (std::pow(2, ps - 1) + 1);
  std::cout << "Using prescale factor " << ps_factor << std::endl;

  std::string rootfile = "ROOTfiles/" + mode + "_replay_production_" + std::to_string(RunNumber) +
                         "_" + std::to_string(nevents) + ".root ";

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
    rootfile = "ROOTfiles_online/" + mode + "_replay_production_" + std::to_string(RunNumber) +
               "_" + std::to_string(nevents) + ".root ";

    if (!gSystem->AccessPathName(rootfile.c_str())) {
      TFile file(rootfile.c_str());
      if (file.IsZombie()) {
        std::cout << rootfile << " is a zombie.\n";
        std::cout << " Did your replay finish?  Check that the it is done before running this "
                     "script.\n";
      } else {
        found_good_file = true;
      }
    }
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
  // SHMS Scaler tree
  ROOT::RDataFrame d_sh("TSP", rootfile);

  // Select HMS singles only
  auto dHMS = d.Filter("fEvtHdr.fEvtType == 2");

  // Good track cuts
  auto dGoodTrack = dHMS.Filter(goodTrack);

  // PID cuts
  auto dEl = dGoodTrack.Filter(eCut);

  // Data frame index
  std::vector<std::pair<std::string, RDFInfo>> dfs = {{"raw", {dHMS, "HMS"}},
                                                      {"tracked", {dGoodTrack, "Cuts: tracking"}},
                                                      {"identified", {dEl, "Cuts: tracking+PID"}}};

  // =========================================================================================
  // Histograms
  // =========================================================================================
  // H.cal.etottracknorm
  using HistoMap = std::map<std::string, Histo1DProxy>;
  HistoMap hHcalEP;
  HistoMap hHcerNphe;
  HistoMap hHdp;
  HistoMap hHth;
  HistoMap hHph;
  HistoMap hHy;
  for (auto& kval : dfs) {
    std::string name{kval.first};
    RDFInfo&    df_info{kval.second};
    // Calorimeter
    hHcalEP[name] = df_info.df.Histo1D({("H.cal.etottracknorm_" + name).c_str(),
                                        (df_info.title + ";HMS E/P;counts").c_str(), 200, -.5, 2.},
                                       "H.cal.etottracknorm");
    // H.cer.npeSum
    hHcerNphe[name] =
        df_info.df.Histo1D({("H.cer.npeSum_" + name).c_str(),
                            (df_info.title + "; HMS Cer #phe; counts ").c_str(), 200, -1, 15},
                           " H.cer.npeSum ");
    // H.gtr.dp
    hHdp[name] = df_info.df.Histo1D({("H.gtr.dp_" + name).c_str(),
                                     (df_info.title + ";#deltap [%];counts").c_str(), 200, -30, 40},
                                    "H.gtr.dp");
    // H.gtr.th
    hHth[name] =
        df_info.df.Histo1D({("H.gtr.th_" + name).c_str(),
                            (df_info.title + ";#theta_{HMS};counts ").c_str(), 200, -0.1, 0.1},
                           "H.gtr.th");
    // H.gtr.ph
    hHph[name] =
        df_info.df.Histo1D({("H.gtr.ph_" + name).c_str(),
                            (df_info.title + ";#phi_{HMS};counts ").c_str(), 200, -0.1, 0.1},
                           "H.gtr.ph");
    // H.gtr.y
    hHy[name] =
        df_info.df.Histo1D({("H.gtr.y_" + name).c_str(),
                            (df_info.title + ";ytar_{HMS};counts ").c_str(), 200, -0.1, 0.1},
                           "H.gtr.y");
  }

  // scalers
  auto total_charge  = d_sh.Max("P.BCM4B.scalerChargeCut");
  auto time_1MHz_cut = d_sh.Max("P.1MHz.scalerTimeCut");

  auto count_raw        = dHMS.Count();
  auto count_tracked    = dGoodTrack.Count();
  auto count_identified = dEl.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------
  double good_total_charge = *total_charge / 1000.0; // mC
  double good_time         = *time_1MHz_cut;         // s

  std::map<std::string, double> counts = {
      {"count_identified", (*count_identified)},
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
    std::ifstream input_file("db2/hms_run_count_list.json");
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
    std::cout << " " << kv.first << ": " << kv.second << "\n";
    if (kv.first.find("count") != std::string::npos) {
      if (kv.first.find("pscorr") == std::string::npos) {
        std::cout << " --> yield (counts / mC): " << kv.second / good_total_charge << " +- "
                  << sqrt(kv.second) / good_total_charge << "\n";
      } else {
        std::cout << " --> yield (counts / mC): " << kv.second / good_total_charge << " +- "
                  << (1 / sqrt(kv.second / ps_factor)) * kv.second / good_total_charge << "\n";
      }
    }
    jruns[run_str][kv.first] = kv.second;
  }

  jruns[run_str]["charge bcm4b 2u cut"] = good_total_charge;
  jruns[run_str]["time 1MHz 2u cut"]    = good_time;
  jruns[run_str]["hms ps factor"]       = ps_factor;

  if (update) {
    std::cout << "Updating db2/jpsi_run_count_list.json with shms counts\n";
    std::ofstream json_output_file("db2/hms_run_count_list.json");
    json_output_file << std::setw(4) << jruns << "\n";
  }
  // =====================================================================================
  // Display
  // =====================================================================================
  // This is a naked pointer that we 'leak' on purpose so the connection stays alive
  // as long as the root session is running
  auto ddisplay    = new hallc::MonitoringDisplay(RunNumber);
  auto TrackingHdp = ddisplay->CreateDisplayPlot(
      "Tracking", "H.gtr.dp",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHdp["raw"]->SetLineColor(kGreen + 2);
        hHdp["raw"]->SetLineColor(kGreen + 2);
        hHdp["raw"]->SetLineWidth(2);
        hHdp["raw"]->DrawClone();
        hHdp["tracked"]->SetLineColor(kMagenta + 2);
        hHdp["tracked"]->SetLineWidth(2);
        hHdp["tracked"]->DrawClone("same");
        hHdp["identified"]->SetLineColor(kBlue + 1);
        hHdp["identified"]->SetLineWidth(2);
        hHdp["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingHth = ddisplay->CreateDisplayPlot(
      "Tracking", "H.gtr.th",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHth["raw"]->SetLineColor(kGreen + 2);
        hHth["raw"]->SetLineWidth(2);
        hHth["raw"]->DrawClone();
        hHth["tracked"]->SetLineColor(kMagenta + 2);
        hHth["tracked"]->SetLineWidth(2);
        hHth["tracked"]->DrawClone("same");
        hHth["identified"]->SetLineColor(kBlue + 1);
        hHth["identified"]->SetLineWidth(2);
        hHth["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingHph = ddisplay->CreateDisplayPlot(
      "Tracking", "H.gtr.ph",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHph["raw"]->SetLineColor(kGreen + 2);
        hHph["raw"]->SetLineWidth(2);
        hHph["raw"]->DrawClone();
        hHph["tracked"]->SetLineColor(kMagenta + 2);
        hHph["tracked"]->SetLineWidth(2);
        hHph["tracked"]->DrawClone("same");
        hHph["identified"]->SetLineColor(kBlue + 1);
        hHph["identified"]->SetLineWidth(2);
        hHph["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingHy = ddisplay->CreateDisplayPlot(
      "Tracking", "H.gtr.y",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHy["raw"]->SetLineColor(kGreen + 2);
        hHy["raw"]->SetLineWidth(2);
        hHy["raw"]->DrawClone();
        hHy["tracked"]->SetLineColor(kMagenta + 2);
        hHy["tracked"]->SetLineWidth(2);
        hHy["tracked"]->DrawClone("same");
        hHy["identified"]->SetLineColor(kBlue + 1);
        hHy["identified"]->SetLineWidth(2);
        hHy["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDcalEP = ddisplay->CreateDisplayPlot(
      "PID", "H.cal.etottracknorm",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHcalEP["raw"]->SetLineColor(kGreen + 2);
        hHcalEP["raw"]->SetLineWidth(2);
        hHcalEP["raw"]->DrawClone();
        hHcalEP["tracked"]->SetLineColor(kMagenta + 2);
        hHcalEP["tracked"]->SetLineWidth(2);
        hHcalEP["tracked"]->DrawClone("same");
        hHcalEP["identified"]->SetLineColor(kBlue + 1);
        hHcalEP["identified"]->SetLineWidth(2);
        hHcalEP["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDcerNphe = ddisplay->CreateDisplayPlot(
      "PID", "H.cer.nphe",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHcerNphe["raw"]->SetLineColor(kGreen + 2);
        hHcerNphe["raw"]->SetLineWidth(2);
        hHcerNphe["raw"]->DrawClone();
        hHcerNphe["tracked"]->SetLineColor(kMagenta + 2);
        hHcerNphe["tracked"]->SetLineWidth(2);
        hHcerNphe["tracked"]->DrawClone("same");
        hHcerNphe["identified"]->SetLineColor(kBlue + 1);
        hHcerNphe["identified"]->SetLineWidth(2);
        hHcerNphe["identified"]->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/good_hms_counter/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
