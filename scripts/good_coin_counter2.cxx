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

#include "monitor/DetectorDisplay.h"
#include "monitor/DisplayPlots.h"
#include "monitor/MonitoringDisplay.h"
//#include "monitor/ExperimentMonitor.h"
//#include "scandalizer/PostProcessors.h"
R__LOAD_LIBRARY(libScandalizer.so)

using RDFNode = decltype(ROOT::RDataFrame{0}.Filter(""));
using Histo1DProxy =
    decltype(ROOT::RDataFrame{0}.Histo1D(ROOT::RDF::TH1DModel{"", "", 128u, 0., 0.}, ""));

struct RDFInfo {
  RDFNode&          df;
  const std::string title;
  RDFInfo(RDFNode& df, std::string_view title) : df{df}, title{title} {}
};

constexpr const double M_P     = .938272;
constexpr const double M_P2    = M_P * M_P;
constexpr const double M_pion  = 0.139;
constexpr const double M_pion2 = M_pion * M_pion;
constexpr const double M_e     = .000511;

// =================================================================================
// Cuts
// =================================================================================
std::string goodTrackSHMS = "P.gtr.dp > -10 && P.gtr.dp < 22";
std::string goodTrackHMS  = "H.gtr.dp > -10 && H.gtr.dp < 10 ";

std::string piCutSHMS =
    "P.aero.npeSum > 1.0 && P.cal.eprtracknorm < 0.2 && P.cal.etottracknorm<1.0";

std::string eCutHMS = "H.cal.etottracknorm > 0.80 && H.cal.etottracknorm < 2. && "
                      "H.cer.npeSum > 1.";

std::string epiCut = "P.aero.npeSum > 1.0 && P.cal.eprtracknorm < 0.2 && "
                     "H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 && "
                     "H.cal.etottracknorm < 2.0 && P.cal.etottracknorm<1.0";

using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// =================================================================================
// reconstruction
// =================================================================================
auto p_pion = [](double px, double py, double pz) {
  return Pvec4D{px * 0.996, py * 0.996, pz * 0.996, M_e};
};
auto p_electron = [](double px, double py, double pz) {
  return Pvec4D{px * 0.994, py * 0.994, pz * 0.994, M_e};
};
auto t = [](const double Egamma, Pvec4D& jpsi) {
  Pvec4D beam{0, 0, Egamma, 0};
  return (beam - jpsi).M2();
};

bool root_file_exists(std::string rootfile) {
  bool found_good_file = false;
  if (!gSystem->AccessPathName(rootfile.c_str())) {
    TFile file(rootfile.c_str());
    if (file.IsZombie()) {
      std::cout << rootfile << " is a zombie.\n";
      std::cout
          << " Did your replay finish?  Check that the it is done before running this script.\n";
      return false;
      // return;
    } else {
      std::cout << " using : " << rootfile << "\n";
      return true;
    }
  }
  return false;
}

void good_coin_counter2(int RunNumber = 7146, int nevents = -1, int prompt = 0, int update = 1,
                        int default_count_goal = 30000, int redo_timing = 0) {

  // ===============================================================================================
  // Initialization
  // ===============================================================================================
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

  bool found_good_file = false;

  std::string rootfile =
      fmt::format("full_online/coin_replay_production_{}_{}.root", RunNumber, nevents);
  found_good_file = root_file_exists(rootfile.c_str());
  if (!found_good_file) {
    rootfile =
        fmt::format("ROOTfiles_volatile/coin_replay_production_{}_{}.root", RunNumber, nevents);
    found_good_file = root_file_exists(rootfile.c_str());
  }
  if (!found_good_file) {
    rootfile = fmt::format("ROOTfiles_csv/coin_replay_production_{}_{}.root", RunNumber, nevents);
    found_good_file = root_file_exists(rootfile.c_str());
  }
  if (!found_good_file) {
    rootfile = fmt::format("ROOTfiles/coin_replay_production_{}_{}.root", RunNumber, nevents);
    found_good_file = root_file_exists(rootfile.c_str());
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
  // int N_scaler_events = *(d_sh.Count());

  auto d_coin = d.Filter("fEvtHdr.fEvtType == 4");

  // Good track cuts
  auto dHMSGoodTrack  = d_coin.Filter(goodTrackHMS);
  auto dSHMSGoodTrack = d_coin.Filter(goodTrackSHMS);
  auto dCOINGoodTrack = dHMSGoodTrack.Filter(goodTrackSHMS)
                            .Define("p_electron", p_electron, {"H.gtr.px", "H.gtr.py", "H.gtr.pz"})
                            .Define("p_pion", p_pion, {"P.gtr.px", "P.gtr.py", "P.gtr.pz"});
  // PID cuts
  auto dHMSEl  = dHMSGoodTrack.Filter(eCutHMS);
  auto dSHMSEl = dSHMSGoodTrack.Filter(piCutSHMS);
  auto dCOINEl = dCOINGoodTrack.Filter(eCutHMS + " && " + piCutSHMS)
                     .Filter(
                         [=](double npe, double dp) {
                           double p_track = P0_shms * (100.0 + dp) / 100.0;
                           // no cerenkov cut needed when momentum is below 2.8 GeV/c
                           if (p_track < 2.8) {
                             return true;
                           }
                           return npe > 1.0;
                         },
                         {"P.hgcer.npeSum", "P.gtr.dp"});

  // Timing cuts
  // Find the timing peak
  // Find the coin peak
  double coin_peak_center = 0;
  if (redo_timing) {
    auto h_coin_time =
        dCOINEl.Histo1D({"coin_time", "coin_time", 8000, 0, 1000}, "CTime.ePositronCoinTime_ROC2");
    h_coin_time->DrawClone();
    int coin_peak_bin = h_coin_time->GetMaximumBin();
    coin_peak_center  = h_coin_time->GetBinCenter(coin_peak_bin);
    std::cout << "COINCIDENCE time peak found at: " << coin_peak_center << std::endl;
  } else {
    coin_peak_center = 43.0; // run 7240-7241
    std::cout << "COINCIDENCE time peak: using pre-calculated value at: " << coin_peak_center
              << std::endl;
    ;
  }
  // timing cut lambdas
  // TODO: evaluate timing cut and offset for random background
  auto timing_cut = [=](double coin_time) { return std::abs(coin_time - coin_peak_center) < 2.; };
  // anti-timing set to 5x width of regular
  auto anti_timing_cut = [=](double coin_time) {
    return std::abs(coin_time - coin_peak_center - 28.) < 10.;
  };

  // timing counts
  auto dHMSElInTime  = dHMSEl.Filter(timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dHMSElRandom  = dHMSEl.Filter(anti_timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dSHMSElInTime = dSHMSEl.Filter(timing_cut, {"CTime.ePositronCoinTime_ROC2"});
  auto dSHMSElRandom = dSHMSEl.Filter(anti_timing_cut, {"CTime.ePositronCoinTime_ROC2"});

  auto dCOINElInTime = dCOINEl.Filter(timing_cut, {"CTime.ePiCoinTime_ROC2"});
  auto dCOINElRandom = dCOINEl.Filter(anti_timing_cut, {"CTime.ePiCoinTime_ROC2"});

  // Output root file
  auto out_file =
      new TFile(fmt::format("monitoring/{}/good_csv_counter.root", RunNumber).c_str(), "UPDATE");
  out_file->cd();

  // =========================================================================================
  // Histograms
  // =========================================================================================
  // 2D correlations
  auto hTheta2DNoCuts = d_coin.Histo2D(
      {"theta2D", "No cuts;#theta_{SHMS};#theta_{HMS};#counts", 50, -.1, .1, 50, -.1, .1},
      "P.gtr.th", "H.gtr.th");
  auto hTheta2DTracking = dCOINGoodTrack.Histo2D(
      {"theta2D", "Cuts: tracking;#theta_{SHMS};#theta_{HMS};#counts", 50, -.1, .1, 50, -.1, .1},
      "P.gtr.th", "H.gtr.th");
  auto hTheta2DPID =
      dCOINEl.Histo2D({"theta2D", "Cuts: tracking+PID;#theta_{SHMS};#theta_{HMS};#counts", 50, -.1,
                       .1, 50, -.1, .1},
                      "P.gtr.th", "H.gtr.th");
  auto hTheta2DTiming =
      dCOINElInTime.Histo2D({"theta2D", "Cuts: tracking+PID;#theta_{SHMS};#theta_{HMS};#counts", 50,
                             -.1, .1, 50, -.1, .1},
                            "P.gtr.th", "H.gtr.th");
  // timing
  auto hCoinTimeNoCuts =
      d_coin.Histo1D({"coin_time.NoCuts", "No Cuts;coin_time;counts", 8000, 0, 1000},
                     "CTime.ePositronCoinTime_ROC2");
  auto hCoinTimeTracking = dCOINGoodTrack.Histo1D(
      {"coin_time.Tracking", "Cuts: Tracking;coin_time;counts", 8000, 0, 1000},
      "CTime.ePositronCoinTime_ROC2");
  auto hCoinTimePID =
      dCOINEl.Histo1D({"coin_time.PID", "Cuts: Tracking+PID;coin_time;counts", 8000, 0, 1000},
                      "CTime.ePositronCoinTime_ROC2");
  auto hCoinTimeTiming = dCOINElInTime.Histo1D(
      {"coin_time.Timing", "Cuts: Tracking+PID+Timing;coin_time;counts", 8000, 0, 1000},
      "CTime.ePositronCoinTime_ROC2");

  auto hRandCoinTimePID = dCOINElRandom.Histo1D(
      {"rand_coin_time.PID", "Cuts: Tracking+PID;coin_time;counts", 8000, 0, 1000},
      "CTime.ePositronCoinTime_ROC2");

  // P.gtr.dp
  auto hPdpNoCuts =
      d_coin.Histo1D({"P.gtr.dp.NoCuts", "No Cuts;#deltap [%];counts", 200, -30, 40}, "P.gtr.dp");
  auto hPdpTracking = dSHMSGoodTrack.Histo1D(
      {"P.gtr.dp.Tracking", "Cuts: Tracking;#deltap [%];counts", 200, -30, 40}, "P.gtr.dp");
  auto hPdpPID = dSHMSEl.Histo1D(
      {"P.gtr.dp.PID", "Cuts: Tracking+PID;#deltap [%];counts", 200, -30, 40}, "P.gtr.dp");
  auto hPdpTiming = dSHMSElInTime.Histo1D(
      {"P.gtr.dp.Timing", "Cuts: Tracking+PID+Timing;#deltap [%];counts", 200, -30, 40},
      "P.gtr.dp");
  // P.gtr.th
  auto hPthNoCuts = d_coin.Histo1D(
      {"P.gtr.th.NoCuts", "No Cuts;#theta_{SHMS};counts", 200, -0.1, 0.1}, "P.gtr.th");
  auto hPthTracking = dSHMSGoodTrack.Histo1D(
      {"P.gtr.th.Tracking", "Cuts: Tracking;#theta_{SHMS};counts", 200, -0.1, 0.1}, "P.gtr.th");
  auto hPthPID = dSHMSEl.Histo1D(
      {"P.gtr.th.PID", "Cuts: Tracking+PID;#theta_{SHMS};counts", 200, -0.1, 0.1}, "P.gtr.th");
  auto hPthTiming = dSHMSElInTime.Histo1D(
      {"P.gtr.th.Timing", "Cuts: Tracking+PID+Timing;#theta_{SHMS};counts", 200, -0.1, 0.1},
      "P.gtr.th");
  // P.gtr.ph
  auto hPphNoCuts =
      d_coin.Histo1D({"P.gtr.ph.NoCuts", "No Cuts;#phi_{SHMS};counts", 200, -0.1, 0.1}, "P.gtr.ph");
  auto hPphTracking = dSHMSGoodTrack.Histo1D(
      {"P.gtr.ph.Tracking", "Cuts: Tracking;#phi_{SHMS};counts", 200, -0.1, 0.1}, "P.gtr.ph");
  auto hPphPID = dSHMSEl.Histo1D(
      {"P.gtr.ph.PID", "Cuts: Tracking+PID;#phi_{SHMS};counts", 200, -0.1, 0.1}, "P.gtr.ph");
  auto hPphTiming = dSHMSElInTime.Histo1D(
      {"P.gtr.ph.Timing", "Cuts: Tracking+PID+Timing;#phi_{SHMS};counts", 200, -0.1, 0.1},
      "P.gtr.ph");
  // P.gtr.y
  auto hPyNoCuts =
      d_coin.Histo1D({"P.gtr.y.NoCuts", "No Cuts;ytar;counts", 200, -10., 10.}, "P.gtr.y");
  auto hPyTracking = dSHMSGoodTrack.Histo1D(
      {"P.gtr.y.Tracking", "Cuts: Tracking;ytar;counts", 200, -10., 10.}, "P.gtr.y");
  auto hPyPID =
      dSHMSEl.Histo1D({"P.gtr.y.PID", "Cuts: Tracking+PID;ytar;counts", 200, -10., 10.}, "P.gtr.y");
  auto hPyTiming = dSHMSElInTime.Histo1D(
      {"P.gtr.y.Timing", "Cuts: Tracking+PID+Timing;ytar;counts", 200, -10., 10.}, "P.gtr.y");
  // P.cal.etottracknorm
  auto hPcalEPNoCuts =
      d_coin.Histo1D({"P.cal.etottracknorm.NoCuts", "No Cuts;SHMS E/P;counts", 200, -.5, 1.5},
                     "P.cal.etottracknorm");
  auto hPcalEPTracking = dSHMSGoodTrack.Histo1D(
      {"P.cal.etottracknorm.Tracking", "Cuts: Tracking;SHMS E/P;counts", 200, -.5, 1.5},
      "P.cal.etottracknorm");
  auto hPcalEPPID = dSHMSEl.Histo1D(
      {"P.cal.etottracknorm.PID", "Cuts: Tracking+PID;SHMS E/P;counts", 200, -.5, 1.5},
      "P.cal.etottracknorm");
  auto hPcalEPAll = dCOINElInTime.Histo1D(
      {"P.cal.etottracknorm.All", "Cuts: Tracking+PID+Coincidence;SHMS E/P;counts", 200, -.5, 1.5},
      "P.cal.etottracknorm");
  // P.ngcer.npeSum
  auto hPcerNpheNoCuts = d_coin.Histo1D(
      {"P.ngcer.npeSum.NoCuts", "No Cuts;SHMS NGC #phe;counts", 200, -5, 76}, "P.ngcer.npeSum");
  auto hPcerNpheTracking = dSHMSGoodTrack.Histo1D(
      {"P.ngcer.npeSum.Tracking", "Cuts: Tracking;SHMS NGC #phe;counts", 200, -5, 76},
      "P.ngcer.npeSum");
  auto hPcerNphePID = dSHMSEl.Histo1D(
      {"P.ngcer.npeSum.PID", "Cuts: Tracking+PID;SHMS NGC #phe;counts", 200, -5, 76},
      "P.ngcer.npeSum");
  auto hPcerNpheAll = dCOINElInTime.Histo1D(
      {"P.ngcer.npeSum.All", "Cuts: Tracking+PID+Coincidence;SHMS NGC #phe;counts", 200, -5, 76},
      "P.ngcer.npeSum");
  // P.hgcer.npeSum
  auto hPhgcerNpheNoCuts = d_coin.Histo1D(
      {"P.hgcer.npeSum.NoCuts", "No Cuts;SHMS HGC #phe;counts", 200, -5, 76}, "P.hgcer.npeSum");
  auto hPhgcerNpheTracking = dSHMSGoodTrack.Histo1D(
      {"P.hgcer.npeSum.Tracking", "Cuts: Tracking;SHMS HGC #phe;counts", 200, -5, 76},
      "P.hgcer.npeSum");
  auto hPhgcerNphePID = dSHMSEl.Histo1D(
      {"P.hgcer.npeSum.PID", "Cuts: Tracking+PID;SHMS HGC #phe;counts", 200, -5, 76},
      "P.hgcer.npeSum");
  auto hPhgcerNpheAll = dCOINElInTime.Histo1D(
      {"P.hgcer.npeSum.All", "Cuts: Tracking+PID+Coincidence;SHMS HGC #phe;counts", 200, -5, 76},
      "P.hgcer.npeSum");
  // H.cal.etottracknorm
  auto hHcalEPNoCuts =
      d_coin.Histo1D({"H.cal.etottracknorm.NoCuts", "No Cuts;HMS E/P;counts", 200, -.5, 1.5},
                     "H.cal.etottracknorm");
  auto hHcalEPTracking = dHMSGoodTrack.Histo1D(
      {"H.cal.etottracknorm.Tracking", "Cuts: Tracking;HMS E/P;counts", 200, -.5, 1.5},
      "H.cal.etottracknorm");
  auto hHcalEPPID = dHMSEl.Histo1D(
      {"H.cal.etottracknorm.PID", "Cuts: Tracking+PID;HMS E/P;counts", 200, -.5, 1.5},
      "H.cal.etottracknorm");
  auto hHcalEPAll = dCOINElInTime.Histo1D(
      {"H.cal.etottracknorm.All", "Cuts: Tracking+PID+Coincidence;HMS E/P;counts", 200, -.5, 1.5},
      "H.cal.etottracknorm");
  // H.cer.npeSum
  auto hHcerNpheNoCuts = d_coin.Histo1D(
      {"H.cer.npeSum.NoCuts", "No Cuts;HMS Cer #phe;counts", 200, -1, 15}, "H.cer.npeSum");
  auto hHcerNpheTracking = dHMSGoodTrack.Histo1D(
      {"H.cer.npeSum.Tracking", "Cuts: Tracking;HMS Cer #phe;counts", 200, -1, 15}, "H.cer.npeSum");
  auto hHcerNphePID = dHMSEl.Histo1D(
      {"H.cer.npeSum.PID", "Cuts: Tracking+PID+Coincidence;HMS Cer #phe;counts", 200, -1, 15},
      "H.cer.npeSum");
  auto hHcerNpheAll = dCOINElInTime.Histo1D(
      {"H.cer.npeSum.PID", "Cuts: Tracking+PID+Coincidence;HMS Cer #phe;counts", 200, -1, 15},
      "H.cer.npeSum");
  // H.gtr.dp
  auto hHdpNoCuts =
      d_coin.Histo1D({"H.gtr.dp.NoCuts", "No Cuts;#deltap [%];counts", 200, -30, 40}, "H.gtr.dp");
  auto hHdpTracking = dHMSGoodTrack.Histo1D(
      {"H.gtr.dp.Tracking", "Cuts: Tracking;#deltap [%];counts", 200, -30, 40}, "H.gtr.dp");
  auto hHdpPID = dHMSEl.Histo1D(
      {"H.gtr.dp.PID", "Cuts: Tracking+PID;#deltap [%];counts", 200, -30, 40}, "H.gtr.dp");
  auto hHdpTiming = dHMSElInTime.Histo1D(
      {"H.gtr.dp.Timing", "Cuts: Tracking+PID+Timing;#deltap [%];counts", 200, -30, 40},
      "H.gtr.dp");
  // H.gtr.th
  auto hHthNoCuts = d_coin.Histo1D(
      {"H.gtr.th.NoCuts", "No Cuts;#theta_{HMS};counts", 200, -0.1, 0.1}, "H.gtr.th");
  auto hHthTracking = dHMSGoodTrack.Histo1D(
      {"H.gtr.th.Tracking", "Cuts: Tracking;#theta_{HMS};counts", 200, -0.1, 0.1}, "H.gtr.th");
  auto hHthPID = dHMSEl.Histo1D(
      {"H.gtr.th.PID", "Cuts: Tracking+PID;#theta_{HMS};counts", 200, -0.1, 0.1}, "H.gtr.th");
  auto hHthTiming = dHMSElInTime.Histo1D(
      {"H.gtr.th.Timing", "Cuts: Tracking+PID+Timing;#theta_{HMS};counts", 200, -0.1, 0.1},
      "H.gtr.th");
  // H.gtr.ph
  auto hHphNoCuts =
      d_coin.Histo1D({"H.gtr.ph.NoCuts", "No Cuts;#phi_{HMS};counts", 200, -0.1, 0.1}, "H.gtr.ph");
  auto hHphTracking = dHMSGoodTrack.Histo1D(
      {"H.gtr.ph.Tracking", "Cuts: Tracking;#phi_{HMS};counts", 200, -0.1, 0.1}, "H.gtr.ph");
  auto hHphPID = dHMSEl.Histo1D(
      {"H.gtr.ph.PID", "Cuts: Tracking+PID;#phi_{HMS};counts", 200, -0.1, 0.1}, "H.gtr.ph");
  auto hHphTiming = dHMSElInTime.Histo1D(
      {"H.gtr.ph.Timing", "Cuts: Tracking+PID+Timing;#phi_{HMS};counts", 200, -0.1, 0.1},
      "H.gtr.ph");
  // H.gtr.y
  auto hHyNoCuts =
      d_coin.Histo1D({"H.gtr.y.NoCuts", "No Cuts;ytar;counts", 200, -10., 10.}, "H.gtr.y");
  auto hHyTracking = dHMSGoodTrack.Histo1D(
      {"H.gtr.y.Tracking", "Cuts: Tracking;ytar;counts", 200, -10., 10.}, "H.gtr.y");
  auto hHyPID =
      dHMSEl.Histo1D({"H.gtr.y.PID", "Cuts: Tracking+PID;ytar;counts", 200, -10., 10.}, "H.gtr.y");
  auto hHyTiming = dHMSElInTime.Histo1D(
      {"H.gtr.y.Timing", "Cuts: Tracking+PID+Timing;ytar;counts", 200, -10., 10.}, "H.gtr.y");

  // scalers
  auto total_charge        = d_sh.Max("P.BCM4B.scalerChargeCut");
  auto shms_el_real_scaler = d_sh.Max("P.pEL_REAL.scaler");
  auto hms_el_real_scaler  = d_sh.Max("P.hEL_REAL.scaler");
  auto time_1MHz           = d_sh.Max("P.1MHz.scalerTime");
  auto time_1MHz_cut       = d_sh.Max("P.1MHz.scalerTimeCut");

  auto yield_all = d.Count();
  // 5 timing cut widths worth of random backgrounds
  auto yield_coin          = d_coin.Count();
  auto yield_HMSGoodTrack  = dHMSGoodTrack.Count();
  auto yield_SHMSGoodTrack = dSHMSGoodTrack.Count();
  auto yield_COINGoodTrack = dCOINGoodTrack.Count();
  auto yield_HMSEl         = dHMSEl.Count();
  auto yield_SHMSEl        = dSHMSEl.Count();
  auto yield_COINEl        = dCOINEl.Count();
  auto yield_HMSElInTime   = dHMSElInTime.Count();
  auto yield_HMSElRandom   = dHMSElRandom.Count();
  auto yield_SHMSElInTime  = dSHMSElInTime.Count();
  auto yield_SHMSElRandom  = dSHMSElRandom.Count();
  auto yield_COINElInTime  = dCOINElInTime.Count();
  auto yield_COINElRandom  = dCOINElRandom.Count();
  auto yield_coin_raw      = dCOINElInTime.Count();
  auto yield_coin_random   = dCOINElRandom.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------
  auto n_coin_good  = *yield_coin_raw - *yield_coin_random / 5.;
  auto n_HMSElGood  = *yield_HMSElInTime - *yield_HMSElRandom / 5;
  auto n_SHMSElGood = *yield_SHMSElInTime - *yield_SHMSElRandom / 5;
  auto n_COINElGood = *yield_COINElInTime - *yield_COINElRandom / 5;

  double good_total_charge = *total_charge / 1000.0; // mC
  double good_time         = *time_1MHz_cut;         // s

  std::map<std::string, double> counts = {
      {"shms_raw_yield_coin", (*yield_SHMSGoodTrack) / (good_total_charge)},
      {"hms_raw_yield_coin", (*yield_HMSGoodTrack) / (good_total_charge)},
      {"coin_raw_yield", (*yield_COINGoodTrack) / (good_total_charge)},
      {"shms_e_yield_coin", (*yield_SHMSEl) / (good_total_charge)},
      {"hms_e_yield_coin", (*yield_HMSEl) / (good_total_charge)},
      {"coin_ee_yield", (*yield_COINEl) / (good_total_charge)},
      {"hms_intime", (*yield_HMSElInTime) / (good_total_charge)},
      {"shms_intime", (*yield_SHMSElInTime) / (good_total_charge)},
      {"coin_intime", (*yield_COINElInTime) / (good_total_charge)},
      {"coin_intime count", (*yield_COINElInTime) },
      {"hms_random", (*yield_HMSElRandom) / (good_total_charge) / 5.},
      {"shms_random", (*yield_SHMSElRandom) / (good_total_charge) / 5.},
      {"coin_random", (*yield_COINElRandom) / (good_total_charge) / 5.},
      {"coin_random count", (*yield_COINElRandom)  / 5.},
      {"coin_bg_corrected count", (*yield_COINElInTime) - (*yield_COINElRandom)  / 5.},
      {"hms_e_good", (n_HMSElGood) / (good_total_charge)},
      {"shms_e_good", (n_SHMSElGood) / (good_total_charge)},
      {"coin_e_good", (n_COINElGood) / (good_total_charge)},
      {"coin count", n_COINElGood},
      {"coin yield", *yield_coin_raw / (good_total_charge)},
      {"coin random background count", *yield_coin_random / 5.},
      {"coin random background", *yield_coin_random / (good_total_charge) / 5.},
      {"coin Good event count", n_coin_good},
      {"coin Good event yield", n_coin_good / (good_total_charge)},
      {"good_total_charge", good_total_charge}};

  // Update counts list
  json jruns;
  {
    std::ifstream input_file("db2/csv_run_count_list.json");
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

  if (1) {
    std::cout << "Updating db2/coin_run_count_list.json with shms counts\n";
    std::ofstream json_output_file("db2/csv_run_count_list.json");
    json_output_file << std::setw(4) << jruns << "\n";
  }

  int count_goal = default_count_goal;

  // if (prompt) {
  //  std::cout << "----------------------------------------------------------\n";
  //  std::cout << "Reference the run plan for this setting found on the wiki\n"
  //               "       https://hallcweb.jlab.org/wiki/index.php/CSV_Fall_2018_Run_Plan\n";
  //  std::cout << "----------------------------------------------------------\n";
  //  std::cout << "Please enter **total count** goal for this setting. \n";
  //  std::cout << "   Desired count [default=30000] : ";
  //  // std::cin >> count_goal ;
  //  // int number = 0;
  //  if (std::cin.peek() == '\n') { // check if next character is newline
  //    // count_goal = 30000; //and assign the default
  //  } else if (!(std::cin >> count_goal)) { // be sure to handle invalid input
  //    std::cout << "Invalid input.\n";
  //    // error handling
  //  }
  //  std::cout << "\n";
  //}

  double n_seconds      = double(*time_1MHz_cut);
  int    nev_tot        = (*yield_all);
  double time_remaining = (count_goal * n_seconds) - (n_seconds);
  double charge_goal    = count_goal * (good_total_charge) / counts["coin_bg_corrected count"];
  double goal_Nevents   = double(nev_tot) * charge_goal / good_total_charge;

  std::cout << " GOOD COUNTS : " << counts["coin_bg_corrected count"] << "\n";

  std::cout << "----------------------------------------------------------\n";
  std::cout << " N events to reach goal  : " << goal_Nevents / 1000000.0 << "M events\n";
  std::cout << " Charge   to reach goal  : " << charge_goal << " mC\n";

  if (update) {
    std::string cmd = "caput hcRunPlanChargeGoal " + std::to_string(charge_goal) + " &> /dev/null ";
    system(cmd.c_str());

    cmd = "caput hcRunPlanNTrigEventsGoal " + std::to_string(goal_Nevents) + " &> /dev/null ";
    system(cmd.c_str());

    cmd = "caput hcRunPlanCountGoal " + std::to_string(count_goal) + " &> /dev/null ";
    system(cmd.c_str());
  }

  // =====================================================================================
  // Display
  // =====================================================================================
  // This is a naked pointer that we 'leak' on purpose so the connection stays alive
  // as long as the root session is running
  auto ddisplay    = new hallc::MonitoringDisplay(RunNumber);
  auto TrackingPdp = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.dp",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPdpNoCuts->SetLineColor(kGreen + 2);
        hPdpNoCuts->SetLineWidth(2);
        hPdpNoCuts->GetYaxis()->SetRangeUser(.8, hPdpNoCuts->GetMaximum() * 1.2);
        hPdpNoCuts->DrawClone();
        hPdpTracking->SetLineColor(kMagenta + 2);
        hPdpTracking->SetLineWidth(2);
        hPdpTracking->DrawClone("same");
        hPdpPID->SetLineColor(kBlue + 1);
        hPdpPID->SetLineWidth(2);
        hPdpPID->DrawClone("same");
        hPdpTiming->SetLineColor(kBlack);
        hPdpTiming->SetLineWidth(2);
        hPdpTiming->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingPth = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.th",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPthNoCuts->SetLineColor(kGreen + 2);
        hPthNoCuts->SetLineWidth(2);
        hPthNoCuts->GetYaxis()->SetRangeUser(.8, hPthNoCuts->GetMaximum() * 1.2);
        hPthNoCuts->DrawClone();
        hPthTracking->SetLineColor(kMagenta + 2);
        hPthTracking->SetLineWidth(2);
        hPthTracking->DrawClone("same");
        hPthPID->SetLineColor(kBlue + 1);
        hPthPID->SetLineWidth(2);
        hPthPID->DrawClone("same");
        hPthTiming->SetLineColor(kBlack);
        hPthTiming->SetLineWidth(2);
        hPthTiming->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingPph = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.ph",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPphNoCuts->SetLineColor(kGreen + 2);
        hPphNoCuts->SetLineWidth(2);
        hPphNoCuts->GetYaxis()->SetRangeUser(.8, hPphNoCuts->GetMaximum() * 1.2);
        hPphNoCuts->DrawClone();
        hPphTracking->SetLineColor(kMagenta + 2);
        hPphTracking->SetLineWidth(2);
        hPphTracking->DrawClone("same");
        hPphPID->SetLineColor(kBlue + 1);
        hPphPID->SetLineWidth(2);
        hPphPID->DrawClone("same");
        hPphTiming->SetLineColor(kBlack);
        hPphTiming->SetLineWidth(2);
        hPphTiming->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingPy = ddisplay->CreateDisplayPlot(
      "Tracking", "P.gtr.y",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPyNoCuts->SetLineColor(kGreen + 2);
        hPyNoCuts->SetLineWidth(2);
        hPyNoCuts->GetYaxis()->SetRangeUser(.8, hPyNoCuts->GetMaximum() * 1.2);
        hPyNoCuts->DrawClone();
        hPyTracking->SetLineColor(kMagenta + 2);
        hPyTracking->SetLineWidth(2);
        hPyTracking->DrawClone("same");
        hPyPID->SetLineColor(kBlue + 1);
        hPyPID->SetLineWidth(2);
        hPyPID->DrawClone("same");
        hPyTiming->SetLineColor(kBlack);
        hPyTiming->SetLineWidth(2);
        hPyTiming->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TrackingHdp = ddisplay->CreateDisplayPlot(
      "Tracking", "H.gtr.dp",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHdpNoCuts->SetLineColor(kGreen + 2);
        hHdpNoCuts->SetLineWidth(2);
        hHdpNoCuts->GetYaxis()->SetRangeUser(.8, hHdpNoCuts->GetMaximum() * 1.2);
        hHdpNoCuts->DrawClone();
        hHdpTracking->SetLineColor(kMagenta + 2);
        hHdpTracking->SetLineWidth(2);
        hHdpTracking->DrawClone("same");
        hHdpPID->SetLineColor(kBlue + 1);
        hHdpPID->SetLineWidth(2);
        hHdpPID->DrawClone("same");
        hHdpTiming->SetLineColor(kBlack);
        hHdpTiming->SetLineWidth(2);
        hHdpTiming->DrawClone("same");
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
        hHthNoCuts->SetLineColor(kGreen + 2);
        hHthNoCuts->SetLineWidth(2);
        hHthNoCuts->GetYaxis()->SetRangeUser(.8, hHthNoCuts->GetMaximum() * 1.2);
        hHthNoCuts->DrawClone();
        hHthTracking->SetLineColor(kMagenta + 2);
        hHthTracking->SetLineWidth(2);
        hHthTracking->DrawClone("same");
        hHthPID->SetLineColor(kBlue + 1);
        hHthPID->SetLineWidth(2);
        hHthPID->DrawClone("same");
        hHthTiming->SetLineColor(kBlack);
        hHthTiming->SetLineWidth(2);
        hHthTiming->DrawClone("same");
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
        hHphNoCuts->SetLineColor(kGreen + 2);
        hHphNoCuts->SetLineWidth(2);
        hHphNoCuts->GetYaxis()->SetRangeUser(.8, hHphNoCuts->GetMaximum() * 1.2);
        hHphNoCuts->DrawClone();
        hHphTracking->SetLineColor(kMagenta + 2);
        hHphTracking->SetLineWidth(2);
        hHphTracking->DrawClone("same");
        hHphPID->SetLineColor(kBlue + 1);
        hHphPID->SetLineWidth(2);
        hHphPID->DrawClone("same");
        hHphTiming->SetLineColor(kBlack);
        hHphTiming->SetLineWidth(2);
        hHphTiming->DrawClone("same");
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
        hHyNoCuts->SetLineColor(kGreen + 2);
        hHyNoCuts->SetLineWidth(2);
        hHyNoCuts->GetYaxis()->SetRangeUser(.8, hHyNoCuts->GetMaximum() * 1.2);
        hHyNoCuts->DrawClone();
        hHyTracking->SetLineColor(kMagenta + 2);
        hHyTracking->SetLineWidth(2);
        hHyTracking->DrawClone("same");
        hHyPID->SetLineColor(kBlue + 1);
        hHyPID->SetLineWidth(2);
        hHyPID->DrawClone("same");
        hHyTiming->SetLineColor(kBlack);
        hHyTiming->SetLineWidth(2);
        hHyTiming->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto TimingCoinTime = ddisplay->CreateDisplayPlot(
      "Timing", "coin_time",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hCoinTimeNoCuts->SetLineColor(kGreen + 2);
        hCoinTimeNoCuts->SetLineWidth(2);
        hCoinTimeNoCuts->GetYaxis()->SetRangeUser(.8, hCoinTimeNoCuts->GetMaximum() * 1.2);
        hCoinTimeNoCuts->DrawClone();
        hCoinTimeTracking->SetLineColor(kMagenta + 2);
        hCoinTimeTracking->SetLineWidth(2);
        hCoinTimeTracking->DrawClone("same");
        hCoinTimePID->SetLineColor(kBlue + 1);
        hCoinTimePID->SetLineWidth(2);
        hCoinTimePID->DrawClone("same");
        hCoinTimeTiming->SetLineColor(kBlack);
        hCoinTimeTiming->SetLineWidth(2);
        hCoinTimeTiming->DrawClone("same");
        hRandCoinTimePID->SetLineColor(kRed);
        hRandCoinTimePID->SetLineWidth(2);
        hRandCoinTimePID->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDHcalEP = ddisplay->CreateDisplayPlot(
      "PID", "H.cal.etottracknorm",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHcalEPNoCuts->SetLineColor(kGreen + 2);
        hHcalEPNoCuts->SetLineWidth(2);
        hHcalEPNoCuts->GetYaxis()->SetRangeUser(.8, hHcalEPNoCuts->GetMaximum() * 1.2);
        hHcalEPNoCuts->DrawClone();
        hHcalEPTracking->SetLineColor(kMagenta + 2);
        hHcalEPTracking->SetLineWidth(2);
        hHcalEPTracking->DrawClone("same");
        hHcalEPPID->SetLineColor(kBlue + 1);
        hHcalEPPID->SetLineWidth(2);
        hHcalEPPID->DrawClone("same");
        hHcalEPAll->SetLineColor(kBlack);
        hHcalEPAll->SetLineWidth(2);
        hHcalEPAll->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDHcerNphe = ddisplay->CreateDisplayPlot(
      "PID", "H.cer.nphe",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hHcerNpheNoCuts->SetLineColor(kGreen + 2);
        hHcerNpheNoCuts->SetLineWidth(2);
        hHcerNpheNoCuts->GetYaxis()->SetRangeUser(.8, hHcerNpheNoCuts->GetMaximum() * 1.2);
        hHcerNpheNoCuts->DrawClone();
        hHcerNpheTracking->SetLineColor(kMagenta + 2);
        hHcerNpheTracking->SetLineWidth(2);
        hHcerNpheTracking->DrawClone("same");
        hHcerNphePID->SetLineColor(kBlue + 1);
        hHcerNphePID->SetLineWidth(2);
        hHcerNphePID->DrawClone("same");
        hHcerNpheAll->SetLineColor(kBlack);
        hHcerNpheAll->SetLineWidth(2);
        hHcerNpheAll->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDPcalEP = ddisplay->CreateDisplayPlot(
      "PID", "P.cal.etottracknorm",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPcalEPNoCuts->SetLineColor(kGreen + 2);
        hPcalEPNoCuts->SetLineWidth(2);
        hPcalEPNoCuts->GetYaxis()->SetRangeUser(.8, hPcalEPNoCuts->GetMaximum() * 1.2);
        hPcalEPNoCuts->DrawClone();
        hPcalEPTracking->SetLineColor(kMagenta + 2);
        hPcalEPTracking->SetLineWidth(2);
        hPcalEPTracking->DrawClone("same");
        hPcalEPPID->SetLineColor(kBlue + 1);
        hPcalEPPID->SetLineWidth(2);
        hPcalEPPID->DrawClone("same");
        hPcalEPAll->SetLineColor(kBlack);
        hPcalEPAll->SetLineWidth(2);
        hPcalEPAll->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDPhgcerNphe = ddisplay->CreateDisplayPlot(
      "PID", "P.hgcer.nphe",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPhgcerNpheNoCuts->SetLineColor(kGreen + 2);
        hPhgcerNpheNoCuts->SetLineWidth(2);
        hPhgcerNpheNoCuts->GetYaxis()->SetRangeUser(.8, hPhgcerNpheNoCuts->GetMaximum() * 1.2);
        hPhgcerNpheNoCuts->DrawClone();
        hPhgcerNpheTracking->SetLineColor(kMagenta + 2);
        hPhgcerNpheTracking->SetLineWidth(2);
        hPhgcerNpheTracking->DrawClone("same");
        hPhgcerNphePID->SetLineColor(kBlue + 1);
        hPhgcerNphePID->SetLineWidth(2);
        hPhgcerNphePID->DrawClone("same");
        hPhgcerNpheAll->SetLineColor(kBlack);
        hPhgcerNpheAll->SetLineWidth(2);
        hPhgcerNpheAll->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  auto PIDPcerNphe = ddisplay->CreateDisplayPlot(
      "PID", "P.ngcer.nphe",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        plt.SetPersist();
        c->SetLogy();
        hPcerNpheNoCuts->SetLineColor(kGreen + 2);
        hPcerNpheNoCuts->SetLineWidth(2);
        hPcerNpheNoCuts->GetYaxis()->SetRangeUser(.8, hPcerNpheNoCuts->GetMaximum() * 1.2);
        hPcerNpheNoCuts->DrawClone();
        hPcerNpheTracking->SetLineColor(kMagenta + 2);
        hPcerNpheTracking->SetLineWidth(2);
        hPcerNpheTracking->DrawClone("same");
        hPcerNphePID->SetLineColor(kBlue + 1);
        hPcerNphePID->SetLineWidth(2);
        hPcerNphePID->DrawClone("same");
        hPcerNpheAll->SetLineColor(kBlack);
        hPcerNpheAll->SetLineWidth(2);
        hPcerNpheAll->DrawClone("same");
        c->BuildLegend();
        return 0;
      },
      [](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/good_coin_counter/";

  out_file->Write();

  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
