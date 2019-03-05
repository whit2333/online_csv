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

// =================================================================================
// Constants
// =================================================================================
constexpr const double M_P  = .938272;
constexpr const double M_P2 = M_P * M_P;
constexpr const double M_J  = 3.096916;
constexpr const double M_J2 = M_J * M_J;
constexpr const double M_e  = .000511;
constexpr const double M_mu = .1056583745;

// =================================================================================
// Cuts -- Looser than the standard analysis cuts
// =================================================================================
std::string goodTrackSHMS = "P.gtr.dp > -15 && P.gtr.dp < 25";
std::string goodTrackHMS  = "H.gtr.dp > -11 && H.gtr.dp < 11";
std::string muCutSHMS     = "P.cal.etottracknorm < 0.4 && P.ngcer.npeSum > -1";
std::string muCutHMS      = "H.cal.etottracknorm < 0.3 && H.cer.npeSum > -1";
std::string coinCut       = "std::abs(coin_time - 42.5) < 2";

// =================================================================================
// Definitions
// =================================================================================
using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// =================================================================================
// J/psi reconstruction --> apply central momentum corrections from Mark
// =================================================================================
auto p_mumin = [](double px, double py, double pz) {
  return Pvec4D{px * 0.996, py * 0.996, pz * 0.996, M_mu};
};
auto p_muplus = [](double px, double py, double pz) {
  return Pvec4D{px * 0.994, py * 0.994, pz * 0.994, M_mu};
};
auto p_jpsi  = [](const Pvec4D& e1, const Pvec4D& e2) { return e1 + e2; };
auto E_gamma = [](const Pvec4D& jpsi) {
  double res =
      (jpsi.M2() - 2. * jpsi.E() * M_P) / (2. * (jpsi.E() - M_P - jpsi.P() * cos(jpsi.Theta())));
  return res;
};
auto abst = [](const double Egamma, Pvec4D& jpsi) {
  Pvec4D beam{0, 0, Egamma, 0};
  return fabs((beam - jpsi).M2());
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

void muskim(int RunNumber = 7146, int nevents = -1) {

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

  std::string coda_type = "COIN";

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
    rootfile = fmt::format("ROOTfiles_jpsi/coin_replay_production_{}_{}.root", RunNumber, nevents);
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

  ROOT::EnableImplicitMT(8);

  //---------------------------------------------------------------------------
  // Detector tree
  ROOT::RDataFrame d("T", rootfile);

  // SHMS Scaler tree
  ROOT::RDataFrame d_sh("TSP", rootfile);

  auto d_coin = d.Filter("fEvtHdr.fEvtType == 4");

  // Good track cuts
  auto dHMSGoodTrack  = d_coin.Filter(goodTrackHMS);
  auto dCOINGoodTrack = dHMSGoodTrack.Filter(goodTrackSHMS)
                            .Define("p_mumin", p_mumin, {"P.gtr.px", "P.gtr.py", "P.gtr.pz"})
                            .Define("p_muplus", p_muplus, {"H.gtr.px", "H.gtr.py", "H.gtr.pz"})
                            .Define("p_jpsi", p_jpsi, {"p_mumin", "p_muplus"})
                            .Define("M_jpsi", "p_jpsi.M()")
                            .Define("E_gamma", E_gamma, {"p_jpsi"})
                            .Define("abst", abst, {"E_gamma", "p_jpsi"})
                            .Define("coin_time", "CTime.ePositronCoinTime_ROC2");

  // PID cuts
  auto dCOINEl = dCOINGoodTrack.Filter(muCutHMS + " && " + muCutSHMS + " && " + coinCut);

  // scalers
  auto   total_charge      = d_sh.Max("P.BCM4B.scalerChargeCut");
  double good_total_charge = *total_charge / 1000.0; // mC

  // ===============================================================================================
  // Create output
  // ===============================================================================================
  std::string ofname = "results/skim/skim_coinmu_" + std::to_string(RunNumber) + ".root";
  dCOINEl.Snapshot("Tjpsi", ofname,
                   {"P.gtr.dp",
                    "P.gtr.th",
                    "P.gtr.ph",
                    "P.gtr.y",
                    "P.gtr.x",
                    "P.cal.eprtrack",
                    "P.cal.eprtracknorm",
                    "P.cal.etot",
                    "P.cal.etotnorm",
                    "P.cal.etottracknorm",
                    "P.cal.etrack",
                    "P.cal.etracknorm",
                    "P.ngcer.npeSum",
                    "H.gtr.dp",
                    "H.gtr.th",
                    "H.gtr.ph",
                    "H.gtr.y",
                    "H.gtr.x",
                    "H.cal.1pr.eplane",
                    "H.cal.2ta.eplane",
                    "H.cal.3ta.eplane",
                    "H.cal.4ta.eplane",
                    "H.cal.eprtrack",
                    "H.cal.eprtracknorm",
                    "H.cal.etot",
                    "H.cal.etotnorm",
                    "H.cal.etottracknorm",
                    "H.cal.etrack",
                    "H.cal.etracknorm",
                    "H.cer.npeSum",
                    "coin_time",
                    "p_mumin",
                    "p_muplus",
                    "p_jpsi",
                    "M_jpsi",
                    "E_gamma",
                    "abst"

                   });

  // also write total charge
  TFile ofile(ofname.c_str(), "update");
  TH1D* ocharge = new TH1D("total_good_charge", "total_good_charge", 1, 0, 1);
  ocharge->SetBinContent(1, good_total_charge);
  ocharge->Write();
  ofile.Close();
}
