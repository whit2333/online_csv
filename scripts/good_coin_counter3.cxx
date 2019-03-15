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
using DFIndex = std::vector<std::pair<std::string, RDFInfo>>;

constexpr const double M_P     = .938272;
constexpr const double M_P2    = M_P * M_P;
constexpr const double M_pion  = 0.139;
constexpr const double M_pion2 = M_pion * M_pion;
constexpr const double M_e     = .000511;

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
auto p_q = [](Pvec4D& pe ) {
  return Pvec4D{0.0,0.0,10.219, M_e}-pe;
};
auto t = [](const double Egamma, Pvec4D& jpsi) {
  Pvec4D beam{0, 0, Egamma, 0};
  return (beam - jpsi).M2();
};
auto z = [](Pvec4D& pq, Pvec4D& ph) {
  return ph.E()/pq.E();
};
auto xbj = [](double Q2,Pvec4D& pq) {
  return Q2/(2.0*0.938*pq.E());
};
auto Q2 = [](Pvec4D& pq) {
  return -1.0*(pq.Dot(pq));
};
auto Wprime2 = [](Pvec4D& pq,Pvec4D& ph) {
  auto Ptot = Pvec4D{0.0,0.0,0.0, M_P} + pq - ph;
  return Ptot.Dot(Ptot);
};
auto W2 = [](Pvec4D& pq) {
  auto Ptot = Pvec4D{0.0,0.0,0.0, M_P} + pq;
  return Ptot.Dot(Ptot);
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

void good_coin_counter3(int RunNumber = 7146, int nevents = -1, int prompt = 0, int update = 1,
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
  std::string hgc_cut = " (p_pion.P() < 2.8 || P.hgcer.npeSum > 1.0) ";//&& P.ngcer.npeSum < 0.5 ";

  // ===============================================================================================
  // Dataframe
  // ===============================================================================================

  ROOT::EnableImplicitMT(24);

  // Detector tree
  ROOT::RDataFrame d("T", rootfile);

  // SHMS Scaler tree
  ROOT::RDataFrame d_sh("TSP", rootfile);

  auto d_coin = d.Filter("fEvtHdr.fEvtType == 4");

  // Good "track" cuts
  auto dHMSGoodTrack  = d_coin.Filter(goodTrackHMS);
  auto dSHMSGoodTrack = d_coin.Filter(goodTrackSHMS);
  auto dCOINGoodTrack = dHMSGoodTrack.Filter(goodTrackSHMS)
                            .Define("p_electron", p_electron, {"H.gtr.px", "H.gtr.py", "H.gtr.pz"})
                            .Define("p_pion", p_pion, {"P.gtr.px", "P.gtr.py", "P.gtr.pz"})
                            .Define("p_pion_HMS", p_pion, {"H.gtr.px", "H.gtr.py", "H.gtr.pz"})
                            .Define("p_q", p_q, {"p_electron"})
                            .Define("z", z, {"p_q","p_pion"})
                            .Define("Q2", Q2, {"p_q"})
                            .Define("xbj", xbj, {"Q2","p_q"})
                            .Define("W2", W2, {"p_q"})
                            .Define("Wp2", Wprime2, {"p_q","p_pion"})
                            .Define("W", "std::sqrt(W2)")
                            .Define("Wp", "std::sqrt(Wp2)")
                            .Define("InvMass","p_electron.Dot(p_pion)")
                            .Define("InvMass_pions","p_pion_HMS.Dot(p_pion)")
                            ;

  // PID cuts
  auto dHMS_electron    = dHMSGoodTrack.Filter(eCutHMS);
  auto dSHMS_pion       = dSHMSGoodTrack.Filter(piCutSHMS);
  auto dCOIN_sidis      = dCOINGoodTrack.Filter(eCutHMS + " && " + piCutSHMS);
  auto dCOIN_sidis_pion = dCOIN_sidis.Filter(hgc_cut);
      //[=](double npe, double dp) {
      //  double p_track = P0_shms * (100.0 + dp) / 100.0;
      //  // no cerenkov cut needed when momentum is below 2.8 GeV/c
      //  if (p_track < 2.8) {
      //    return true;
      //  }
      //  return npe > 1.0;
      //},
      //{"P.hgcer.npeSum", "P.gtr.dp"});

  // --------------------------------------------
  // Timing cuts
  // Find the timing peak
  // Find the coin peak
  double coin_peak_center = 0;
  if (redo_timing) {
    auto h_coin_time = dCOIN_sidis_pion.Histo1D({"coin_time", "coin_time", 8000, 0, 1000},
                                                "CTime.ePositronCoinTime_ROC2");
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
  //auto timing_cut = [=](double coin_time) { return std::abs(coin_time - coin_peak_center) < 2.; };
  //// anti-timing set to 5x width of regular
  //auto anti_timing_cut = [=](double coin_time) {
  //  return std::abs(coin_time - coin_peak_center - 28.) < 10.;
  //};

  std::string timing_cut = " std::abs(CTime.ePiCoinTime_ROC2 -43.0) < 2.0 ";
  std::string anti_timing_cut = " std::abs(CTime.ePiCoinTime_ROC2 -43.0 - 28.0) < 10.0 ";

  // timing counts
  auto dHMS_electron_intime  = dHMS_electron.Filter(timing_cut);//, {"CTime.ePiCoinTime_ROC2"});
  auto dHMS_electron_randbg  = dHMS_electron.Filter(anti_timing_cut);//, {"CTime.ePiCoinTime_ROC2"});
  auto dSHMS_electron_intime = dSHMS_pion.Filter(timing_cut);//, {"CTime.ePiCoinTime_ROC2"});
  auto dSHMS_electron_randbg = dSHMS_pion.Filter(anti_timing_cut);//, {"CTime.ePiCoinTime_ROC2"});

  auto dCOIN_sidis_pion_intime = dCOIN_sidis_pion.Filter(timing_cut);
  auto dCOIN_sidis_pion_randbg = dCOIN_sidis_pion.Filter(anti_timing_cut);

  // Output root file
  auto out_file =
      new TFile(fmt::format("monitoring/{}/good_csv_counter.root", RunNumber).c_str(), "UPDATE");
  out_file->cd();

  DFIndex dfs = {{"1_tracks", {dCOINGoodTrack, "Cuts: + Delta"}},
                 {"2_sidis_K_and_pi", {dCOIN_sidis, "Cuts: + e and hadron"}},
                 {"3_pion_sidis", {dCOIN_sidis_pion, "Cuts: + pion HGC"}},
                 {"4_timing", {dCOIN_sidis_pion_intime, "Cuts: + timing"}},
                 {"5_rand_bg", {dCOIN_sidis_pion_randbg, "Cuts: random bg coin"}}};

  // =========================================================================================
  // Histograms
  // =========================================================================================
  using HistoMap    = std::map<std::string, Histo1DProxy>;
  using HistoMapMap = std::map<std::string, HistoMap>;
  HistoMapMap histos;

  for (auto& kval : dfs) {
    std::string name{kval.first};
    RDFInfo&    df_info{kval.second};
    std::string title{df_info.title};
    auto&       df{df_info.df};

    histos["H.gtr.dp"][name] =
        df.Histo1D({("H.gtr.dp-" + name).c_str(), (title + ";#deltap_{HMS} [%];counts").c_str(),
                    200, -30., 40.},
                   "H.gtr.dp");
    histos["P.gtr.dp"][name] =
        df.Histo1D({("P.gtr.dp-" + name).c_str(), (title + ";#deltap_{SHMS} [%];counts").c_str(),
                    200, -30, 40},
                   "P.gtr.dp");
    histos["H.gtr.ph"][name] = df.Histo1D(
        {("H.gtr.ph-" + name).c_str(), (title + ";#phi_{HMS};counts").c_str(), 200, -0.1, 0.1},
        "H.gtr.ph");
    histos["P.gtr.ph"][name] = df.Histo1D(
        {("P.gtr.ph-" + name).c_str(), (title + ";#phi_{SHMS};counts").c_str(), 200, -0.1, 0.1},
        "P.gtr.ph");
    histos["H.gtr.th"][name] = df.Histo1D(
        {("H.gtr.th-" + name).c_str(), (title + ";#theta_{HMS};counts").c_str(), 200, -0.1, 0.1},
        "H.gtr.th");
    histos["P.gtr.th"][name] = df.Histo1D(
        {("P.gtr.th-" + name).c_str(), (title + ";#theta_{SHMS};counts").c_str(), 200, -0.1, 0.1},
        "P.gtr.th");
    histos["H.gtr.y"][name] = df.Histo1D(
        {("H.gtr.y-" + name).c_str(), (title + ";ytar_{HMS};counts").c_str(), 200, -10, 10},
        "H.gtr.y");
    histos["P.gtr.y"][name] = df.Histo1D(
        {("P.gtr.y-" + name).c_str(), (title + ";ytar_{SHMS};counts").c_str(), 200, -10, 10},
        "P.gtr.y");
    histos["H.cer.npeSum"][name] = df.Histo1D(
        {("H.cer.npeSum-" + name).c_str(), (title + ";HMS Cer #phe;counts").c_str(), 200, -1, 15},
        "H.cer.npeSum");
    histos["P.ngcer.npeSum"][name] =
        df.Histo1D({("P.ngcer.npeSum-" + name).c_str(), (title + ";SHMS NGC #phe;counts").c_str(),
                    200, -5, 75},
                   "P.ngcer.npeSum");
    histos["P.hgcer.npeSum"][name] =
        df.Histo1D({("P.hgcer.npeSum-" + name).c_str(), (title + ";SHMS HGC #phe;counts").c_str(),
                    200, -5, 75},
                   "P.hgcer.npeSum");
    histos["P.aero.npeSum"][name] =
        df.Histo1D({("P.aero.npeSum-" + name).c_str(), (title + ";SHMS aerogel #phe;counts").c_str(),
                    200, -5, 75},
                   "P.aero.npeSum");
    histos["H.cal.etottracknorm"][name] =
        df.Histo1D({("H.cal.etottracknorm-" + name).c_str(), (title + ";HMS E/P;counts").c_str(),
                    200, -.5, 1.5},
                   "H.cal.etottracknorm");
    histos["P.cal.etottracknorm"][name] =
        df.Histo1D({("P.cal.etottracknorm-" + name).c_str(), (title + ";SHMS E/P;counts").c_str(),
                    200, -.5, 1.5},
                   "P.cal.etottracknorm");
    histos["H.cal.eprtracknorm"][name] =
        df.Histo1D({("H.cal.eprtracknorm-" + name).c_str(),
                    (title + ";HMS Preshower;counts").c_str(), 200, -.5, 1.5},
                   "H.cal.eprtracknorm");
    histos["P.cal.eprtracknorm"][name] =
        df.Histo1D({("P.cal.eprtracknorm-" + name).c_str(),
                    (title + ";SHMS Preshower;counts").c_str(), 200, -.5, 1.5},
                   "P.cal.eprtracknorm");

     histos["coin.time"][name] =
        df.Histo1D({("coin_time" + name).c_str(), (title + ";coin time [ns];counts").c_str(),
                    600, 0,200},
                   "CTime.ePiCoinTime_ROC2");
     histos["coin.z"][name] =
        df.Histo1D({("z" + name).c_str(), (title + ";z;counts").c_str(),
                    100, 0,1},
                   "z");
     histos["coin.xbj"][name] =
        df.Histo1D({("xbj" + name).c_str(), (title + ";x;counts").c_str(),
                    100, 0,1},
                   "xbj");
     histos["coin.Q2"][name] =
        df.Histo1D({("Q2" + name).c_str(), (title + ";Q2;counts").c_str(),
                    200, 0.0,8.0},
                   "Q2");
     histos["coin.W2"][name] =
        df.Histo1D({("W2" + name).c_str(), (title + ";W2;counts").c_str(),
                    200, 0.0,8.0},
                   "W2");
     histos["coin.Wp2"][name] =
        df.Histo1D({("Wp2" + name).c_str(), (title + ";Wp2;counts").c_str(),
                    200, 0.0,8.0},
                   "Wp2");
     histos["coin.W"][name] =
        df.Histo1D({("W" + name).c_str(), (title + ";W;counts").c_str(),
                    200, 0.0,8.0},
                   "W");
     histos["coin.Wp"][name] =
        df.Histo1D({("Wp" + name).c_str(), (title + ";Wp;counts").c_str(),
                    200, 0.0,8.0},
                   "Wp");
     histos["coin.InvMass"][name] =
        df.Histo1D({("InvMass" + name).c_str(), (title + ";InvMass ;counts").c_str(),
                    200, 0.0,8.0},
                   "InvMass");
     histos["coin.InvMass_pions"][name] =
        df.Histo1D({("InvMass_pions" + name).c_str(), (title + ";InvMass_pions ;counts").c_str(),
                    200, 0.0,8.0},
                   "InvMass_pions");

    // J/psi invariant mass
    // histos["Jpsi.mass"][name] =
    //    df.Histo1D({("Jpsi_mass-" + name).c_str(), (title + ";M_{J/#psi} [GeV];counts").c_str(),
    //                100, 2.5, 3.5},
    //               "M_jpsi");
    // histos["Jpsi.Egamma"][name] = df.Histo1D(
    //    {("Jpsi_Egamma-" + name).c_str(), (title + ";E_{#gamma} [GeV];counts").c_str(), 100, 8,
    //    11}, "E_gamma");
    // histos["Jpsi.abst"][name] = df.Histo1D(
    //    {("Jpsi_abst-" + name).c_str(), (title + ";|t| [GeV^{2}];counts").c_str(), 100, -.5, 6.},
    //    "abst");
  }

  //// scalers
  // auto total_charge        = d_sh.Max("P.BCM4B.scalerChargeCut");
  // auto shms_el_real_scaler = d_sh.Max("P.pEL_REAL.scaler");
  // auto hms_el_real_scaler  = d_sh.Max("P.hEL_REAL.scaler");
  // auto time_1MHz           = d_sh.Max("P.1MHz.scalerTime");
  // auto time_1MHz_cut       = d_sh.Max("P.1MHz.scalerTimeCut");

  // auto yield_all = d.Count();
  //// 5 timing cut widths worth of random backgrounds
  // auto yield_coin          = d_coin.Count();
  // auto yield_HMSGoodTrack  = dHMSGoodTrack.Count();
  // auto yield_SHMSGoodTrack = dSHMSGoodTrack.Count();
  // auto yield_COINGoodTrack = dCOINGoodTrack.Count();
  // auto yield_HMSEl         = dHMSEl.Count();
  // auto yield_SHMSEl        = dSHMSEl.Count();
  // auto yield_COINEl        = dCOINEl.Count();
  // auto yield_HMSElInTime   = dHMSElInTime.Count();
  // auto yield_HMSElRandom   = dHMSElRandom.Count();
  // auto yield_SHMSElInTime  = dSHMSElInTime.Count();
  // auto yield_SHMSElRandom  = dSHMSElRandom.Count();
  // auto yield_COINElInTime  = dCOINElInTime.Count();
  // auto yield_COINElRandom  = dCOINElRandom.Count();
  // auto yield_coin_raw      = dCOINElInTime.Count();
  // auto yield_coin_random   = dCOINElRandom.Count();

  // -------------------------------------
  // End lazy eval
  // -------------------------------------
  // auto n_coin_good  = *yield_coin_raw - *yield_coin_random / 5.;
  // auto n_HMSElGood  = *yield_HMSElInTime - *yield_HMSElRandom / 5;
  // auto n_SHMSElGood = *yield_SHMSElInTime - *yield_SHMSElRandom / 5;
  // auto n_COINElGood = *yield_COINElInTime - *yield_COINElRandom / 5;
  ////// end of lazy eval
  ////auto n_jpsi = d_jpsi.Count();
  ////std::cout << "Found " << *n_jpsi << " J/psi candidates!\n";

  auto ddisplay = new hallc::MonitoringDisplay(RunNumber);
  for (auto& hskval : histos) {
    std::string full_name = hskval.first;
    auto&       hs        = hskval.second;
    // get display names
    std::string section_name = full_name.substr(0, full_name.find("."));
    if (section_name == "H") {
      section_name = "HMS";
    } else if (section_name == "P") {
      section_name = "SHMS";
    }
    std::string plot_name = full_name;
    // add display figure
    // std::cout << "Creating figure: " << section_name << " " << full_name << std::endl;
    auto figure = ddisplay->CreateDisplayPlot(
        section_name, plot_name,
        [&](hallc::DisplayPlot& plt) {
          auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
          plt.SetPersist();
          c->SetLogy();
          int              idx    = 0;
          std::vector<int> colors = {kGreen + 2, kMagenta + 2, kBlue + 1, kBlack, kRed + 1};
          for (auto& h : hs) {
            if (idx == 0) {
              h.second->SetTitle(plot_name.c_str());
            }
            h.second->GetYaxis()->SetRangeUser(.8, h.second->GetMaximum() * 1.5);
            h.second->SetLineColor(colors[idx]);
            h.second->SetLineWidth(2);
            h.second->DrawClone((idx == 0) ? "hist" : "histsame");
            h.second->Write();
            idx += 1;
          }
          c->BuildLegend();
          return 0;
        },
        [](hallc::DisplayPlot& plt) { return 0; });
  }
  ddisplay->_data._root_folder = "/csv_coin/";

  out_file->Write();
  ddisplay->InitAll();
  ddisplay->UpdateAll();

  // double good_total_charge = *total_charge / 1000.0; // mC
  // double good_time         = *time_1MHz_cut;         // s

  // std::map<std::string, double> counts = {
  //    {"shms_raw_yield_coin", (*yield_SHMSGoodTrack) / (good_total_charge)},
  //    {"hms_raw_yield_coin", (*yield_HMSGoodTrack) / (good_total_charge)},
  //    {"coin_raw_yield", (*yield_COINGoodTrack) / (good_total_charge)},
  //    {"shms_e_yield_coin", (*yield_SHMSEl) / (good_total_charge)},
  //    {"hms_e_yield_coin", (*yield_HMSEl) / (good_total_charge)},
  //    {"coin_ee_yield", (*yield_COINEl) / (good_total_charge)},
  //    {"hms_intime", (*yield_HMSElInTime) / (good_total_charge)},
  //    {"shms_intime", (*yield_SHMSElInTime) / (good_total_charge)},
  //    {"coin_intime", (*yield_COINElInTime) / (good_total_charge)},
  //    {"coin_intime count", (*yield_COINElInTime)},
  //    {"hms_random", (*yield_HMSElRandom) / (good_total_charge) / 5.},
  //    {"shms_random", (*yield_SHMSElRandom) / (good_total_charge) / 5.},
  //    {"coin_random", (*yield_COINElRandom) / (good_total_charge) / 5.},
  //    {"coin_random count", (*yield_COINElRandom) / 5.},
  //    {"coin_bg_corrected count", (*yield_COINElInTime) - (*yield_COINElRandom) / 5.},
  //    {"hms_e_good", (n_HMSElGood) / (good_total_charge)},
  //    {"shms_e_good", (n_SHMSElGood) / (good_total_charge)},
  //    {"coin_e_good", (n_COINElGood) / (good_total_charge)},
  //    {"coin count", n_COINElGood},
  //    {"coin yield", *yield_coin_raw / (good_total_charge)},
  //    {"coin random background count", *yield_coin_random / 5.},
  //    {"coin random background", *yield_coin_random / (good_total_charge) / 5.},
  //    {"coin Good event count", n_coin_good},
  //    {"coin Good event yield", n_coin_good / (good_total_charge)},
  //    {"good_total_charge", good_total_charge}};

  //// Update counts list
  // json jruns;
  //{
  //  std::ifstream input_file("db2/csv_run_count_list.json");
  //  try {
  //    input_file >> jruns;
  //  } catch (json::parse_error) {
  //    std::cerr << "error: json file is incomplete or has broken syntax.\n";
  //    std::quick_exit(-127);
  //  }
  //}
  // std::string run_str = std::to_string(RunNumber);
  // std::cout << "----------------------------------------------------------" << std::endl;
  // for (const auto& kv : counts) {
  //  std::cout << " " << kv.first;
  //  if (kv.first.find("yield") != std::string::npos) {
  //    std::cout << " (counts / mC)";
  //  }
  //  std::cout << ": " << kv.second;
  //  std::cout << "\n";
  //  jruns[run_str][kv.first] = kv.second;
  //}

  // jruns[run_str]["charge bcm4b 2u cut"] = good_total_charge;
  // jruns[run_str]["time 1MHz 2u cut"]    = good_time;

  // if (1) {
  //  std::cout << "Updating db2/coin_run_count_list.json with shms counts\n";
  //  std::ofstream json_output_file("db2/csv_run_count_list.json");
  //  json_output_file << std::setw(4) << jruns << "\n";
  //}

  // int count_goal = default_count_goal;

  //// if (prompt) {
  ////  std::cout << "----------------------------------------------------------\n";
  ////  std::cout << "Reference the run plan for this setting found on the wiki\n"
  ////               "       https://hallcweb.jlab.org/wiki/index.php/CSV_Fall_2018_Run_Plan\n";
  ////  std::cout << "----------------------------------------------------------\n";
  ////  std::cout << "Please enter **total count** goal for this setting. \n";
  ////  std::cout << "   Desired count [default=30000] : ";
  ////  // std::cin >> count_goal ;
  ////  // int number = 0;
  ////  if (std::cin.peek() == '\n') { // check if next character is newline
  ////    // count_goal = 30000; //and assign the default
  ////  } else if (!(std::cin >> count_goal)) { // be sure to handle invalid input
  ////    std::cout << "Invalid input.\n";
  ////    // error handling
  ////  }
  ////  std::cout << "\n";
  ////}

  // double n_seconds      = double(*time_1MHz_cut);
  // int    nev_tot        = (*yield_all);
  // double time_remaining = (count_goal * n_seconds) - (n_seconds);
  // double charge_goal    = count_goal * (good_total_charge) / counts["coin_bg_corrected count"];
  // double goal_Nevents   = double(nev_tot) * charge_goal / good_total_charge;

  // std::cout << " GOOD COUNTS : " << counts["coin_bg_corrected count"] << "\n";

  // std::cout << "----------------------------------------------------------\n";
  // std::cout << " N events to reach goal  : " << goal_Nevents / 1000000.0 << "M events\n";
  // std::cout << " Charge   to reach goal  : " << charge_goal << " mC\n";

  // if (update) {
  //  std::string cmd = "caput hcRunPlanChargeGoal " + std::to_string(charge_goal) + " &> /dev/null
  //  "; system(cmd.c_str());

  //  cmd = "caput hcRunPlanNTrigEventsGoal " + std::to_string(goal_Nevents) + " &> /dev/null ";
  //  system(cmd.c_str());

  //  cmd = "caput hcRunPlanCountGoal " + std::to_string(count_goal) + " &> /dev/null ";
  //  system(cmd.c_str());
  //}

  //// =====================================================================================
  //// Display
  //// =====================================================================================
  //// This is a naked pointer that we 'leak' on purpose so the connection stays alive
  //// as long as the root session is running
  // auto ddisplay    = new hallc::MonitoringDisplay(RunNumber);
  // auto TrackingPdp = ddisplay->CreateDisplayPlot(
  //    "Tracking", "P.gtr.dp",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hPdpNoCuts->SetLineColor(kGreen + 2);
  //      hPdpNoCuts->SetLineWidth(2);
  //      hPdpNoCuts->GetYaxis()->SetRangeUser(.8, hPdpNoCuts->GetMaximum() * 1.2);
  //      hPdpNoCuts->DrawClone();
  //      hPdpTracking->SetLineColor(kMagenta + 2);
  //      hPdpTracking->SetLineWidth(2);
  //      hPdpTracking->DrawClone("same");
  //      hPdpPID->SetLineColor(kBlue + 1);
  //      hPdpPID->SetLineWidth(2);
  //      hPdpPID->DrawClone("same");
  //      hPdpTiming->SetLineColor(kBlack);
  //      hPdpTiming->SetLineWidth(2);
  //      hPdpTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TrackingPth = ddisplay->CreateDisplayPlot(
  //    "Tracking", "P.gtr.th",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hPthNoCuts->SetLineColor(kGreen + 2);
  //      hPthNoCuts->SetLineWidth(2);
  //      hPthNoCuts->GetYaxis()->SetRangeUser(.8, hPthNoCuts->GetMaximum() * 1.2);
  //      hPthNoCuts->DrawClone();
  //      hPthTracking->SetLineColor(kMagenta + 2);
  //      hPthTracking->SetLineWidth(2);
  //      hPthTracking->DrawClone("same");
  //      hPthPID->SetLineColor(kBlue + 1);
  //      hPthPID->SetLineWidth(2);
  //      hPthPID->DrawClone("same");
  //      hPthTiming->SetLineColor(kBlack);
  //      hPthTiming->SetLineWidth(2);
  //      hPthTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TrackingPph = ddisplay->CreateDisplayPlot(
  //    "Tracking", "P.gtr.ph",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hPphNoCuts->SetLineColor(kGreen + 2);
  //      hPphNoCuts->SetLineWidth(2);
  //      hPphNoCuts->GetYaxis()->SetRangeUser(.8, hPphNoCuts->GetMaximum() * 1.2);
  //      hPphNoCuts->DrawClone();
  //      hPphTracking->SetLineColor(kMagenta + 2);
  //      hPphTracking->SetLineWidth(2);
  //      hPphTracking->DrawClone("same");
  //      hPphPID->SetLineColor(kBlue + 1);
  //      hPphPID->SetLineWidth(2);
  //      hPphPID->DrawClone("same");
  //      hPphTiming->SetLineColor(kBlack);
  //      hPphTiming->SetLineWidth(2);
  //      hPphTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TrackingPy = ddisplay->CreateDisplayPlot(
  //    "Tracking", "P.gtr.y",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hPyNoCuts->SetLineColor(kGreen + 2);
  //      hPyNoCuts->SetLineWidth(2);
  //      hPyNoCuts->GetYaxis()->SetRangeUser(.8, hPyNoCuts->GetMaximum() * 1.2);
  //      hPyNoCuts->DrawClone();
  //      hPyTracking->SetLineColor(kMagenta + 2);
  //      hPyTracking->SetLineWidth(2);
  //      hPyTracking->DrawClone("same");
  //      hPyPID->SetLineColor(kBlue + 1);
  //      hPyPID->SetLineWidth(2);
  //      hPyPID->DrawClone("same");
  //      hPyTiming->SetLineColor(kBlack);
  //      hPyTiming->SetLineWidth(2);
  //      hPyTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TrackingHdp = ddisplay->CreateDisplayPlot(
  //    "Tracking", "H.gtr.dp",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hHdpNoCuts->SetLineColor(kGreen + 2);
  //      hHdpNoCuts->SetLineWidth(2);
  //      hHdpNoCuts->GetYaxis()->SetRangeUser(.8, hHdpNoCuts->GetMaximum() * 1.2);
  //      hHdpNoCuts->DrawClone();
  //      hHdpTracking->SetLineColor(kMagenta + 2);
  //      hHdpTracking->SetLineWidth(2);
  //      hHdpTracking->DrawClone("same");
  //      hHdpPID->SetLineColor(kBlue + 1);
  //      hHdpPID->SetLineWidth(2);
  //      hHdpPID->DrawClone("same");
  //      hHdpTiming->SetLineColor(kBlack);
  //      hHdpTiming->SetLineWidth(2);
  //      hHdpTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TrackingHth = ddisplay->CreateDisplayPlot(
  //    "Tracking", "H.gtr.th",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hHthNoCuts->SetLineColor(kGreen + 2);
  //      hHthNoCuts->SetLineWidth(2);
  //      hHthNoCuts->GetYaxis()->SetRangeUser(.8, hHthNoCuts->GetMaximum() * 1.2);
  //      hHthNoCuts->DrawClone();
  //      hHthTracking->SetLineColor(kMagenta + 2);
  //      hHthTracking->SetLineWidth(2);
  //      hHthTracking->DrawClone("same");
  //      hHthPID->SetLineColor(kBlue + 1);
  //      hHthPID->SetLineWidth(2);
  //      hHthPID->DrawClone("same");
  //      hHthTiming->SetLineColor(kBlack);
  //      hHthTiming->SetLineWidth(2);
  //      hHthTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TrackingHph = ddisplay->CreateDisplayPlot(
  //    "Tracking", "H.gtr.ph",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hHphNoCuts->SetLineColor(kGreen + 2);
  //      hHphNoCuts->SetLineWidth(2);
  //      hHphNoCuts->GetYaxis()->SetRangeUser(.8, hHphNoCuts->GetMaximum() * 1.2);
  //      hHphNoCuts->DrawClone();
  //      hHphTracking->SetLineColor(kMagenta + 2);
  //      hHphTracking->SetLineWidth(2);
  //      hHphTracking->DrawClone("same");
  //      hHphPID->SetLineColor(kBlue + 1);
  //      hHphPID->SetLineWidth(2);
  //      hHphPID->DrawClone("same");
  //      hHphTiming->SetLineColor(kBlack);
  //      hHphTiming->SetLineWidth(2);
  //      hHphTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TrackingHy = ddisplay->CreateDisplayPlot(
  //    "Tracking", "H.gtr.y",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hHyNoCuts->SetLineColor(kGreen + 2);
  //      hHyNoCuts->SetLineWidth(2);
  //      hHyNoCuts->GetYaxis()->SetRangeUser(.8, hHyNoCuts->GetMaximum() * 1.2);
  //      hHyNoCuts->DrawClone();
  //      hHyTracking->SetLineColor(kMagenta + 2);
  //      hHyTracking->SetLineWidth(2);
  //      hHyTracking->DrawClone("same");
  //      hHyPID->SetLineColor(kBlue + 1);
  //      hHyPID->SetLineWidth(2);
  //      hHyPID->DrawClone("same");
  //      hHyTiming->SetLineColor(kBlack);
  //      hHyTiming->SetLineWidth(2);
  //      hHyTiming->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto TimingCoinTime = ddisplay->CreateDisplayPlot(
  //    "Timing", "coin_time",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hCoinTimeNoCuts->SetLineColor(kGreen + 2);
  //      hCoinTimeNoCuts->SetLineWidth(2);
  //      hCoinTimeNoCuts->GetYaxis()->SetRangeUser(.8, hCoinTimeNoCuts->GetMaximum() * 1.2);
  //      hCoinTimeNoCuts->DrawClone();
  //      hCoinTimeTracking->SetLineColor(kMagenta + 2);
  //      hCoinTimeTracking->SetLineWidth(2);
  //      hCoinTimeTracking->DrawClone("same");
  //      hCoinTimePID->SetLineColor(kBlue + 1);
  //      hCoinTimePID->SetLineWidth(2);
  //      hCoinTimePID->DrawClone("same");
  //      hCoinTimeTiming->SetLineColor(kBlack);
  //      hCoinTimeTiming->SetLineWidth(2);
  //      hCoinTimeTiming->DrawClone("same");
  //      hRandCoinTimePID->SetLineColor(kRed);
  //      hRandCoinTimePID->SetLineWidth(2);
  //      hRandCoinTimePID->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto PIDHcalEP = ddisplay->CreateDisplayPlot(
  //    "PID", "H.cal.etottracknorm",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hHcalEPNoCuts->SetLineColor(kGreen + 2);
  //      hHcalEPNoCuts->SetLineWidth(2);
  //      hHcalEPNoCuts->GetYaxis()->SetRangeUser(.8, hHcalEPNoCuts->GetMaximum() * 1.2);
  //      hHcalEPNoCuts->DrawClone();
  //      hHcalEPTracking->SetLineColor(kMagenta + 2);
  //      hHcalEPTracking->SetLineWidth(2);
  //      hHcalEPTracking->DrawClone("same");
  //      hHcalEPPID->SetLineColor(kBlue + 1);
  //      hHcalEPPID->SetLineWidth(2);
  //      hHcalEPPID->DrawClone("same");
  //      hHcalEPAll->SetLineColor(kBlack);
  //      hHcalEPAll->SetLineWidth(2);
  //      hHcalEPAll->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto PIDHcerNphe = ddisplay->CreateDisplayPlot(
  //    "PID", "H.cer.nphe",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hHcerNpheNoCuts->SetLineColor(kGreen + 2);
  //      hHcerNpheNoCuts->SetLineWidth(2);
  //      hHcerNpheNoCuts->GetYaxis()->SetRangeUser(.8, hHcerNpheNoCuts->GetMaximum() * 1.2);
  //      hHcerNpheNoCuts->DrawClone();
  //      hHcerNpheTracking->SetLineColor(kMagenta + 2);
  //      hHcerNpheTracking->SetLineWidth(2);
  //      hHcerNpheTracking->DrawClone("same");
  //      hHcerNphePID->SetLineColor(kBlue + 1);
  //      hHcerNphePID->SetLineWidth(2);
  //      hHcerNphePID->DrawClone("same");
  //      hHcerNpheAll->SetLineColor(kBlack);
  //      hHcerNpheAll->SetLineWidth(2);
  //      hHcerNpheAll->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto PIDPcalEP = ddisplay->CreateDisplayPlot(
  //    "PID", "P.cal.etottracknorm",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hPcalEPNoCuts->SetLineColor(kGreen + 2);
  //      hPcalEPNoCuts->SetLineWidth(2);
  //      hPcalEPNoCuts->GetYaxis()->SetRangeUser(.8, hPcalEPNoCuts->GetMaximum() * 1.2);
  //      hPcalEPNoCuts->DrawClone();
  //      hPcalEPTracking->SetLineColor(kMagenta + 2);
  //      hPcalEPTracking->SetLineWidth(2);
  //      hPcalEPTracking->DrawClone("same");
  //      hPcalEPPID->SetLineColor(kBlue + 1);
  //      hPcalEPPID->SetLineWidth(2);
  //      hPcalEPPID->DrawClone("same");
  //      hPcalEPAll->SetLineColor(kBlack);
  //      hPcalEPAll->SetLineWidth(2);
  //      hPcalEPAll->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto PIDPhgcerNphe = ddisplay->CreateDisplayPlot(
  //    "PID", "P.hgcer.nphe",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hPhgcerNpheNoCuts->SetLineColor(kGreen + 2);
  //      hPhgcerNpheNoCuts->SetLineWidth(2);
  //      hPhgcerNpheNoCuts->GetYaxis()->SetRangeUser(.8, hPhgcerNpheNoCuts->GetMaximum() * 1.2);
  //      hPhgcerNpheNoCuts->DrawClone();
  //      hPhgcerNpheTracking->SetLineColor(kMagenta + 2);
  //      hPhgcerNpheTracking->SetLineWidth(2);
  //      hPhgcerNpheTracking->DrawClone("same");
  //      hPhgcerNphePID->SetLineColor(kBlue + 1);
  //      hPhgcerNphePID->SetLineWidth(2);
  //      hPhgcerNphePID->DrawClone("same");
  //      hPhgcerNpheAll->SetLineColor(kBlack);
  //      hPhgcerNpheAll->SetLineWidth(2);
  //      hPhgcerNpheAll->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // auto PIDPcerNphe = ddisplay->CreateDisplayPlot(
  //    "PID", "P.ngcer.nphe",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      plt.SetPersist();
  //      c->SetLogy();
  //      hPcerNpheNoCuts->SetLineColor(kGreen + 2);
  //      hPcerNpheNoCuts->SetLineWidth(2);
  //      hPcerNpheNoCuts->GetYaxis()->SetRangeUser(.8, hPcerNpheNoCuts->GetMaximum() * 1.2);
  //      hPcerNpheNoCuts->DrawClone();
  //      hPcerNpheTracking->SetLineColor(kMagenta + 2);
  //      hPcerNpheTracking->SetLineWidth(2);
  //      hPcerNpheTracking->DrawClone("same");
  //      hPcerNphePID->SetLineColor(kBlue + 1);
  //      hPcerNphePID->SetLineWidth(2);
  //      hPcerNphePID->DrawClone("same");
  //      hPcerNpheAll->SetLineColor(kBlack);
  //      hPcerNpheAll->SetLineWidth(2);
  //      hPcerNpheAll->DrawClone("same");
  //      c->BuildLegend();
  //      return 0;
  //    },
  //    [](hallc::DisplayPlot& plt) { return 0; });
  // ddisplay->_data._root_folder = "/good_coin_counter/";

  // out_file->Write();

  // ddisplay->InitAll();
  // ddisplay->UpdateAll();
}
