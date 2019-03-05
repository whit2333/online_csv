#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TCanvas.h"
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libGenVector.so)

#include "monitor/DetectorDisplay.h"
#include "monitor/DisplayPlots.h"
#include "monitor/MonitoringDisplay.h"
//#include "monitor/ExperimentMonitor.h"
//#include "scandalizer/PostProcessors.h"
R__LOAD_LIBRARY(libScandalizer.so)

#include <fmt/core.h>
R__LOAD_LIBRARY(libfmt.so)

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
std::string track_cut = "P_gtr_dp > -10 && P_gtr_dp < 22"
                        "&& H_gtr_dp > -8 && H_gtr_dp < 8 ";
std::string cal_cut = "P_cal_etottracknorm > .08 && P_cal_etottracknorm < 0.25 "
                      "&& H_cal_etottracknorm < 0.1 && P_cal_eprtracknorm < .07"
                      "&& H_cal_etottracknorm < 0.1 && P_cal_eprtracknorm < .07"

    ;
// No cut on HMS E/p for now
// std::string cher_cut = "P_ngcer_npeSum > .5 && H_cer_npeSum > -1.";
std::string cher_cut = "P_ngcer_npeSum > -1.5 && H_cer_npeSum > 0.01 "
                       " && H_cal_1pr_eplanenorm < 0.04 "
                       " && H_cal_2ta_eplanenorm < 0.04 "
                       " && H_cal_3ta_eplanenorm < 0.04 "
                       " && H_cal_4ta_eplanenorm < 0.04 ";

std::string jpsi_cut = "M_jpsi > 3.07 && M_jpsi < 3.105";

// =================================================================================
// Definitions
// =================================================================================
using Pvec4D = ROOT::Math::PxPyPzMVector;

void muons(const std::string& setting = "phase1") {

  const std::string skim_file = fmt::format("results/skim/skim_coinmu_{}.root", setting);
  std::cout << "Looking for mu+mu- coincidence event in: " << skim_file << "\n";

  // ===============================================================================================
  // Dataframe
  // ===============================================================================================

  ROOT::EnableImplicitMT(24);

  //---------------------------------------------------------------------------
  // Detector tree
  ROOT::RDataFrame d0("Tjpsi", skim_file);
  auto             d = d0.Define("H_cal_1pr_eplanenorm", "H_cal_1pr_eplane/p_muplus.P()")
               .Define("H_cal_2ta_eplanenorm", "H_cal_2ta_eplane/p_muplus.P()")
               .Define("H_cal_3ta_eplanenorm", "H_cal_3ta_eplane/p_muplus.P()")
               .Define("H_cal_4ta_eplanenorm", "H_cal_4ta_eplane/p_muplus.P()");
  // cuts
  auto d_track = d.Filter(track_cut);
  auto d_cal   = d_track.Filter(cal_cut);
  auto d_cer   = d_cal.Filter(cher_cut);
  auto d_jpsi  = d_cer.Filter(jpsi_cut);

  // Data frame index
  std::vector<std::pair<std::string, RDFInfo>> dfs = {{"1. delta", {d_track, "Cuts: + Delta"}},
                                                      {"2. cal", {d_cal, "Cuts: + E/P"}},
                                                      {"3. cer", {d_cer, "Cuts: + Cherenkov"}},
                                                      {"4. jpsi", {d_jpsi, "Cuts: + J/psi"}}};

  // =========================================================================================
  // Histograms
  // =========================================================================================
  using HistoMap    = std::map<std::string, Histo1DProxy>;
  using HistoMapMap = std::map<std::string, HistoMap>;
#if 0
  HistoMap hP_dp;
  HistoMap hH_dp;
  HistoMap hP_ph;
  HistoMap hH_ph;
  HistoMap hP_th;
  HistoMap hH_th;
  HistoMap hP_y;
  HistoMap hH_y;
  HistoMap hP_cal_ep;
  HistoMap hP_cal_pre;
  HistoMap hH_cal_ep;
  HistoMap hP_cer;
  HistoMap hH_cer;
  HistoMap hP_jpsi_mass;
  HistoMap hH_jpsi_mass;
  HistoMap hP_jpsi_E;
  HistoMap hH_jpsi_E;
  HistoMap hP_jpsi_abst;
  HistoMap hH_jpsi_abst;
#endif
  HistoMapMap histos;

  for (auto& kval : dfs) {
    std::string name{kval.first};
    RDFInfo&    df_info{kval.second};
    std::string title{df_info.title};
    auto&       df{df_info.df};

    histos["H.gtr.dp"][name] =
        df.Histo1D({("H.gtr.dp-" + name).c_str(), (title + ";#deltap_{HMS} [%];counts").c_str(),
                    200, -30., 40.},
                   "H_gtr_dp");
    histos["P.gtr.dp"][name] =
        df.Histo1D({("P.gtr.dp-" + name).c_str(), (title + ";#deltap_{SHMS} [%];counts").c_str(),
                    200, -30, 40},
                   "P_gtr_dp");
    histos["H.gtr.ph"][name] = df.Histo1D(
        {("H.gtr.ph-" + name).c_str(), (title + ";#phi_{HMS};counts").c_str(), 200, -0.1, 0.1},
        "H_gtr_ph");
    histos["P.gtr.ph"][name] = df.Histo1D(
        {("P.gtr.ph-" + name).c_str(), (title + ";#phi_{SHMS};counts").c_str(), 200, -0.1, 0.1},
        "P_gtr_ph");
    histos["H.gtr.th"][name] = df.Histo1D(
        {("H.gtr.th-" + name).c_str(), (title + ";#theta_{HMS};counts").c_str(), 200, -0.1, 0.1},
        "H_gtr_th");
    histos["P.gtr.th"][name] = df.Histo1D(
        {("P.gtr.th-" + name).c_str(), (title + ";#theta_{SHMS};counts").c_str(), 200, -0.1, 0.1},
        "P_gtr_th");
    histos["H.gtr.y"][name] = df.Histo1D(
        {("H.gtr.y-" + name).c_str(), (title + ";ytar_{HMS};counts").c_str(), 200, -10, 10},
        "H_gtr_y");
    histos["P.gtr.y"][name] = df.Histo1D(
        {("P.gtr.y-" + name).c_str(), (title + ";ytar_{SHMS};counts").c_str(), 200, -10, 10},
        "P_gtr_y");
    histos["H.cer.npeSum"][name] = df.Histo1D(
        {("H.cer.npeSum-" + name).c_str(), (title + ";HMS Cer #phe;counts").c_str(), 200, -1, 15},
        "H_cer_npeSum");
    histos["P.cer.npeSum"][name] =
        df.Histo1D({("P.ngcer.npeSum-" + name).c_str(), (title + ";SHMS NGC #phe;counts").c_str(),
                    200, -5, 75},
                   "P_ngcer_npeSum");
    histos["H.cal.etottracknorm"][name] =
        df.Histo1D({("H.cal.etottracknorm-" + name).c_str(), (title + ";HMS E/P;counts").c_str(),
                    200, -.5, 1.5},
                   "H_cal_etottracknorm");
    histos["P.cal.etottracknorm"][name] =
        df.Histo1D({("P.cal.etottracknorm-" + name).c_str(), (title + ";SHMS E/P;counts").c_str(),
                    200, -.5, 1.5},
                   "P_cal_etottracknorm");
    histos["H.cal.eprtracknorm"][name] =
        df.Histo1D({("H.cal.eprtracknorm-" + name).c_str(),
                    (title + ";HMS Preshower;counts").c_str(), 200, -.5, 1.5},
                   "H_cal_eprtracknorm");
    histos["H.cal.1pr.eplane"][name] =
        df.Histo1D({("H_cal_1pr_eplane-" + name).c_str(),
                    (title + ";HMS cal layer 1;counts").c_str(), 200, -.5, 1.5},
                   "H_cal_1pr_eplane");
    histos["H.cal.2ta.eplane"][name] =
        df.Histo1D({("H_cal_2ta_eplane-" + name).c_str(),
                    (title + ";HMS cal layer 2;counts").c_str(), 200, -.5, 1.5},
                   "H_cal_2ta_eplane");
    histos["H.cal.3ta.eplane"][name] =
        df.Histo1D({("H_cal_3ta_eplane-" + name).c_str(),
                    (title + ";HMS cal layer 3;counts").c_str(), 200, -.5, 1.5},
                   "H_cal_3ta_eplane");
    histos["H.cal.4ta.eplane"][name] =
        df.Histo1D({("H_cal_4ta_eplane-" + name).c_str(),
                    (title + ";HMS cal layer 4;counts").c_str(), 200, -.5, 1.5},
                   "H_cal_4ta_eplane");

    histos["H.cal.1pr.eplanenorm"][name] =
        df.Histo1D({("H_cal_1pr_eplanenorm-" + name).c_str(),
                    (title + ";HMS cal layer 1;counts").c_str(), 200, -.5, 0.5},
                   "H_cal_1pr_eplanenorm");
    histos["H.cal.2ta.eplanenorm"][name] =
        df.Histo1D({("H_cal_2ta_eplanenorm-" + name).c_str(),
                    (title + ";HMS cal layer 2;counts").c_str(), 200, -.5, 0.5},
                   "H_cal_2ta_eplanenorm");
    histos["H.cal.3ta.eplanenorm"][name] =
        df.Histo1D({("H_cal_3ta_eplanenorm-" + name).c_str(),
                    (title + ";HMS cal layer 3;counts").c_str(), 200, -.5, 0.5},
                   "H_cal_3ta_eplanenorm");
    histos["H.cal.4ta.eplanenorm"][name] =
        df.Histo1D({("H_cal_4ta_eplanenorm-" + name).c_str(),
                    (title + ";HMS cal layer 4;counts").c_str(), 200, -.5, 0.5},
                   "H_cal_4ta_eplanenorm");
    histos["P.cal.eprtracknorm"][name] =
        df.Histo1D({("P.cal.eprtracknorm-" + name).c_str(),
                    (title + ";SHMS Preshower;counts").c_str(), 200, -.5, 1.5},
                   "P_cal_eprtracknorm");
    // J/psi invariant mass
    histos["Jpsi.mass"][name] =
        df.Histo1D({("Jpsi_mass-" + name).c_str(), (title + ";M_{J/#psi} [GeV];counts").c_str(),
                    100, 2.5, 3.5},
                   "M_jpsi");
    histos["Jpsi.Egamma"][name] = df.Histo1D(
        {("Jpsi_Egamma-" + name).c_str(), (title + ";E_{#gamma} [GeV];counts").c_str(), 100, 8, 11},
        "E_gamma");
    histos["Jpsi.abst"][name] = df.Histo1D(
        {("Jpsi_abst-" + name).c_str(), (title + ";|t| [GeV^{2}];counts").c_str(), 100, -.5, 6.},
        "abst");
  }

  // end of lazy eval
  auto n_jpsi = d_jpsi.Count();
  std::cout << "Found " << *n_jpsi << " J/psi candidates!\n";

  // =====================================================================================
  // Display
  // =====================================================================================
  // This is a naked pointer that we 'leak' on purpose so the connection stays
  // alive as long as the root session is running
  auto ddisplay = new hallc::MonitoringDisplay(-11);
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
            idx += 1;
          }
          c->BuildLegend();
          return 0;
        },
        [](hallc::DisplayPlot& plt) { return 0; });
  }
  ddisplay->_data._root_folder = "/muons-" + setting + "/";

  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
