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

#include "THaPostProcess.h"
#include "hcana/HallC_Data.h"
R__LOAD_LIBRARY(libHallC.so)

#include "monitor/DetectorDisplay.h"
#include "monitor/DisplayPlots.h"
#include "monitor/MonitoringDisplay.h"
//#include "monitor/ExperimentMonitor.h"
//#include "scandalizer/PostProcessors.h"
R__LOAD_LIBRARY(libScandalizer.so)


using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// VecOps::RVec is like std::vector with some extra bells and whistles
using inters   = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;
using nlohmann::json;

void muon_pion_separation(int RunNumber = 7146, int nevents = -1, int prompt = 0, int update = 1) {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list_shms.json");
    try {
      json_input_file >> j;
    } catch (json::parse_error) {
      std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }

  auto runnum_str = std::to_string(RunNumber);
  if (j.find(runnum_str) == j.end()) {
    std::cout << "Run " << RunNumber << " not found in db2/run_list_shms.json\n";
    std::cout << "Check that run number and replay exists. \n";
    std::cout << "If problem persists please contact Whit (717-341-1080)\n";
  }
  double P0_shms_setting = j[runnum_str]["spectrometers"]["shms_momentum"].get<double>();
  double P0_shms         = std::abs(P0_shms_setting);

  //double P0_hms_setting = j[runnum_str]["spectrometers"]["hms_momentum"].get<double>();
  //double P0_hms         = std::abs(P0_hms_setting);

  // int ps2 = 0;
  // if (j[runnum_str].find("daq") != j[runnum_str].end()) {
  //  ps2 = j[runnum_str]["daq"]["ps2"].get<int>();
  //} else {
  //  std::cout << " using default ps2 = 0 \n";
  //}
  ////  The way the input rates are prescaled follows:
  ////       input-rate/(2^{val - 1} + 1)
  // double singles_ps_value = std::pow(2.0, ps2);
  // std::cout << "prescale value " << singles_ps_value << "\n";

  std::string coda_type = "SHMS";
  std::string rootfile  = "ROOTfiles/";
  rootfile += std::string("shms_replay_production_");
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

  Pvec4D Pbeam(0, 0, 10.598, 0.000511);

  //---------------------------------------------------------------------------
  // Detector tree
  ROOT::RDataFrame d("T", rootfile);

  std::string hpdelta = "P.gtr.dp > -10 && P.gtr.dp < 20";

  double pion_threshold = 5.2;
  double muon_threshold = 4.7;

  std::cout << " Pion threshold : " << pion_threshold << "\n";
  std::cout << " Muon threshold : " << muon_threshold << "\n";
  std::cout << " P0 shms        : " << P0_shms << "\n";


  auto d0 = d.Filter(hpdelta).Filter(
      [=](double dp) {
        double p_track = P0_shms * (100.0 + dp) / 100.0;
        return (p_track > muon_threshold) && (p_track < pion_threshold);
      },
      {"P.gtr.dp"}).Filter("P.hgcer.npeSum > 3.0");

  auto d_muons = d0.Filter("P.ngcer.npeSum > 3.0");
  auto d_pions = d0.Filter("P.ngcer.npeSum < 0.5");

  // Apply the electron cuts
  //auto d_spec_cuts           = d0.Filter(hpdelta).Filter(epiCut);
  //auto d_spec_cuts_coin_only = d_spec_cuts.Filter("fEvtHdr.fEvtType == 6");
  //auto d_spec_cuts_shms_only = d_spec_cuts.Filter("fEvtHdr.fEvtType == 1");

  auto h_EOverP_0 = d_muons.Histo1D<doublers>( {"shms_e_EoverP_0", "SHMS total shower muons; SHMS E/P", 200, 0.001, 1.3}, "P.cal.etottracknorm");
  auto h_EOverP_2 = d_pions.Histo1D<doublers>( {"shms_e_EoverP_2", "SHMS total shower pions; SHMS E/P", 200, 0.001, 1.3}, "P.cal.etottracknorm");

  //auto h_h_EOverP_0 = d_muons.Histo1D<doublers>( {"hms_e_EoverP_0", "HMS total shower muons; HMS E/P", 200, 0.001, 1.3}, "H.cal.etottracknorm");
  //auto h_h_EOverP_2 = d_pions.Histo1D<doublers>( {"hms_e_EoverP_2", "HMS total shower pions; HMS E/P", 200, 0.001, 1.3}, "H.cal.etottracknorm");

  auto h_EprOverP_0 = d_muons.Histo1D<doublers>( {"shms_e_EoverP_0", "SHMS pre-shower muons; SHMS E/P", 200, 0.001, 1.3}, "P.cal.eprtracknorm");
  auto h_EprOverP_2 = d_pions.Histo1D<doublers>( {"shms_e_EoverP_2", "SHMS pre-shower pions; SHMS E/P", 200, 0.001, 1.3}, "P.cal.eprtracknorm");

  //auto h_h_EprOverP_0 = d_muons.Histo1D<doublers>( {"hms_e_EoverP_0", "HMS pre-shower muons; HMS E/P", 200, 0.001, 1.3}, "H.cal.eprtracknorm");
  //auto h_h_EprOverP_2 = d_pions.Histo1D<doublers>( {"hms_e_EoverP_2", "HMS pre-shower pions; HMS E/P", 200, 0.001, 1.3}, "H.cal.eprtracknorm");

  auto h_ngc_0 = d_muons.Histo1D<doublers>( {"shms_ngcer", "SHMS NGC muons; npe", 100, 0, 30}, "P.ngcer.npeSum");
  auto h_ngc_1 = d_pions.Histo1D<doublers>( {"shms_ngcer2","SHMS NGC pions; npe", 100, 0, 30}, "P.ngcer.npeSum");

  auto h_hgc_0 = d_muons.Histo1D<doublers>( {"shms_hgcer", "SHMS HGC muons; npe", 100, 0, 30}, "P.hgcer.npeSum");
  auto h_hgc_1 = d_pions.Histo1D<doublers>({"shms_hgcer2", "SHMS HGC pions; npe", 100, 0, 30}, "P.hgcer.npeSum");

  //auto h_h_cer_0 = d_muons.Histo1D<doublers>( {"hms_cer", "HMS cer muons; npe", 100, 0, 30}, "H.cer.npeSum");
  //auto h_h_cer_1 = d_pions.Histo1D<doublers>( {"hms_cer2","HMS cer pions; npe", 100, 0, 30}, "H.cer.npeSum");

  auto muon_counts = d_muons.Count();
  auto pion_counts = d_pions.Count();
  std::cout << *muon_counts << " muons\n";
  std::cout << *pion_counts << " pions\n";
  // -----------------------------------------------------------
  //
  TCanvas* c    = nullptr;
  int      b1   = 0;
  int      b2   = 0;
  double   hmax = 0.0;
  THStack* hs   = nullptr;
  // TLatex latex;

  gSystem->mkdir("results/good_shms_counter", true);

  auto ddisplay      = new hallc::MonitoringDisplay(RunNumber);
  auto waveform_plot = ddisplay->CreateDisplayPlot(
      "muon_pion_sep", "calorimeter",
      [&](hallc::DisplayPlot& plt) {
        c = plt.SetCanvas( new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        //c = plt._plot_data._canvas; // new TCanvas("c1", "c1", 1600, 1200);
        c->Divide(2, 2);
        c->cd(1);

        // This call starts the loop over the data.
        // A DrawCopy is used so that the histogram is not deleted at the end of
        // scope, and thus stays visible on the canvas.
        h_EOverP_0->DrawCopy();
        h_EOverP_2->SetLineColor(2);
        h_EOverP_2->DrawCopy("same");


        c->cd(2);
        // h_beta_0->DrawCopy();

        h_EprOverP_0->DrawCopy();
        h_EprOverP_2->SetLineColor(2);
        h_EprOverP_2->DrawCopy("same");

        //c->cd(3);
        //h_h_EOverP_0->DrawCopy();
        //h_h_EOverP_2->SetLineColor(2);
        //h_h_EOverP_2->DrawCopy("same");


        //c->cd(4);
        //// h_beta_0->DrawCopy();
        //h_h_EprOverP_0->DrawCopy();
        //h_h_EprOverP_2->SetLineColor(2);
        //h_h_EprOverP_2->DrawCopy("same");
        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });

  auto plot2 = ddisplay->CreateDisplayPlot(
      "muon_pion_sep", "cherenkovs",
      [&](hallc::DisplayPlot& plt) {
        c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        c->Divide(2, 2);
        c->cd(1);
        // gPad->SetLogy(true);
        h_ngc_0->DrawCopy();
        h_ngc_1->SetLineColor(2);
        h_ngc_1->DrawCopy("same");
        gPad->BuildLegend();

        //c->cd(2);
        //h_hgc_0->DrawCopy();
        //h_hgc_1->SetLineColor(2);
        //h_hgc_1->DrawCopy("same");
        //gPad->BuildLegend();

        //c->cd(3);
        //// gPad->SetLogy(true);
        //h_h_cer_0->DrawCopy();
        //h_h_cer_1->SetLineColor(2);
        //h_h_cer_1->DrawCopy("same");
        //gPad->BuildLegend();

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });

  ddisplay->_data._root_folder = "/whit/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();

  // gSystem->mkdir("results/muon_pion_separation", true);
  // c->SaveAs((std::string("results/muon_pion_separation/c1_") + std::to_string(RunNumber) +
  // ".pdf").c_str()); c->SaveAs((std::string("results/muon_pion_separation/c1_") +
  // std::to_string(RunNumber) + ".png").c_str());
}
