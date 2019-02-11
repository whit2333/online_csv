#include "nlohmann/json.hpp"
#include <chrono>
#include <cmath>
#include <cstdlib>
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

#include "THaPostProcess.h"
#include "hcana/HallC_Data.h"
R__LOAD_LIBRARY(libHallC.so)

#include "monitor/DisplayPlots.h"
#include "monitor/DetectorDisplay.h"
#include "monitor/MonitoringDisplay.h"
//#include "monitor/ExperimentMonitor.h"
//#include "scandalizer/PostProcessors.h"
R__LOAD_LIBRARY(libScandalizer.so)

#include "THStack.h"
using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// VecOps::RVec is like std::vector with some extra bells and whistles
using inters   = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;
using nlohmann::json;

void plot_waveforms(int RunNumber = 7111, int nevents = -1) {

  //using nlohmann::json;
  //json j;
  //{
  //  std::ifstream json_input_file("db2/test.json");
  //  try {
  //    json_input_file >> j;
  //  } catch (json::parse_error) {
  //    std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
  //    std::quick_exit(-127);
  //  }
  //}

  //auto runnum_str = std::to_string(RunNumber);
  //if (j.find(runnum_str) == j.end()) {
  //  std::cout << "Run " << RunNumber << " not found in ddb2/run_list.json\n";
  //  std::cout << "Check that run number and replay exists. \n";
  //  std::cout << "If problem persists please contact Whit (717-341-1080)\n";
  //  std::cout << "In the meantime use: good_coin_counter_old.cxx \n";
  //}
  //double P0_shms_setting = j[runnum_str]["spectrometers"]["shms_momentum"].get<double>();
  //double P0_shms         = std::abs(P0_shms_setting);

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
      std::cout << " using : " << rootfile << "\n";
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
        std::cout << " using : " << rootfile << "\n";
      }
    }
  }
  if (!found_good_file) {
    std::cout << " Error: suitable root file not found\n";
    return;
  }
  // new TBrowser;
  //

  // ROOT::EnableImplicitMT(24);
  // Detector tree
  ROOT::RDataFrame d("T", rootfile);
  // HMS Scaler tree
  //ROOT::RDataFrame d_sh("TSH", rootfile);

  //auto bcm4b_charge   = d_sh.Max("H.BCM4B.scalerChargeCut");
  //auto el_real_scaler = d_sh.Max("H.hEL_REAL.scaler");
  //auto time_1MHz      = d_sh.Max("H.1MHz.scalerTime");
  //auto time_1MHz_cut  = d_sh.Max("H.1MHz.scalerTimeCut");
  //auto total_charge   = bcm4b_charge;

  std::string hpdelta = "P.gtr.dp > -10 && P.gtr.dp < 20 && "
                        "H.gtr.dp > -10 && H.gtr.dp < 10";

  auto ddisplay = new hallc::MonitoringDisplay(RunNumber);
  // auto exp_monitor = new hallc::ExperimentMonitor(RunNumber,"derp");
  auto waveform_plot = ddisplay->CreateDisplayPlot(
      "test", "A",
      [&](hallc::DisplayPlot& plt) {
        plt._plot_data._canvas  = new TCanvas(plt.GetName().c_str(), plt.GetName().c_str());
        plt._plot_data._graphs1 = {new TGraph(), new TGraph(), new TGraph(), new TGraph()};
        for(int iwf = 0 ; iwf<4 ; iwf++) {
            plt._plot_data._graphs1[iwf]->SetPoint(0, 0, iwf);
        }
        plt._plot_data._canvas->cd();
        plt._plot_data._graphs1[0]->Draw("al");
        plt._plot_data._graphs1[1]->SetLineColor(2);
        plt._plot_data._graphs1[1]->Draw("l");
        plt._plot_data._graphs1[2]->SetLineColor(4);
        plt._plot_data._graphs1[2]->Draw("l");
        plt._plot_data._graphs1[3]->SetLineColor(8);
        plt._plot_data._graphs1[3]->Draw("l");
        return 0;
      },
      [&](hallc::DisplayPlot& plt) {

plt._plot_data._graphs1[0]->Set(0);
plt._plot_data._graphs1[1]->Set(0);
plt._plot_data._graphs1[2]->Set(0);
plt._plot_data._graphs1[3]->Set(0);
plt._plot_data._graphs1[0]->Set(1);
plt._plot_data._graphs1[1]->Set(1);
plt._plot_data._graphs1[2]->Set(1);
plt._plot_data._graphs1[3]->Set(1);

        // do nothing
        plt._plot_data._canvas->cd();
        plt._plot_data._canvas->Clear();
        plt._plot_data._graphs1[0]->Draw("al");
        plt._plot_data._graphs1[1]->SetLineColor(2);
        plt._plot_data._graphs1[1]->Draw("l");
        plt._plot_data._graphs1[2]->SetLineColor(4);
        plt._plot_data._graphs1[2]->Draw("l");
        plt._plot_data._graphs1[3]->SetLineColor(8);
        plt._plot_data._graphs1[3]->Draw("l");
        return 0;
      });
  ddisplay->_data._root_folder = "/derp/";

  ddisplay->InitAll();

  auto d1 = d.Range(2005,0);

  d1.Foreach(
    [&](std::vector<hallc::data::PulseWaveForm>& wf_data) {
      ddisplay->Process();
      static int wf_count = 0;
      int iwf = 0;
        std::cout << wf_data.size() << " wfs\n";
        if(wf_data.size()  > 0 ) {
          for (hallc::data::PulseWaveForm& pwf : wf_data) {
            int isamp = 0;
            for (auto samp : pwf._buffer) {
              //std::cout <<  isamp << " sample\n";
              waveform_plot->_plot_data._graphs1[iwf]->SetPoint(isamp, isamp, samp);
              isamp++;
            }
            iwf++;
            if(iwf>3) break;
          }
          wf_count++;
        }
        if( wf_count >= 1) {
          ddisplay->UpdateAll();
          wf_count = 0;
          std::this_thread::sleep_for(std::chrono::seconds(2));
        }
        return 0;
      },
      {"P_ngcer_waveforms"});
}
