#include <fmt/core.h>
#include <fmt/ostream.h>
R__LOAD_LIBRARY(libfmt.so)

#ifdef __cpp_lib_filesystem
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

#include "nlohmann/json.hpp"

#include "TGraph.h"
#include "TObject.h"

#include "THaPostProcess.h"
#include "hcana/HallC_Data.h"
R__LOAD_LIBRARY(libHallC.so)

#include "monitor/DetectorDisplay.h"
#include "monitor/DisplayPlots.h"
#include "monitor/MonitoringDisplay.h"
//#include "monitor/ExperimentMonitor.h"
//#include "scandalizer/PostProcessors.h"
R__LOAD_LIBRARY(libScandalizer.so)

#include "TBufferJSON.h"

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

void csv_ratios(int start_run = 0) {

  using nlohmann::json;
  json j2;
  {
    std::ifstream json_input_file("db2/csv_run_count_list.json");
    json_input_file >> j2;
  }
  double total_charge = 0.000;

  //// 35/36-1
  // std::vector<int> positive_runs = {7699, 7700, 7701, 7702, 7703};
  // std::vector<int> negative_runs = {7671, 7672};

  //// 35/36-3
  // std::vector<int> positive_runs = { 7706,7707,7708,7709,7710 };
  // std::vector<int> negative_runs = {7679,7680,7681,7682};

  // 35/36-5
  int              setting_pair  = 35365;
  std::vector<int> positive_runs = {7716, 7717, 7718, 7719, 7720, 7721, 7722};
  std::vector<int> negative_runs = {7684, 7685, 7686, 7687, 7688, 7689, 7690, 7691, 7692};

  // Dummy
  // std::vector<int> kine_0_runs       = { 7242, 7243, 7244};

  // total histograms
  TH1D* h_z_Total_positive = nullptr;
  TH1D* h_z_Total_negative = nullptr;

  double total_charge_positive = 0.00001;
  double total_charge_negative = 0.00001;

  for (auto run_number : positive_runs) {
    std::string file_name = fmt::format("monitoring/{}/good_csv_counter.root", run_number);
    if (!root_file_exists(file_name)) {
      continue;
    }
    if (j2[std::to_string(run_number)].find("good_total_charge") !=
        j2[std::to_string(run_number)].end()) {
      total_charge_positive += j2[std::to_string(run_number)]["good_total_charge"].get<double>();
    }

    auto hist_file = new TFile(file_name.c_str(), "READ");

    TH1D* z_run = (TH1D*)hist_file->Get("z4_timing");
    if (!h_z_Total_positive) {
      h_z_Total_positive = (TH1D*)z_run->Clone("z");
    } else {
      h_z_Total_positive->Add((TH1D*)z_run->Clone());
    }
  }

  for (auto run_number : negative_runs) {
    std::string file_name = fmt::format("monitoring/{}/good_csv_counter.root", run_number);
    if (!root_file_exists(file_name)) {
      continue;
    }
    if (j2[std::to_string(run_number)].find("good_total_charge") !=
        j2[std::to_string(run_number)].end()) {
      total_charge_negative += j2[std::to_string(run_number)]["good_total_charge"].get<double>();
    }

    auto hist_file = new TFile(file_name.c_str(), "READ");

    TH1D* z_run = (TH1D*)hist_file->Get("z4_timing");
    if (!h_z_Total_negative) {
      h_z_Total_negative = (TH1D*)z_run->Clone("z");
    } else {
      h_z_Total_negative->Add((TH1D*)z_run->Clone());
    }
  }

  h_z_Total_positive->Scale(1.0 / total_charge_positive);
  h_z_Total_negative->Scale(1.0 / total_charge_negative);

  auto h_z_ratio = (TH1D*)h_z_Total_negative->Clone("ratio");
  h_z_ratio->Divide(h_z_Total_positive);

  // ------------------------------
  //
  TGraphErrors* gr = nullptr;
  TMultiGraph*  mg = nullptr;
  TCanvas*      c  = nullptr;

  auto ddisplay = new hallc::MonitoringDisplay(setting_pair);

  ddisplay->CreateDisplayPlot("z_plots", "z_yield_positive",
                              [&](hallc::DisplayPlot& plt) {
                                auto c = plt.SetCanvas(
                                    new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
                                plt.SetPersist();
                                c->cd();
                                // hTotal->SetTitle(fmt::format("with {}
                                // C",total_charge/1000).c_str());
                                // h_z_Total_positive->SetTitle(
                                //    fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 *
                                //    (2.60943) / 66.0).c_str());
                                h_z_Total_positive->Draw("e1");
                                h_z_Total_positive->Draw("hist same");
                                c->SaveAs("results/csv_ratios/c1-2.png");
                                c->SaveAs("results/csv_ratios/c1-2.pdf");
                                return 0;
                              },
                              [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot("z_plots", "z_yield_negative",
                              [&](hallc::DisplayPlot& plt) {
                                auto c = plt.SetCanvas(
                                    new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
                                plt.SetPersist();
                                c->cd();
                                // hTotal->SetTitle(fmt::format("with {}
                                // C",total_charge/1000).c_str());
                                // h_z_Total_negative->SetTitle(
                                //    fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 *
                                //    (2.60943) / 66.0).c_str());
                                h_z_Total_negative->Draw("e1");
                                h_z_Total_negative->Draw("hist same");
                                c->SaveAs("results/csv_ratios/c2-2.png");
                                c->SaveAs("results/csv_ratios/c2-2.pdf");
                                return 0;
                              },
                              [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot("z_plots", "z_ratio",
                              [&](hallc::DisplayPlot& plt) {
                                auto c = plt.SetCanvas(
                                    new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
                                plt.SetPersist();
                                c->cd();
                                // hTotal->SetTitle(fmt::format("with {}
                                // C",total_charge/1000).c_str());
                                // h_z_Total_negative->SetTitle(
                                //    fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 *
                                //    (2.60943) / 66.0).c_str());
                                h_z_ratio->Draw("e1");
                                h_z_ratio->Draw("hist same");
                                c->SaveAs("results/csv_ratios/ratio-1.png");
                                c->SaveAs("results/csv_ratios/ratio-1.pdf");
                                return 0;
                              },
                              [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/ratios/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
