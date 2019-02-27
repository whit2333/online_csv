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

void plot_jpsi_masses_phase3(int start_run = 0) {

  auto mcFile          = new TFile("monitoring/simulation/mc-histos-phase3.root");
  auto mc_yield_Egamma = (TH1D*)mcFile->Get("hSumEgamma_00");
  auto mc_yield_Abst   = (TH1D*)mcFile->Get("hSumAbst_00");
  //  auto mc_yield_Egamma_30 = (TH1D*)mcFile->Get("hSumEgamma_30");

  using nlohmann::json;
  json j2;
  {
    std::ifstream json_input_file("db2/jpsi_run_count_list.json");
    json_input_file >> j2;
  }
  double total_charge = 0.000;

  std::vector<int> kine_0_runs = {
      7308, 7312, 7313, 7314, 7315, 7316, 7317, 7318, 7319, 7320, 7322, 7324, 7325, 7326, 7327,
      7328, 7329, 7330, 7331, 7332, 7333, 7334, 7335, 7336, 7337, 7338, 7339, 7340, 7341, 7342,
      7343, 7344, 7345, 7346, 7347, 7348, 7349, 7350, 7351, 7352, 7353, 7354, 7355, 7356, 7357,
      7358, 7359, 7360, 7361, 7362, 7363, 7364, 7365, 7366, 7367, 7368, 7369, 7370, 7371, 7372,
      7373, 7374, 7375, 7376, 7377, 7378, 7379, 7380, 7381, 7384, 7385, 7386, 7387, 7388, 7389,
      7390, 7391, 7392, 7393, 7394, 7395, 7396, 7397, 7398, 7399, 7400, 7401, 7402, 7403, 7404};

  // total histograms
  TH1D* hTotal       = nullptr;
  TH1D* hEgammaTotal = nullptr;
  TH1D* hAbstTotal   = nullptr;

  for (auto run_number : kine_0_runs) {

    std::string file_name = fmt::format("monitoring/{}/good_jpsi_counter.root", run_number);

    if (!root_file_exists(file_name)) {
      continue;
    }
    auto hist_file = new TFile(file_name.c_str(), "READ");

    TH1D* hnew  = (TH1D*)hist_file->Get("JpsiMassAfterTiming");
    TH1D* hnew2 = (TH1D*)hist_file->Get("JpsiEgammaFreeAfterCuts");
    TH1D* hnew3 = (TH1D*)hist_file->Get("JpsiAbstAfterCuts");

    if (!hnew) {
      std::cout << "JpsiMassAfterTiming not found \n";
      continue;
    }
    if (!hnew2) {
      std::cout << "JpsiEgammaFreeAfterCuts not found \n";
      continue;
    }
    if (!hnew3) {
      std::cout << "JpsiAbstAfterCuts not found \n";
      continue;
    }
    if (!hTotal) {
      hTotal = (TH1D*)hnew->Clone("total");
      std::cout << " asdfaksdfkj\n";
    } else {
      hTotal->Add((TH1D*)hnew->Clone());
    }
    if (!hEgammaTotal) {
      hEgammaTotal = (TH1D*)hnew2->Clone("Egamma_total");
    } else {
      hEgammaTotal->Add(hnew2);
    }
    if (!hAbstTotal) {
      hAbstTotal = (TH1D*)hnew3->Clone("Abst_total");
    } else {
      hAbstTotal->Add(hnew3);
    }

    if (j2[std::to_string(run_number)].find("good_total_charge") !=
        j2[std::to_string(run_number)].end()) {
      total_charge += j2[std::to_string(run_number)]["good_total_charge"].get<double>();
    }
  }

  std::cout << "Total charge: " << total_charge / 1000.0 << " C\n";
  bool make_yields = true;
  if (make_yields) {
    hTotal->Scale(1000.0 / total_charge);
    hEgammaTotal->Scale(1000.0 / total_charge);
    hAbstTotal->Scale(1000.0 / total_charge);

    hTotal->GetYaxis()->SetTitle("Yield [counts/C]");
    hEgammaTotal->GetYaxis()->SetTitle("Yield [counts/C]");
    hAbstTotal->GetYaxis()->SetTitle("Yield [counts/C]");
  }

  // ------------------------------
  //
  TGraphErrors* gr = nullptr;
  TMultiGraph*  mg = nullptr;
  TCanvas*      c  = nullptr;

  auto ddisplay = new hallc::MonitoringDisplay(-1);
  ddisplay->CreateDisplayPlot(
      "total", "inv_mass_phase3",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        // hTotal->SetTitle(fmt::format("with {} C",total_charge/1000).c_str());
        hTotal->SetTitle(
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str());
        hTotal->Draw("e1");
        hTotal->Draw("hist same");
        c->SaveAs("results/plot_jpsi_masses/c1-3.png");
        c->SaveAs("results/plot_jpsi_masses/c1-3.pdf");
        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot(
      "total", "E_gamma_phase3",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        hEgammaTotal->SetTitle(
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str());
        hEgammaTotal->SetMarkerStyle(20);
        hEgammaTotal->SetMarkerSize(1.7);
        hEgammaTotal->Draw("e1");

        mc_yield_Egamma->SetLineColor(kRed + 1);
        mc_yield_Egamma->Draw("hist same");
#if 0
        mc_yield_Egamma_30->SetLineColor(kBlue + 1);
        mc_yield_Egamma_30->Draw("hist same");
#endif

#if 0
        mc_50_yield_Egamma->SetLineColor(kBlue + 1);
        mc_50_yield_Egamma->Draw("hist same");
#endif
        TLegend* leg = new TLegend(0.15, 0.6, 0.52, 0.8);

        leg->AddEntry(
            hEgammaTotal,
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str(),
            "lep");
        leg->AddEntry(mc_yield_Egamma, "t-channel production only");
#if 0
        leg->AddEntry(mc_yield_Egamma_30, "+ 3\% pentaquark coupling");
#endif
        leg->Draw();
        c->SaveAs("results/plot_jpsi_masses/c2-3.png");
        c->SaveAs("results/plot_jpsi_masses/c2-3.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot(
      "total", "Abs_t_phase3",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        hAbstTotal->SetMarkerStyle(20);
        hAbstTotal->SetMarkerSize(1.7);
        hAbstTotal->Draw("e1");

        mc_yield_Abst->SetLineColor(kRed + 1);
        mc_yield_Abst->Draw("hist same");
        TLegend* leg = new TLegend(0.55, 0.6, 0.85, 0.8);

        leg->AddEntry(
            hEgammaTotal,
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (1.67) / 66.0).c_str(), "lep");
        leg->AddEntry(mc_yield_Abst, "t-channel production only");
        leg->Draw();

        c->SaveAs("results/plot_jpsi_masses/c3-3.png");
        c->SaveAs("results/plot_jpsi_masses/c3-3.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/progress/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
