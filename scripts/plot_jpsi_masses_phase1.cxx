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

void plot_jpsi_masses_phase1(int start_run = 0) {

  auto mcFile          = new TFile("monitoring/simulation/mc-histos-phase1.root");
  auto mc_yield_Egamma = (TH1D*)mcFile->Get("hSumEgamma_00");
  auto mc_yield_Abst   = (TH1D*)mcFile->Get("hSumAbst_00");

  using nlohmann::json;
  json j2;
  {
    std::ifstream json_input_file("db2/jpsi_run_count_list.json");
    json_input_file >> j2;
  }
  double total_charge = 0.000;

  std::vector<int> kine_0_runs_old = {7225, 7226, 7228, 7229, 7230, 7233, 7234, 7235, 7236, 7237,
                                      7238, 7239, 7240, 7241, 7245, 7246, 7247, 7248, 7249, 7264,
                                      7265, 7266, 7267, 7268, 7269, 7270, 7271, 7272, 7273};
  std::vector<int> kine_0_runs     = {7225, 7226, 7228, 7229, 7230, 7231, 7232, 7233, 7234, 7235,
                                  7236, 7237, 7238, 7239, 7240, 7241, 7245, 7246, 7247, 7248,
                                  7249, 7264, 7265, 7266, 7267, 7268, 7269, 7270, 7271, 7272,
                                  7273, 7411, 7412, 7413, 7414, 7415, 7417, 7418, 7419, 7420,
                                  7421, 7422, 7423, 7424, 7425, 7426};

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
      "total", "inv_mass_phase1",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        // hTotal->SetTitle(fmt::format("with {} C",total_charge/1000).c_str());
        //hTotal->SetTitle(
        //    fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str());
        hTotal->Draw("e1");
        hTotal->Draw("hist same");
        c->SaveAs("results/plot_jpsi_masses/c1-1.png");
        c->SaveAs("results/plot_jpsi_masses/c1-1.pdf");
        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot(
      "total", "E_gamma_phase1",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        //hEgammaTotal->SetTitle(
        //    fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str());
        hEgammaTotal->SetMarkerStyle(20);
        hEgammaTotal->SetMarkerSize(1.7);
        hEgammaTotal->Draw("e1");

        mc_yield_Egamma->SetLineColor(kRed + 1);
        mc_yield_Egamma->Draw("hist same");

        //mc_yield_Egamma_30->SetLineColor(kBlue + 1);
        //mc_yield_Egamma_30->Draw("hist same");

        //mc_50_yield_Egamma->SetLineColor(kBlue + 1);
        //mc_50_yield_Egamma->Draw("hist same");
        
        TLegend* leg = new TLegend(0.15, 0.6, 0.52, 0.8);

        //leg->AddEntry(
        //    hEgammaTotal,
        //    fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str(),
        //    "lep");
        leg->AddEntry(mc_yield_Egamma, "t-channel production only");
        //        leg->AddEntry(mc_yield_Egamma_30, "+ 3\% pentaquark coupling");
        leg->Draw();
        c->SaveAs("results/plot_jpsi_masses/c2-1.png");
        c->SaveAs("results/plot_jpsi_masses/c2-1.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot(
      "total", "Abs_t_phase1",
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

        //leg->AddEntry(
        //    hEgammaTotal,
        //    fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (1.67) / 66.0).c_str(), "lep");
        leg->AddEntry(mc_yield_Abst, "t-channel production only");
        leg->Draw();

        c->SaveAs("results/plot_jpsi_masses/c3-1.png");
        c->SaveAs("results/plot_jpsi_masses/c3-1.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/progress/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
