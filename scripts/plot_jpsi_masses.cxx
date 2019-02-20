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

void plot_jpsi_masses(int start_run = 0) {

  auto mcFile          = new TFile("monitoring/simulation/pcsim.106-kin-2-0.52p32m.c00.root");
  auto mc_yield_Egamma = (TH1D*)mcFile->Get("hsum_Egamma_rc");
  auto mc_yield_t      = (TH1D*)mcFile->Get("hsum_fabs(t_rc)");

  auto mcFile_40          = new TFile("monitoring/simulation/pcsim.106-kin-2-0.52p32m.c40.root");
  auto mc_40_yield_Egamma = (TH1D*)mcFile_40->Get("hsum_Egamma_rc");
  auto mc_40_yield_t      = (TH1D*)mcFile_40->Get("hsum_fabs(t_rc)");

  auto mcFile_50          = new TFile("monitoring/simulation/pcsim.106-kin-2-0.52p32m.c50.root");
  auto mc_50_yield_Egamma = (TH1D*)mcFile_50->Get("hsum_Egamma_rc");
  auto mc_50_yield_t      = (TH1D*)mcFile_50->Get("hsum_fabs(t_rc)");

  using nlohmann::json;
  json j2;
  {
    std::ifstream json_input_file("db2/jpsi_run_count_list.json");
    json_input_file >> j2;
  }
  double total_charge = 0.000;

  //  std::vector<int> kine_0_runs = {7226, 7229, 7230, 7233, 7234, 7235, 7236, 7237, 7238,
  //  7239, 7240, 7241, 7245, 7246, 7247, 7248, 7249, 7264, 7265, 7266, 7267, 7268, 7269 };
  // std::vector<int> kine_0_runs = {7264, 7265, 7266, 7267, 7268, 7269, 7270, 7271, 7272, 7273};
  std::vector<int> kine_0_runs = {7225, 7226, 7228, 7229, 7230, 7233, 7234, 7235, 7236, 7237,
                                  7238, 7239, 7240, 7241, 7245, 7246, 7247, 7248, 7249, 7264,
                                  7265, 7266, 7267, 7268, 7269, 7270, 7271, 7272, 7273};

  // Dummy
  // std::vector<int> kine_0_runs       = { 7242, 7243, 7244};

  // total histograms
  TH1D* hTotal           = nullptr;
  TH1D* hEgammaTotal     = nullptr;
  TH1D* hAbstTotal       = nullptr;
  TH1D* hAbstTotal_lowE  = nullptr;
  TH1D* hAbstTotal_highE = nullptr;

  json j;
  json j_pe;
  json j_t;
  json j_t_lowE;
  json j_t_highE;
  for (auto run_number : kine_0_runs) {

    {
      std::ifstream json_input_file(
          fmt::format("monitoring/{}/Invariant_Mass_after_timing.json", run_number));
      try {
        json_input_file >> j;
      } catch (json::parse_error) {
        std::cerr << "error: json file, "
                  << fmt::format("monitoring/{}/Invariant_Mass_after_timing.json", run_number)
                  << " is incomplete or has broken syntax.\n";
        // std::quick_exit(-127);
        continue;
      }
    }
    {
      std::ifstream json_input_file(
          fmt::format("monitoring/{}/Photon_Energy_after_cuts.json", run_number));
      try {
        json_input_file >> j_pe;
      } catch (json::parse_error) {
        std::cerr << "error: json file, "
                  << fmt::format("monitoring/{}/Photon_Energy2_after_cuts.json", run_number)
                  << " is incomplete or has broken syntax.\n";
        // std::quick_exit(-127);
        continue;
      }
    }
    {
      std::ifstream json_input_file(fmt::format("monitoring/{}/Abst_after_cuts.json", run_number));
      try {
        json_input_file >> j_t;
      } catch (json::parse_error) {
        std::cerr << "error: json file, "
                  << fmt::format("monitoring/{}/Abst_after_cuts.json", run_number)
                  << " is incomplete or has broken syntax.\n";
        // std::quick_exit(-127);
        continue;
      }
    }
    {
      std::ifstream json_input_file(
          fmt::format("monitoring/{}/Abst_after_cuts_lowE.json", run_number));
      try {
        json_input_file >> j_t_lowE;
      } catch (json::parse_error) {
        std::cerr << "error: json file, "
                  << fmt::format("monitoring/{}/Abst_after_cuts_lowE.json", run_number)
                  << " is incomplete or has broken syntax.\n";
        // std::quick_exit(-127);
        continue;
      }
    }
    {
      std::ifstream json_input_file(
          fmt::format("monitoring/{}/Abst_after_cuts_highE.json", run_number));
      try {
        json_input_file >> j_t_highE;
      } catch (json::parse_error) {
        std::cerr << "error: json file, "
                  << fmt::format("monitoring/{}/Abst_after_cuts_highE.json", run_number)
                  << " is incomplete or has broken syntax.\n";
        // std::quick_exit(-127);
        continue;
      }
    }

    TH1D* hnew  = nullptr;
    TH1D* hnew2 = nullptr;
    TH1D* hnew3 = nullptr;

    TBufferJSON::FromJSON(hnew, j.dump().c_str());
    TBufferJSON::FromJSON(hnew2, j_pe.dump().c_str());
    TBufferJSON::FromJSON(hnew3, j_t.dump().c_str());

    TH1D* hnew3_lowE  = nullptr;
    TH1D* hnew3_highE = nullptr;
    TBufferJSON::FromJSON(hnew3_lowE, j_t_lowE.dump().c_str());
    TBufferJSON::FromJSON(hnew3_highE, j_t_highE.dump().c_str());

    if (hnew && hnew2 && hnew3 && hnew3_lowE && hnew3_highE) {

      if (!hTotal) {
        hTotal = (TH1D*)hnew->Clone("total");
      } else {
        hTotal->Add(hnew);
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
      if (!hAbstTotal_lowE) {
        hAbstTotal_lowE = (TH1D*)hnew3_lowE->Clone("Abst_total_lowE");
      } else {
        hAbstTotal_lowE->Add(hnew3_lowE);
      }
      if (!hAbstTotal_highE) {
        hAbstTotal_highE = (TH1D*)hnew3_highE->Clone("Abst_total_highE");
      } else {
        hAbstTotal_highE->Add(hnew3_highE);
      }

      if (j2[std::to_string(run_number)].find("good_total_charge") !=
          j2[std::to_string(run_number)].end()) {
        total_charge += j2[std::to_string(run_number)]["good_total_charge"].get<double>();
        // fmt::print("/ {:>10.1f} ", total_charge);
      }
    }
  }

  std::cout << "Total charge: " << total_charge / 1000.0 << " C\n";
  bool make_yields = true;
  if (make_yields) {
    hTotal->Scale(1000.0 / total_charge);
    hEgammaTotal->Scale(1000.0 / total_charge);
    hAbstTotal->Scale(1000.0 / total_charge);
    hAbstTotal_lowE->Scale(1000.0 / total_charge);
    hAbstTotal_highE->Scale(1000.0 / total_charge);

    hTotal->GetYaxis()->SetTitle("Yield [counts/C]");
    hEgammaTotal->GetYaxis()->SetTitle("Yield [counts/C]");
    hAbstTotal->GetYaxis()->SetTitle("Yield [counts/C]");
    hAbstTotal_lowE->GetYaxis()->SetTitle("Yield [counts/C]");
    hAbstTotal_highE->GetYaxis()->SetTitle("Yield [counts/C]");
  }

  // ------------------------------
  //
  TGraphErrors* gr = nullptr;
  TMultiGraph*  mg = nullptr;
  TCanvas*      c  = nullptr;

  auto ddisplay = new hallc::MonitoringDisplay(-1);
  ddisplay->CreateDisplayPlot(
      "total", "inv_mass",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        // hTotal->SetTitle(fmt::format("with {} C",total_charge/1000).c_str());
        hTotal->SetTitle(
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str());
        hTotal->Draw("e1");
        hTotal->Draw("hist same");
        c->SaveAs("results/plot_jpsi_masses/c1.png");
        c->SaveAs("results/plot_jpsi_masses/c1.pdf");
        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot(
      "total", "E_gamma",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        hEgammaTotal->SetTitle(
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str());
        hEgammaTotal->SetMarkerStyle(20);
        hEgammaTotal->SetMarkerSize(1.7);
        hEgammaTotal->Draw("e1");

        mc_yield_Egamma->SetLineColor(2);
        mc_yield_Egamma->Draw("hist same");

        // mc_40_yield_Egamma->SetLineColor(kMagenta+2);
        // mc_40_yield_Egamma->Draw("hist same");

        mc_50_yield_Egamma->SetLineColor(kBlue + 1);
        mc_50_yield_Egamma->Draw("hist same");
        TLegend* leg = new TLegend(0.15, 0.6, 0.52, 0.8);

        leg->AddEntry(
            hEgammaTotal,
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (2.60943) / 66.0).c_str(),
            "lep");
        leg->AddEntry(mc_yield_Egamma, "t-channel production only");
        // leg->AddEntry(mc_40_yield_Egamma, "+ 4\% coupling");
        leg->AddEntry(mc_50_yield_Egamma, "+ 5\% pentaquark coupling");
        leg->Draw();
        c->SaveAs("results/plot_jpsi_masses/c2.png");
        c->SaveAs("results/plot_jpsi_masses/c2.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->CreateDisplayPlot(
      "total", "Abs_t",
      [&](hallc::DisplayPlot& plt) {
        auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
        // plt.SetPersist();
        c->cd();
        hAbstTotal->SetMarkerStyle(20);
        hAbstTotal->SetMarkerSize(1.7);
        hAbstTotal->Draw("e1");

        mc_yield_t->SetLineColor(2);
        mc_yield_t->Draw("hist same");

        // mc_40_yield_t->SetLineColor(kMagenta+2);
        // mc_40_yield_t->Draw("hist same");

        mc_50_yield_t->SetLineColor(kBlue + 1);
        mc_50_yield_t->Draw("hist same");
        TLegend* leg = new TLegend(0.55, 0.6, 0.85, 0.8);

        leg->AddEntry(
            hEgammaTotal,
            fmt::format("{:3.1f}\% of scheduled beam time ", 100.0 * (1.67) / 66.0).c_str(), "lep");
        leg->AddEntry(mc_yield_t, "t-channel production only");
        // leg->AddEntry(mc_40_yield_t, "+ 4\% coupling");
        leg->AddEntry(mc_50_yield_t, "+ 5\% pentaquark coupling");
        leg->Draw();

        c->SaveAs("results/plot_jpsi_masses/c3.png");
        c->SaveAs("results/plot_jpsi_masses/c3.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  // ddisplay->CreateDisplayPlot(
  //    "total", "t_fits",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas(new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      c->Divide(2, 2);
  //      plt.SetPersist();
  //      c->cd(1);
  //      hAbstTotal->SetTitle(
  //          fmt::format("{:3.1f}\% of scheduled beam time ", 100 * (total_charge / 1000) / 66.0)
  //              .c_str());
  //      hAbstTotal->Fit("expo", "", "", 0.55, 1.5);
  //      hAbstTotal->Draw("e1");
  //      mc_yield_t->SetLineColor(2);
  //      mc_yield_t->Draw("hist same");
  //      c->cd(3);
  //      hAbstTotal_lowE->Fit("expo", "", "", 0.55, 1.5);
  //      hAbstTotal_lowE->Draw("e1");
  //      c->cd(4);
  //      hAbstTotal_highE->SetLineColor(4);
  //      hAbstTotal_highE->Fit("expo", "", "", 0.55, 1.5);
  //      hAbstTotal_highE->Draw("e1");
  //      return 0;
  //    },
  //    [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/progress/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();
}
