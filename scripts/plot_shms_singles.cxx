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

void plot_shms_singles(int start_run = 0) {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/jpsi_run_count_list.json");
    try {
      json_input_file >> j;
    } catch (json::parse_error) {
      std::cerr
          << "error: json file, db2/run_count_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }
  json j2;
  {
    std::ifstream json_input_file("db2/run_list_shms.json");
    try {
      json_input_file >> j2;
    } catch (json::parse_error) {
      std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }

  std::vector<double> runs;
  std::vector<double> hms_yields;
  std::vector<double> hms_yield_uncs;
  std::vector<double> run_current;

  std::vector<double> runs2;
  std::vector<double> hms_yields2;
  std::vector<double> hms_yield_uncs2;
  std::vector<double> run_current2;

  std::vector<double> runs3;
  std::vector<double> yield_ratios;
  std::vector<double> yield_ratio_uncs;
  std::vector<double> run_current3;

  std::vector<double> radin_runs;
  std::vector<double> radin_hms_yields;
  std::vector<double> radin_hms_yield_uncs;
  std::vector<double> radin_run_current;

  std::vector<double> radin_runs2;
  std::vector<double> radin_hms_yields2;
  std::vector<double> radin_hms_yield_uncs2;
  std::vector<double> radin_run_current2;

  std::vector<double> radin_runs3;
  std::vector<double> radin_yield_ratios;
  std::vector<double> radin_yield_ratio_uncs;
  std::vector<double> radin_run_current3;

  std::vector<double> target_ids;

  double scaler_reduce_factor = 0.4;

  for (json::iterator it = j.begin(); it != j.end(); ++it) {

    if (j2.find(it.key()) == j2.end()) {
      std::cout << it.key() << " not found \n";
      continue;
    }
    if (j2[it.key()]["target"]["target_id"].get<int>() != 2) {
      continue;
    }

    // Select shms polarity
    if (j2[it.key()]["spectrometers"]["shms_momentum"].get<double>() > 0.0) {
      continue;
    }
    bool radiator_in = (j2[it.key()]["radiator"]["radiator_status"].get<std::string>() == "IN");

    int    run_number  = std::stoi(it.key());
    int    run_number2 = 0;
    double yield       = 0.0;
    double yield2      = 0.0;

    if (run_number < start_run) {
      continue;
    }

    // charge bcm4b 2u cut
    // time 1MHz 2u cut
    double beam_current = 0.0;
    if (it.value().find("time 1MHz 2u cut") != it.value().end()) {
      if (it.value().find("charge bcm4b 2u cut") != it.value().end()) {
        beam_current = it.value()["charge bcm4b 2u cut"].get<double>() /
                       it.value()["time 1MHz 2u cut"].get<double>();
      }
    }

    if (it.value().find("shms e yield") != it.value().end()) {
      try {
        yield = it.value()["shms e yield"].get<double>();
        if (radiator_in) {
          radin_runs.push_back(double(run_number));
          radin_hms_yields.push_back(yield);
          radin_run_current.push_back(beam_current);
        } else {
          runs.push_back(double(run_number));
          run_current.push_back(beam_current);
          hms_yields.push_back(yield);
          if (it.value().find("shms e yield unc.") != it.value().end()) {
            hms_yield_uncs.push_back(it.value()["shms e yield unc."].get<double>());
          } else {
            hms_yield_uncs.push_back(0.0);
          }
        }
      } catch (std::domain_error) {
        continue;
        // you suck
      }
    }
    // Scaler yields
    if (it.value().find("shms scaler yield") != it.value().end()) {
      run_number2 = std::stoi(it.key());
      yield2      = it.value()["shms scaler yield"].get<double>();
      if (radiator_in) {
        radin_runs2.push_back(double(run_number2));
        radin_hms_yields2.push_back(scaler_reduce_factor * yield2);
        radin_run_current2.push_back(beam_current);
      } else {
        run_current2.push_back(beam_current);
        runs2.push_back(double(run_number2));
        hms_yields2.push_back(scaler_reduce_factor * yield2);
        if (it.value().find("shms scaler yield unc.") != it.value().end()) {
          hms_yield_uncs2.push_back(scaler_reduce_factor *
                                    it.value()["shms scaler yield unc."].get<double>());
        } else {
          hms_yield_uncs2.push_back(0.0);
        }
      }
    }
    if (it.value().find("shms scaler yield") != it.value().end()) {
      if (it.value().find("shms e yield") != it.value().end()) {
        double yield_ratio = yield / yield2;
        if (radiator_in) {
          radin_runs3.push_back(double(run_number2));
          radin_yield_ratios.push_back(yield_ratio);
          radin_run_current3.push_back(beam_current);
        } else {
          runs3.push_back(double(run_number2));
          yield_ratios.push_back(yield_ratio);
          run_current3.push_back(beam_current);
          if ((it.value().find("shms e yield unc.") != it.value().end()) &&
              (it.value().find("shms scaler yield unc.") != it.value().end())) {
            double y1 = it.value()["shms e yield unc."].get<double>();
            double y2 = it.value()["shms scaler yield unc."].get<double>();
            yield_ratio_uncs.push_back(
                yield_ratio * std::sqrt(y1 * y1 / (yield * yield) + y2 * y2 / (yield2 * yield2)));
          } else {
            yield_ratio_uncs.push_back(0.0);
          }
        }
      }
    }
  }

  // ------------------------------
  //
  TGraphErrors* gr = nullptr;
  TMultiGraph*  mg = nullptr;
  TCanvas*      c  = nullptr;

  auto ddisplay      = new hallc::MonitoringDisplay(-1);
  auto waveform_plot = ddisplay->CreateDisplayPlot(
      "test", "A",
      [&](hallc::DisplayPlot& plt) {
        plt._plot_data._canvas = new TCanvas(plt.GetName().c_str(), plt.GetName().c_str());
        // = {new TGraph(), new TGraph(), new TGraph(), new TGraph()};
        // ------------------------------
        //
        c = plt._plot_data._canvas; // new TCanvas("c1", "c1", 1600, 1200);
        c->Divide(1, 3);

        c->cd(1);
        mg = new TMultiGraph();
        mg->SetTitle("; run number ; counts/mC");

        std::vector<double> zeros(runs.size());
        gr = new TGraphErrors(runs.size(), runs.data(), hms_yields.data(), zeros.data(),
                              hms_yield_uncs.data());
        gr->SetMarkerStyle(20);
        mg->Add(gr, "p");

        TGraph* gr1 = new TGraph(radin_runs.size(), radin_runs.data(), radin_hms_yields.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(4);
        gr1->SetMarkerSize(1.8);
        mg->Add(gr1, "p");

        std::vector<double> zeros2(runs2.size());
        gr = new TGraphErrors(runs2.size(), runs2.data(), hms_yields2.data(), zeros2.data(),
                              hms_yield_uncs2.data());
        gr->SetMarkerStyle(23);
        gr->SetMarkerColor(2);
        mg->Add(gr, "p");

        gr1 = new TGraph(radin_runs2.size(), radin_runs2.data(), radin_hms_yields2.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(2);
        gr1->SetMarkerSize(1.8);
        mg->Add(gr1, "p");

        mg->Draw("a");

        c->cd(2);
        mg = new TMultiGraph();
        mg->SetTitle("; run number ; ratio counter/scaler");
        std::vector<double> zeros3(runs3.size());
        gr = new TGraphErrors(runs3.size(), runs3.data(), yield_ratios.data(), zeros3.data(),
                              yield_ratio_uncs.data());
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.6);
        mg->Add(gr, "p");

        gr1 = new TGraph(radin_runs3.size(), radin_runs3.data(), radin_yield_ratios.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(2);
        gr1->SetMarkerSize(1.6);
        mg->Add(gr1, "p");
        mg->Draw("a");


        c->cd(3);
        mg = new TMultiGraph();
        mg->SetTitle("; run number ; beam current");
        gr1 = new TGraph(radin_runs3.size(), radin_runs3.data(), radin_run_current3.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(8);
        gr1->SetMarkerSize(1.6);
        mg->Add(gr1, "p");

        gr1 = new TGraph(runs3.size(), runs3.data(), run_current3.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(4);
        gr1->SetMarkerSize(1.6);
        mg->Add(gr1, "p");


        // gr->Draw("alp");
        // gr->GetYaxis()->SetRangeUser(0.3, 0.45);
        mg->Draw("a");

        c->SaveAs("results/hms_singles_vs_run_number.png");
        c->SaveAs("results/hms_singles_vs_run_number.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  auto plot2 = ddisplay->CreateDisplayPlot(
      "test", "vs_beam_current",
      [&](hallc::DisplayPlot& plt) {
        plt._plot_data._canvas = new TCanvas(plt.GetName().c_str(), plt.GetName().c_str());
        // = {new TGraph(), new TGraph(), new TGraph(), new TGraph()};
        // ------------------------------
        //
        c = plt._plot_data._canvas; // new TCanvas("c1", "c1", 1600, 1200);
        c->Divide(1, 2);

        c->cd(1);
        mg = new TMultiGraph();
        mg->SetTitle("; run number ; counts/mC");

        std::vector<double> zeros(runs.size());
        gr = new TGraphErrors(runs.size(), run_current.data(), hms_yields.data(), zeros.data(),
                              hms_yield_uncs.data());
        gr->SetMarkerStyle(20);
        mg->Add(gr, "p");

        TGraph* gr1 = new TGraph(radin_runs.size(), radin_run_current.data(), radin_hms_yields.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(4);
        gr1->SetMarkerSize(1.8);
        mg->Add(gr1, "p");

        std::vector<double> zeros2(runs2.size());
        gr = new TGraphErrors(runs2.size(), run_current2.data(), hms_yields2.data(), zeros2.data(),
                              hms_yield_uncs2.data());
        gr->SetMarkerStyle(23);
        gr->SetMarkerColor(2);
        mg->Add(gr, "p");

        gr1 = new TGraph(radin_runs2.size(), radin_run_current2.data(), radin_hms_yields2.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(2);
        gr1->SetMarkerSize(1.8);
        mg->Add(gr1, "p");

        mg->Draw("a");
        mg->GetXaxis()->SetLimits(0,0.05);

        c->cd(2);
        mg = new TMultiGraph();
        mg->SetTitle("; run number ; ratio counter/scaler");
        std::vector<double> zeros3(runs3.size());
        gr = new TGraphErrors(runs3.size(), run_current3.data(), yield_ratios.data(), zeros3.data(),
                              yield_ratio_uncs.data());
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.6);
        mg->Add(gr, "p");

        gr1 = new TGraph(radin_runs3.size(), radin_run_current3.data(), radin_yield_ratios.data());
        gr1->SetMarkerStyle(20);
        gr1->SetMarkerColor(2);
        gr1->SetMarkerSize(1.6);
        mg->Add(gr1, "p");

        // gr->Draw("alp");
        // gr->GetYaxis()->SetRangeUser(0.3, 0.45);
        //mg->Draw("a");

        //c->cd(3);
        //mg = new TMultiGraph();
        //mg->SetTitle("; beam current ; ");

        //gr1 = new TGraph(runs3.size(), run_current3.data(), run_current3.data());
        //gr1->SetMarkerStyle(20);
        //gr1->SetMarkerColor(4);
        //gr1->SetMarkerSize(1.6);
        //mg->Add(gr1, "p");

        //gr1 = new TGraph(radin_runs3.size(), radin_run_current3.data(), radin_run_current3.data());
        //gr1->SetMarkerStyle(20);
        //gr1->SetMarkerColor(8);
        //gr1->SetMarkerSize(1.6);
        //mg->Add(gr1, "p");


        // gr->Draw("alp");
        // gr->GetYaxis()->SetRangeUser(0.3, 0.45);
        mg->Draw("a");
        mg->GetXaxis()->SetLimits(0,0.05);

        c->SaveAs("results/hms_singles_vs_run_number.png");
        c->SaveAs("results/hms_singles_vs_run_number.pdf");

        return 0;
      },
      [&](hallc::DisplayPlot& plt) { return 0; });
  ddisplay->_data._root_folder = "/derp/";
  ddisplay->InitAll();
  ddisplay->UpdateAll();

  std::this_thread::sleep_for(std::chrono::seconds(200));
}
