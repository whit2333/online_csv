#ifdef __cpp_lib_filesystem
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

#include "nlohmann/json.hpp"

#include "TObject.h"
#include "TGraph.h"


void plot_hms_singles(int start_run = 0){

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_count_list.json");
    try {
      json_input_file >> j;
    } catch(json::parse_error) {
      std::cerr << "error: json file, db2/run_count_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }
  json j2;
  {
    std::ifstream json_input_file("db2/run_list.json");
    try {
      json_input_file >> j2;
    } catch(json::parse_error)  {
      std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }


  std::vector<double> runs;
  std::vector<double> hms_yields;
  std::vector<double> hms_yield_uncs;

  std::vector<double> runs2;
  std::vector<double> hms_yields2;
  std::vector<double> hms_yield_uncs2;

  std::vector<double> runs3;
  std::vector<double> yield_ratios;
  std::vector<double> yield_ratio_uncs;

  std::vector<double> target_ids;

  TGraphErrors* gr_hms_count  = new TGraphErrors();
  TGraphErrors* gr_hms_scaler = new TGraphErrors();
  TGraphErrors* gr_yield_ratio = new TGraphErrors();


  double scaler_reduce_factor = 0.4;

  for (json::iterator it = j.begin(); it != j.end(); ++it) {

    if (j2.find(it.key()) == j2.end()) {
      continue;
    }
    if (j2[it.key()]["target"]["target_id"].get<int>() != 3) {
      continue;
    }

    int run_number = std::stoi(it.key());
    int run_number2 = 0;
    double yield      = 0.0;
    double yield2      =0.0;

    if(run_number < start_run ) {
      continue;
    }

    if (it.value().find("hms e yield") != it.value().end()) {
      try {
      yield      = it.value()["hms e yield"].get<double>();
      runs.push_back(double(run_number));
      hms_yields.push_back(yield);
      if (it.value().find("hms e yield unc.") != it.value().end()) {
        hms_yield_uncs.push_back(it.value()["hms e yield unc."].get<double>());
      } else {
        hms_yield_uncs.push_back(0.0);
      }
      } catch(std::domain_error ) {
        continue;
        //you suck
      }
    }
    if (it.value().find("hms scaler yield") != it.value().end()) {
      run_number2 = std::stoi(it.key());
      yield2      = it.value()["hms scaler yield"].get<double>();
      runs2.push_back(double(run_number2));
      hms_yields2.push_back(scaler_reduce_factor*yield2);
      if (it.value().find("hms scaler yield unc.") != it.value().end()) {
        hms_yield_uncs2.push_back(scaler_reduce_factor*it.value()["hms scaler yield unc."].get<double>());
      } else {
        hms_yield_uncs2.push_back(0.0);
      }
    }
    if (it.value().find("hms scaler yield") != it.value().end()) {
      if (it.value().find("hms e yield") != it.value().end()) {
        runs3.push_back(double(run_number2));
        double yield_ratio = yield / yield2;
        yield_ratios.push_back(yield_ratio);
        if ((it.value().find("hms e yield unc.") != it.value().end()) &&
            (it.value().find("hms scaler yield unc.") != it.value().end())) {
          double y1 = it.value()["hms e yield unc."].get<double>();
          double y2 = it.value()["hms scaler yield unc."].get<double>();
          yield_ratio_uncs.push_back(yield_ratio*std::sqrt(y1 * y1/(yield*yield) + y2 * y2/(yield2*yield2)));
        } else {
          yield_ratio_uncs.push_back(0.0);
        }
      }
    }
  }

  // ------------------------------
  //
  TGraphErrors* gr = nullptr;
  TMultiGraph*  mg = nullptr;
  TCanvas*      c  = nullptr;

  // ------------------------------
  //
  c  = new TCanvas("c1","c1",1600,1200);
  c->Divide(1, 2);

  c->cd(1);
  mg = new TMultiGraph();
  mg->SetTitle("; run number ; counts/mC");

  std::vector<double> zeros(runs.size());
  gr = new TGraphErrors(runs.size(), runs.data(), hms_yields.data(),zeros.data(),hms_yield_uncs.data());
  gr->SetMarkerStyle(20);
  mg->Add(gr,"p");

  std::vector<double> zeros2(runs2.size());
  gr = new TGraphErrors(runs2.size(), runs2.data(), hms_yields2.data(),zeros2.data(),hms_yield_uncs2.data());
  gr->SetMarkerStyle(23);
  gr->SetMarkerColor(2);
  mg->Add(gr,"p");

  mg->Draw("a");

  c->cd(2);
  std::vector<double> zeros3(runs3.size());
  gr = new TGraphErrors(runs3.size(), runs3.data(), yield_ratios.data(),zeros3.data(),yield_ratio_uncs.data());
  gr->SetTitle("; run number ; ps4*T2 + T6 counts / scaler");
  gr->SetMarkerStyle(20);
  //gr->Draw("alp");
  //gr->GetYaxis()->SetRangeUser(0.3, 0.45);
  gr->Draw("alp");

  c->SaveAs("results/hms_singles_vs_run_number.png");
  c->SaveAs("results/hms_singles_vs_run_number.pdf");
}
