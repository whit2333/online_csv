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


void plot_hms_singles(){

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_count_list.json");
    json_input_file >> j;
  }

  std::cout << j.dump(2);

  std::cout << " runs : ";
  std::vector<double> runs;
  std::vector<double> hms_yields;
  std::vector<double> runs2;
  std::vector<double> hms_yields2;

  std::vector<double> runs3;
  std::vector<double> yield_ratios;

  for (json::iterator it = j.begin(); it != j.end(); ++it) {

    int run_number = 0;
    int run_number2 = 0;
    double yield      = 0.0;
    double yield2      =0.0;
    if(it.value().find("hms e yield") != it.value().end()) {
      run_number = std::stoi(it.key());
      yield      = it.value()["hms e yield"].get<double>();
      runs.push_back(double(run_number));
      hms_yields.push_back(yield);
    }
    if(it.value().find("hms scaler yield") != it.value().end()) {
      run_number2 = std::stoi(it.key());
      yield2      = it.value()["hms scaler yield"].get<double>();
      runs2.push_back(double(run_number2));
      hms_yields2.push_back(yield2);
    }
    if(it.value().find("hms scaler yield") != it.value().end()) {
      if(it.value().find("hms e yield") != it.value().end()) {
        runs3.push_back(double(run_number2));
        yield_ratios.push_back(yield/yield2);
      }
    }
  }

  // ------------------------------
  //
  TGraph*      gr = nullptr;
  TMultiGraph* mg = nullptr;
  TCanvas*     c  = nullptr;

  // ------------------------------
  //
  c  = new TCanvas();
  c->Divide(1,2);

  c->cd(1);
  mg = new TMultiGraph();

  gr = new TGraph(runs.size(), runs.data(), hms_yields.data());
  gr->SetMarkerStyle(20);
  mg->Add(gr,"p");

  gr = new TGraph(runs2.size(), runs2.data(), hms_yields2.data());
  gr->SetMarkerStyle(23);
  gr->SetMarkerColor(2);
  mg->Add(gr,"p");

  mg->Draw("a");

  c->cd(2);
  gr = new TGraph(runs3.size(), runs3.data(), yield_ratios.data());
  gr->SetMarkerStyle(20);
  gr->Draw("alp");
}
