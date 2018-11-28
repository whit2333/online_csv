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

#include <sstream>
#include <iomanip>
#include <ctime>

void plot_charge_vs_time(int start_run = 0){

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
  std::vector<double> time1;
  std::vector<double> daq_triggers;
  std::vector<double> daq_triggers2;

  std::vector<double> time2;
  std::vector<double> runs2;
  std::vector<double> charges;
  std::vector<double> charges2;

  std::vector<double> hms_yields2;
  std::vector<double> hms_yield_uncs2;

  std::vector<double> runs3;
  std::vector<double> yield_ratios;
  std::vector<double> yield_ratio_uncs;

  std::vector<double> target_ids;

  TGraph* gr_triggers = new TGraph();
  TGraph* gr_charge   = new TGraph();

  TGraphErrors* gr_yield_ratio = new TGraphErrors();

  double total_accumulated_charge2 = 0.0;
  double total_accumulated_charge = 0.0;
  double total_accumulated_triggers = 0.0;
  double total_accumulated_triggers2 = 0.0;

  for (json::iterator it = j.begin(); it != j.end(); ++it) {

    int run_number = std::stoi(it.key());
    int run_number2 = 0;
    double yield    = 0.0;
    double yield2   = 0.0;


    if (j2.find(it.key()) == j2.end()) {
      continue;
    }
    //if (j2[it.key()]["target"]["target_id"].get<int>() != 3) {
    //  continue;
    //}

    if(run_number < start_run ) {
      continue;
    }

    double P0_shms_setting =
    j2[it.key()]["spectrometers"]["shms_momentum"].get<double>();

    double sign = 1.0;
    if( P0_shms_setting < 0 ) {
      sign = -1.0;
    }

    std::string start_time = j2[it.key()]["start_time"].get<std::string>();
    std::cout << start_time << "\n";

    std::tm a_tm = {};
    std::stringstream time_date_ss(start_time );
    time_date_ss >> std::get_time(&a_tm, "%H:%M:%S %m/%d/%y EST");


    if (it.value().find("total trigger events") != it.value().end()) {
      try {
        int run_triggers =  it.value()["pion bg sub. counts"].get<int>();
        total_accumulated_triggers += double(run_triggers);
        total_accumulated_triggers2 += sign*double(run_triggers);
        runs.push_back(double(run_number));
        time1.push_back(double(a_tm.tm_mday)+ double(a_tm.tm_hour)/24.0 + double(a_tm.tm_min)/60.0/24.0);
        daq_triggers.push_back(total_accumulated_triggers);
        daq_triggers2.push_back(total_accumulated_triggers2);
      } catch(std::domain_error ) {
        std::cout << " you suck\n";
        //you suck
        continue;
      }
    }
    if (it.value().find("charge bcm4b 2u cut") != it.value().end()) {
      try {
        int run_charge =  it.value()["charge bcm4b 2u cut"].get<double>();
        total_accumulated_charge += double(run_charge/1000.0);
        total_accumulated_charge2 += sign*double(run_charge/1000.0);
        runs2.push_back(double(run_number));
        double atime = a_tm.tm_mday+ double(a_tm.tm_hour)/24.0 + double(a_tm.tm_min)/60.0/24.0;
        //std::cout << atime << "\n";
        //std::cout << "d " <<  a_tm.tm_mday << "\n";
        //std::cout << "h " <<  a_tm.tm_hour << "\n";
        //std::cout << "m " <<  a_tm.tm_min << "\n";
        time2.push_back(atime);
        charges.push_back(total_accumulated_charge);
        charges2.push_back(total_accumulated_charge2);
      } catch(std::domain_error ) {
        std::cout << " you suck\n";
        continue;
      }
    }
    //if (it.value().find("hms e yield") != it.value().end()) {
    //  try {
    //  yield      = it.value()["hms e yield"].get<double>();
    //  runs.push_back(double(run_number));
    //  hms_yields.push_back(yield);
    //  if (it.value().find("hms e yield unc.") != it.value().end()) {
    //    hms_yield_uncs.push_back(it.value()["hms e yield unc."].get<double>());
    //  } else {
    //    hms_yield_uncs.push_back(0.0);
    //  }
    //  } catch(std::domain_error ) {
    //    continue;
    //    //you suck
    //  }
    //}
    //if (it.value().find("hms scaler yield") != it.value().end()) {
    //  run_number2 = std::stoi(it.key());
    //  yield2      = it.value()["hms scaler yield"].get<double>();
    //  runs2.push_back(double(run_number2));
    //  hms_yields2.push_back(scaler_reduce_factor*yield2);
    //  if (it.value().find("hms scaler yield unc.") != it.value().end()) {
    //    hms_yield_uncs2.push_back(scaler_reduce_factor*it.value()["hms scaler yield unc."].get<double>());
    //  } else {
    //    hms_yield_uncs2.push_back(0.0);
    //  }
    //}
    //if (it.value().find("hms scaler yield") != it.value().end()) {
    //  if (it.value().find("hms e yield") != it.value().end()) {
    //    runs3.push_back(double(run_number2));
    //    double yield_ratio = yield / yield2;
    //    yield_ratios.push_back(yield_ratio);
    //    if ((it.value().find("hms e yield unc.") != it.value().end()) &&
    //        (it.value().find("hms scaler yield unc.") != it.value().end())) {
    //      double y1 = it.value()["hms e yield unc."].get<double>();
    //      double y2 = it.value()["hms scaler yield unc."].get<double>();
    //      yield_ratio_uncs.push_back(yield_ratio*std::sqrt(y1 * y1/(yield*yield) + y2 * y2/(yield2*yield2)));
    //    } else {
    //      yield_ratio_uncs.push_back(0.0);
    //    }
    //  }
    //}
  }

  //// ------------------------------
  ////
  TGraphErrors* gr = nullptr;
  TMultiGraph*  mg = nullptr;
  TCanvas*      c  = nullptr;

  // ---------------------------------------------------------
  //
  // ---------------------------------------------------------
  c  = new TCanvas("chargeVsTime","charge",1200,800);
  c->Divide(1, 2);

  c->cd(1);
  mg = new TMultiGraph();
  mg->SetTitle("; run number ; triggers");
  //std::vector<double> zeros(runs.size());
  gr = new TGraphErrors(runs.size(), runs.data(), daq_triggers.data());
  gr->SetMarkerStyle(20);
  mg->Add(gr,"p");
  mg->Draw("a");

  c->cd(2);
  mg = new TMultiGraph();
  mg->SetTitle("; run number ; charge");
  gr = new TGraphErrors(runs2.size(), runs2.data(), charges.data());
  gr->SetMarkerStyle(23);
  gr->SetMarkerColor(2);
  mg->Add(gr,"p");

  //std::vector<double> zeros2(runs2.size());
  //gr = new TGraphErrors(runs2.size(), runs2.data(), hms_yields2.data(),zeros2.data(),hms_yield_uncs2.data());
  //mg->Add(gr,"p");
  mg->Draw("a");

  c->SaveAs("results/plot_charge_vs_time.png");
  c->SaveAs("results/plot_charge_vs_time.pdf");

  // ---------------------------------------------------------
  //
  // ---------------------------------------------------------
  c  = new TCanvas("chargeVsTime2","charge2",1200,800);
  c->Divide(1, 2);

  c->cd(1);
  mg = new TMultiGraph();
  mg->SetTitle("; days ; triggers");

  //std::vector<double> zeros(runs.size());
  
  gr = new TGraphErrors(time1.size(), time1.data(), daq_triggers.data());
  gr->SetMarkerStyle(20);
  mg->Add(gr,"p");
  mg->Draw("a");

  c->cd(2);
  mg = new TMultiGraph();
  mg->SetTitle("; days ; charge");

  gr = new TGraphErrors(time2.size(), time2.data(), charges.data());
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(2);
  mg->Add(gr,"p");

  mg->Draw("a");

  c->SaveAs("results/plot_charge_vs_time2.png");
  c->SaveAs("results/plot_charge_vs_time2.pdf");
  
  // ---------------------------------------------------------
  //
  // ---------------------------------------------------------
  c  = new TCanvas("chargeVsTime3","charge3",1200,800);
  c->Divide(1, 2);

  c->cd(1);
  mg = new TMultiGraph();
  mg->SetTitle("; days ; triggers");

  gr = new TGraphErrors(time1.size(), time1.data(), daq_triggers2.data());
  gr->SetMarkerStyle(20);
  mg->Add(gr,"p");
  mg->Draw("a");

  c->cd(2);
  mg = new TMultiGraph();
  mg->SetTitle("; days ; charge");

  gr = new TGraphErrors(time2.size(), time2.data(), charges2.data());
  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(2);
  mg->Add(gr,"p");

  mg->Draw("a");

  c->SaveAs("results/plot_charge_vs_time3.png");
  c->SaveAs("results/plot_charge_vs_time3.pdf");
}
