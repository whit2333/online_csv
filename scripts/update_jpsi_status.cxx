#include "nlohmann/json.hpp"
#include <cmath>
#include <iostream>

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

#include "monitor/DetectorDisplay.h"
#include "monitor/DisplayPlots.h"
#include "monitor/MonitoringDisplay.h"
R__LOAD_LIBRARY(libScandalizer.so)

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TSystem.h"

void update_jpsi_status() {

  std::cout << "Analyzing J/psi run info...\n";

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list_coin.json");
    json_input_file >> j;
  }
  json j2;
  {
    std::ifstream json_input_file("db2/jpsi_run_count_list.json");
    json_input_file >> j2;
  }
  json j3;
  {
    std::ifstream json_input_file("db2/run_list_extra.json");
    try {
      json_input_file >> j3;
    } catch (json::parse_error) {
      std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }
  json jsum;
  {
    std::ifstream json_input_file{"db2/jpsi_settings.json"};
    try {
      json_input_file >> jsum;
    } catch (json::parse_error) {
      {
        std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
        std::quick_exit(-127);
      }
    }
  }

  // set all counters to zero
  for (json::iterator iset = jsum.begin(); iset != jsum.end(); ++iset) {
    auto& sumjs = iset.value();

    sumjs["total_charge"]     = 0.;
    sumjs["jpsi_count"]       = 0.;
    sumjs["status"]           = "0%";
    std::vector<int> run_list = {};
    sumjs["good_run_list"]    = run_list;
  }

  // loop over all runs and add info to the summary where relevant
  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    int  run   = std::stoi(it.key());
    auto runjs = it.value();

    const std::string radiator = runjs["radiator"]["radiator_status"];
    const std::string target   = runjs["target"]["target_label"];
    const double      p_hms    = runjs["spectrometers"]["hms_momentum"].get<double>();
    const double      th_hms   = runjs["spectrometers"]["hms_angle"].get<double>();
    const double      p_shms   = runjs["spectrometers"]["shms_momentum"].get<double>();
    const double      th_shms  = runjs["spectrometers"]["shms_angle"].get<double>();
    if (j2[it.key()].find("good_total_charge") == j2[it.key()].end()) {
      continue;
    }
    // skip "empty " runs
    const double charge = j2[it.key()]["good_total_charge"].get<double>();
    if (charge < 1.) {
      continue;
    }
    // skip junk runs
    if (j3.count(it.key()) > 0 && j3[it.key()].find("comment") != j3[it.key()].end()) {
      std::string comment = j3[it.key()]["comment"];
      if (comment.find("junk") != std::string::npos || comment.find("Junk") != std::string::npos ||
          comment.find("JUNK") != std::string::npos) {
        continue;
      }
    }
    // loop over all settings
    for (json::iterator iset = jsum.begin(); iset != jsum.end(); ++iset) {
      auto& sumjs = iset.value();

      if (sumjs["radiator"] == radiator && sumjs["target"] == target &&
          fabs(sumjs["hms_angle"].get<double>() - th_hms) < .015 &&
          fabs(sumjs["hms_momentum"].get<double>() - p_hms) < .015 &&
          fabs(sumjs["shms_angle"].get<double>() - th_shms) < .015 &&
          fabs(sumjs["shms_momentum"].get<double>() - p_shms) < .015) {
        sumjs["total_charge"] = sumjs["total_charge"].get<double>() + charge;
        sumjs["jpsi_count"]   = sumjs["jpsi_count"].get<double>() +
                              j2[it.key()]["J/psi Good event count"].get<double>();
        sumjs["status"] = std::to_string(sumjs["total_charge"].get<double>() /
                                         sumjs["charge_goal"].get<double>() * 100) +
                          "%";
        auto run_list = sumjs["good_run_list"].get<std::vector<int>>();
        run_list.push_back(run);
        sumjs["good_run_list"] = run_list;

      } else {
        ; // do nothing
      }
    }
  }

  {
    std::ofstream json_output_file{"db2/jpsi_status.json"};
    std::ifstream json_input_file{"db2/jpsi_settings.json"};
    json_output_file << std::setw(4) << jsum << "\n";
  }
  std::cout << "Updated status written to db2/jpsi_status.json\n";
}
