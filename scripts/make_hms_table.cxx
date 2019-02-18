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

#include "TObject.h"

void make_hms_table() {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list_coin.json");
    json_input_file >> j;
    // also merge in single-arm runs
    json j_single_arm;
    {
      std::ifstream json_input_file("db2/run_list_hms.json");
      json_input_file >> j_single_arm;
      j.merge_patch(j_single_arm);
    }
  }
  json j2;
  {
    std::ifstream json_input_file("db2/hms_run_count_list.json");
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

  // std::cout << j.dump(2);
  auto print_header = []() {
    // print header
    std::cout << "\n";
    fmt::print(" {:<5} ", "Run");
    fmt::print(" {:^5} ", "Target");
    fmt::print(" {:^7} ", "radiator");
    fmt::print(" {:^7} ", "P_hms");
    fmt::print(" {:^7} ", "th_hms");
    fmt::print(" {:^20} ", "start time");
    fmt::print(" {:^19} ", "end time");
    fmt::print(" {:>9} ", "count");
    fmt::print(" {:<7} ", "/ charge");
    fmt::print(" {:^24} ", "corrected yield");
    fmt::print(" {:<9} ", "triggers");
    fmt::print(" {:<} ", "comment");
    std::cout << "\n";
  };

  std::cout << " runs : \n";

  std::string old_target          = "";
  std::string old_radiator_status = "";

  double p_hms  = 0.0;
  double th_hms = 0.0;

  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    auto runjs = it.value();

    std::string target_lab = runjs["target"]["target_label"].get<std::string>();
    if (target_lab != old_target) {
      print_header();
    }
    std::string radiator_status = runjs["radiator"]["radiator_status"].get<std::string>();
    if (radiator_status != old_radiator_status) {
      print_header();
    }

    if ((p_hms != runjs["spectrometers"]["hms_momentum"].get<double>()) ||
        (th_hms != runjs["spectrometers"]["hms_angle"].get<double>())) {
      print_header();
    }

    p_hms  = runjs["spectrometers"]["hms_momentum"].get<double>();
    th_hms = runjs["spectrometers"]["hms_angle"].get<double>();

    old_target          = target_lab;
    old_radiator_status = radiator_status;

    fmt::print(" {:<5} ", std::stoi(it.key()));
    fmt::print(" {:^5} ", target_lab);
    fmt::print(" {:^7} ", radiator_status);
    fmt::print(" {:>7.3f} ", runjs["spectrometers"]["hms_momentum"].get<double>());
    fmt::print(" {:>7.2f} ", runjs["spectrometers"]["hms_angle"].get<double>());
    if (runjs["run_info"].find("start_time") != runjs["run_info"].end()) {
      fmt::print(" {:^21} ", runjs["run_info"]["start_time"].get<std::string>());
    } else {
      fmt::print(" {:21} ", "");
    }
    if (runjs["run_info"].find("stop_time") != runjs["run_info"].end()) {
      fmt::print(" {:^21} ", runjs["run_info"]["stop_time"].get<std::string>());
    } else {
      fmt::print(" {:21} ", "");
    }

    double total_charge = 0.0001;
    if (runjs.find("good_total_charge") != runjs.end()) {
      total_charge = runjs["good_total_charge"].get<double>() / 1000.0;
    }
    if (j2.count(it.key()) != 0) {
      try {

        const double count     = j2[it.key()]["count_e"].get<double>();
        const double ps_factor = j2[it.key()]["ps_factor"];
        fmt::print(" {:>9.1f} ", count);
        if (j2[it.key()].find("good_total_charge") != j2[it.key()].end()) {
          total_charge = j2[it.key()]["good_total_charge"].get<double>();
          fmt::print("/ {:<6.1f} ", total_charge);
        } else {
          fmt::print(" {:>11} ", "");
        }
        fmt::print(" {:>9.5f} ", (total_charge > 0) ? count * ps_factor / total_charge : 0);
        fmt::print(" +- {:>8.5f} ",
                   (total_charge > 0) ? std::sqrt(count) * ps_factor / total_charge : 0);
        int n_events = j2[it.key()]["count_raw"].get<int>();
        fmt::print(" {:>9d} ", n_events);
      } catch (std::domain_error) {
        // you suck
        // std::cout << "you suck\n";
      }
      // double live_time = j2[it.key()]["live time"].get<int>();
    } else {
      fmt::print(" {:>9} ", "");
      fmt::print(" {:>9} ", "");
      fmt::print(" {:>10} ", "");
    }
    std::string comment;
    if (j3.count(it.key()) != 0) {
      if (j3[it.key()].find("comment") != j3[it.key()].end()) {
        try {
          comment = j3[it.key()]["comment"].get<std::string>();
        } catch (std::domain_error) {
          // you suck
        }
      }
    }
    fmt::print(" {:<} ", comment);
    std::cout << "\n";
  }
  print_header();
}
