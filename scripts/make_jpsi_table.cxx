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

void make_jpsi_table() {

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
  json j_hms;
  {
    std::ifstream json_input_file("db2/hms_run_count_list.json");
    json_input_file >> j_hms;
  }
  json j_shms;
  {
    std::ifstream json_input_file("db2/shms_run_count_list.json");
    json_input_file >> j_shms;
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
  std::ofstream runs_text("jpsi_runs.txt");

  // std::cout << j.dump(2);
  auto print_header = []() {
    // print header
    std::cout << "\n";
    fmt::print(" {:<5} ", "Run");
    fmt::print(" {:^5} ", "Target");
    fmt::print(" {:^7} ", "radiator");
    fmt::print(" {:^7} ", "P_hms");
    fmt::print(" {:^7} ", "th_hms");
    fmt::print(" {:^7} ", "P_shms");
    fmt::print(" {:^7} ", "th_shms");
    fmt::print(" {:^21} ", "start time");
    fmt::print(" {:^21} ", "end time");
    fmt::print(" {:>9} ", "count");
    fmt::print(" {:>9} ", "charge");
    fmt::print(" {:>9} ", "yield");
    fmt::print(" {:>9} ", "triggers");
    fmt::print(" {:<} ", "comment");
    std::cout << "\n";
  };

  std::cout << " runs : \n";

  std::string old_target          = "";
  std::string old_radiator_status = "";

  double p_hms   = 0.0;
  double th_hms  = 0.0;
  double p_shms  = 0.0;
  double th_shms = 0.0;

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
        (th_hms != runjs["spectrometers"]["hms_angle"].get<double>()) ||
        (p_shms != runjs["spectrometers"]["shms_momentum"].get<double>()) ||
        (th_shms != runjs["spectrometers"]["shms_angle"].get<double>())) {
      print_header();
    }

    p_hms   = runjs["spectrometers"]["hms_momentum"].get<double>();
    th_hms  = runjs["spectrometers"]["hms_angle"].get<double>();
    p_shms  = runjs["spectrometers"]["shms_momentum"].get<double>();
    th_shms = runjs["spectrometers"]["shms_angle"].get<double>();

    old_target          = target_lab;
    old_radiator_status = radiator_status;

    runs_text << std::stoi(it.key()) << "\n";

    fmt::print(" {:<5} ", std::stoi(it.key()));
    fmt::print(" {:^5} ", target_lab);
    fmt::print(" {:^7} ", radiator_status);
    fmt::print(" {:>7.3f} ", runjs["spectrometers"]["hms_momentum"].get<double>());
    fmt::print(" {:>7.2f} ", runjs["spectrometers"]["hms_angle"].get<double>());
    fmt::print(" {:>7.3f} ", runjs["spectrometers"]["shms_momentum"].get<double>());
    fmt::print(" {:>7.2f} ", runjs["spectrometers"]["shms_angle"].get<double>());
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
      // fmt::print("/ {:>10.1f} ", total_charge);
    }
    // else {
    //  fmt::print(" {:>11} ", "");
    //}
    if (j_hms.count(it.key()) != 0) {
      auto   rl_hms = j_hms[it.key()];
      double n_hms  = rl_hms["count_e"].get<double>();
      double hms_yield =
          n_hms * rl_hms["ps_factor"].get<double>() / rl_hms["good_total_charge"].get<double>();
      double hms_unc = (n_hms > 0) ? hms_yield / sqrt(n_hms) : 0;
      if (j2.count(it.key()) != 0) {
        try {

          double shms_yield = j2[it.key()]["J/psi Good event count"].get<double>();
          // double pi_yield = j2[it.key()]["pion bg sub. counts"].get<double>();
          fmt::print(" {:>9.1f} ", shms_yield);
          if (j2[it.key()].find("good_total_charge") != j2[it.key()].end()) {
            total_charge = j2[it.key()]["good_total_charge"].get<double>();
            // total_charge = runjs["total_charge"].get<double>() / 1000.0;
            fmt::print("/ {:<6.1f} ", total_charge);
          } else {
            fmt::print(" {:>11} ", "");
          }
          fmt::print(" {:>9.5f} ", shms_yield / total_charge);
          fmt::print(" +- {:>8.5f} ", std::sqrt(shms_yield) / total_charge);
          int n_events = j2[it.key()]["total trigger events"].get<int>();
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
