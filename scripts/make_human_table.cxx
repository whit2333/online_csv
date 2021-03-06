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

void make_human_table() {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list_shms.json");
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
    } catch(json::parse_error)  {
      std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }
  std::ofstream runs_text("csv_runs.txt");

  // std::cout << j.dump(2);
  auto print_header = [](){
  // print header
    std::cout << "\n";
    fmt::print(" {:<5} ", "Run");
  fmt::print(" {:^5} ", "Target");
  fmt::print(" {:^7} ", "P_hms");
  fmt::print(" {:^7} ", "th_hms");
  fmt::print(" {:^7} ", "P_shms");
  fmt::print(" {:^7} ", "th_shms");
  fmt::print(" {:^21} ", "start time");
  fmt::print(" {:^21} ", "end time");
  fmt::print(" {:>9} ", "pi_count");
  fmt::print(" {:>9} ", "charge");
  fmt::print(" {:>9} ", "yield");
  fmt::print(" {:>9} ", "triggers");
  fmt::print(" {:<} ", "comment");
  std::cout << "\n";
  };

  std::cout << " runs : \n";

  std::string old_target = "";

  double p_hms   = 0.0;
  double th_hms  = 0.0;
  double p_shms  = 0.0;
  double th_shms = 0.0;

  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    auto runjs = it.value();

    std::string target_lab = runjs["target"]["target_label"].get<std::string>();
    if(target_lab != old_target) { 
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
    th_shms  = runjs["spectrometers"]["shms_angle"].get<double>();


    old_target = target_lab;

    runs_text << std::stoi(it.key()) << "\n";

    fmt::print(" {:<5} ", std::stoi(it.key()));
    fmt::print(" {:^5} ", target_lab);
    fmt::print(" {:>7.3f} ", runjs["spectrometers"]["hms_momentum"].get<double>());
    fmt::print(" {:>7.2f} ", runjs["spectrometers"]["hms_angle"].get<double>());
    fmt::print(" {:>7.3f} ", runjs["spectrometers"]["shms_momentum"].get<double>());
    fmt::print(" {:>7.2f} ", runjs["spectrometers"]["shms_angle"].get<double>());
    if (runjs.find("start_time") != runjs.end()) {
      fmt::print(" {:^21} ", runjs["start_time"].get<std::string>());
    } else {
      fmt::print(" {:21} ", "");
    }
    if (runjs.find("end_time") != runjs.end()) {
      fmt::print(" {:^21} ", runjs["end_time"].get<std::string>());
    }else{
      fmt::print(" {:21} ", "");
    }

    double total_charge = 0.0001;
    if (runjs.find("total_charge") != runjs.end()) {
      total_charge = runjs["total_charge"].get<double>() / 1000.0;
      //fmt::print("/ {:>10.1f} ", total_charge);
    }
    //else {
    //  fmt::print(" {:>11} ", "");
    //}
    if (j2.count(it.key()) != 0) {
      try {
        double pi_yield = j2[it.key()]["pion bg sub. counts"].get<double>();
        fmt::print(" {:>9.1f} ", pi_yield);
    if (runjs.find("total_charge") != runjs.end()) {
      //total_charge = runjs["total_charge"].get<double>() / 1000.0;
      fmt::print("/ {:<6.1f} ", total_charge);
    } else {
      fmt::print(" {:>11} ", "");
    }
        fmt::print(" {:>9.1f} ", pi_yield / total_charge);
        int n_events = j2[it.key()]["total trigger events"].get<int>();
        fmt::print(" {:>9d} ", n_events);
      } catch(std::domain_error ) {
        //you suck
      }
      //double live_time = j2[it.key()]["live time"].get<int>();
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
        } catch(std::domain_error ) {
          //you suck
        }
      }
    }
    fmt::print(" {:<} ", comment);
    std::cout << "\n";
  }
  print_header(); 


}
