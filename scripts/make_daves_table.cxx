#include <algorithm>
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

void make_daves_table() {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list.json");
    json_input_file >> j;
  }
  json j2;
  {
    std::ifstream json_input_file("db2/run_count_list.json");
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

  // std::cout << j.dump(2);


  for (json::iterator it = j.begin(); it != j.end(); ++it) {

    auto runjs = it.value();
    // Dave's header
    // Run Type  Target P_HMS    Theta_HMS	 P_SHMS	   Theta_SHMS	 BCM4b_charge  #Trig6(k) SHMS_trkeff  HMS_trkeff CPU-Live-time	  #good-coin	  Current    Comment
    // Get the kinematic setting, eg 13-5
    std::string kin_setting = "";
    if (j3.count(it.key()) != 0) {
      if (j3[it.key()].find("kinematic") != j3[it.key()].end()) {
        try {
          kin_setting = j3[it.key()]["kinematic"].get<std::string>();
        } catch(std::domain_error ) {
          //you suck
          kin_setting = "";
        }
      }
    }

    //todo  : safety try/catch
    int run_number = std::stoi(it.key());
    //fmt::print(" {:<5} ", std::stoi(it.key()));
    
    std::string target_label = runjs["target"]["target_label"].get<std::string>();
    //fmt::print(" {:^5} ", runjs["target"]["target_label"].get<std::string>());

    double P_hms   = runjs["spectrometers"]["hms_momentum"].get<double>();
    double th_hms  = runjs["spectrometers"]["hms_angle"].get<double>();
    double P_shms  = runjs["spectrometers"]["shms_momentum"].get<double>();
    double th_shms = runjs["spectrometers"]["shms_angle"].get<double>();

    std::string start_time = "";
    std::string end_time   = "";
    if (runjs.find("start_time") != runjs.end()) {
      // todo: add try/catch
      start_time = runjs["start_time"].get<std::string>();
    }
    if (runjs.find("end_time") != runjs.end()) {
      // todo: add try/catch
      end_time = runjs["end_time"].get<std::string>();
    }

    double beam_current = runjs["spectrometers"]["shms_angle"].get<double>();

    double total_charge = 0.00000001;
    if (runjs.find("total_charge") != runjs.end()) {
      // todo: add try/catch
      total_charge = runjs["total_charge"].get<double>() / 1000.0;
      //fmt::print(" {:>11.1f} ", total_charge);
    }

    double pi_yield = 0.0;
    double pi_count = 0.0;
    int    n_events = 0;
    double time_2uA_cut = 0;
    double charge_2uA_cut = 0;
    double current_2uA_cut = 0;
    if (j2.count(it.key()) != 0) {
      if (j2[it.key()].count("time 1MHz 2u cut") != 0) {
        time_2uA_cut = j2[it.key()]["time 1MHz 2u cut"].get<double>();
      }
      if (j2[it.key()].count("charge bcm4b 2u cut") != 0) {
        charge_2uA_cut = j2[it.key()]["charge bcm4b 2u cut"].get<double>();
      }
      current_2uA_cut =  charge_2uA_cut/time_2uA_cut;
      try {
        pi_count = j2[it.key()]["pion bg sub. counts"].get<double>();
        pi_yield = pi_count/total_charge;
        //fmt::print(" {:>9.1f} ", pi_yield);
        //fmt::print(" {:>9.1f} ", pi_yield / total_charge);
        n_events = j2[it.key()]["total trigger events"].get<int>();
        //fmt::print(" {:>9d} ", n_events);
      } catch(std::domain_error ) {
        //you suck
      }
      //double live_time = j2[it.key()]["live time"].get<int>();
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

    fmt::print("{:<4} ", kin_setting);
    fmt::print("{:>4d} ", run_number);
    fmt::print("{:^6} ", target_label);
    fmt::print("{:>7.3f} ", P_hms);
    fmt::print("{:<9.3f} ", th_hms);
    fmt::print("{:>7.3f} ", P_shms);
    fmt::print("{:<10.3f} ", th_shms);
    fmt::print("{:>7.2f} ", total_charge);
    fmt::print("{:>8d} ", n_events);
    fmt::print("{:^12} ", "SHMS_trkeff");
    fmt::print("{:^12} ", "HMS_trkeff");
    fmt::print("{:^12} ", "CPU-Live-time");
    fmt::print("{:^11} ", pi_count);
    fmt::print("{:^12.1f} ",current_2uA_cut);
    fmt::print("{:<20} ", comment.substr(0,std::min((unsigned long)20,comment.length())));
    //fmt::print("{:^21} ", start_time);
    //fmt::print("{:^21} ", end_time);

    std::cout << "\n";
  }

  // Run Type  Target P_HMS Theta_HMS P_SHMS Theta_SHMS BCM4b_charge #Trig6(k) SHMS_trkeff  HMS_trkeff CPU-Live-time	  #good-coin	  Current    Comment
  // print header
  fmt::print("{:<4} ", "Kin");
  fmt::print("{:<4} ", "Run");
  fmt::print("{:^6} ", "Target");
  fmt::print("{:>7} ", "P_HMS");
  fmt::print("{:<9} ", "Theta_HMS");
  fmt::print("{:>7} ", "P_SHMS");
  fmt::print("{:<10} ", "Theta_SHMS");
  fmt::print("{:>7} ", "charge");
  fmt::print("{:>8} ", "Trig6(k)");
  fmt::print("{:^12} ", "SHMS_trkeff");
  fmt::print("{:^12} ", "HMS_trkeff");
  fmt::print("{:^12} ", "CPU-Live-time");
  fmt::print("{:^11} ", "#good-coin");
  fmt::print("{:^12} ", "Current");
  fmt::print("{:<20} ",   "comment");
  //fmt::print("{:^21} ", "start time");
  //fmt::print("{:^21} ", "end time");
  //fmt::print("{:>9} ", "yield");
  //fmt::print("{:>9} ", "triggers");
  std::cout << "\n";

}
