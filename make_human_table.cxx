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
    std::ifstream json_input_file("db2/run_list.json");
    json_input_file >> j;
  }
  json j2;
  {
    std::ifstream json_input_file("db2/run_count_list.json");
    json_input_file >> j2;
  }

  // std::cout << j.dump(2);

  std::cout << " runs : \n";
  std::vector<int> runs;
  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    auto runjs = it.value();
    runs.push_back(std::stoi(it.key()));
    fmt::print(" {:<5} ", std::stoi(it.key()));
    fmt::print(" {:^5} ", runjs["target"]["target_label"].get<std::string>());
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
      fmt::print(" {:>11.1f} ", total_charge);
    }
    if (j2.count(it.key()) != 0) {
      double pi_yield = j2[it.key()]["pion bg sub. counts"].get<double>();
      fmt::print(" {:>9.0f} ", pi_yield);
      fmt::print(" {:>9.1f} ", pi_yield / total_charge);

      int n_events = j2[it.key()]["total trigger events"].get<int>();
      fmt::print(" {:>9d} ", n_events);

      //double live_time = j2[it.key()]["live time"].get<int>();
    }
    std::cout << "\n";
  }

  // print header
  fmt::print(" {:<5} ", "Run");
  fmt::print(" {:^5} ", "Target");
  fmt::print(" {:^7} ", "P_hms");
  fmt::print(" {:^7} ", "th_hms");
  fmt::print(" {:^7} ", "P_shms");
  fmt::print(" {:^7} ", "th_shms");
  fmt::print(" {:^21} ", "start time");
  fmt::print(" {:^21} ", "end time");
  fmt::print(" {:>11} ", "charge");
  fmt::print(" {:>9} ", "pi_count");
  fmt::print(" {:>9} ", "yield");
  fmt::print(" {:>9} ", "triggers");
  std::cout << "\n";

}
