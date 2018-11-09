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

void make_human_table(){

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

  //std::cout << j.dump(2);

  std::cout << " runs : \n";
  std::vector<int> runs;
  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    auto runjs = it.value();
    runs.push_back(std::stoi(it.key()));
    fmt::print(" {:>6} ",std::stoi(it.key()));
    //fmt::print(" \"{:^20}\n ",runjs["target"]["target_name"].get<std::string>());
    fmt::print(" {:^7} ",runjs["target"]["target_label"].get<std::string>());
    // runjs["target"]["target_name"]
    // << runjs["target"]["target_id"] 

    fmt::print(" {:>5.2f} ",runjs["spectrometers"]["hms_momentum"].get<double>());
    fmt::print(" {:>5.2f} ",runjs["spectrometers"]["hms_angle"].get<double>());
    fmt::print(" {:>5.2f} ",runjs["spectrometers"]["shms_momentum"].get<double>());
    fmt::print(" {:>5.2f} ",runjs["spectrometers"]["shms_angle"].get<double>());
    if(runjs.find("start_time") != runjs.end()) {
      fmt::print(" {} ",runjs["start_time"].get<std::string>());
      fmt::print(" {} ",runjs["end_time"].get<std::string>());
    }

    double total_charge = 0.0001;
    if(runjs.find("total_charge") != runjs.end()) {
      total_charge = runjs["total_charge"].get<double>()/1000.0;
      fmt::print(" {:>8.1f} ",total_charge);
    }
    if( j2.count(it.key()) != 0 ){
      //std::cout << j2[it.key()] << std::endl;
      double pi_yield = j2[it.key()]["pion bg sub. counts"].get<double>();
      fmt::print(" {:>5.1f} ",pi_yield);
      fmt::print(" {:>5.1f} ",pi_yield/total_charge);
      int n_events = j2[it.key()]["total trigger events"].get<int>();
      fmt::print(" {:>9d} ",n_events);
      
    }
    
    std::cout << "\n";

  }

}
