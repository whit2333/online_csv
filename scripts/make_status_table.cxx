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

void make_status_table() {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/jpsi_status.json");
    json_input_file >> j;
  }

  // std::cout << j.dump(2);
  auto print_header = []() {
    // print header
    std::cout << "\n";
    fmt::print(" {:>27} ", "Setting");
    fmt::print(" {:^7} ", "P_hms");
    fmt::print(" {:^7} ", "th_hms");
    fmt::print(" {:^7} ", "P_shms");
    fmt::print(" {:^7} ", "th_shms");
    fmt::print(" {:^7} ", "count");
    fmt::print(" {:^18} ", "yield");
    fmt::print(" {:>6} ", "charge");
    fmt::print(" {:<6} ", "/ goal");
    fmt::print(" {:<7} ", "status");
    std::cout << "\n";
  };

  std::cout << " runs : \n";

  std::string old_target          = "";
  std::string old_radiator_status = "";

  print_header();
  for (auto setting : j) {
    const std::string name      = setting["name"];
    const double      p_hms     = setting["hms_momentum"].get<double>();
    const double      th_hms    = setting["hms_angle"].get<double>();
    const double      p_shms    = setting["shms_momentum"].get<double>();
    const double      th_shms   = setting["shms_angle"].get<double>();
    const double      count     = setting["jpsi_count"].get<double>();
    const double      charge    = setting["total_charge"].get<double>();
    const double      goal      = setting["charge_goal"].get<double>();
    const double      yield     = (charge > 0) ? count / charge : 0;
    const double      yield_err = (charge > 0) ? sqrt(count) / charge : 0;
    const double      status    = charge / goal * 100.;

    fmt::print(" {:>27} ", name);
    fmt::print(" {:>7.3f} ", p_hms);
    fmt::print(" {:>7.2f} ", th_hms);
    fmt::print(" {:>7.3f} ", p_shms);
    fmt::print(" {:>7.2f} ", th_shms);
    fmt::print(" {:>7} ", count);
    fmt::print(" {:>8.2f} +- {:<8.2f}", yield, yield_err);
    fmt::print(" {:>8.2f} / {:<8.2f}", charge, goal);
    fmt::print(" {:>5.1f}%", status);
    std::cout << "\n";
  }
  std::cout << std::endl;
}
