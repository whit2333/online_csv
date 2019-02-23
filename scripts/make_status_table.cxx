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

void print_hline() {
  std::cout << "-----------------------------------------------------------------------------------"
               "----------------------------------------";
}
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
    fmt::print(" {:>12} ", "Setting");
    fmt::print(" {:^7} ", "target");
    fmt::print(" {:^5} ", "rad");
    fmt::print(" {:>6} ", "P_hms");
    fmt::print(" {:>7} ", "th_hms");
    fmt::print(" {:>6} ", "P_shms");
    fmt::print(" {:>7} ", "th_shms");
    fmt::print(" {:>5} ", "n_runs");
    fmt::print(" {:>7} ", "J/psi");
    fmt::print(" {:^17} ", "yield");
    fmt::print(" {:>6} / {:<5} ", "charge", "goal");
    fmt::print(" {:>7} ", "status");
    std::cout << "\n";
  };

  std::string old_name     = "";
  double      total_charge = 0;
  double      total_goal   = 0;

  print_header();
  print_hline();
  for (auto setting : j) {
    const std::string name      = setting["name"];
    const std::string target    = setting["target"];
    const std::string radiator  = setting["radiator"];
    const double      p_hms     = setting["hms_momentum"].get<double>();
    const double      th_hms    = setting["hms_angle"].get<double>();
    const double      p_shms    = setting["shms_momentum"].get<double>();
    const double      th_shms   = setting["shms_angle"].get<double>();
    const int         n_runs    = setting["good_run_list"].size();
    const int         count     = setting["jpsi_count"].get<double>();
    const double      charge    = setting["total_charge"].get<double>();
    const double      goal      = setting["charge_goal"].get<double>();
    const double      yield     = (charge > 0) ? count / charge : 0;
    const double      yield_err = (charge > 0) ? sqrt(count) / charge : 0;
    const double      status    = charge / goal * 100.;

    if (name != old_name) {
      std::cout << "\n";
    }
    old_name = name;
    total_charge += charge;
    total_goal += goal;

    fmt::print(" {:>12} ", name);
    fmt::print(" {:^7} ", target);
    fmt::print(" {:^5} ", radiator);
    fmt::print(" {:>6.3f} ", p_hms);
    fmt::print(" {:>7.2f} ", th_hms);
    fmt::print(" {:>6.3f} ", p_shms);
    fmt::print(" {:>7.2f} ", th_shms);
    fmt::print(" {:>5} ", n_runs);
    fmt::print(" {:>7} ", count);
    fmt::print(" {:>8.4f} +- {:<8.4f} ", yield, yield_err);
    fmt::print(" {:>6.0f} / {:<5.0f} ", charge, goal);
    fmt::print(" {:>6.1f}%", status);
    std::cout << "\n";
  }
  print_hline();
  std::cout << "\n";
  fmt::print(" {:>97}  {:>5.0f} / {:<5.0f}  {:>6.1f}%", "TOTAL GOOD CHARGE", total_charge,
             total_goal, total_charge / total_goal * 100.);
  std::cout << "\n";

  std::cout << std::endl;
}
