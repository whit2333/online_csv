#ifdef __cpp_lib_filesystem
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#include <fstream>
#include "nlohmann/json.hpp"

#include "TObject.h"


void creat_standard_kinematics(){

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list_update.json");
    json_input_file >> j;
  }
  int runstart,runend;
  int i;
  std::ofstream out_file;
  std::cout<<"What range do you want? 1.custom 2.spring runs"<<"\n";
  std::cin>>i;
  switch (i)
  {
    case 1:
      std::cout<<"What range do you want? eg. 7593 7693"<<"\n";
      std::cin>>runstart>>runend;
      if (runstart>runend){std::cout<<"wrong input";
        std::quick_exit();}
      out_file.open("custom_standard.kinematics");
      break;
    case 2:
      runstart = 7593;
      runend = 7830;
      out_file.open("spring_standard.kinematics");
      break;
    default:
      std::cout<<"wrong input"<<"\n";
      std::quick_exit();
      break;
  }

  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    if (std::stoi(it.key())>=runstart && std::stoi(it.key()) <=runend){
    out_file << it.key() << "\n";
    out_file << "gpbeam = " << it.value()["run_info"]["beam_energy"].get<double>()<<"\n";
    out_file << "gtargmass = " << it.value()["target"]["target_mass_amu"].get<double>()<<"\n";
    out_file << "htheta_lab = " << it.value()["spectrometers"]["hms_angle"].get<double>()<<"\n";
    out_file << "ptheta_lab = " << it.value()["spectrometers"]["shms_angle"].get<double>()<<"\n";
    out_file << "hpcentral = " << std::abs(it.value()["spectrometers"]["hms_momentum"].get<double>())<<"\n";
    out_file << "ppcentral = " << std::abs(it.value()["spectrometers"]["shms_momentum"].get<double>())<<"\n";
    out_file << "ppartmass = " << 0.93827231<<"\n";
    out_file << "hpartmass = " << 0.0005109<<"\n";
    
    
    out_file<<"\n";
    
    }
  }
}
