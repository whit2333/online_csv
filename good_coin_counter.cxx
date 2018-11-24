#include "nlohmann/json.hpp"
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libGenVector.so)

#include "THStack.h"
using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// VecOps::RVec is like std::vector with some extra bells and whistles
using inters = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;
using nlohmann::json;

void good_coin_counter(int RunNumber = 6018, int nevents = -1, int prompt = 1, int update = 1, int default_count_goal = 30000) {

  using nlohmann::json;
  json j;
  {
    std::ifstream json_input_file("db2/run_list.json");
    try {
      json_input_file >> j;
    } catch(json::parse_error)  {
      std::cerr << "error: json file, db2/run_list.json, is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }

  auto runnum_str = std::to_string(RunNumber);
  if (j.find(runnum_str) == j.end()) {
    std::cout << "Run " << RunNumber << " not found in ddb2/run_list.json\n";
    std::cout << "Check that run number and replay exists. \n";
    std::cout << "If problem persists please contact Whit (717-341-1080)\n";
    std::cout << "In the meantime use: good_coin_counter_old.cxx \n";
  }
  double P0_shms_setting =
      j[runnum_str]["spectrometers"]["shms_momentum"].get<double>();
  double P0_shms = std::abs(P0_shms_setting);

  std::string coda_type = "COIN";
  std::string rootfile = "ROOTfiles_csv/";
  rootfile += std::string("coin_replay_production_");
  rootfile +=
      std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";

  //{
  //  TFile file(rootfile.c_str());
  //  if (file.IsZombie()) {
  //    std::cout << " Did your replay finish?  Check that the it is done before running this script.\n";
  //    return;
  //  }
  //}
  bool found_good_file = false; 
  if(!gSystem->AccessPathName(rootfile.c_str())) {
    TFile file(rootfile.c_str());
    if (file.IsZombie()) {
      std::cout << rootfile << " is a zombie.\n";
      std::cout << " Did your replay finish?  Check that the it is done before running this script.\n";
      //return;
    } else {
      std::cout <<  " using : " << rootfile << "\n";
      found_good_file = true;
    }
  }
  if(!found_good_file) {
    rootfile = "ROOTfiles_online/";
    rootfile += std::string("coin_replay_production_");
    rootfile += std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";

    if(!gSystem->AccessPathName(rootfile.c_str())) {
      TFile file(rootfile.c_str());
      if (file.IsZombie()) {
        std::cout << rootfile << " is a zombie.\n";
        std::cout << " Did your replay finish?  Check that the it is done before running this script.\n";
      } else {
        found_good_file = true;
        std::cout <<  " using : " << rootfile << "\n";
      }
    }
  }
  if(!found_good_file) {
    std::cout << " Error: suitable root file not found\n";
    return;
  }
  // new TBrowser;
  //

  ROOT::EnableImplicitMT(24);

  Pvec4D Pbeam(0, 0, 10.598, 0.000511);

  // Detector tree 
  ROOT::RDataFrame d("T", rootfile);
  // HMS Scaler tree
  ROOT::RDataFrame d_sh("TSH", rootfile);
  
  auto bcm4b_charge        = d_sh.Max("H.BCM4B.scalerChargeCut");
  auto el_real_scaler      = d_sh.Max("H.hEL_REAL.scaler");
  auto time_1MHz           = d_sh.Max("H.1MHz.scalerTime");
  auto time_1MHz_cut       = d_sh.Max("H.1MHz.scalerTimeCut");
  auto total_charge        = bcm4b_charge;

  std::string hpdelta = "P.gtr.dp > -10 && P.gtr.dp < 20 && "
                        "H.gtr.dp > -10 && H.gtr.dp < 10";
  std::string epiCut = "P.aero.npeSum > 1.0 && P.cal.eprtracknorm < 0.2 && "
                       "H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 && "
                       "H.cal.etottracknorm < 2.0 && P.cal.etottracknorm<1.0";
  // && H.cal.eprtracknorm  > 0.2

  auto c_n_events_total = d.Count();

  auto c_n_events_coin =
      d.Filter("(fEvtHdr.fEvtType == 6) || (fEvtHdr.fEvtType == 4)").Count();

  auto d0 = d.Define("hms_e_EoverP",
                     [](doublers &EOverP) {
                       return EOverP[EOverP > 0.5 && EOverP < 1.8];
                     },
                     {"H.cal.etottracknorm"})
                .Define("hms_e_EoverP_nGood",
                        [](doublers &EOverP) { return (int)EOverP.size(); },
                        {"hms_e_EoverP"})
                .Define("shms_pi_preshower",
                        [](doublers &EOverP) { return EOverP[EOverP < 0.3]; },
                        {"P.cal.eprtracknorm"})
                .Define("shms_pi_preshower_nGood",
                        [](doublers &EOverP) { return (int)EOverP.size(); },
                        {"shms_pi_preshower"});

  auto d2 = d0.Filter([](double npe) { return npe > 1.0; }, {"H.cer.npeSum"});
  auto d3 = d2.Filter([](int n_good) { return n_good > 0; },
                      {"shms_pi_preshower_nGood"});

  auto d_spec_cuts = d2.Filter(hpdelta).Filter(epiCut);

  // Find the coin peak
  auto h_coin_time = d_spec_cuts.Histo1D(
      {"coin_time", "coin_time", 800, 0, 100}, "CTime.ePiCoinTime_ROC2");
  int coin_peak_bin = h_coin_time->GetMaximumBin();
  double coin_peak_center = h_coin_time->GetBinCenter(coin_peak_bin);
  // std::cout << " coin peak at " << coin_peak_center << " ns\n";

  // Cut on the coin time
  auto d_coin_time_cut = d_spec_cuts.Filter(
      [=](double coin_time) {
        return std::abs(coin_time - coin_peak_center) < 2.0;
      },
      {"CTime.ePiCoinTime_ROC2"});
  auto h_coin_time2 = d_coin_time_cut.Histo1D(
      {"coin_time2", "coin_time2", 800, 0, 100}, "CTime.ePiCoinTime_ROC2");

  auto d_random = d_spec_cuts.Filter(
      [=](double coin_time) {
        return std::abs(coin_time - coin_peak_center - 28) < 10.0;
      },
      {"CTime.ePiCoinTime_ROC2"});
  auto h_coin_time3 = d_random.Histo1D(
      {"coin_time2", "coin_time2", 800, 0, 100}, "CTime.ePiCoinTime_ROC2");

  auto d2_2 = d_coin_time_cut.Filter("fEvtHdr.fEvtType == 2");
  auto d2_coin = d_coin_time_cut.Filter(
      "(fEvtHdr.fEvtType == 6) || (fEvtHdr.fEvtType == 4)");

  // Histograms lazy evaluated -> Nothing happens until they are used below (eg
  // Draw() is called)
  auto h_EOverP_0 = d0.Histo1D<doublers>(
      {"hms_e_EoverP_0", "SHMS total shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.etottracknorm");
  auto h_EOverP_2 = d2.Histo1D<doublers>(
      {"hms_e_EoverP_2", "SHMS total shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.etottracknorm");
  auto h_EOverP_3 = d3.Histo1D<doublers>(
      {"hms_e_EoverP_3", "SHMS total shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.etottracknorm");

  auto h_EprOverP_0 = d0.Histo1D<doublers>(
      {"hms_e_EoverP_0", "SHMS pre-shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.eprtracknorm");
  auto h_EprOverP_2 = d2.Histo1D<doublers>(
      {"hms_e_EoverP_2", "SHMS pre-shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.eprtracknorm");
  auto h_EprOverP_3 = d3.Histo1D<doublers>(
      {"hms_e_EoverP_3", "SHMS pre-shower; SHMS E/P", 100, 0.05, 1.8},
      "P.cal.eprtracknorm");

  // auto h_EOverP_1       = d0.Histo1D<doublers>({"hms_e_EoverP_1","E/P cut;
  // HMS E/P",100,0.05,1.8},"hms_e_EoverP"); auto h_EOverP_nGood_0 =
  // d0.Histo1D({"h_EOverP_nGood_0","h_EOverP_nGood_0",10,0,10},"P.cal.ntracks");
  // auto h_EOverP_nGood_1 =
  // d0.Histo1D({"h_EOverP_nGood_1","h_EOverP_nGood_0",10,0,10},"hms_e_EoverP_nGood");

  auto h_event_type =
      d0.Histo1D({"event_type", "event_type", 10, 0, 10}, "fEvtHdr.fEvtType");


  auto d_hgc_cut = d2_coin.Filter(
      [=](double npe, double dp) {
        double p_track = P0_shms * (100.0 + dp) / 100.0;
        // no cerenkov cut needed when momentum is below 2.8 GeV/c
        if (p_track < 2.8) {
          return true;
        }
        return npe > 1.0;
      },
      {"P.hgcer.npeSum", "P.gtr.dp"});

  auto hms_electron_counts = d2.Count();
  auto hms_electron_counts2 = d2_2.Count();
  auto coin_counts = d2_coin.Count();

  auto rand_bg_count = d_random.Count();

  auto pion_count = d_hgc_cut.Count();

  double random_bg    = double(*rand_bg_count) / 5.0;
  double hms_e_yield  = (*hms_electron_counts) / (*total_charge);
  double hms_e_yield2 = (*hms_electron_counts2) / (*total_charge);
  double coin_yield   = (*coin_counts) / (*total_charge);

  double pi_K_ratio =
      double(*pion_count) / (double(*coin_counts) - double(*pion_count));
  double kaon_counts = (double(*coin_counts) - double(*pion_count));

  double pion_corrected =
      double(double(*pion_count) - random_bg * double(*pion_count) / (double(*coin_counts)));

  std::cout << " ----------------------------------------------    \n";
  std::cout << " # of good coin  = " << int(pion_corrected) << "    \n";
  std::cout << "    out of  " << *c_n_events_total << " total triggers\n";
  std::cout << "        and " << *c_n_events_coin << " coin triggers\n";
  std::cout << "  coin yield: " << coin_yield << " events/mC \n";

  json jruns;
  {
    std::ifstream input_file("db2/run_count_list.json");
    try {
      input_file >> jruns;
    } catch(json::parse_error)  {
      std::cerr << "error: json file is incomplete or has broken syntax.\n";
      std::quick_exit(-127);
    }
  }


  // std::cout << jruns << std::endl;;
  // Open json db file to get current values (if run already exists, replace)
  // fill the run values
  json j_current_run;
  std::string run_str = std::to_string(RunNumber);
  j_current_run["total trigger events"] = int(*c_n_events_total);
  j_current_run["coin trigger events"] = int(*c_n_events_coin);
  j_current_run["pi+K counts"] = int(*coin_counts);
  j_current_run["random background"] = double(random_bg);
  j_current_run["pion counts"] = int(*pion_count);
  j_current_run["pi/K ratio"] = pi_K_ratio;
  j_current_run["kaon counts"] = kaon_counts;
  j_current_run["pion bg sub. counts"] = pion_corrected;
  j_current_run["kaon bg sub. counts"] = double(
    double(kaon_counts) - random_bg * kaon_counts / (double(*coin_counts)));
  // std::cout << std::setw(4) << j_current_run << "\n";

  jruns[run_str] = nlohmann::json::parse(j_current_run.dump());
  // jruns[run_str] = nlohmann::json::parse(j_current_run.dump());
  if( update ) {
    std::ofstream json_output_file("db2/run_count_list.json");
    json_output_file << std::setw(4) << jruns << "\n";
  }

  int count_goal = default_count_goal;

  if (prompt) {
    std::cout << "----------------------------------------------------------\n";
    std::cout << "Reference the run plan for this setting found on the wiki\n"
    "       https://hallcweb.jlab.org/wiki/index.php/CSV_Fall_2018_Run_Plan\n";
    std::cout << "----------------------------------------------------------\n";
    std::cout << "Please enter **total count** goal for this setting. \n";
    std::cout << "   Desired count [default=30000] : ";
    // std::cin >> count_goal ;
    // int number = 0;
    if (std::cin.peek() == '\n') { // check if next character is newline
      // count_goal = 30000; //and assign the default
    } else if (!(std::cin >> count_goal)) { // be sure to handle invalid input
      std::cout << "Invalid input.\n";
      // error handling
    }
    std::cout << "\n";
  }

  //std::cout << "count goal : " << count_goal << '\n';
  //std::cout << "time       : " << (*time_1MHz_cut) / 60.0 << "\n";
  //std::cout << "charge     : " << (*total_charge) << "\n";
  double n_seconds      = double(*time_1MHz_cut);
  int    nev_tot        = (*c_n_events_total);
  double goal_Nevents   = (count_goal / pion_corrected) * nev_tot;
  double time_remaining = (count_goal * n_seconds) / pion_corrected - (n_seconds);
  double charge_goal    = count_goal * (*total_charge) / (pion_corrected);

  std::cout << "----------------------------------------------------------\n";
  std::cout << " N events to reach goal  : " << goal_Nevents / 1000000 << "M events\n";
  std::cout << " Charge   to reach goal  : " << charge_goal / 1000.0 << " mC\n";


  if( update ) {
    std::string cmd =
    "caput hcRunPlanChargeGoal " + std::to_string(charge_goal / 1000.0) + " &> /dev/null ";
    system(cmd.c_str());

    cmd = "caput hcRunPlanNTrigEventsGoal " + std::to_string(goal_Nevents) + " &> /dev/null ";
    system(cmd.c_str());

    cmd = "caput hcRunPlanCountGoal " + std::to_string(count_goal) + " &> /dev/null ";
    system(cmd.c_str());
  }

  // std::cout << " pions+kaons   : " << *coin_counts << "\n";
  // std::cout << " pions         : " << *pion_count << "\n";
  // std::cout << " random        : " << random_bg << "\n";
  // std::cout << " ratio         : "  << pi_K_ratio << " \n";

  // std::cout << " counts  : " << *hms_electron_counts << "\n";
  // std::cout << " charge  : " << *total_charge << " uC\n";
  // std::cout << " yield   : " << (*hms_electron_counts)  << " cnts, " <<
  // hms_e_yield  << " cnts/uC\n"; std::cout << " singles : " <<
  // (*hms_electron_counts2) << " cnts, " << hms_e_yield2 << " cnts/uC\n";
  // std::cout << " coin    : " << (*coin_counts) << " cnts, " << coin_yield <<
  // " cnts/uC\n";

  // auto s_dc_x_fp           = d2.Histo1D({"s_dc_x_fp ","xy fp; x",100,-50,50},
  // "P.dc.x_fp"); auto s_dc_y_fp           = d2.Histo1D({"s_dc_y_fp ","xy fp;
  // y",100,-50,50}, "P.dc.y_fp"); auto s_dc_xp_fp          =
  // d2.Histo1D({"s_dc_xp_fp","xy fp; xp",100,-50,50},"P.dc.xp_fp"); auto
  // s_dc_yp_fp          = d2.Histo1D({"s_dc_yp_fp","xy fp;
  // xp",100,-50,50},"P.dc.yp_fp");

  // -----------------------------------------------------------
  //
  TCanvas *c = nullptr;
  int b1 = 0;
  int b2 = 0;
  double hmax = 0.0;
  THStack *hs = nullptr;
  // TLatex latex;

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  c->Divide(2, 2);
  c->cd(1);

  // This call starts the loop over the data.
  // A DrawCopy is used so that the histogram is not deleted at the end of
  // scope, and thus stays visible on the canvas.
  h_EOverP_0->DrawCopy();

  c->cd(2);
  h_EOverP_2->DrawCopy();

  c->cd(3);
  h_event_type->DrawCopy();

  c->cd(4);
  h_coin_time->DrawCopy();
  h_coin_time2->SetLineColor(4);
  h_coin_time2->DrawCopy("same");
  h_coin_time3->SetLineColor(2);
  h_coin_time3->DrawCopy("same");

  gPad->BuildLegend();

  gSystem->mkdir("results/df_example", true);
  c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) +
             ".pdf")
                .c_str());
  c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) +
             ".png")
                .c_str());

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  hs = new THStack("SHMS_cal", "SHMS calorimeter; E/P");

  h_EOverP_0->SetLineColor(1);
  h_EOverP_2->SetLineColor(4);
  h_EOverP_3->SetLineColor(2);
  h_EprOverP_0->SetLineColor(8);
  h_EprOverP_2->SetLineColor(6);
  h_EprOverP_3->SetLineColor(9);

  hs->Add((TH1 *)h_EOverP_0->Clone());
  hs->Add((TH1 *)h_EOverP_2->Clone());
  hs->Add((TH1 *)h_EOverP_3->Clone());

  hs->Add((TH1 *)h_EprOverP_0->Clone());
  hs->Add((TH1 *)h_EprOverP_2->Clone());
  hs->Add((TH1 *)h_EprOverP_3->Clone());

  hs->Draw("nostack");
  gPad->SetLogy(true);
  gPad->BuildLegend();

  c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
             ".pdf")
                .c_str());
  c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
             ".png")
                .c_str());
}
