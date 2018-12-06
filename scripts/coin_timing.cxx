#include <iostream>
#include <cmath>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
//#include "ROOT/TCanvas.hxx"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libGenVector.so)

#include "THStack.h"
#include "TBufferJSON.h"

#include "nlohmann/json.hpp"
#include "THcParmList.h"
R__LOAD_LIBRARY(libHallA.so)
R__LOAD_LIBRARY(libdc.so)
R__LOAD_LIBRARY(libHallC.so)

// fmt - string formatting library
#include "fmt/core.h"
#include "fmt/ostream.h"
R__LOAD_LIBRARY(libfmt.so)

using Pvec3D   = ROOT::Math::XYZVector;
using Pvec4D   = ROOT::Math::PxPyPzMVector;
using inters = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;
using FVec     = std::vector<float>;

void coin_timing(int RunNumber = 5893, const char* codatype = "COIN", int nevents = 100000) {

  std::string coda_type = codatype;
  //std::string coda_type = "SHMS";
  //std::string coda_type = "COIN";

  std::string hallc_replay_dir = "./";

  std::string rootfile    = std::string("ROOTfiles/coin_replay_production_");
  if(coda_type == "SHMS") {
    rootfile         = std::string("ROOTfiles/shms_replay_production_all_");// + std::to_string(RunNumber) + "_50000.root";
  }

  rootfile  += std::to_string(RunNumber) + "_" + std::to_string(nevents) +".root";

  std::string run_list_json  = "DBASE/run_list.json";

  // //std::string(std::getenv("DB_DIR"))
  //if (setenv("DB_DIR", hallc_replay_dir.c_str(), 1)) {
  //  std::cout << "Failed to set env var DB_DIR\n";
  //  std::exit(EXIT_FAILURE);
  //}

  std::string db_filename = hallc_replay_dir+"DBASE/"+coda_type+"/standard.database";

  //THcParmList* hc_parms = new THcParmList();
  //hc_parms->Define("gen_run_number", "Run Number", RunNumber);
  //hc_parms->AddString("g_ctp_database_filename", db_filename.c_str() );
  //hc_parms->Load(hc_parms->GetString("g_ctp_database_filename"), RunNumber);
  //hc_parms->Load(hc_parms->GetString("g_ctp_parm_filename"));
  //hc_parms->Load(hc_parms->GetString("g_ctp_kinematics_filename"), RunNumber);

  ////auto j = nlohmann::json::parse(hc_parms->PrintJSON(RunNumber));
  ////auto htheta_lab = hc_parms->Find("htheta_lab");
  ////if(htheta_lab){
  ////  htheta_lab->Print();
  ////}
  //auto ptheta_lab = hc_parms->Find("ptheta_lab");
  //ptheta_lab->Print();
  ////auto hpcentral = hc_parms->Find("hpcentral");
  ////hpcentral->Print();
  //auto ppcentral = hc_parms->Find("ppcentral");
  //ppcentral->Print();
  ////std::cout << j.dump()  << "\n";

  //double P0_SHMS = hcana::json::FindVarValueOr(hc_parms,"ppcentral",0.0) ;
  //std::cout <<  "SHMS P0 = " << P0_SHMS << "\n";

  //double THETA_SHMS = hcana::json::FindVarValueOr(hc_parms,"ptheta_lab",0.0) ;

  //nlohmann::json j2;
  //{
  //  ifstream  in_run_file(run_list_json);
  //  in_run_file >> j2;
  //}
  //if( j2.find(std::to_string(RunNumber)) == j2.end() ) {
  //}
  //j2[std::to_string(RunNumber)] = j[std::to_string(RunNumber)];

  //{
  //  //std::cout << j.dump(2) << "\n";
  //  // write prettified JSON to another file
  //  std::ofstream o(run_list_json);
  //  o << std::setw(4) << j2 << std::endl;
  //}

  //std::cout << j2.dump() << "\n";

  ROOT::EnableImplicitMT(30);

  Pvec4D  Pbeam(0,0,10.6,0.000511);

  ROOT::RDataFrame d("T",rootfile);

  auto d0 = d
  .Define("shms_e_EoverP",
          [](doublers& EOverP) {
            return EOverP[EOverP > 0.8 && EOverP < 1.8];
          }, {"P.cal.etottracknorm"})
  .Define("shms_e_EoverP_nGood",
          [](doublers& EOverP) {
            return (int)EOverP.size();
          }, {"shms_e_EoverP"})
  .Define("t1_t4_ROC1", 
          [](doublers& pos, doublers& neg) {
            doublers res;
            for(const auto& p :pos) {
              for(const auto& n : neg) {
                res.push_back(p-n);
              }
            }
            return res;
          }, {"T.coin.hTRIG1_ROC1_tdcTime","T.coin.hTRIG4_ROC1_tdcTime"})
  .Define("t1_t4_ROC2", 
          [](doublers& pos, doublers& neg) {
            doublers res;
            for(const auto& p :pos) {
              for(const auto& n : neg) {
                res.push_back(p-n);
              }
            }
            return res;
          }, {"T.coin.hTRIG1_ROC2_tdcTime","T.coin.hTRIG4_ROC2_tdcTime"})
  .Define("H_hod_1x_time_diff",
          [](doublers& pos, doublers& neg) {
            doublers res;
            for(const auto& p :pos) {
              for(const auto& n : neg) {
                res.push_back(p-n);
              }
            }
            return res;
          },{"H.hod.1x.posTdcTime","H.hod.1x.negTdcTime"})
  .Define("H_hod_1y_time_diff",[](doublers& pos, doublers& neg) {
    doublers res;
    for(const auto& p :pos) {
      for(const auto& n : neg) {
        res.push_back(p-n);
      }
    }
    return res;
  }, {"H.hod.1y.posTdcTime","H.hod.1y.negTdcTime"})
  .Define("H_hod_2x_time_diff",[](doublers& pos, doublers& neg) {
    doublers res;
    for(const auto& p :pos) {
      for(const auto& n : neg) {
        res.push_back(p-n);
      }
    }
    return res;
  }, {"H.hod.2x.posTdcTime","H.hod.2x.negTdcTime"})
  .Define("H_hod_2y_time_diff",[](doublers& pos, doublers& neg) {
    doublers res;
    for(const auto& p :pos) {
      for(const auto& n : neg) {
        res.push_back(p-n);
      }
    }
    return res;
  },
  {"H.hod.2y.posTdcTime","H.hod.2y.negTdcTime"})
  .Define("test_pulse",[](doublers& pos, doublers& neg) {
    std::vector<double> res;
    int i = 0; 
    for(const auto& p :pos) {
      if(i < neg.size()) {
        res.push_back(neg.at(i));
      }
      i++;
    }
    return res;
  }, {"H.hod.1x.posTdcTime", "H.hod.1x.posAdcPulseAmp"})
  .Define("test_time",[](doublers& pos, doublers& neg) {
    std::vector<double> res;
    int i = 0; 
    for(const auto& n :neg) {
      if(i < pos.size()) {
        res.push_back(pos.at(i));
      }
      i++;
    }
    return res;
  }, {"H.hod.1x.posTdcTime", "H.hod.1x.posAdcPulseAmp"})
  //.Define("test_adc1",[](doublers& pos) {
  //  doublers res;
  //  for(const auto& n :neg) {
  //    if(i < pos.size()) {
  //      res.push_back(n);
  //    }
  //    i++;
  //  }
  //  return res;
  //}, {"H.hod.1x.posTdcTime", "H.hod.1x.posAdcPulseAmp"})

  //.Define("H_hod_2y_time_t1t4_diff",
  //        [](double E, doublers& neg) {
  //  doublers res;
  //  for(const auto& n : neg) {
  //    res.push_back(E-n);
  //  }
  //  return res;
  //}, {"H_hod_2y_time_diff","t1_t4_ROC1"})
  ;
    
  auto d2 = d0.Filter(
    [](int n_good) {
        if( n_good >0 ) {
          return true;
        }
      return  false;
    },{"shms_e_EoverP_nGood"});


  gStyle->SetOptStat(1111);

  auto h_Trig1_roc1 =
      d2.Histo1D({"h_Trig1_roc1", "h_Trig1_roc1", 200, -600, 800}, "T.coin.hTRIG1_ROC1_tdcTime");
  auto h_Trig4_roc1 =
      d2.Histo1D({"h_Trig4_roc1", "h_Trig4_roc1", 200, -600, 800}, "T.coin.hTRIG4_ROC1_tdcTime");

  auto h_Trig1_m_Trig4_roc1 = d2.Histo1D(
      {"h_Trig1_m_Trig4_roc1", "ROC1;ROC1: trig1 - trig4", 200, -600, 800}, "t1_t4_ROC1");
  auto h_Trig1_roc2 =
      d2.Histo1D({"h_Trig1_roc2", "h_Trig1_roc2", 200, -600, 800}, "T.coin.hTRIG1_ROC2_tdcTime");
  auto h_Trig4_roc2 =
      d2.Histo1D({"h_Trig4_roc1", "h_Trig4_roc1", 200, -600, 800}, "T.coin.hTRIG4_ROC2_tdcTime");

  auto h_Trig1_m_Trig4_roc2 = d2.Histo1D(
      {"h_Trig1_m_Trig4_roc2", "ROC2;ROC2: trig1 - trig4", 200, -600, 800}, "t1_t4_ROC2");

  auto h_Hhod1x_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod1x_neg_tdc",";",    200, -3600, 2800},"H.hod.1x.negTdcTime");
  auto h_Hhod1y_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod1y_neg_tdc",";",    200, -3600, 2800},"H.hod.1y.negTdcTime");
  auto h_Hhod2x_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod2x_neg_tdc",";",    200, -3600, 2800},"H.hod.2x.negTdcTime");
  auto h_Hhod2y_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod2y_neg_tdc",";",    200, -3600, 2800},"H.hod.2y.negTdcTime");

  auto h_Hhod1x_pos_tdc  = d2.Histo1D<doublers>({"h_Hhod1x_pos_tdc",";",    200, -3600, 2800},"H.hod.1x.posTdcTime");
  auto h_Hhod1y_pos_tdc  = d2.Histo1D<doublers>({"h_Hhod1y_pos_tdc",";",    200, -3600, 2800},"H.hod.1y.posTdcTime");
  auto h_Hhod2x_pos_tdc  = d2.Histo1D<doublers>({"h_Hhod2x_pos_tdc",";",    200, -3600, 2800},"H.hod.2x.posTdcTime");
  auto h_Hhod2y_pos_tdc  = d2.Histo1D<doublers>({"h_Hhod2y_pos_tdc",";",    200, -3600, 2800},"H.hod.2y.posTdcTime");

  auto h_Hhod1x_tdc_diff  = d2.Histo1D<doublers>({"h_Hhod1x_tdc_diff",";",    200, -3600, 2800},"H_hod_1x_time_diff");
  auto h_Hhod1y_tdc_diff  = d2.Histo1D<doublers>({"h_Hhod1y_tdc_diff",";",    200, -3600, 2800},"H_hod_1y_time_diff");
  auto h_Hhod2x_tdc_diff  = d2.Histo1D<doublers>({"h_Hhod2x_tdc_diff",";",    200, -3600, 2800},"H_hod_2x_time_diff");
  auto h_Hhod2y_tdc_diff  = d2.Histo1D<doublers>({"h_Hhod2y_tdc_diff",";",    200, -3600, 2800},"H_hod_2y_time_diff");

  auto h_Hhod2y_adc_tdc_diff =
      d2.Histo2D({"h_Hhod2y_tdc_diff", ";x;y", 100, -3600, 3600, 100, 0, 2000},
                                     "test_time", "test_pulse");

  //auto h_Hhod1x_tdc_diff  = d0.Histo1D<doublers>({"h_Hhod1x_tdc_diff",";",    200, -3600, 2800},"H.hod.1y.negTdcTime");
  //auto h_Hhod2_negtdc  = d0.Histo1D<doublers>({"h_Hhod2x_negtdc",";",    200, -3600, 2800},"H.hod.2x.negTdcTime");

  //auto h_EOverP_1       = d0.Histo1D<doublers>({"hms_e_EoverP_1","E/P cut; HMS E/P",100,0.05,1.8},"shms_e_EoverP");
  //auto h_EOverP_nGood_0 = d0.Histo1D({"h_EOverP_nGood_0","h_EOverP_nGood_0",10,0,10},"P.cal.ntracks");
  //auto h_EOverP_nGood_1 = d0.Histo1D({"h_EOverP_nGood_1","h_EOverP_nGood_0",10,0,10},"shms_e_EoverP_nGood");

  //auto s_HGC_xy_pos      = d3.Fill<double,double>(TH2D("s_HGC_xy_pos ","xy fp; x;y",120,-60,60,120,-60,60), {"hgc_xpos","hgc_ypos"});

  //auto h_Hhod1_tdc_diff = d0.Fill(
  //    TH2D("h_Hhod1_tdc_diff ", "xy fp; x;y", 120, -6000, 6000, 120, -6000, 6000),
  //    {"t1_t4_ROC1", "H_hod_1x_time_diff"});

  // -----------------------------------------------------------
  //
  TCanvas* c = nullptr;
  int b1 = 0;
  int b2 = 0;
  double hmax = 0.0;
  TLatex latex;
  TH1* h1 = nullptr;
  TH2* h2 = nullptr;

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  auto hs_trigs = new THStack("hs_trigs","trig1 - trig4");

  h1 = (TH1*)h_Trig1_roc1->Clone();
  hs_trigs->Add(h1);
  h1 = (TH1*)h_Trig4_roc1->Clone();
  h1->SetLineColor(4);
  hs_trigs->Add(h1);
  h1 = (TH1*)h_Trig1_roc2->Clone();
  h1->SetLineColor(2);
  hs_trigs->Add(h1);
  h1 = (TH1*)h_Trig4_roc2->Clone();
  h1->SetLineColor(8);
  hs_trigs->Add(h1);

  gPad->SetLogy(true);
  hs_trigs->Draw("nostack");


  c->cd(2);
  gPad->SetLogy(true);
  auto hs_Hhod1 = new THStack("hs_Hhod1","H hod; tdcTime");

  h1 = (TH1*)h_Hhod1x_neg_tdc->Clone();
  h1->SetLineColor(1);
  hs_Hhod1->Add(h1);

  h1 = (TH1*)h_Hhod1y_neg_tdc->Clone();
  h1->SetLineColor(2);
  hs_Hhod1->Add(h1);

  h1 = (TH1*)h_Hhod2x_neg_tdc->Clone();
  h1->SetLineColor(4);
  hs_Hhod1->Add(h1);

  h1 = (TH1*)h_Hhod2y_neg_tdc->Clone();
  h1->SetLineColor(8);
  hs_Hhod1->Add(h1);

  hs_Hhod1->Draw("nostack");

  c->cd(3);
  h_Hhod2y_adc_tdc_diff->DrawCopy("colz");

  c->cd(4);
  gPad->SetLogy(true);
  auto hs_Hhod_diff = new THStack("hs_Hhod_diff","H hod; tdcTime");

  h1 = (TH1*)h_Hhod1x_tdc_diff->Clone();
  h1->SetLineColor(1);
  hs_Hhod_diff->Add(h1);

  h1 = (TH1*)h_Hhod1y_tdc_diff->Clone();
  h1->SetLineColor(2);
  hs_Hhod_diff->Add(h1);

  h1 = (TH1*)h_Hhod2x_tdc_diff->Clone();
  h1->SetLineColor(4);
  hs_Hhod_diff->Add(h1);

  h1 = (TH1*)h_Hhod2y_tdc_diff->Clone();
  h1->SetLineColor(8);
  hs_Hhod_diff->Add(h1);

  hs_Hhod_diff->Draw("nostack");

}
