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

void rf_timing(int RunNumber = 6195, const char* codatype = "COIN", int nevents = -1) {

  std::string coda_type        = codatype;
  std::string hallc_replay_dir = "./";
  std::string rootfile         = std::string("ROOTfiles_online/coin_replay_production_");
  if (coda_type == "SHMS") {
    rootfile = std::string("ROOTfiles_online/shms_replay_production_all_");
  }
  rootfile += std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";
  std::string run_list_json = "DBASE/run_list.json";

  std::string db_filename = hallc_replay_dir + "DBASE/" + coda_type + "/standard.database";

  ROOT::EnableImplicitMT(30);

  Pvec4D Pbeam(0, 0, 10.6, 0.000511);

  ROOT::RDataFrame d("T", rootfile);

  std::string hpdelta = "P.gtr.dp > -10 && P.gtr.dp < 20 && "
                        "H.gtr.dp > -10 && H.gtr.dp < 10";
  std::string epiCut = "P.aero.npeSum > 1.0 && P.cal.eprtracknorm < 0.2 && "
                       "H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 && "
                       "H.cal.etottracknorm < 2.0 && P.cal.etottracknorm<1.0";
  // && H.cal.eprtracknorm  > 0.2

  //auto d0 = d
  auto h_coin_time = d.Histo1D(
      {"coin_time", "coin_time", 1000, -10, 100}, "CTime.ePiCoinTime_ROC2");
  auto h_rf_time = d.Histo1D(
      {"h_rf_time", "rf_time; rf_time", 1000, -100, 100}, "T.coin.pRF_tdcTime");

  auto d1 = d.Define("rf_minus_fp_time","T.coin.pRF_tdcTime - P.hod.starttime");
  auto h_time_diff = d1.Histo1D( {"h_rf_time", "rf_time; rf_time", 1000, -10, 10}, "rf_minus_fp_time");
  auto h_time_diff2 = d1.Filter("P.hgcer.npeSum>1.0").Histo1D(
      {"h_rf_time", "rf_time; rf_time", 1000, -10, 10}, "rf_minus_fp_time");


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
  h_coin_time->DrawCopy();

  c->cd(2);
  h_rf_time->DrawCopy();

  c->cd(3);
  h_time_diff->DrawCopy();
h_time_diff2->SetLineColor(2);
  h_time_diff2->DrawCopy("same");



}
