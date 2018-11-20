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

void coin_timing3(int RunNumber = 6335, int nevents = -1) {
  std::string coda_type = "COIN";

  std::string rootfile    = std::string("ROOTfiles_csv/coin_replay_production_");
  rootfile  += std::to_string(RunNumber) + "_" + std::to_string(nevents) +".root";

  std::string db_filename = "DBASE/"+coda_type+"/standard.database";
  std::cout << rootfile << "\n";

  TFile* f = new TFile(rootfile.c_str()); 
  new TBrowser();

  ROOT::EnableImplicitMT(4);
  ROOT::RDataFrame d("T",rootfile);

  // Beam energy
  Pvec4D  Pbeam(0,0,10.6,0.000511);

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
          }, {"T.coin.hTRIG1_ROC2_tdcTime","T.coin.hTRIG4_ROC2_tdcTime"});

  auto d2 = d0.Filter("(fEvtHdr.fEvtType == 6) || (fEvtHdr.fEvtType == 4)");
  //.Filter(
  //  [](int n_good) {
  //      if( n_good >0 ) {
  //        return true;
  //      }
  //    return  false;
  //  },{"shms_e_EoverP_nGood"});


  gStyle->SetOptStat(1111);

  auto dderp3 =
      d2.Define("time_fix_maybe", "T.coin.pTRIG1_ROC2_tdcTime  - T.coin.p1X_tdcTimeRaw")
          .Define("time_fix_maybe2", "T.coin.p1X_tdcTimeRaw - T.coin.h1X_tdcTimeRaw")

          .Define("delta_s1xs_1", "(T.coin.p1X_tdcTimeRaw - T.coin.pTRIG1_ROC2_tdcTimeRaw) - "
                                  "(T.coin.h1X_tdcTimeRaw-T.coin.pTRIG4_ROC1_tdcTimeRaw )")
          .Define("delta_s1xs_2", "(T.coin.p1X_tdcTimeRaw - T.coin.pTRIG4_ROC2_tdcTimeRaw) - "
                                  "(T.coin.h1X_tdcTimeRaw-T.coin.pTRIG1_ROC1_tdcTimeRaw )")
          .Define("delta_s1xs_3", "T.coin.p1X_tdcTimeRaw - T.coin.h1X_tdcTimeRaw");

  auto h_test2 = dderp3.Histo1D({"h_test2", "(p1x-T1_R2) - (h1x-T4_R1) RAW", 800, -400, 400}, "delta_s1xs_1" );
  auto h_test3 = dderp3.Histo1D({"h_test3", "(p1x-T4_R2) - (h1x-T1_R1) RAW", 800, -400, 400 }, "delta_s1xs_2" );
  auto h_test4 = dderp3.Histo1D({"h_test3", "time_fix_maybe", 800, -400, 400 }, "delta_s1xs_3" );

  auto h_test  = dderp3.Histo2D({"h_test", ";CTime.ePiCoinTime_ROC2;(p1x-T1_R2) - (h1x-T4_R1) RAW", 100, 35, 55, 100, -20, 60}, "CTime.ePiCoinTime_ROC2","delta_s1xs_1" );

  auto h_coin_time1 = d2.Histo1D( {"coin_time", "coin_time", 1000, -20,120}, "CTime.ePiCoinTime_ROC1");
  auto h_coin_time2 = d2.Histo1D( {"coin_time", "coin_time", 1000, -20,120}, "CTime.ePiCoinTime_ROC2");

  auto h_coin_time_test = dderp3
  .Define("coin_time_2", "CTime.ePiCoinTime_ROC1 - T.coin.h1X_tdcTime")
  .Histo1D( {"h_coin_time_test", "coin_time_test", 2000, -400, 400}, "coin_time_2");

  auto h_test_fixed = d2
      .Define("fix_try1","T.coin.hTRIG1_ROC1_tdcTime ")
      .Histo1D({"h_Trig1_roc1", "h_Trig1_roc1", 200, 0, 800}, "fix_try1");

  auto h_Trig1_roc1 =
      d2.Histo1D({"h_Trig1_roc1", "h_Trig1_roc1", 200, 0, 800}, "T.coin.hTRIG1_ROC1_tdcTime");
  auto h_Trig4_roc1 =
      d2.Histo1D({"h_Trig4_roc1", "h_Trig4_roc1", 200, 0, 800}, "T.coin.hTRIG4_ROC1_tdcTime");

  auto h_Trig1_m_Trig4_roc1 = d2.Histo1D(
      {"h_Trig1_m_Trig4_roc1", "ROC1;ROC1: trig1 - trig4", 200, 0, 800}, "t1_t4_ROC1");
  auto h_Trig1_roc2 =
      d2.Histo1D({"h_Trig1_roc2", "h_Trig1_roc2", 200, 0, 800}, "T.coin.hTRIG1_ROC2_tdcTime");
  auto h_Trig4_roc2 =
      d2.Histo1D({"h_Trig4_roc1", "h_Trig4_roc1", 200, 0, 800}, "T.coin.hTRIG4_ROC2_tdcTime");

  auto h_Trig1_m_Trig4_roc2 = d2.Histo1D(
      {"h_Trig1_m_Trig4_roc2", "ROC2;ROC2: trig1 - trig4", 200, 0, 800}, "t1_t4_ROC2");

  auto P_hod_goodstarttime = d2.Histo1D(
      {"P_hod_goodstarttime", "P_hod_goodstarttime;", 200, -5, 5}, "P.hod.goodstarttime");
  auto H_hod_goodstarttime = d2.Histo1D(
      {"H_hod_goodstarttime", "H_hod_goodstarttime;", 200, -5, 5}, "H.hod.goodstarttime");

  auto h_Hhod1x_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod1x_neg_tdc",";",    200, -200, 200},"P.hod.1x.GoodNegTdcTimeUnCorr");
  auto h_Hhod1y_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod1y_neg_tdc",";",    200, -200, 200},"P.hod.1y.GoodNegTdcTimeUnCorr");
  auto h_Hhod2x_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod2x_neg_tdc",";",    200, -200, 200},"P.hod.2x.GoodNegTdcTimeUnCorr");
  auto h_Hhod2y_neg_tdc  = d2.Histo1D<doublers>({"h_Hhod2y_neg_tdc",";",    200, -200, 200},"P.hod.2y.GoodNegTdcTimeUnCorr");

  auto d3 =  d2.Filter("Ndata.P.hod.2x.GoodNegTdcTimeUnCorr > 0 ")
  .Define("hod1x1",[](const doublers& t) {return t.at(0);},{"P.hod.2x.GoodNegTdcTimeUnCorr"})
  //.Define("coin_t1r1_1",[](const doublers& t) {return t.at(0);},{"T.coin.hTRIG1_ROC2_tdcTime"})
  ;
  auto h_Hhod2y_adc_tdc_diff = d3.Histo2D<doublers>({"h_Hhod2y_tdc_diff", ";x;y", 100, 200, 300, 100, -200, 200},
                                     "T.coin.hTRIG1_ROC2_tdcTime", "hod1x1");

  //TH2F* hS1x = new TH2F("hS1x","hS1x", 50, 270, 300, 50, 20, 80);
  //TH2F* hS2x = new TH2F("hS2x","hS2x", 50, 270, 300, 50, 20, 80);
  //TH2F* hS1y = new TH2F("hS1y","hS1y", 50, 270, 300, 50, 20, 80);
  //TH2F* hS2y = new TH2F("hS2y","hS2y", 50, 270, 300, 50, 20, 80);
  //TH2F* hS1xp = new TH2F("hS1xp","hS1xp", 50, 270, 300, 50, 20, 80);
  //TH2F* hS2xp = new TH2F("hS2xp","hS2xp", 50, 270, 300, 50, 20, 80);
  //TH2F* hS1yp = new TH2F("hS1yp","hS1yp", 50, 270, 300, 50, 20, 80);
  //TH2F* hS2yp = new TH2F("hS2yp","hS2yp", 50, 270, 300, 50, 20, 80);

  //TH2F* pS1xp = new TH2F("pS1xp","pS1xp", 50, 270, 300, 50, 20, 80);
  //TH2F* pS2xp = new TH2F("pS2xp","pS2xp", 50, 270, 300, 50, 20, 80);
  //TH2F* pS1yp = new TH2F("pS1yp","pS1yp", 50, 270, 300, 50, 20, 80);
  //TH2F* pS2yp = new TH2F("pS2yp","pS2yp", 50, 270, 300, 50, 20, 80);

  //TH2F* pS1xn = new TH2F("pS1xn","pS1xn", 50, 270, 300, 50, 20, 80);
  //TH2F* pS2xn = new TH2F("pS2xn","pS2xn", 50, 270, 300, 50, 20, 80);
  //TH2F* pS1yn = new TH2F("pS1yn","pS1yn", 50, 270, 300, 50, 20, 80);
  //TH2F* pS2yn = new TH2F("pS2yn","pS2yn", 50, 270, 300, 50, 20, 80);


  //d2.Foreach([=](double t,const doublers& hx1,const doublers& hy1,const doublers& hx2 ,const doublers& hy2,
  //               const doublers& px1,const doublers& py1,const doublers& px2 ,const doublers& py2) {
  //  for(const auto& h0 : hx1) {
  //    hS1x->Fill(t, h0);
  //  }
  //  for(const auto& h0 : hy1) {
  //    hS1y->Fill(t, h0);
  //  }
  //  for(const auto& h0 : hx2) {
  //    hS2x->Fill(t, h0);
  //  }
  //  for(const auto& h0 : hy2) {
  //    hS2y->Fill(t, h0);
  //  }
  //  for(const auto& h0 : px1) {
  //    pS1xn->Fill(t, h0);
  //  }
  //  for(const auto& h0 : py1) {
  //    pS1yn->Fill(t, h0);
  //  }
  //  for(const auto& h0 : px2) {
  //    pS2xn->Fill(t, h0);
  //  }
  //  for(const auto& h0 : py2) {
  //    pS2yn->Fill(t, h0);
  //  }
  //},{"T.coin.hTRIG1_ROC2_tdcTime",
  //  "P.hod.1x.GoodNegTdcTimeUnCorr",
  //  "P.hod.1y.GoodNegTdcTimeUnCorr",
  //  "P.hod.2x.GoodNegTdcTimeUnCorr",
  //  "P.hod.2y.GoodNegTdcTimeUnCorr",
  //  "H.hod.1x.GoodNegTdcTimeUnCorr",
  //  "H.hod.1y.GoodNegTdcTimeUnCorr",
  //  "H.hod.2x.GoodNegTdcTimeUnCorr",
  //  "H.hod.2y.GoodNegTdcTimeUnCorr"
  //});

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


  //c = new TCanvas();
  //c->Divide(4,1);
  //c->cd(1);
  //hS1x->Draw("colz");
  //c->cd(2);
  //hS1y->Draw("colz");
  //c->cd(3);
  //hS2x->Draw("colz");
  //c->cd(4);
  //hS2y->Draw("colz");

  //c = new TCanvas();
  //c->Divide(4,1);
  //c->cd(1);
  //pS1xn->Draw("colz");
  //c->cd(2);
  //pS1yn->Draw("colz");
  //c->cd(3);
  //pS2xn->Draw("colz");
  //c->cd(4);
  //pS2yn->Draw("colz");

  // ---------------------------------------------------------
  //
  c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  //auto hs_trigs = new THStack("hs_trigs","trig1 - trig4");

  //h1 = (TH1*)h_Trig1_roc1->Clone();
  //hs_trigs->Add(h1);

  //h1 = (TH1*)h_Trig4_roc1->Clone();
  //h1->SetLineColor(4);
  //hs_trigs->Add(h1);

  //h1 = (TH1*)h_Trig1_roc2->Clone();
  //h1->SetLineColor(2);
  //hs_trigs->Add(h1);

  //h1 = (TH1*)h_Trig4_roc2->Clone();
  //h1->SetLineColor(8);
  //hs_trigs->Add(h1);

  //gPad->SetLogy(true);
  //hs_trigs->Draw("nostack");
  h_test_fixed->DrawCopy();

  c->cd(2);
  h_test->DrawCopy("colz");
  //h23->Draw("colz");
  //h_Hhod2y_adc_tdc_diff->DrawCopy("colz");
  //H_hod_goodstarttime->SetLineColor(2);
  //H_hod_goodstarttime->DrawCopy("same");
  //auto hs_Hhod1 = new THStack("hs_Hhod1","H hod; tdcTime");

  //h1 = (TH1*)h_Hhod1x_neg_tdc->Clone();
  //h1->SetLineColor(1);
  //hs_Hhod1->Add(h1);

  //h1 = (TH1*)h_Hhod1y_neg_tdc->Clone();
  //h1->SetLineColor(2);
  //hs_Hhod1->Add(h1);

  //h1 = (TH1*)h_Hhod2x_neg_tdc->Clone();
  //h1->SetLineColor(4);
  //hs_Hhod1->Add(h1);

  //h1 = (TH1*)h_Hhod2y_neg_tdc->Clone();
  //h1->SetLineColor(8);
  //hs_Hhod1->Add(h1);

  //hs_Hhod1->Draw("nostack");

  c->cd(3);
  auto hs_trigs_test2 = new THStack("hs_trigs_test","trigstest");

  hs_trigs_test2->Add((TH1*)h_test2->Clone());

  h_test3->SetLineColor(2);
  hs_trigs_test2->Add((TH1*)h_test3->Clone());

  //h_test4->SetLineColor(4);
  //hs_trigs_test2->Add((TH1*)h_test4->Clone());

  hs_trigs_test2->Draw("nostack");
  gPad->BuildLegend(0.7,0.7,0.9,0.9);

  c->cd(4);
  auto hs_trigs_test = new THStack("hs_trigs_test","trigstest");

  hs_trigs_test->Add((TH1*)h_coin_time1->Clone());

  h_coin_time2->SetLineColor(4);
  hs_trigs_test->Add((TH1*)h_coin_time2->Clone());

  //h_coin_time_test->SetLineColor(2);
  //hs_trigs_test->Add((TH1*)h_coin_time_test->Clone());


  hs_trigs_test->Draw("nostack");


}
