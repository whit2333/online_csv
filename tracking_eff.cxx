#include "nlohmann/json.hpp"
#include <cmath>
#include <iostream>

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

//#include "THcParmList.h"
// R__LOAD_LIBRARY(libPodd.so)
// R__LOAD_LIBRARY(libHallA.so)
// R__LOAD_LIBRARY(libdc.so)
// R__LOAD_LIBRARY(libHallC.so)

// fmt - string formatting library
//#include "fmt/core.h"
//#include "fmt/ostream.h"
// R__LOAD_LIBRARY(libfmt.so)

using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;

// VecOps::RVec is like std::vector with some extra bells and whistles
using inters = ROOT::VecOps::RVec<int>;
using doublers = ROOT::VecOps::RVec<double>;
using floaters = ROOT::VecOps::RVec<float>;
using shorters = ROOT::VecOps::RVec<short>;
using nlohmann::json;

auto get_trk_eff_df = [] (auto d, int n1, int n2) {
  ROOT::EnableThreadSafety();
  auto ret =  d
      //.Filter(std::string("fEvtHdr.fEvtNum >= ") + std::to_string(n1))
      //.Filter(std::string("fEvtHdr.fEvtNum <  ") + std::to_string(n2))
      .Define("shmsDC1Planes_large", "(P.dc.1x1.nhit + P.dc.1u2.nhit + P.dc.1u1.nhit + "
                                     "P.dc.1v1.nhit + P.dc.1x2.nhit + P.dc.1v2.nhit) > 20")
      .Define("shmsDC2Planes_large", "(P.dc.2x1.nhit + P.dc.2u2.nhit + P.dc.2u1.nhit + "
                                     "P.dc.2v1.nhit + P.dc.2x2.nhit + P.dc.2v2.nhit) > 20")
      .Define("shmsDCany_large", "shmsDC1Planes_large || shmsDC2Planes_large")
      .Define("pcut_cer_ng_elec", "P.ngcer.npeSum > 0.5")
      .Define("pcut_cer_ng_pi", "P.ngcer.npeSum <= 0.5")
      .Define("pcut_cer_hg_elec", "P.hgcer.npeSum > 0.5")
      .Define("pcut_cer_hg_pi", "P.hgcer.npeSum <= 0.5")
      .Define("pcut_cer_pi_both", "pcut_cer_ng_pi && pcut_cer_hg_pi")
      .Define("pcut_cal_elec", "P.cal.etracknorm > 0.6 && P.cal.etracknorm < 1.6")
      .Define("pcut_cal_pi", "P.cal.etracknorm <= 0.6 && P.cal.etracknorm > 0.")
      .Define("pcut_elec_all", "pcut_cer_ng_elec && pcut_cer_hg_elec && pcut_cal_elec")
      .Define("pcut_pi_all", "pcut_cer_ng_pi && pcut_cer_hg_pi && pcut_cal_pi")
      .Define("shmsScinGood", "P.hod.goodscinhit == 1")
      .Define("shmsGoodBetanotrk", "(P.hod.betanotrack > 0.5) && (P.hod.betanotrack < 1.4)")
      .Define("shmsScinShould", "shmsScinGood && shmsGoodBetanotrk && !shmsDCany_large")
      // shmsScinShoulde     shmsScinShould &&  P.cal.etotnorm > 0.6 && P.cal.etotnorm < 1.6 &&
      // P.hgcer.npeSum > 0.5 #shmsScinShouldh     shmsScinShould && P.cal.etotnorm <= 0.6 &&
      //P.cal.etotnorm > 0. && P.hgcer.npeSum < 0.2 shmsScinShouldh     shmsScinShould &&
      // P.cal.etotnorm <= 0.6 && P.cal.etotnorm > 0.
      .Define("shmsScinDid", "shmsScinShould && P.dc.ntrack > 0")
      // shmsScinDid         shmsScinShould && P.dc.ntrack > 0
      // shmsScinDide        shmsScinShoulde && P.dc.ntrack > 0
      ;
  return ret;
};

void tracking_eff(int RunNumber = 6244, int nevents = -1) {

  ROOT::EnableThreadSafety();
  ROOT::EnableImplicitMT();
  ROOT::Experimental::TTaskGroup tg;

  std::string coda_type = "COIN";
  std::string rootfile  = "ROOTfiles_csv/";
  rootfile += std::string("coin_replay_production_");
  rootfile += std::to_string(RunNumber) + "_" + std::to_string(nevents) + ".root";

  auto   file = new TFile(rootfile.c_str(), "READ");
  TTree* t    = nullptr;
  new TBrowser;


  Pvec4D Pbeam(0, 0, 10.598, 0.000511);


  TGraphErrors* gr = new TGraphErrors();
  TCanvas*      c  = new TCanvas("c", "x hist");


  int n_entries_0 = 0;

  for (int iplot = 0; iplot < 500000; iplot++) {

    file->ReadKeys();
    if (t) {
      delete t;
      t = nullptr;
    }
    file->GetObject("T", t);
    if (!t) {
      t = (TTree*)gROOT->FindObject("T");
    }
    if (!t) {
      std::cout << " error : no tree\n";
      std::quick_exit(-127);
    }

    int n_entries_1 = t->GetEntries();

    if ((n_entries_1 - n_entries_0) < 1000) {
      std::cout <<  n_entries_1 << " 1  entries \n";
      std::cout <<  n_entries_0 << " 0  entries \n";
      std::cout << "sleeping...\n";
      gSystem->Sleep(10);
      continue;
    }
    std::cout << "tree : " << t << "\n";

    ROOT::RDataFrame d(*t);
    // auto             d0       = d;
    std::cout << n_entries_0 << " to " << n_entries_1 << "\n";
    auto d0       = get_trk_eff_df(d, n_entries_0, n_entries_1);
    //auto d0       = get_trk_eff_df(d.Range(n_entries_0, n_entries_1), n_entries_0, n_entries_1);
    auto d_should = d0.Filter("shmsScinShould");
    // auto             d_did    = d0.Filter("shmsScinDid");

    auto h = d_should.Histo1D({"hshmsScinDid", "hshmsScinDid", 2, 0, 2}, "shmsScinDid");
    //if (iplot == 0) {
      h.OnPartialResult(1000, [c, gr](TH1D& h_) {
        c->cd();
        int    n1     = h_.GetBinContent(1);
        int    n2     = h_.GetBinContent(2);
        int    ipoint = gr->GetN();
        if(n1+n2 >100) {
        double val    = double(n2) / double(n1 + n2);
        std::cout << val << "\n";
        gr->SetPoint(ipoint, ipoint, val);
        gr->SetPointError(ipoint, 0, std::sqrt(double(n1)) / double(n1 + n2));
        gr->Draw("al");
        h_.Reset();
        c->Update();
        }
      });
    //}
    //std::cout << " asdf\n";

    c->cd();
    //h->Clone(); // event loop runs here, this `Draw` is executed after the event loop is
    int    n1     = h->GetBinContent(1);
    int    n2     = h->GetBinContent(2);
    int    ipoint = gr->GetN();
    double val    = double(n2) / double(n1 + n2);
    if (!std::isnan(val)) {
      std::cout << val << "\n";
      gr->SetPoint(ipoint, ipoint, val);
      gr->SetPointError(ipoint, 0, std::sqrt(double(n1)) / double(n1 + n2));
      h->Reset();
      n_entries_0 = n_entries_1;
    }
    // finished
    gr->Draw("al"); // event loop runs here, this `Draw` is executed after the event loop is finished
    c->Update();
    std::cout << " 2\n";
  }

  c->cd();
  gr->Draw("al"); // event loop runs here, this `Draw` is executed after the event loop is finished
  c->Update();

  //# Add cuts to stricter times and apply them to drift distances
  // pcut_cer_ng_elec    	  P.ngcer.npeSum > 0.5
  // pcut_cer_ng_pi      	  P.ngcer.npeSum <= 0.5
  // pcut_cer_hg_elec	          P.hgcer.npeSum > 0.5
  // pcut_cer_hg_pi	  	  P.hgcer.npeSum <= 0.5
  // pcut_cer_pi_both       	  pcut_cer_ng_pi && pcut_cer_hg_pi
  // pcut_cer_elec_both     	  pcut_cer_ng_elec && pcut_cer_hg_elec
  // pcut_cal_elec          	  P.cal.etracknorm > 0.6 && P.cal.etracknorm < 1.6
  // pcut_cal_pi            	  P.cal.etracknorm <= 0.6 && P.cal.etracknorm > 0.
  // pcut_elec_all          	  pcut_cer_ng_elec && pcut_cer_hg_elec && pcut_cal_elec
  // pcut_pi_all            	  pcut_cer_ng_pi && pcut_cer_hg_pi && pcut_cal_pi
  //
  // shmsScinGood        P.hod.goodscinhit == 1
  // shmsGoodBetanotrk   P.hod.betanotrack > 0.5 && P.hod.betanotrack < 1.4
  //
  // shmsScinShould      shmsScinGood && shmsGoodBetanotrk && !shmsDCany_large
  // shmsScinShoulde     shmsScinShould &&  P.cal.etotnorm > 0.6 && P.cal.etotnorm < 1.6 &&
  // P.hgcer.npeSum > 0.5 #shmsScinShouldh     shmsScinShould && P.cal.etotnorm <= 0.6 &&
  // P.cal.etotnorm > 0. && P.hgcer.npeSum < 0.2 shmsScinShouldh     shmsScinShould &&
  // P.cal.etotnorm <= 0.6 && P.cal.etotnorm > 0. shmsScinDid         shmsScinShould &&
  // P.dc.ntrack > 0 shmsScinDide        shmsScinShoulde && P.dc.ntrack > 0

  // HMS Scaler tree
  // ROOT::RDataFrame d_sh("TSH", rootfile);

  // int N_scaler_events = *(d_sh.Count());

  // int N_total = *c_n_events_total;

  // std::vector<std::array<int,2>> cut_ranges;

  // d_sh.Foreach([&](auto e,  auto y){
  //}, {"evNumber", "H.BCM4B.scalerCharge"})

  // std::string hpdelta = "H.gtr.dp > -10 && H.gtr.dp < 10";
  // std::string epiCut = "H.cer.npeSum > 1.0 && H.cal.etottracknorm > 0.6 && "
  //                     "H.cal.etottracknorm < 2.0 && P.cal.etottracknorm<1.0";
  //// && H.cal.eprtracknorm  > 0.2

  // auto c_n_events_hms = d.Filter("(fEvtHdr.fEvtType == 6) || (fEvtHdr.fEvtType == 2)").Count();

  // auto d0 = d.Define("hms_e_EoverP",
  //                   [](doublers& EOverP) { return EOverP[EOverP > 0.5 && EOverP < 1.8]; },
  //                   {"H.cal.etottracknorm"})
  //              .Define("hms_e_EoverP_nGood", [](doublers& EOverP) { return (int)EOverP.size();
  //              },
  //                      {"hms_e_EoverP"});

  // auto d0_coin_only = d.Filter("fEvtHdr.fEvtType == 6");
  // auto d0_hms_only  = d.Filter("fEvtHdr.fEvtType == 2");

  //// Apply the electron cuts
  // auto d_spec_cuts = d0.Filter(hpdelta).Filter(epiCut);
  // auto d_spec_cuts_coin_only = d_spec_cuts.Filter("fEvtHdr.fEvtType == 6");
  // auto d_spec_cuts_hms_only  = d_spec_cuts.Filter("fEvtHdr.fEvtType == 2");

  // auto h_EOverP_0 = d0.Histo1D<doublers>(
  //    {"hms_e_EoverP_0", "HMS total shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.etottracknorm");
  // auto h_EOverP_2 = d_spec_cuts.Histo1D<doublers>(
  //    {"hms_e_EoverP_2", "HMS total shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.etottracknorm");

  // auto h_EprOverP_0 = d0.Histo1D<doublers>(
  //    {"hms_e_EoverP_0", "HMS pre-shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.eprtracknorm");
  // auto h_EprOverP_2 = d_spec_cuts.Histo1D<doublers>(
  //    {"hms_e_EoverP_2", "HMS pre-shower; HMS E/P", 100, 0.05, 1.8}, "H.cal.eprtracknorm");

  // auto h_event_type   = d0.Histo1D({"EvtType", "EvtType", 10, 0, 10}, "fEvtHdr.fEvtType");
  // auto h_event_type_2 = d_spec_cuts.Histo1D({"EvtType2", "EvtType", 10, 0, 10},
  // "fEvtHdr.fEvtType");

  // auto h_beta_0 = d_spec_cuts.Histo1D({"h_beta_0", "h_beta_0", 150, 0, 1.5}, "H.gtr.beta");

  // auto bcm4b_charge        = d_sh.Max("H.BCM4B.scalerChargeCut");
  // auto el_real_scaler      = d_sh.Max("H.hEL_REAL.scaler");
  // auto time_1MHz           = d_sh.Max("H.1MHz.scalerTime");
  ////auto hTRIG1_ROC1_npassed = d_sh.Max("H.hTRIG1_ROC1.npassed");
  ////auto H_hTRIG1_scaler     = d_sh.Max("H.hTRIG1.scaler");

  ////{(hTRIG1_ROC1.npassed / H.hTRIG1.scaler)*100.0:%3.4f}
  ////H.hEL_REAL.scaler/H.1MHz.scalerTime)/1000
  // auto total_charge = bcm4b_charge;

  // auto c_e_yield_raw  = d.Count();
  // auto c_e_yield      = d_spec_cuts.Count();
  // auto c_T2_yield_raw = d0_hms_only.Count();
  // auto c_T6_yield_raw = d0_coin_only.Count();
  // auto c_T2_yield     = d_spec_cuts_hms_only.Count();
  // auto c_T6_yield     = d_spec_cuts_coin_only.Count();

  //// -------------------------------------
  //// End lazy eval
  //// -------------------------------------

  // double good_total_charge = *bcm4b_charge / 1000.0; // mC
  ////double hms_live_time = double(*hTRIG1_ROC1_npassed) / double(*H_hTRIG1_scaler);

  // double hms_scaler_yield     = ((*el_real_scaler) / good_total_charge);
  // double hms_scaler_yield_unc = (std::sqrt(*el_real_scaler) / good_total_charge);

  // double hms_e_yield      = (*c_e_yield_raw) / (*total_charge);
  // double hms_e_yield2     = (*c_e_yield) / (*total_charge);

  // double ps_cor_hms_e_yield  = hms_e_yield *singles_ps_value;
  // double ps_cor_hms_e_yield2 = hms_e_yield2*singles_ps_value;

  // std::cout << " good_total_charge " << good_total_charge << " \n";
  // std::cout << " hms_scaler_yield " << hms_scaler_yield << " \n";

  //// Update counts list
  // json jruns;
  //{
  //  std::ifstream input_file("db2/run_count_list.json");
  //  input_file >> jruns;
  //}
  // std::string run_str                = std::to_string(RunNumber);
  // jruns[run_str]["hms e raw counts"] = int(*c_T2_yield_raw + *c_T6_yield_raw);
  // jruns[run_str]["hms e counts"]     = int(*c_T2_yield + *c_T6_yield);
  // jruns[run_str]["ps cor. hms e raw counts"] = int((*c_T2_yield_raw)*singles_ps_value +
  // (*c_T6_yield_raw)); jruns[run_str]["ps cor. hms e counts"]     =
  // int((*c_T2_yield_raw)*singles_ps_value + (*c_T6_yield_raw)); jruns[run_str]["charge bcm4b 2u
  // cut"]     = good_total_charge; jruns[run_str]["hms e raw yield"] =
  // double((*c_T2_yield_raw)*singles_ps_value + (*c_T6_yield_raw))/good_total_charge;
  // jruns[run_str]["hms e yield"]     = double((*c_T2_yield)*singles_ps_value +
  // (*c_T6_yield))/good_total_charge; jruns[run_str]["hms e yield unc."]     =
  // double(std::sqrt(*c_T2_yield)*singles_ps_value + std::sqrt(*c_T6_yield))/good_total_charge;
  // jruns[run_str]["hms ps4 factor"]        = singles_ps_value;
  // jruns[run_str]["hms scaler yield"]      = hms_scaler_yield;
  // jruns[run_str]["hms scaler yield unc."] = hms_scaler_yield_unc;
  ////jruns[run_str]["hms live time"]         = hms_live_time;

  // std::ofstream json_output_file("db2/run_count_list.json");
  // json_output_file << std::setw(4) << jruns << "\n";

  // std::cout << " ----------------------------------------------    \n";
  // std::cout << " # of good coin  = "
  //          << ps_cor_hms_e_yield2
  //          << "    \n";
  // std::cout << " ----------------------------------------------    \n";
  // std::cout << " of  " << *c_n_events_total << " total triggers\n";
  // std::cout << " and " << *c_n_events_hms << " hms triggers\n";
  //// std::cout << " pions+kaons   : " << *coin_counts << "\n";
  //// std::cout << " pions         : " << *pion_count << "\n";
  //// std::cout << " random        : " << random_bg << "\n";
  //// std::cout << " ratio         : "  << pi_K_ratio << " \n";

  //// std::cout << " counts  : " << *hms_electron_counts << "\n";
  //// std::cout << " charge  : " << *total_charge << " uC\n";
  //// std::cout << " yield   : " << (*hms_electron_counts)  << " cnts, " <<
  //// hms_e_yield  << " cnts/uC\n"; std::cout << " singles : " <<
  //// (*hms_electron_counts2) << " cnts, " << hms_e_yield2 << " cnts/uC\n";
  //// std::cout << " coin    : " << (*coin_counts) << " cnts, " << coin_yield <<
  //// " cnts/uC\n";

  //// auto s_dc_x_fp           = d2.Histo1D({"s_dc_x_fp ","xy fp; x",100,-50,50},
  //// "P.dc.x_fp"); auto s_dc_y_fp           = d2.Histo1D({"s_dc_y_fp ","xy fp;
  //// y",100,-50,50}, "P.dc.y_fp"); auto s_dc_xp_fp          =
  //// d2.Histo1D({"s_dc_xp_fp","xy fp; xp",100,-50,50},"P.dc.xp_fp"); auto
  //// s_dc_yp_fp          = d2.Histo1D({"s_dc_yp_fp","xy fp;
  //// xp",100,-50,50},"P.dc.yp_fp");

  //// -----------------------------------------------------------
  ////
  // TCanvas *c = nullptr;
  // int b1 = 0;
  // int b2 = 0;
  // double hmax = 0.0;
  // THStack *hs = nullptr;
  //// TLatex latex;

  //// ---------------------------------------------------------
  ////
  // c = new TCanvas();
  // c->Divide(2, 2);
  // c->cd(1);

  //// This call starts the loop over the data.
  //// A DrawCopy is used so that the histogram is not deleted at the end of
  //// scope, and thus stays visible on the canvas.
  // h_EOverP_0->DrawCopy();

  // c->cd(2);
  // h_beta_0->DrawCopy();

  // c->cd(3);
  // gPad->SetLogy(true);
  // h_event_type->DrawCopy();
  // h_event_type_2->SetLineColor(2);
  // h_event_type_2->DrawCopy("same");

  // c->cd(4);
  ////h_coin_time->DrawCopy();
  ////h_coin_time2->SetLineColor(4);
  ////h_coin_time2->DrawCopy("same");
  ////h_coin_time3->SetLineColor(2);
  ////h_coin_time3->DrawCopy("same");

  // gPad->BuildLegend();

  // gSystem->mkdir("results/df_example", true);
  // c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) +
  //           ".pdf")
  //              .c_str());
  // c->SaveAs((std::string("results/df_example/c1_") + std::to_string(RunNumber) +
  //           ".png")
  //              .c_str());

  //// ---------------------------------------------------------
  ////
  // c = new TCanvas();
  // hs = new THStack("SHMS_cal", "SHMS calorimeter; E/P");

  // h_EOverP_0->SetLineColor(1);
  // h_EOverP_2->SetLineColor(4);
  // h_EprOverP_0->SetLineColor(8);
  // h_EprOverP_2->SetLineColor(6);

  // hs->Add((TH1 *)h_EOverP_0->Clone());
  // hs->Add((TH1 *)h_EOverP_2->Clone());

  // hs->Add((TH1 *)h_EprOverP_0->Clone());
  // hs->Add((TH1 *)h_EprOverP_2->Clone());

  // hs->Draw("nostack");
  // gPad->SetLogy(true);
  // gPad->BuildLegend();

  // c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
  //           ".pdf")
  //              .c_str());
  // c->SaveAs((std::string("results/df_example/c2_") + std::to_string(RunNumber) +
  //           ".png")
  //              .c_str());
}
