#include <fmt/core.h>
#include <fmt/ostream.h>
R__LOAD_LIBRARY(libfmt.so)

#include <chrono>
#include <iostream>

using namespace std;

#include "spdlog/spdlog.h"

//#include "THaPostProcess.h"
//
R__LOAD_LIBRARY(libHallC.so)
#include "hcana/HallC_Data.h"
//
//R__LOAD_LIBRARY(libScandalizer.so)
//#include "monitor/DetectorDisplay.h"
//#include "monitor/DisplayPlots.h"
//#include "monitor/MonitoringDisplay.h"
//#include "scandalizer/PostProcessors.h"
//#include "scandalizer/ScriptHelpers.h"
//
//#include "THaPostProcess.h"
//#include "monitor/ExperimentMonitor.h"
//#include "scandalizer/PostProcessors.h"
#include "THcAnalyzer.h"
#include "THaCut.h"
#include "THcGlobals.h"
#include "THcHallCSpectrometer.h"
#include "THcDetectorMap.h"
#include "THcCherenkov.h"
#include "THcDC.h"
#include "THcHodoscope.h"
#include "THcParmList.h"
#include "THaGoldenTrack.h"
#include "THcHodoEff.h"
#include "THcScalerEvtHandler.h"
#include "THcShower.h"
#include "THaReactionPoint.h"
#include "THcExtTarCor.h"
#include "THcRasteredBeam.h"
#include "THcRun.h"
#include "THcCoinTime.h"
#include "THcConfigEvtHandler.h"
#include "THcTrigDet.h"
#include "THcTrigApp.h"
#include "THcSecondaryKine.h"
#include "THcAerogel.h"
#include "THcPrimaryKine.h"
void replay_production_coin(Int_t RunNumber = 7224, Int_t MaxEvent = 50000) {

  //hallc::helper::script_requires_hcana();

  spdlog::set_level(spdlog::level::trace);
  spdlog::flush_every(std::chrono::seconds(5));

  // Get RunNumber and MaxEvent if not provided.
  if (RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    cin >> RunNumber;
    if (RunNumber <= 0)
      return;
  }
  if (MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if (MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      throw 0;
    }
  }

  // Create file name patterns.
  const char*     RunFileNamePattern = "coin_all_%05d.dat";
  vector<TString> pathList;
  pathList.push_back(".");
  pathList.push_back("./DATA/raw");
  pathList.push_back("./raw");
  pathList.push_back("./raw.copiedtotape");
  pathList.push_back("./raw/../raw.copiedtotape");
  pathList.push_back("./cache");

  // const char* RunFileNamePattern = "raw/coin_all_%05d.dat";
  //const char* ROOTFileNamePattern = "/lcrc/project/jlab/data/hallc/jpsi-007/replay/data/coin_replay_production_%d_%d.root";
  const char* ROOTFileNamePattern = "ROOTfiles/coin_replay_production_%d_%d.root";

  // Load global parameters
  gHcParms->Define("gen_run_number", "Run Number", RunNumber);
  gHcParms->AddString("g_ctp_database_filename", "DBASE/COIN/standard.database");
  gHcParms->Load(gHcParms->GetString("g_ctp_database_filename"), RunNumber);
  gHcParms->Load(gHcParms->GetString("g_ctp_parm_filename"));
  gHcParms->Load(gHcParms->GetString("g_ctp_kinematics_filename"), RunNumber);
  // Load params for COIN trigger configuration
  gHcParms->Load(gHcParms->GetString("g_ctp_trigdet_filename"));
  // Load fadc debug parameters
  gHcParms->Load("PARAM/HMS/GEN/h_fadc_debug.param");
  gHcParms->Load("PARAM/SHMS/GEN/p_fadc_debug.param");
  // Load the Hall C detector map
  gHcDetectorMap = new THcDetectorMap();
  gHcDetectorMap->Load(gHcParms->GetString("g_ctp_map_filename"));

  //=:=:=:=
  // SHMS
  //=:=:=:=

  // Set up the equipment to be analyzed.
  THcHallCSpectrometer* SHMS = new THcHallCSpectrometer("P", "SHMS");
  SHMS->SetEvtType(1);
  SHMS->AddEvtType(4);
  SHMS->AddEvtType(5);
  SHMS->AddEvtType(6);
  SHMS->AddEvtType(7);
  gHaApps->Add(SHMS);
  // Add Noble Gas Cherenkov to SHMS apparatus
  THcCherenkov* pngcer = new THcCherenkov("ngcer", "Noble Gas Cherenkov");
  SHMS->AddDetector(pngcer);
  // Add drift chambers to SHMS apparatus
  THcDC* pdc = new THcDC("dc", "Drift Chambers");
  SHMS->AddDetector(pdc);
  // Add hodoscope to SHMS apparatus
  THcHodoscope* phod = new THcHodoscope("hod", "Hodoscope");
  SHMS->AddDetector(phod);
  // Add Heavy Gas Cherenkov to SHMS apparatus
  THcCherenkov* phgcer = new THcCherenkov("hgcer", "Heavy Gas Cherenkov");
  SHMS->AddDetector(phgcer);
  // Add Aerogel Cherenkov to SHMS apparatus
  THcAerogel* paero = new THcAerogel("aero", "Aerogel");
  SHMS->AddDetector(paero);
  // Add calorimeter to SHMS apparatus
  THcShower* pcal = new THcShower("cal", "Calorimeter");
  SHMS->AddDetector(pcal);

  // Add rastered beam apparatus
  THaApparatus* pbeam = new THcRasteredBeam("P.rb", "Rastered Beamline");
  gHaApps->Add(pbeam);
  // Add physics modules
  // Calculate reaction point
  THaReactionPoint* prp = new THaReactionPoint("P.react", "SHMS reaction point", "P", "P.rb");
  gHaPhysics->Add(prp);
  // Calculate extended target corrections
  THcExtTarCor* pext =
      new THcExtTarCor("P.extcor", "HMS extended target corrections", "P", "P.react");
  gHaPhysics->Add(pext);
  // Calculate golden track quantites
  THaGoldenTrack* pgtr = new THaGoldenTrack("P.gtr", "SHMS Golden Track", "P");
  gHaPhysics->Add(pgtr);
  // Calculate the hodoscope efficiencies
  THcHodoEff* peff = new THcHodoEff("phodeff", "SHMS hodo efficiency", "P.hod");
  gHaPhysics->Add(peff);

  // Add event handler for scaler events
  THcScalerEvtHandler* pscaler = new THcScalerEvtHandler("P", "Hall C scaler event type 1");
  pscaler->AddEvtType(1);
  pscaler->AddEvtType(4);
  pscaler->AddEvtType(5);
  pscaler->AddEvtType(6);
  pscaler->AddEvtType(7);
  pscaler->AddEvtType(129);
  pscaler->SetDelayedType(129);
  pscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(pscaler);

  //=:=:=
  // HMS
  //=:=:=

  // Set up the equipment to be analyzed.
  THcHallCSpectrometer* HMS = new THcHallCSpectrometer("H", "HMS");
  HMS->SetEvtType(2);
  HMS->AddEvtType(4);
  HMS->AddEvtType(5);
  HMS->AddEvtType(6);
  HMS->AddEvtType(7);
  gHaApps->Add(HMS);
  // Add drift chambers to HMS apparatus
  THcDC* hdc = new THcDC("dc", "Drift Chambers");
  HMS->AddDetector(hdc);
  // Add hodoscope to HMS apparatus
  THcHodoscope* hhod = new THcHodoscope("hod", "Hodoscope");
  HMS->AddDetector(hhod);
  // Add Cherenkov to HMS apparatus
  THcCherenkov* hcer = new THcCherenkov("cer", "Heavy Gas Cherenkov");
  HMS->AddDetector(hcer);
  // Add Aerogel Cherenkov to HMS apparatus
  // THcAerogel* haero = new THcAerogel("aero", "Aerogel");
  // HMS->AddDetector(haero);
  // Add calorimeter to HMS apparatus
  THcShower* hcal = new THcShower("cal", "Calorimeter");
  HMS->AddDetector(hcal);

  // Add rastered beam apparatus
  THaApparatus* hbeam = new THcRasteredBeam("H.rb", "Rastered Beamline");
  gHaApps->Add(hbeam);
  // Add physics modules
  // Calculate reaction point
  THaReactionPoint* hrp = new THaReactionPoint("H.react", "HMS reaction point", "H", "H.rb");
  gHaPhysics->Add(hrp);
  // Calculate extended target corrections
  THcExtTarCor* hext =
      new THcExtTarCor("H.extcor", "HMS extended target corrections", "H", "H.react");
  gHaPhysics->Add(hext);
  // Calculate golden track quantities
  THaGoldenTrack* hgtr = new THaGoldenTrack("H.gtr", "HMS Golden Track", "H");
  gHaPhysics->Add(hgtr);
  // Calculate the hodoscope efficiencies
  THcHodoEff* heff = new THcHodoEff("hhodeff", "HMS hodo efficiency", "H.hod");
  gHaPhysics->Add(heff);

  // Add event handler for scaler events
  THcScalerEvtHandler* hscaler = new THcScalerEvtHandler("H", "Hall C scaler event type 4");
  hscaler->AddEvtType(2);
  hscaler->AddEvtType(4);
  hscaler->AddEvtType(5);
  hscaler->AddEvtType(6);
  hscaler->AddEvtType(7);
  hscaler->AddEvtType(129);
  hscaler->SetDelayedType(129);
  hscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(hscaler);

  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
  // Kinematics Modules
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

  // Add Physics Module to calculate primary (scattered electrons) beam kinematics
  THcPrimaryKine* hkin_primary =
      new THcPrimaryKine("H.kin.primary", "HMS Single Arm Kinematics", "H", "H.rb");
  gHaPhysics->Add(hkin_primary);
  // Add Physics Module to calculate secondary (scattered hadrons) beam kinematics
  THcSecondaryKine* pkin_secondary =
      new THcSecondaryKine("P.kin.secondary", "SHMS Single Arm Kinematics", "P", "H.kin.primary");
  gHaPhysics->Add(pkin_secondary);

  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=
  // Global Objects & Event Handlers
  //=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=:=

  // Add trigger apparatus
  THaApparatus* TRG = new THcTrigApp("T", "TRG");
  gHaApps->Add(TRG);
  // Add trigger detector to trigger apparatus
  THcTrigDet* coin = new THcTrigDet("coin", "Coincidence Trigger Information");
  // Suppress missing reference time warnings for these event types
  coin->SetEvtType(1);
  coin->AddEvtType(2);
  TRG->AddDetector(coin);
  // Add event handler for prestart event 125.
  THcConfigEvtHandler* ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  gHaEvtHandlers->Add(ev125);
  // Add event handler for EPICS events
  THaEpicsEvtHandler* hcepics = new THaEpicsEvtHandler("epics", "HC EPICS event type 180");
  gHaEvtHandlers->Add(hcepics);

  // Add coin physics module
  THcCoinTime* coinTime =
      new THcCoinTime("CTime", "Coincidende Time Determination", "P", "H", "T.coin");
  gHaPhysics->Add(coinTime);

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Acpparatus's and PhysicsModules,
  // and executes the output routines.
  THcAnalyzer* analyzer = new THcAnalyzer;

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  THcRun* run = new THcRun(pathList, Form(RunFileNamePattern, RunNumber));

  // Set to read in Hall C run database parameters
  run->SetRunParamClass("THcRunParameters");

  // Eventually need to learn to skip over, or properly analyze the pedestal events
  run->SetEventRange(1, MaxEvent); // Physics Event number, does not include scaler or control events.
  run->SetNscan(1);
  run->SetDataRequired(0x7);
  run->Print();

  // Define the analysis parameters
  TString ROOTFileName = Form(ROOTFileNamePattern, RunNumber, RunNumber, MaxEvent);
  analyzer->SetCountMode(2); // 0 = counter is # of physics triggers
                             // 1 = counter is # of all decode reads
                             // 2 = counter is event number

  //auto ddisplay      = new hallc::MonitoringDisplay({"cdaql1.jlab.org",9091},RunNumber);
  //ddisplay->_data._root_folder = "/replay_production_coin/";

  //auto hms_event_display  = new hallc::event_display::HMSEventDisplay(ddisplay);
  //hms_event_display->_hod = hhod;
  //hms_event_display->_run_number = RunNumber;

  //auto shms_event_display  = new hallc::event_display::SHMSEventDisplay(ddisplay);
  //shms_event_display->_hod = phod;
  //shms_event_display->_run_number = RunNumber;

  ////ddisplay->InitAll();
  ////ddisplay->UpdateAll();

  //analyzer->AddPostProcess(hms_event_display);
  //analyzer->AddPostProcess(shms_event_display);

  analyzer->SetEvent(event);
  // Set EPICS event type
  analyzer->SetEpicsEvtType(180);
  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap.dat");
  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());
  // Define DEF-file+
  // analyzer->SetOdefFile("UTIL_SIDIS/DEF-files/coin_production_sidis.def");
  analyzer->SetOdefFile("DEF-files/COIN/PRODUCTION/coin_production_hElec_pPion_csv.def"); //most skimed
  //analyzer->SetOdefFile("DEF-files/COIN/PRODUCTION/coin_production_hElec_pPion.def");  //for tracking efficiency
  //analyzer->SetOdefFile("DEF-files/COIN/PRODUCTION/coin_production_ep.def");  //everything
  // Define cuts file
  // analyzer->SetCutFile("UTIL_SIDIS/DEF-files/coin_production_sidis_cuts.def");  // optional
  analyzer->SetCutFile("DEF-files/COIN/PRODUCTION/CUTS/coin_production_cuts.def"); // optional
  // File to record accounting information for cuts
  analyzer->SetSummaryFile(Form("/lcrc/project/jlab/data/hallc/csv/replay/report/COIN/summary_production_%d_%d.report",
                                RunNumber, MaxEvent)); // optional
  // Start the actual analysis.
  analyzer->Process(run);
  // Create report file from template
  analyzer->PrintReport("TEMPLATES/COIN/PRODUCTION/coin_production.template",
                        Form("/lcrc/project/jlab/data/hallc/csv/replay/report/COIN/replay_coin_production_%d_%d.report",
                             RunNumber, MaxEvent)); // optional

  //ddisplay->CreateDisplayPlot(
  //    "hms", "tar",
  //    [&](hallc::DisplayPlot& plt) {
  //      auto c = plt.SetCanvas( new TCanvas(plt.GetName().c_str(), plt.GetName().c_str()));
  //      //c = plt._plot_data._canvas; // new TCanvas("c1", "c1", 1600, 1200);
  //      c->Divide(2, 1);
  //      c->cd(1);
  //      TH2* h2 = (TH2*) gROOT->FindObject("hgtrxp_vs_pgtrxp");
  //      h2->Draw("colz");

  //      c->cd(2);
  //      // h_beta_0->DrawCopy();

  //      //h_EprOverP_0->DrawCopy();
  //      //h_EprOverP_2->SetLineColor(2);
  //      //h_EprOverP_2->DrawCopy("same");
  //      return 0;
  //    },
  //    [&](hallc::DisplayPlot& plt) { return 0; });
  //ddisplay->UpdateAll();
}
