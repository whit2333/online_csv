// ----------------------------------------------
R__LOAD_LIBRARY(libsimple_epics.so)
#include "simple_epics/PVList.h"

R__LOAD_LIBRARY(libScandalizer.so)
#include "scandalizer/PostProcessors.h"
#include "monitor/DetectorDisplay.h"
#include "monitor/ExperimentMonitor.h"
#include "monitor/EventDisplays.h"

void detector_monitor(Int_t RunNumber = 6216, Int_t MaxEvent = 500000) {

  spdlog::set_level(spdlog::level::warn);
  spdlog::flush_every(std::chrono::seconds(5));

  // Get RunNumber and MaxEvent if not provided.
  if(RunNumber == 0) {
    cout << "Enter a Run Number (-1 to exit): ";
    if( RunNumber<=0 ) return;
  }
  if(MaxEvent == 0) {
    cout << "\nNumber of Events to analyze: ";
    cin >> MaxEvent;
    if(MaxEvent == 0) {
      cerr << "...Invalid entry\n";
      exit;
    }
  }

  // Create file name patterns.
  const char* RunFileNamePattern = "coin_all_%05d.dat";
  vector<TString> pathList;
  pathList.push_back(".");
  pathList.push_back("./raw");
  pathList.push_back("./raw/../raw.copiedtotape");
  pathList.push_back("./cache");

  //const char* RunFileNamePattern = "raw/coin_all_%05d.dat";
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

  // ========
  //  SHMS 
  // ========

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

  auto trkEff = new hcana::TrackingEfficiency("P.trackEfficiency", "SHMS tracking Efficiency", "P.hod");
  gHaPhysics->Add(trkEff);
  // Add physics modules
  // Calculate reaction point
  THaReactionPoint* prp = new THaReactionPoint("P.react", "SHMS reaction point", "P", "P.rb");
  gHaPhysics->Add(prp);
  // Calculate extended target corrections
  THcExtTarCor* pext = new THcExtTarCor("P.extcor", "HMS extended target corrections", "P", "P.react");
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

  // ========
  //  HMS 
  // ========

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
  THcExtTarCor* hext = new THcExtTarCor("H.extcor", "HMS extended target corrections", "H", "H.react");
  gHaPhysics->Add(hext);
  // Calculate golden track quantities
  THaGoldenTrack* hgtr = new THaGoldenTrack("H.gtr", "HMS Golden Track", "H");
  gHaPhysics->Add(hgtr);
  // Calculate the hodoscope efficiencies
  THcHodoEff* heff = new THcHodoEff("hhodeff", "HMS hodo efficiency", "H.hod");
  gHaPhysics->Add(heff);

  // Add event handler for scaler events
  THcScalerEvtHandler *hscaler = new THcScalerEvtHandler("H", "Hall C scaler event type 4");  
  hscaler->AddEvtType(2);
  hscaler->AddEvtType(4);
  hscaler->AddEvtType(5);
  hscaler->AddEvtType(6);
  hscaler->AddEvtType(7);
  hscaler->AddEvtType(129);
  hscaler->SetDelayedType(129);
  hscaler->SetUseFirstEvent(kTRUE);
  gHaEvtHandlers->Add(hscaler);

  // =================================
  //  Kinematics Modules
  // =================================

  // Add Physics Module to calculate primary (scattered electrons) beam kinematics
  THcPrimaryKine* hkin_primary = new THcPrimaryKine("H.kin.primary", "HMS Single Arm Kinematics", "H", "H.rb");
  gHaPhysics->Add(hkin_primary);
  // Add Physics Module to calculate secondary (scattered hadrons) beam kinematics
  THcSecondaryKine* pkin_secondary = new THcSecondaryKine("P.kin.secondary", "SHMS Single Arm Kinematics", "P", "H.kin.primary");
  gHaPhysics->Add(pkin_secondary);
  
  // =================================
  //  Global Objects & Event Handlers
  // =================================

  // Add trigger apparatus
  THaApparatus* TRG = new THcTrigApp("T", "TRG");
  gHaApps->Add(TRG);
  // Add trigger detector to trigger apparatus
  THcTrigDet* coin = new THcTrigDet("coin", "Coincidence Trigger Information");
  // Suppress missing reference time warnings for these event types
  coin->SetEvtType(1);
  coin->AddEvtType(2);
  TRG->AddDetector(coin); 

  // Add helicity detector to grigger apparatus
  THcHelicity* helicity = new THcHelicity("helicity","Helicity Detector");
  TRG->AddDetector(helicity);

  // Add event handler for prestart event 125.
  THcConfigEvtHandler* ev125 = new THcConfigEvtHandler("HC", "Config Event type 125");
  gHaEvtHandlers->Add(ev125);
  // Add event handler for EPICS events
  THaEpicsEvtHandler* hcepics = new THaEpicsEvtHandler("epics", "HC EPICS event type 180");
  gHaEvtHandlers->Add(hcepics);

  // Add coin physics module 
  THcCoinTime* coinTime = new THcCoinTime("CTime", "Coincidende Time Determination", "P", "H", "T.coin");
  gHaPhysics->Add(coinTime);
 
  // -----------------------------------------------------------
  // Scandalizer 
  // -----------------------------------------------------------
  //
  hcana::Scandalizer* analyzer = new hcana::Scandalizer;
  //analyzer->_skip_events = 100;
  analyzer->SetCodaVersion(2);
  //analyzer->EnableBenchmarks(true);

  // The following analyzes the first 2000 events (for pedestals, is required) 
  // then  repeatedly skips 3000 events and processes 1000.
  //auto pp0 = new hallc::scandalizer::SkipPeriodicAfterPedestal();
  auto pp0 = new hallc::scandalizer::SkipAfterPedestal();
  pp0->_analyzer = analyzer;

  auto ddisplay = new hallc::MonitoringDisplay(RunNumber);
  //auto exp_monitor = new hallc::ExperimentMonitor(RunNumber);
  ddisplay->_data._root_folder = "/detectors/";

  auto plt1 = ddisplay->CreateDisplayPlot(
      "shms/dc/", "SHMS_DC_Occupancy",
      [&](hallc::DisplayPlot& plt) {
        plt._plot_data._canvas = new TCanvas(plt.GetName().c_str(), plt.GetName().c_str());
        plt._plot_data._hists1 = {new TH1F("h_dc1_occ", "h_dc1_occ", 30, 0, 30),
                                  new TH1F("h_dc2_occ", "h_dc2_occ", 30, 0, 30)};
        plt._plot_data._canvas->Divide(2, 1);
        plt._plot_data._canvas->cd(1);
        plt._plot_data._hists1[0]->Draw();
        plt._plot_data._canvas->cd(2);
        plt._plot_data._hists1[1]->SetLineColor(2);
        plt._plot_data._hists1[1]->Draw();
        return 0;
      },
      [&](hallc::DisplayPlot& plt) {
        int           shmsDC1Planes_nhits = 0;
        int           shmsDC2Planes_nhits = 0;
        for (int ip = 0; ip < 6; ip++) {
          shmsDC1Planes_nhits += pdc->GetPlane(ip)->GetNHits();
        }
        for (int ip = 6; ip < 12; ip++) {
          shmsDC2Planes_nhits += pdc->GetPlane(ip)->GetNHits();
        }
        plt._plot_data._hists1[0]->Fill(shmsDC1Planes_nhits);
        plt._plot_data._hists1[1]->Fill(shmsDC2Planes_nhits);
        return 0;
      });

  auto hms_mon_pp                = new hallc::ExperimentMonitorPostProcess(ddisplay);
  hms_mon_pp->_analyzer          = analyzer;
  hms_mon_pp->_hod               = hhod;
  hms_mon_pp->_hgcer             = hcer;
  hms_mon_pp->_dc                = hdc;
  hms_mon_pp->_spectrometer_name = "HMS";

  auto shms_mon_pp                = new hallc::ExperimentMonitorPostProcess(ddisplay);
  shms_mon_pp->_analyzer          = analyzer;
  shms_mon_pp->_hod               = phod;
  shms_mon_pp->_hgcer             = phgcer;
  shms_mon_pp->_dc                = pdc;
  shms_mon_pp->_spectrometer_name = "SHMS";

  //hallc::DisplayPostProcess * display_pp = new hallc::DisplayPostProcess(ddisplay);
  //analyzer->AddPostProcess(display_pp);

  analyzer->AddPostProcess(pp0);
  analyzer->AddPostProcess(hms_mon_pp);
  analyzer->AddPostProcess(shms_mon_pp);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  THcRun* run = new THcRun( pathList, Form(RunFileNamePattern, RunNumber) );

  // Set to read in Hall C run database parameters
  run->SetRunParamClass("THcRunParameters");
  
  // Eventually need to learn to skip over, or properly analyze the pedestal events
  run->SetEventRange(1, MaxEvent); // Physics Event number, does not include scaler or control events.
  run->SetNscan(1);
  run->SetDataRequired(0x7);
  //run->Print();

  // Define the analysis parameters
  TString ROOTFileName = Form(ROOTFileNamePattern, RunNumber, MaxEvent);
  analyzer->SetCountMode(2);  // 0 = counter is # of physics triggers
                              // 1 = counter is # of all decode reads
                              // 2 = counter is event number

  analyzer->SetEvent(event);
  // Set EPICS event type
  analyzer->SetEpicsEvtType(180);
  // Define crate map
  analyzer->SetCrateMapFileName("MAPS/db_cratemap.dat");
  // Define output ROOT file
  analyzer->SetOutFile(ROOTFileName.Data());
  // Define DEF-file+
  analyzer->SetOdefFile("online_monitor/scandalizer_monitor.def");
  // Define cuts file
  analyzer->SetCutFile("UTIL_SIDIS/DEF-files/coin_production_sidis_cuts.def");  // optional
  // File to record accounting information for cuts
  analyzer->SetSummaryFile(Form("REPORT_OUTPUT/COIN/PRODUCTION/summary_production_%d_%d.report", RunNumber, MaxEvent));  // optional
  // Start the actual analysis.
  analyzer->Process(run);
  // Create report file from template
  analyzer->PrintReport("TEMPLATES/COIN/PRODUCTION/coin_production.template",
  			Form("REPORT_OUTPUT/COIN/PRODUCTION/replay_coin_production_%d_%d.report", RunNumber, MaxEvent));  // optional

  //gSystem->mkdir(fmt::format("monitoring/{}",RunNumber).c_str(),true);
  //std::string plt1_output = fmt::format("monitoring/{}/plot_{}.json",RunNumber,plt1->_id);
  //plt1->_canvas->SaveAs(plt1_output.c_str());
  //delete analyzer;
  //
  //display_pp->_run_number = RunNumber;
  //display_pp->Close();
}
