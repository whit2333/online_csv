// ----------------------------------------------

R__LOAD_LIBRARY(libScandalizer.so)
#include "monitor/DisplayServer.h"



void display_server(){


  hallc::DisplayServer srv;

  srv.Run();

}
