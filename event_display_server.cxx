// ----------------------------------------------

R__LOAD_LIBRARY(libScandalizer.so)
#include "monitor/DisplayServer.h"

void event_display_server(){


  hallc::DisplayServer srv;
  srv.StartSocketServer();

  srv.Run();

}
