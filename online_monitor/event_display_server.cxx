// ----------------------------------------------

R__LOAD_LIBRARY(libScandalizer.so)
#include "monitor/DisplayServer.h"

void event_display_server(){


  hallc::DisplayServer srv(8889,"cdaql1.jlab.org");
  srv.StartSocketServer(9091);

  srv.Run();

}
