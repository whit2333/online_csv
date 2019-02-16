// ----------------------------------------------

R__LOAD_LIBRARY(libScandalizer.so)
#include "monitor/DisplayServer.h"

void monitor_display_server(){


  hallc::DisplayServer srv(9990,"cdaql1.jlab.org");
  srv.StartSocketServer(9890);

  srv.Run();

}
