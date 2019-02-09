// ----------------------------------------------

R__LOAD_LIBRARY(libScandalizer.so)
#include "monitor/DisplayServer.h"

void monitor_display_server(){


  hallc::MonitorDisplayServer srv;
  srv.StartSocketServer(9090);

  srv.Run();

}
