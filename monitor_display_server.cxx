// ----------------------------------------------

R__LOAD_LIBRARY(libScandalizer.so)
#include "monitor/DisplayServer.h"

void monitor_display_server(){


  hallc::DisplayServer srv(8888,"129.57.168.41");
  srv.StartSocketServer(9090);

  srv.Run();

}
