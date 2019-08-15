<Request>
  <Email email="shuojia@jlab.org" request="false" job="true"/>
  <Project name="c-csv"/>
  <Track  name="analysis"/>
  <Name   name="kingfallcalib"/>
  <OS     name="centos7"/>
  <Memory space="3" unit="GB"/>
  <List name="runs">
  6068
  6073
  6088
  6091
  6045
  6548
  6550
  6111
  6115
  6124
  6136
  6486
  6489
  6490
  6496
  6245
  6248
  6263
  6270
  6451
  6459
  6465
  6473
  6290
  6292
  6306
  6348
  6359
  6367
  6375
  6378
  6524
  6527
  
  </List>
  <ForEach list="runs">
    <Job>
      <!-- Properties overridden here -->
      <Command><![CDATA[
/bin/bash <<EOF
source /home/shuojia/.bashrc
cd /group/c-csv/shuo/hallc_replay_sidis_fall18/


hcana -q
"mycalibration/shms_cal_calib/pcal_calib.cpp+(\"coin_replay_production_${runs}_-1\")" 

EOF
        ]]></Command>
    </Job>
  </ForEach> 

</Request>

