<Request>
  <Email email="shuojia@jlab.org" request="false" job="true"/>
  <Project name="c-csv"/>
  <Track  name="analysis"/>
  <Name   name="kingfallcalib"/>
  <OS     name="centos7"/>
  <Memory space="3" unit="GB"/>
  <List name="runs">
  6068
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

