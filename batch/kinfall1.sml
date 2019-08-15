<Request>
  <Email email="shuojia@jlab.org" request="false" job="true"/>
  <Project name="c-csv"/>
  <Track  name="reconstruction"/>
  <Name   name="kingfalladd"/>
  <OS     name="centos7"/>
  <DiskSpace space = "50" unit = "GB"/>
  <Memory space="3" unit="GB"/>
  <List name="runs">
  6137
  6138
  6497
  6498
  6499
  6271
  6272
  6273
  6348
  6349
  6350
  6351
  6378
  6379
  6380
  6381
  
  </List>
  <ForEach list="runs">
    <Job>
      <Input src="mss:/mss/hallc/spring17/raw/coin_all_0${runs}.dat" dest="coin_all_0${runs}.dat"/>
      <!-- Properties overridden here -->
      <Command><![CDATA[
/bin/bash <<EOF
source /home/shuojia/.bashrc
cd /group/c-csv/shuo/hallc_replay_sidis_fall18/


hcana -q "SCRIPTS/COIN/PRODUCTION/replay_production_coin_hElec_pProt.C(${runs},-1)" 

EOF
        ]]></Command>
    </Job>
  </ForEach> 

</Request>

