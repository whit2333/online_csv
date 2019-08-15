<Request>
  <Email email="shuojia@jlab.org" request="false" job="true"/>
  <Project name="c-csv"/>
  <Track  name="reconstruction"/>
  <Name   name="kingfalladd"/>
  <OS     name="centos7"/>
  <DiskSpace space = "50" unit = "GB"/>
  <Memory space="3" unit="GB"/>
  <List name="runs">
  6073
  </List>
  <ForEach list="runs">
    <Job>
      <Input src="mss:/mss/hallc/spring17/raw/coin_all_0${runs}.dat" dest="coin_all_0${runs}.dat"/>
      <!-- Properties overridden here -->
      <Command><![CDATA[
/bin/bash <<EOF
source /home/shuojia/.bashrc
cd /group/c-csv/shuo/hallc_replay_sidis_fall18/


hcana -q
"SCRIPTS/COIN/PRODUCTION/replay_production_coin_hElec_pProt.C(${runs},100000)" 
hcana -q
"SCRIPTS/COIN/PRODUCTION/replay_production_coin_hElec_pProt.C(${runs},50000)" 
hcana -q
"SCRIPTS/COIN/PRODUCTION/replay_production_coin_hElec_pProt.C(${runs},30000)" 
hcana -q
"SCRIPTS/COIN/PRODUCTION/replay_production_coin_hElec_pProt.C(${runs},10000)" 

EOF
        ]]></Command>
    </Job>
  </ForEach> 

</Request>

