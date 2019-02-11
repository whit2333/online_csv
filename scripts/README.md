# scripts directory

## Data Quality  and Monitoring 

<dl>
<dt>
<%= markdown("
[good_coin_counter.cxx](good_coin_counter.cxx)
  ") %>
</dt>
<dd>  The coin (bean) counter for CSV sidis.</dd>

<dt>[good_hms_counter.cxx](good_hms_counter.cxx)</dt>
<dd>   The good HMS electron singles yield counter.</dd>
<dd>   Also looks at scalers for monitoring</dd>
<dt>[plot_hms_singles.cxx](plot_hms_singles.cxx)</dt>
<dd>   Plots the results of running `good_hms_counter.cxx` for a range of runs.</dd>

<dt>[make_human_table.cxx](make_human_table.cxx)</dt>
<dd>   Prints run list in a nice table</dd>

<dt>[make_daves_table.cxx](make_daves_table.cxx)</dt>
<dd>   Prints Dave's table</dd>

<dt>[plot_charge_vs_time.cxx](plot_charge_vs_time.cxx)</dt>
<dd>   Plots the experiment charge/yields vs time or run number.</dd>

<dt>[plot_waveforms.cxx](plot_waveforms.cxx)</dt>
<dd>   Plots waveforms on the display server.</dd>
</dl>

## Replay Scripts

`replay_production_coin_sidis.cxx`

## Examples and Templates

`json_load.cxx`
:   Shows how to use nlohmann json library.
`mongo_test.cxx`
:   Basic test for using mongodb



`bsa0.cxx`
`hms_s1_times.cxx`
`hms_trk_eff.cxx`
`rf_timing.cxx`
`scaler_sync.cxx`
`tracking_eff.cxx`

