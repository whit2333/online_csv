CSV online files
================

## Quick links
* [`scripts` directory](scripts/README.md)
* [`online_monitor` directory](online_monitor/README.md)

## COIN/SHMS/HMS replay

```
./bin/hc_shms_replay -r run_number -n Nevents
./bin/hc_shms_replay -h
./bin/hc_coin_replay -r run_number -n Nevents
```

### Tips

#### Running 10 replays in parallel using gnu parallel
```
parallel -N1 -j10 ./bin/hc_shms_replay -r {}  ::: $(seq 7100 7142)
```

###
```
camonitor ibcm1 | stdbuf -o0 awk '{print $4}' | feedgnuplot  --stream --lines  \
--points --title "Test plot" --xlen 100  --unset grid --terminal 'dumb'
```


## Quick Start

[![asciicast](https://asciinema.org/a/222886.svg)](https://asciinema.org/a/222886)

```
git clone https://github.com/whit2333/online_csv.git
git clone https://github.com/whit2333/hallc_replay_sidis_fall18.git
cd online_csv
make_hallc_replay_symlinks -b ../hallc_replay_sidis_fall18
make_hallc_replay -c
```
Now  `online_csv`  is a replay directory that uses the "standard" 
`hallc_replay`. It can be renamed (e.g., "official_offline_replay") and it 
still has https://github.com/whit2333/online_csv.git as the upstream. Also, 
changes made to the symlinked files will be in the replay. Good changes should 
be pushed upstream for this repo too.


## Usage

The main files in the repository are found in the `db2` directory.
* `db2/run_list.json`
* `db2/run_count_list.json`
* `db2/run_list_extra.json`
* `db2/auto_standard.kinematics`

Using `hallc_tools`'s `make_hallc_replay_symlinks` and `make_hallc_replay` 
commands a replay directory can be created from these files.

## Scandalizer 

The `scandalizer` is a modified analyzer that has numerous benefits.

The current `scandalizer` is run through the replay script 
`online_bin/scandalizer` which runs a replay script,
`online_monitor/scandalizer_test.cxx` (for now).


```
online_bin/scandalizer -r <run_number> -n <num>
OPTIONS:
            -r,--run           Required run number
            -n,--n-events      Required number of eventsrun number
            -c,--container     Run using singularity container 
            -C,--calibrate     Replay N-calibrate events and run calibrations before full replay
            --N-calibrate      Number of events to replay for calibration (Default: 50000)
```


## A Beginners Guide (to the basics)

### Environment

Advice: use bash and tmux. 

If you are on the farm and never have used bash. A good setup can be configured 
by doing the following: (this will overwrite existing files)
```
cp /group/c-csv/local/stow/Templates/template.bashrc  $HOME/.bashrc
cp /group/c-csv/local/stow/Templates/template.bash_aliases  $HOME/.bash_aliases
cp /group/c-csv/local/stow/Templates/inputrc  $HOME/.inputrc
cp /group/c-csv/local/stow/Templates/template.tmux.conf  $HOME/.tmux.conf
cp -r /group/c-csv/local/stow/Templates/template.tmux  $HOME/.tmux
```
If you want to change the colors in `.tmux.conf` run
```
/group/c-csv/local/stow/Templates/print_colors
```
and change the colour number on the line that contains
```
set-option -g status-bg colour25
```

Type `bash` to start a bash session and `tmux` to start a new tmux session.

Type `module avail` to see all the software modules and `module list` to see 
those that are loaded.  See [this page for more 
details](https://hallcweb.jlab.org/wiki/index.php/CSV_software#Working_on_the_farm).


### Getting started as a replay directory

```
git clone git@github.com:whit2333/online_csv.git
cd online_csv
ls
```
Then  (using hallc_tools)

```
make_hallc_replay_symlinks -c
```

## cdaq notes

cdaql1 129.57.168.41

