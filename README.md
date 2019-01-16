CSV online files
================

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


## A Beginners Guide

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
