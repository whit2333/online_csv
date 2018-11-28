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

## Guide

```
git clone
cd online_csv
ls
```

```
make_hallc_replay_symlinks -c
```
