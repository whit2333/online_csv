#!/usr/bin/env python3

import argparse
import json
import os
from hallc.error import HallCError

outdir='results/skim/'

class InvalidSettingError(HallCError):
    def __init(self, setting):
        self.message = 'Unknown setting: ' + setting

def skim(setting, db, script):
    run_list = db[setting]['good_run_list']
    print('Creating skim for setting ', setting)
    for run in run_list:
        print('Skimming run {}'.format(run))
        os.system('root -b -q "{}({})"'.format(script, run))
    if len(run_list) > 0:
        if_names = ['{}/skim_coinel_{}.root'.format(outdir, run) for run in run_list]
        odat_name = '{}/skim_coinel_{}.root'.format(outdir, setting)
        print('Combining runs into: ', odat_name)
        os.system('hadd {} {}'.format(odat_name, ' '.join(if_names)))
    else:
        print('No runs taken for ', setting)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Create skimmed root files from good runs.')
    parser.add_argument('-s', '--setting',
            default='all',
            help='ID of the setting we want to skim (e.g. "phase1"), (def.: "all")',
            dest='setting')
    parser.add_argument('-c', '--command',
            default='scripts/skim.cxx',
            help='Skim script to be used (def.: "scripts/skim.cxx")',
            dest='script')
    parser.add_argument('-d', '--database',
            default='db2/jpsi_status.json',
            help='J/psi-007 status database (default: db2/jpsi_status.json)',
            dest='db')
    args = parser.parse_args()

    print('J/psi-007 SKIM')
    print('Creating e+e- coincidence skim')

    print('Reading configuration from: ', args.db)
    with open(args.db, 'r') as dbfile:
        dbstring = dbfile.readlines()
        db = json.loads(''.join(dbstring))

        if args.setting is 'all':
            for setting in db:
                skim(setting, db, args.script)
        else:
            if args.setting not in db:
                raise InvalidSettingError(args.setting)
            skim(args.setting, db, args.script)
