#!/usr/bin/env python3

import argparse
import json
import os
from hallc.error import HallCError

outdir='results/skim'
script={
    'el': 'scripts/skim.cxx+',
    'mu': 'scripts/muskim.cxx+'}

class InvalidSettingError(HallCError):
    def __init__(self, setting):
        self.message = 'Unknown setting: ' + setting
class InvalidParticleError(HallCError):
    def __init__(self, part):
        self.message = 'Unknown particle: ' + part

def skim(setting, db, script, particle, force_update=False):
    run_list = db[setting]['good_run_list']
    print('Creating skim for setting ', setting)
    skims = []
    for run in run_list:
        run_working = '{}/skim_coin{}_{}.root'.format(outdir, particle, run)
        run_full = '{}/{}/skim_coin{}_{}.root'.format(outdir, setting, particle, run)
        skims.append(run_full)
        if not force_update and os.path.isfile('{}/{}/skim_coin{}_{}.root'.format(outdir, setting, particle, run)):
            print('Skipping run', run, ' (file already present)')
        else:
            print('Skimming run {}'.format(run))
            os.system('root -b -q "{}({})"'.format(script, run))
            os.system('mv {} {}'.format(run_working, run_full))
    if len(run_list) > 0:
        odat_name = '{}/skim_coin{}_{}.root'.format(outdir, particle, setting)
        print('Combining runs into: ', odat_name)
        os.system('hadd -f {} {}'.format(odat_name, ' '.join(skims)))
    else:
        print('No runs taken for ', setting)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description='Create skimmed root files from good runs.')
    parser.add_argument('-s', '--setting',
            default='all',
            help='ID of the setting we want to skim (e.g. "phase1"), (def.: "all")',
            dest='setting')
    parser.add_argument('-p', '--particle',
            default='el',
            help='Particle to look for ("el" or "mu"), (def.: "el")',
            dest='part')
    parser.add_argument('-d', '--database',
            default='db2/jpsi_status.json',
            help='J/psi-007 status database (default: db2/jpsi_status.json)',
            dest='db')
    parser.add_argument('-f', '--force-update',
            action='store_true',
            help='Force update mode: also re-skim previously skimmed files.',
            dest='force_update')
    args = parser.parse_args()

    print('J/psi-007 SKIM')

    if args.part not in script:
        raise InvalidParticleError(args.part)
    if args.part == 'el':
        print('Creating e+e- coincidence skim')
    elif args.part == 'mu':
        print('Creating mu+mu- coincidence skim')
    else:
        print('Creating skim for', args.part)
        

    print('Reading configuration from: ', args.db)
    with open(args.db, 'r') as dbfile:
        dbstring = dbfile.readlines()
        db = json.loads(''.join(dbstring))

        if args.setting == 'all':
            for setting in db:
                skim(setting, db, script[args.part], args.part, args.force_update)
        else:
            if args.setting not in db:
                raise InvalidSettingError(args.setting)
            skim(args.setting, db, script[args.part], args.part, args.force_update)
