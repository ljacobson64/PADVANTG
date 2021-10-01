#!/usr/bin/python3

import os

def main():
    # Get list of tallies
    with open('tally_list.txt', 'r') as r:
        tally_list = [int(x) for x in r.read().split()]
    direc_list = ['forward'] + ['adjoint-%u' % (x) for x in tally_list]

    # Read template input file
    with open('advantg_template.txt', 'r') as r: template = str(r.read())

    # Run ADVANTG in each directory
    for direc in direc_list:
        os.system('rm -rf ' + direc)
        os.system('mkdir -pv ' + direc)
        os.chdir(direc)
        if direc.startswith('forward'): vals = ['true', 'false', '']
        else:
            tally = int(direc.split('-')[1])
            vals = ['false', 'true', 'mcnp_tallies%s%u' % (' ' * 16, tally)]
        with open('advantg.inp', 'w') as w: w.write(template.format(*vals))
        os.system('cp -pv ../mcnp.i .')
        os.system('bash ../../run_advantg_once.sh')
        os.chdir('..')

if __name__ == '__main__': main()
