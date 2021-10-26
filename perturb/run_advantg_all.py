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
        #os.system('rm -rf ' + direc)
        os.system('mkdir -pv ' + direc)
        os.chdir(direc)

        # Write ADVANTG input file
        if direc.startswith('forward'): vals = ['true', 'false', '']
        else:
            tally = int(direc.split('-')[1])
            vals = ['false', 'true', 'mcnp_tallies%s%u' % (' ' * 16, tally)]
        with open('advantg.inp', 'w') as w: w.write(template.format(*vals))

        if direc.startswith('forward'):
            # Write MCNP input file
            with open('../mcnp.i', 'r') as r: lines = r.readlines()
            s = ''
            for line in lines:
                tokens = line.split()
                if len(tokens) == 0:
                    s += line
                    continue
                if tokens[0].lower() == 'read':
                    with open('../%s' % (tokens[3]), 'r') as r: s += r.read()
                    continue
                s += line
            with open('mcnp.i', 'w') as w: w.write(s)
        else:
            # Copy MCNP input and output from forward/ directory
            os.system('ln -snfv ../forward/mcnp.i .')
            os.system('ln -snfv ../forward/model  .')

        # Run ADVANTG
        os.system('bash ../../perturb/run_advantg_once.sh')
        os.chdir('..')

if __name__ == '__main__': main()
