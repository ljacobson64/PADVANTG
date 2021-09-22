#!/usr/bin/python2

import os
import sys

def run_advantg():
    advantg_dir = '/home/lucas/PADVANTG/ADVANTG-3.0.3/lib/python2.7/site-packages'
    exnihilo_dir = '/opt/software_native/exnihilo-5.4.0/python'
    sys.path = [x for x in sys.path if 'python3' not in x]
    sys.path.insert(0, advantg_dir)
    os.environ['PYTHONPATH'] = exnihilo_dir
    from advantg.__main__ import main as advantg_main

    #os.system('rm -rf adj_solution fwd_solution model output custom_output')
    os.system('rm -rf output*')

    advantg_main(['-f', '-v', 'advantg.inp'])

if __name__ == '__main__': run_advantg()
