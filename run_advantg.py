#!/usr/bin/python2

import os
import sys

def run_advantg():
    advantg_dir = ('/home/lucas/git/PADVANTG/build/ADVANTG-3.0.3/'
                   'lib/python2.7/site-packages')
    exnihilo_dir = '/home/lucas/git/PADVANTG/build/exnihilo-5.4.0/python'
    sys.path.insert(0, advantg_dir)
    os.environ['PYTHONPATH'] = exnihilo_dir
    from advantg.__main__ import main as advantg_main

    #os.system('rm -rf adj_solution fwd_solution model output text_files')
    os.system('rm -rf output*')

    advantg_main(['-f', '-v', 'advantg.inp'])

if __name__ == '__main__': run_advantg()
