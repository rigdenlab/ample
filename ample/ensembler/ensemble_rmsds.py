#!/usr/bin/env ccp4-python
import os
import shutil
import numpy as np
import iotbx.pdb
from pyjob import Job

def ensemble_rmsds(ensemble, gesamt_exe):
    """Quick attempt to add PHASER RMSD commands to ensemble
    Currently just for testing"""
    
    def parse_rmsd_pair(fh):
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if 'RMSD             :' in line:
                rmsd = float(line.split()[2])
                break
        assert rmsd
        return [rmsd / 2.0] * 2

    def parse_rmsd_multi(fh):
        mlist = []
        reading = 0
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if '____' in line and reading == 2:
                break
            if '----' in line and reading == 1:
                reading = 2
                continue
            if '(o) pairwise r.m.s.d. (consensus r.m.s.d. on diagonal):' in line:
                reading = 1
                continue
            if reading == 2:
                mlist.append([float(r) for r in line.split()[1:]])
        return list(np.array(mlist).diagonal())

    nmodels = iotbx.pdb.pdb_input(ensemble).construct_hierarchy().models_size()
    nmodels = 2
     
    cmd = '{} '.format(gesamt_exe)
    for i in range(nmodels):
        cmd += '{} -s "/{}/" '.format(ensemble, i+1)
     
    script = './gesamt.sh'
    logfile = 'gesamt.log'
    with open(script, 'w') as w:
        w.write('#!/bin/bash\n{}\n\n'.format(cmd))
    os.chmod(script, 0o777)
     
    j = Job('local')
    j.submit(script, log=logfile)
    j.wait(interval=1)
    
    assert nmodels > 1
    with open(logfile) as fh:
        if nmodels > 2:
            rmsds = parse_rmsd_multi(fh)
        else:
            rmsds = parse_rmsd_pair(fh)
    os.unlink(script)
    os.unlink(logfile)
    return rmsds

def add_phaser_rmsds(ensemble, gesamt_exe):
    rmsds = ensemble_rmsds(ensemble, gesamt_exe)
    with open(ensemble) as fh:
        lines = fh.readlines()
    ensemble_bak = ensemble + '.bak'
    shutil.move(ensemble, ensemble_bak)
    plines = ['REMARK PHASER ENSEMBLE MODEL {} RMS {:.4}\n'.format(i, rmsds[i]) for i in range(len(rmsds))]
    lines = plines + lines
    with open(ensemble, 'w') as fh:
        fh.writelines(lines)
    return
