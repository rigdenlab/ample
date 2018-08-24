#!/usr/bin/env ccp4-python
import os
import shutil
import numpy as np
import iotbx.pdb
from pyjob import Job


def parse_rmsd_sheaf(fh, num_models):
    reading = -1
    model_idx = 0
    rmsd_matrix = np.zeros([num_models, num_models])
    for line in fh:
        if line.startswith(' ===== CROSS-RMSDs') or reading == 0:
            # find start of RMSDS and skip blank line
            reading += 1
            continue
        if reading == 1:
            fields = line.strip().split('|')
            model_idx = int(fields[0])
            rmsd_txt = fields[2].strip()
            # poke into distance matrix
            rmsds = [ float(r) for r in rmsd_txt.split() ]
            for j in range(len(rmsds)):
                if j == model_idx - 1:
                    continue
                rmsd_matrix[model_idx-1][j] = rmsds[j]
            if model_idx == num_models:
                reading = -1

    # Get consensus from rmsds
    return [np.sqrt(np.sum(row**2 / len(row))) for row in rmsd_matrix]


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


def ensemble_rmsds(ensemble, gesamt_exe, mode='sheaf'):
    """Quick attempt to add PHASER RMSD commands to ensemble
    Currently just for testing"""

    num_models = iotbx.pdb.pdb_input(ensemble).construct_hierarchy().models_size()
     
    cmd = '{} '.format(gesamt_exe)
    for i in range(num_models):
        cmd += '{} -s "/{}/" '.format(ensemble, i+1)
    if mode == 'sheaf':
        cmd += ' -sheaf-x'
     
    script = './gesamt.sh'
    logfile = 'gesamt.log'
    with open(script, 'w') as w:
        w.write('#!/bin/bash\n{}\n\n'.format(cmd))
    os.chmod(script, 0o777)
     
    j = Job('local')
    j.submit(script, log=logfile)
    j.wait(interval=1)
    
    assert num_models > 1
    with open(logfile) as fh:
        if mode == 'sheaf':
            rmsds = parse_rmsd_sheaf(fh, num_models)
        else:
            if num_models > 2:
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
    os.unlink(ensemble_bak)
    return

if __name__ == '__main__':
    add_phaser_rmsds('/Users/jmht/Desktop/c3_t3_r3_polyala.pdb', 'gesamt')
