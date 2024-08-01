#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="mmd_test_fs"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt

dd = pickle.load(open('mus_domesticus_all_snps_all_pops_no_relatives_dd.bpkl', 'rb'))
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [16, 10]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns, polarized = False)
fs.to_file('IRA_FRA_test2.fs')