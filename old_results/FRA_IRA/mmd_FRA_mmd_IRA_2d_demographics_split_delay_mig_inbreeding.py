#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="IRA_FRA_split_delay_mig_inbreeding"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=10
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1-5
from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt
import numpy as np
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

import sys, os
print('Script running\n')
if 'SLURM_SUBMIT_DIR' in os.environ:
    sys.path.insert(0, os.environ['SLURM_SUBMIT_DIR'])
process_ii = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1)) - 1
print(process_ii)

def split_delay_mig(params, ns, pts):
    """
    params = (nu1,nu2,Tpre,Tmig,m12,m21)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration after some time has passed post split.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Tpre: Time in the past after split but before migration (in units of 2*Na generations) 
    Tmig: Time in the past after migration starts (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2 (2*Na*m21)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,Tpre,Tmig,m12,m21, F1, F2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Tpre, nu1, nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tmig, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1, F2), (2, 2))
    return fs

def split_delay_mig_demography(fs, ns, pts):
    demo_model = split_delay_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    lower = [1e-2, 1e-2, 1e-4, 1e-3, 1e-3, 1e-3, 1e-5, 1e-5, 0]
    upper = [3, 6, 1, 1, 3, 1, 1, 1, 1]
    try:
        fid = open(f'demo_results/IRA_FRA_split_delay_mig_inbreeding_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/IRA_FRA_split_delay_mig_inbreeding_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 1, upper_bound = upper, lower_bound = lower)
        print('Beginning split_delay_mig_inbreeding optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished split_delay_mig_inbreeding optimization ' + str(i) + '*'*20)
    fid.close()
def main():
    data_fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

    
    split_delay_mig_demography(data_fs, ns, pts_l)
if __name__ == '__main__':
    main()