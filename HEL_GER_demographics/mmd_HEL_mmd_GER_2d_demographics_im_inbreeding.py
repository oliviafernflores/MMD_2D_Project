#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="FRA_GER_IM_inbreeding"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
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

def IM(params, ns, pts):
    """
    ns = (n1,n2)
    params = (s,nu1,nu2,T,m12,m21, F)

    Isolation-with-migration model with exponential pop growth.

    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    s,nu1,nu2,T,m12,m21, F1, F2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1, F2), (2, 2))
    return fs

def im_demography(fs, ns, pts):
    demo_model = IM
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [0.1, 1, 1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    lower = [1e-3, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 0]
    upper = [1, 3, 10, 1, 3, 1, 1, 1, 1]
    try:
        fid = open(f'demo_results/GER_HEL_im_inbreeding_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/GER_HEL_im_inbreeding_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning im_inbreeding optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished im_inbreeding optimization ' + str(i) + '*'*20)
    fid.close()
def main():
    data_fs = dadi.Spectrum.from_file('GER_HEL_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

    im_demography(data_fs, ns, pts_l)
if __name__ == '__main__':
    main()