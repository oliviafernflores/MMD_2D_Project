#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="all_mouse_pops_syn_fs"
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

import sys, os
print('Script running\n')
if 'SLURM_SUBMIT_DIR' in os.environ:
    sys.path.insert(0, os.environ['SLURM_SUBMIT_DIR'])
process_ii = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1)) - 1
print(process_ii)

def snm_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.snm_2d
    model = demo_model([], ns, pts)
    plt.figure(figsize = (8, 6))
    dadi.Plotting.plot_2d_comp_multinom(model, fs, vmin = 1, resid_range = 3, pop_ids = ['Mmd_IRA', 'Mmd_FRA'], show = False)
    plt.savefig('IRA_FRA_snm.png')
def im_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.IM
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [0.1, 1, 1, 0.1, 0.1, 0.1, 0.1]
    lower = [1e-3, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 0]
    upper = [1, 3, 3, 1, 1, 1, 1]
    try:
        fid = open(f'demo_results/IRA_FRA_im_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/IRA_FRA_im_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished optimization ' + str(i) + '*'*20)
    fid.close()

def main():
    data_fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

    snm_demography(data_fs, ns, 100)
    im_demography(data_fs, ns, pts_l)
if __name__ == '__main__':
    main()