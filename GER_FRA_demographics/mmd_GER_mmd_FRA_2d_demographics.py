#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="IM_pre_FRA_GER_demographics"
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
    dadi.Plotting.plot_2d_comp_multinom(model, fs, show = False)
    plt.savefig('FRA_GER_snm.png')
def im_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.IM
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [0.1, 1, 1, 0.1, 0.1, 0.1, 0.1]
    lower = [1e-3, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 0]
    upper = [1, 3, 3, 1, 1, 1, 1]
    try:
        fid = open(f'demo_results/FRA_GER_im_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_im_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning im optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished im optimization ' + str(i) + '*'*20)
    fid.close()
def im_pre_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.IM_pre
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [1, 0.1, 0.1, 1, 1, 0.1, 3, 0.1, 0.1]
    lower = [1e-2, 1e-3, 1e-3, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 0]
    upper = [3, 1, 1, 3, 3, 1, 10, 1, 1]
    try:
        fid = open(f'demo_results/FRA_GER_im_pre_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_im_pre_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning im_pre optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished im_pre optimization ' + str(i) + '*'*20)
    fid.close()
def bottlegrowth_demography(fs,  ns, pts):
    demo_model = dadi.Demographics2D.bottlegrowth_2d
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [10, 1, 0.1, 0.1]
    lower = [1e-2, 1e-2, 1e-3, 0]
    upper = [100, 3, 1, 1]
    try:
        fid = open(f'demo_results/FRA_GER_bottlegrowth_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_bottlegrowth_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning bottlegrowth optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished bottlegrowth optimization ' + str(i) + '*'*20)
    fid.close()
def bottlegrowth_split_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.bottlegrowth_split
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [1, 1, 0.1, 0.1, 0.1]
    lower = [1e-2, 1e-2, 1e-3, 1e-3, 0]
    upper = [3, 3, 1, 1, 1]
    try:
        fid = open(f'demo_results/FRA_GER_bottlegrowth_split_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_bottlegrowth_split_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning bottlegrowth_split optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished bottlegrowth_split optimization ' + str(i) + '*'*20)
    fid.close()
def bottlegrowth_split_mig_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.bottlegrowth_split_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [1, 1, 0.1, 0.1, 0.1, 0.1]
    lower = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 0]
    upper = [3, 3, 10, 1, 1, 1]
    try:
        fid = open(f'demo_results/FRA_GER_bottlegrowth_split_mig_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_bottlegrwoth_split_mig_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning bottlegrowth_split_mig optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished bottlegrowth_split_mig optimization ' + str(i) + '*'*20)
    fid.close()
def split_asym_mig_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.split_asym_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [1, 1, 0.1, 0.1, 0.1, 0.1]
    lower = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 0]
    upper = [3, 3, 1, 10, 10, 1]
    try:
        fid = open(f'demo_results/FRA_GER_split_asym_mig_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_split_asym_mig_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning split_asym_mig optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished split_asym_mig optimization ' + str(i) + '*'*20)
    fid.close()
def split_delay_mig_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.split_delay_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [1, 1, 0.1, 0.1, 0.1, 0.1, 0.1]
    lower = [1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 0]
    upper = [3, 3, 1, 1, 1, 1, 1]
    try:
        fid = open(f'demo_results/FRA_GER_split_delay_mig_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_split_delay_mig_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning split_delay_mig optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished split_delay_mig optimization ' + str(i) + '*'*20)
    fid.close()
def split_mig_demography(fs, ns, pts):
    demo_model = dadi.Demographics2D.split_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    params = [1, 1, 0.1, 0.1, 0.1]
    lower = [1e-2, 1e-2, 1e-3, 1e-3, 0]
    upper = [3, 3, 1, 10, 1]
    try:
        fid = open(f'demo_results/FRA_GER_split_mig_demo_fits{process_ii}.txt', 'a')
    except:
        fid = open(f'demo_results/FRA_GER_split_mig_demo_fits{process_ii}.txt', 'w')
    for i in range(20):
        p0 = dadi.Misc.perturb_params(params, fold = 0, upper_bound = upper, lower_bound = lower)
        print('Beginning split_mig optimization ' + str(i) + '*'*20)
        popt, ll_model = dadi.Inference.opt(p0, fs, demo_model_ex, pts, upper_bound = upper, lower_bound = lower)
        model_fs = demo_model_ex(popt, ns, pts)
        theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
        res = [ll_model] + list(popt) + [theta0]
        fid.write('\t'.join([str(ele) for ele in res]) + '\n')
        print('Finished split_mig optimization ' + str(i) + '*'*20)
    fid.close()
def main():
    data_fs = dadi.Spectrum.from_file('FRA_GER_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

    # snm_demography(data_fs, ns, 100)
    # im_demography(data_fs, ns, pts_l)
    im_pre_demography(data_fs, ns, pts_l)
    # bottlegrowth_demography(data_fs, ns, pts_l)
    # bottlegrowth_split_demography(data_fs, ns, pts_l)
    # bottlegrowth_split_mig_demography(data_fs, ns, pts_l)
    # split_asym_mig_demography(data_fs, ns, pts_l)
    # split_delay_mig_demography(data_fs, ns, pts_l)
    # split_mig_demography(data_fs, ns, pts_l)
if __name__ == '__main__':
    main()