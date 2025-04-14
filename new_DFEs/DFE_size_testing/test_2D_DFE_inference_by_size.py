#!/usr/bin/env python
#SBATCH --job-name=make_GO_SFS_1D
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=10
#SBATCH --ntasks=50
#SBATCH --time=24:00:00

import os
import pickle
import numpy as np
import dadi
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend

import matplotlib.pyplot as plt
import pandas as pd
import random


# From demography results summary - im pre with inbreeding 
theta = 9168
theta_ns = 9168 * 2.4

# Inference Set Up
cache_2d = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/FRA_IRA/mmd_FRA_mmd_IRA_2d_cache.bpkl', 'rb'))
cache_1d = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/FRA_IRA/mmd_FRA_mmd_IRA_1d_cache.bpkl', 'rb'))


# biv_lognormal symmetric
dfe_func_2d = cache_2d.integrate
dfe_func_2d = dadi.Numerics.make_anc_state_misid_func(dfe_func_2d)
sele_dist_2d = dadi.DFE.PDFs.biv_lognormal_py
func_args_2d = [sele_dist_2d, theta_ns]
params_2d = [1, 1, 0.5, 0.5]
lower_2d = [-10, 0.1, 0.001, 0.001]
upper_2d = [100, 100, 0.999, 0.5]

# lognormal
dfe_func_1d = cache_1d.integrate
dfe_func_1d = dadi.Numerics.make_anc_state_misid_func(dfe_func_1d)
sele_dist_1d = dadi.DFE.PDFs.lognormal
func_args_1d = [sele_dist_1d, theta_ns]
params_1d = [1, 1, 0.5]
lower_1d = [-10, 0.1, 0]
upper_1d = [100, 100, 0.5]

# pops = ['Mmd_IRA', 'Mmd_FRA', 'Mmd_GER', 'Mmd_HEL', ['Mmd_IRA', 'Mmd_FRA'], ['Mmd_FRA', 'Mmd_GER'], ['Mmd_GER', 'Mmd_HEL']]

sample_sizes = [50000, 25000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 500, 100]

start_dir = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/DFE_size_testing'

for size in sample_sizes:
    # iterates through the sample sizes
    for i in range(10):
        # iterates through the 10 replicates of random SFS for each sample size
        sfs_loc = os.path.join(start_dir, f'{str(size)}_random_samples_SFS')
        sfs_loc = os.path.join(sfs_loc, 'Mmd_FRA_Mmd_IRA_2d')
        sfs_loc = os.path.join(sfs_loc, f'{str(size)}_random_samples_Mmd_FRA_Mmd_IRA_2D_random_SFS_{str(i)}_nsyn_unfolded.fs')
        sfs = dadi.Spectrum.from_file(sfs_loc)
        ns = sfs.sample_sizes
        for j in range(10):
            # runs optimization 10 times for each random SFS
            # need to get the right file name so we know which belong to the same SFS and are just different optimization runs
            fid_2d_dir = f'/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/DFE_size_testing/Mmd_IRA_Mmd_FRA_results/bivariate_lognormal/{str(size)}_random_samples_Mmd_IRA_Mmd_FRA_random_SFS_{str(i)}_nsyn_unfolded_bivariate_lognormal_DFE_inference_{str(j)}.txt'
            fid_1d_dir = f'/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/DFE_size_testing/Mmd_IRA_Mmd_FRA_results/lognormal/{str(size)}_random_samples_Mmd_IRA_Mmd_FRA_random_SFS_{str(i)}_nsyn_unfolded_lognormal_DFE_inference_{str(j)}.txt'
            
            # biv_lognormal
            p0_2d = dadi.Misc.perturb_params(params_2d, fold=1, upper_bound=upper_2d, lower_bound=lower_2d)
            popt_2d, ll_model_2d = dadi.Inference.opt(p0_2d, sfs, dfe_func_2d, pts=None, func_args=func_args_2d, lower_bound=lower_2d, upper_bound=upper_2d,maxeval=600, multinom=False, verbose=0)
            try:
                fid_2d = open(fid_2d_dir, 'a')
            except:
                fid_2d = open(fid_2d_dir, 'w')
            res_2d = [ll_model_2d] + list(popt_2d) + [theta_ns]
            fid_2d.write('\t'.join([str(ele) for ele in res_2d]) + '\n')
            fid_2d.close()
            
            # lognormal
            p0_1d = dadi.Misc.perturb_params(params_1d, fold=1, upper_bound=upper_1d, lower_bound=lower_1d)
            popt_1d, ll_model_1d = dadi.Inference.opt(p0_1d, sfs, dfe_func_1d, pts=None, func_args=func_args_1d, lower_bound=lower_1d, upper_bound=upper_1d,maxeval=600, multinom=False, verbose=0)
            try:
                fid_1d = open(fid_1d_dir, 'a')
            except:
                fid_1d = open(fid_1d_dir, 'w')
            res_1d = [ll_model_1d] + list(popt_1d) + [theta_ns]
            fid_1d.write('\t'.join([str(ele) for ele in res_1d]) + '\n')
            fid_1d.close()
        print(f'Finished DFE inferences for {str(size)}_random_samples_Mmd_IRA_Mmd_FRA_random_SFS_{str(i)}')
    print('*'*20)
    print(f'Finished DFE inferences for all SFS from {str(size)}')