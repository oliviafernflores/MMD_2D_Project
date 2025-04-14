#!/usr/bin/env python
#SBATCH --job-name=DFE_size_testing
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=20
#SBATCH --ntasks=100
#SBATCH --time=48:00:00

import os
import sys
import pickle
import numpy as np
import dadi

# Ensure sample size is passed
if len(sys.argv) != 2:
    print("Usage: sbatch script_name.py <sample_size>")
    sys.exit(1)

sample_size = sys.argv[1]  # Keep as string for folder names

# Constants
theta = 9168
theta_ns = theta * 2.4

# Cache loading
cache_2d = pickle.load(open('/groups/rgutenk/oliviafernflores/IRA_FRA/mmd_FRA_mmd_IRA_2d_cache.bpkl', 'rb'))
cache_1d = pickle.load(open('/groups/rgutenk/oliviafernflores/IRA_FRA/mmd_FRA_mmd_IRA_1d_cache.bpkl', 'rb'))

# 2D model setup
dfe_func_2d = dadi.Numerics.make_anc_state_misid_func(cache_2d.integrate)
sele_dist_2d = dadi.DFE.PDFs.biv_lognormal_py
func_args_2d = [sele_dist_2d, theta_ns]
params_2d = [1, 1, 0.5, 0.5]
lower_2d = [-10, 0.1, 0.001, 0.001]
upper_2d = [100, 100, 0.999, 0.5]

# 1D model setup
dfe_func_1d = dadi.Numerics.make_anc_state_misid_func(cache_1d.integrate)
sele_dist_1d = dadi.DFE.PDFs.lognormal
func_args_1d = [sele_dist_1d, theta_ns]
params_1d = [1, 1, 0.5]
lower_1d = [-10, 0.1, 0]
upper_1d = [100, 100, 0.5]

start_dir = '/groups/rgutenk/oliviafernflores/DFE_size_testing'
results_dir = '/groups/rgutenk/oliviafernflores/DFE_size_testing/Mmd_IRA_Mmd_FRA_results'

for i in range(10):  # 10 replicates
    sfs_path = os.path.join(
        start_dir,
        f'{sample_size}_random_samples_SFS',
        'Mmd_FRA_Mmd_IRA_2d',
        f'{sample_size}_random_samples_Mmd_FRA_Mmd_IRA_2D_random_SFS_{i}_nsyn_unfolded.fs'
    )
    sfs = dadi.Spectrum.from_file(sfs_path)

    # Run optimization 10 times
    for j in range(10):
        # Output paths
        outdir_2d = os.path.join(results_dir, 'bivariate_lognormal', sample_size)
        outdir_1d = os.path.join(results_dir, 'lognormal', sample_size)
        os.makedirs(outdir_2d, exist_ok=True)
        os.makedirs(outdir_1d, exist_ok=True)

        outfile_2d = os.path.join(outdir_2d, f'{sample_size}_random_samples_Mmd_IRA_Mmd_FRA_random_SFS_{i}_nsyn_unfolded_bivariate_lognormal_DFE_inference_{j}.txt')
        outfile_1d = os.path.join(outdir_1d, f'{sample_size}_random_samples_Mmd_IRA_Mmd_FRA_random_SFS_{i}_nsyn_unfolded_lognormal_DFE_inference_{j}.txt')

        # Bivariate lognormal
        p0_2d = dadi.Misc.perturb_params(params_2d, fold=1, upper_bound=upper_2d, lower_bound=lower_2d)
        popt_2d, ll_model_2d = dadi.Inference.opt(p0_2d, sfs, dfe_func_2d, pts=None, func_args=func_args_2d, lower_bound=lower_2d, upper_bound=upper_2d, maxeval=600, multinom=False, verbose=0)
        with open(outfile_2d, 'a') as f:
            res_2d = [ll_model_2d] + list(popt_2d) + [theta_ns]
            f.write('\t'.join(map(str, res_2d)) + '\n')

        # Lognormal
        p0_1d = dadi.Misc.perturb_params(params_1d, fold=1, upper_bound=upper_1d, lower_bound=lower_1d)
        popt_1d, ll_model_1d = dadi.Inference.opt(p0_1d, sfs, dfe_func_1d, pts=None, func_args=func_args_1d, lower_bound=lower_1d, upper_bound=upper_1d, maxeval=600, multinom=False, verbose=0)
        with open(outfile_1d, 'a') as f:
            res_1d = [ll_model_1d] + list(popt_1d) + [theta_ns]
            f.write('\t'.join(map(str, res_1d)) + '\n')

    print(f'Finished DFE inferences for replicate {i} at sample size {sample_size}')

print('*' * 30)
print(f'Finished all DFE inferences for sample size {sample_size}')
