#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="IRA_FRA_mixture_DFE"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1-5
import pickle, random
import numpy as np
import nlopt
import dadi
import dadi.DFE as DFE
import matplotlib.pyplot as plt
from dadi.DFE import *



data_fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')
ns = data_fs.sample_sizes
#synonymous theta from IM_pre_inbreeding demography fit
theta0 = 26681.214011656422
theta_ns = theta0 * 2.4
#load the 1D cache
cache1d = pickle.load(open('mmd_FRA_mmd_IRA_1d_cache.bpkl', 'rb'))
#load the 2D cache
cache2d = pickle.load(open('mmd_FRA_mmd_IRA_2d_cache.bpkl', 'rb'))
#define the DFE function
dfe_func = DFE.mixture
#add the misidentification
dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)
#define the selection distributions
sele_dist1d = DFE.PDFs.lognormal
sele_dist2d = DFE.PDFs.biv_lognormal
#arguments for optimization
func_args = [cache1d, cache2d, sele_dist1d, sele_dist2d, theta_ns]
#starting parameters for DFE inference
params = [2, 2, 0, 0.01, 0.01]
#bounds and fixed parameters for DFE inference
lower_bounds = [1e-2, 1e-2, None, 1e-3, 1e-3]
upper_bounds = [100, 10, None, 1-1e-3, 1]
fixed_params = [None, None, 0, None, None]
#checking if the file exists and opening it
try:
    fid = open('mmd_FRA_mmd_IRA_mixture_DFE.txt','a')
except:
    fid = open('mmd_FRA_mmd_IRA_mixture_DFE.txt','w')
#starting a for loop that will be used to run a bunch of optimizations
for i in range(200):
    #perturb starting pararmeters - this will give you a different starting point each time
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)
    #run the optimization
    popt, ll_model = dadi.Inference.opt(p0, data_fs, dfe_func, pts=None, func_args=func_args, lower_bound=lower_bounds, upper_bound=upper_bounds, fixed_params=fixed_params, maxeval=1000, multinom=False, verbose=30)
    #writing the results to the file
    res = [ll_model] + list(popt) + [theta_ns]
    fid.write('\t'.join([str(ele) for ele in res])+'\n')
fid.close()