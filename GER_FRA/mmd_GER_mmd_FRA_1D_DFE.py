#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="FRA_GER_1d_DFE"
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



data_fs = dadi.Spectrum.from_file('FRA_GER_syn_unfolded.fs')
ns = data_fs.sample_sizes
#synonymous theta from IM_pre_inbreeding demography fit
theta0 = 31366.855164380948
theta_ns = theta0 * 2.4
#load the 1D cache
cache1d = pickle.load(open('mmd_GER_mmd_FRA_1d_cache.bpkl', 'rb'))
#check if there's large negative values in the cache
if (cache1d.spectra<0).sum() > 0:
    print(
        '!!!WARNING!!!\nPotentially large negative values!\nMost negative value is: '+str(cache1d.spectra.min())+
        '\nIf negative values are very negative (<-0.001), rerun with larger values for pts_l'
        )
#define the DFE function
dfe_func = cache1d.integrate
#add the misidentification
dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)
#define the selection distribution
sele_dist1d = DFE.PDFs.lognormal
#arguments for optimization
func_args = [sele_dist1d, theta_ns]
#starting parameters for DFE inference
params = [0.2, 0.2, 0.1]
#bounds for DFE inference
lower_bounds = [1e-2, 1e-2, 1e-3]
upper_bounds = [100, 100, 0.99]
#checking if the file exists and opening it
try:
    fid = open('mmd_GER_mmd_FRA_1D_lognormal_DFE.txt','a')
except:
    fid = open('mmd_GER_mmd_FRA_1D_lognormal_DFE.txt','w')
#starting a for loop that will be used to run a bunch of optimizations
for i in range(200):
    #perturb starting pararmeters - this will give you a different starting point each time
    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)
    #run the optimization
    popt, ll_model = dadi.Inference.opt(p0, data_fs, dfe_func, pts=None, func_args=func_args, lower_bound=lower_bounds, upper_bound=upper_bounds, maxeval=1000, multinom=False, verbose=100)
    #writing the results to the file
    res = [ll_model] + list(popt) + [theta_ns]
    fid.write('\t'.join([str(ele) for ele in res])+'\n')
fid.close()