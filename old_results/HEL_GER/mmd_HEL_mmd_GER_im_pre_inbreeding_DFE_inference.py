import pickle, random
import numpy as np
import nlopt
import dadi
import dadi.DFE as DFE
import matplotlib.pyplot as plt
from dadi.DFE import *

def biv_lognormal_same_mu_sigma(theta_ns, data_fs, dfe_func):
    #define the selection distribution
    sele_dist2d = DFE.PDFs.biv_lognormal
    func_args = [sele_dist2d, theta_ns]
    #starting parameters for DFE inference
    params = [0.2, 0.2, 0.9, 0.1]
    #bounds for DFE inference
    lower_bounds = [-10, 0.01, 1e-1, 1e-2]
    upper_bounds = [10, 100, 0.99, 0.99]

    for x in range(1):
        p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)
        popt, ll_model = dadi.Inference.opt(p0, data_fs, dfe_func, pts=None, func_args=func_args, lower_bound=lower_bounds, upper_bound=upper_bounds, maxeval=6000, multinom=False, verbose=0)
        #get DFE spectrum
        model_fs = dfe_func(popt, None, sele_dist2d, theta_ns, None)
        fid = open('biv_lognormal_dfe_same_mu_sigma.txt','a')
        #write results to a file
        res = [ll_model] + list(popt) + [theta_ns]
        fid.write('\t'.join([str(ele) for ele in res])+'\n')
        fid.close()
        #print the iteration
        print('Finished optimization: {0}'.format(x+1))
def biv_lognormal_independent_mu_sigma(theta_ns, data_fs, dfe_func):
    #define the selection distribution
    sele_dist2d = DFE.PDFs.biv_lognormal
    func_args = [sele_dist2d, theta_ns]
    #starting parameters for DFE inference
    params = [0.2, 0.2, 0.2, 0.2, 0.9, 0.1]
    #bounds for DFE inference
    lower_bounds = [-10, -10, 0.01, 0.01, 1e-1, 1e-2]
    upper_bounds = [10, 10, 100, 100, 0.99, 0.99]

    for x in range(1):
        p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)
        popt, ll_model = dadi.Inference.opt(p0, data_fs, dfe_func, pts=None, func_args=func_args, lower_bound=lower_bounds, upper_bound=upper_bounds, maxeval=6000, multinom=False, verbose=0)
        #get DFE spectrum
        model_fs = dfe_func(popt, None, sele_dist2d, theta_ns, None)
        fid = open('biv_lognormal_dfe_different_mu_sigma.txt','a')
        #write results to a file
        res = [ll_model] + list(popt) + [theta_ns]
        fid.write('\t'.join([str(ele) for ele in res])+'\n')
        fid.close()
        #print the iteration
        print('Finished optimization: {0}'.format(x+1))

def main():
    data_fs = dadi.Spectrum.from_file('GER_HEL_syn_unfolded.fs')
    ns = data_fs.sample_sizes

    #synonymous theta from IM_pre_inbreeding_pop1_only demography fit
    theta0 = 22685.18157942911
    theta_ns = theta0 * 2.4

    pts_l = [max(ns)+140, max(ns)+150, max(ns)+160]

    demo_sel_model = DFE.DemogSelModels.IM_pre
    #params from IM_pre_inbreeding demography fit, no misid and no inbreeding coefficients
    demo_popt = [0.3060088614080329, 0.05367337350992898, 0.14744426867687283, 0.03573416869342573, 0.6340741046187689, 0.045314335100249334, 0.11760025775130044, 0.19431372592724705]
    # DFE.Cache2D
    # cache2d = DFE.Cache2D(popt, ns, demo_sel_model, pts=pts_l, gamma_bounds=[1e-4, 2000], gamma_pts=50)

    cache2d = pickle.load(open('mmd_HEL_mmd_GER_2d_cache.bpkl', 'rb'))

    #check if there's large negative values in the cache
    if (cache2d.spectra<0).sum() > 0:
        print(
            '!!!WARNING!!!\nPotentially large negative values!\nMost negative value is: '+str(cache2d.spectra.min())+
            '\nIf negative values are very negative (<-0.001), rerun with larger values for pts_l'
            )

    #define the DFE function
    dfe_func = cache2d.integrate
    dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)


    biv_lognormal_same_mu_sigma(theta_ns, data_fs, dfe_func)
    biv_lognormal_independent_mu_sigma(theta_ns, data_fs, dfe_func)
if __name__ == '__main__':
    main()