import pickle, random
import numpy as np
import nlopt
import dadi
import dadi.DFE as DFE
import matplotlib.pyplot as plt
from dadi.DFE import *

def get_params(fname):
    fp = open(fname, 'r')
    lines = fp.readlines()
    best_params = lines[0].strip('\n').split('\t')
    popt = best_params[1:-1]
    for i in range(len(popt)):
        popt[i] = float(popt[i])
    return popt
def plot_biv_lognormal_dfe(popt, dfe_func, sele_dist2d, theta_ns, data_fs, figname):
    model_fs = dfe_func(popt, None, sele_dist2d, theta_ns, None)
    # Plot
    fig = plt.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_Poisson(model_fs, data_fs,vmin=1)
    fig.savefig('plots/' + figname + '.png')
def main():
    #loading the cache, creating the extrapolation function, specifying the distribution, and loading our raw data frequency spectrum
    cache2d = pickle.load(open('mmd_HEL_mmd_GER_2d_cache.bpkl', 'rb'))
    dfe_func = cache2d.integrate
    dfe_func = dadi.Numerics.make_anc_state_misid_func(dfe_func)
    sele_dist2d = DFE.PDFs.biv_lognormal
    data_fs = dadi.Spectrum.from_file('GER_HEL_syn_unfolded.fs')
    #synonymous theta from IM_pre_inbreeding demography fit
    theta0 = 26681.214011656422
    theta_ns = theta0 * 2.4
    #getting the best fit parameters for the biv lognormmal with the same mu and sigma
    popt_same_mu_sigma = get_params('biv_lognormal_dfe_same_mu_sigma.txt')
    #plotting the DFE for IRA vs. FRA with the same mu and sigma and biv_lognormal distribution
    plot_biv_lognormal_dfe(popt_same_mu_sigma, dfe_func, sele_dist2d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_biv_lognormal_same_sigma_mu_DFE')
    #getting the best fit parameters for the biv lognormmal with the different mu and sigma
    popt_same_mu_sigma = get_params('biv_lognormal_dfe_different_mu_sigma.txt')
    #plotting the DFE for IRA vs. FRA with different mu and sigma and biv_lognormal distribution
    plot_biv_lognormal_dfe(popt_same_mu_sigma, dfe_func, sele_dist2d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_biv_lognormal_different_sigma_mu_DFE')
if __name__ == '__main__':
    main()