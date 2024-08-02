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
    for line in lines:
        log_likelihood = str(line.strip('\n').split('\t')[0].strip('-'))
        if float(log_likelihood) < float(best_params[0].strip('-')):
            best_params = line.strip('\n').split('\t')
    popt = best_params[1:-1]
    optimal_log_likelihood = best_params[0]
    for i in range(len(popt)):
        popt[i] = float(popt[i])
    return popt, optimal_log_likelihood
def plot_lognormal_dfe(popt, dfe_func, sele_dist1d, theta_ns, data_fs, figname):
    model_fs = dfe_func(popt, None, sele_dist1d, theta_ns, None)
    # Plot
    fig = plt.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_Poisson(model_fs, data_fs, vmin=1)
    fig.savefig('plots/' + figname + '.png')
def plot_biv_lognormal_dfe(popt, dfe_func, sele_dist2d, theta_ns, data_fs, figname):
    model_fs = dfe_func(popt, None, sele_dist2d, theta_ns, None)
    # Plot
    fig = plt.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_Poisson(model_fs, data_fs, vmin=1)
    fig.savefig('plots/' + figname + '.png')
def plot_mixture_simple_dfe(popt, s1, s2, sele_dist1d, sele_dist2d, theta_ns, data_fs, figname):
    dfe_func = dadi.Numerics.make_anc_state_misid_func(DFE.mixture)
    #s1 = cache1d
    #s2 = cache2d
    # [cache1d, cache2d, sele_dist1d, sele_dist2d, theta_ns]
    model_fs = dfe_func(popt,None,s1, s2, sele_dist1d, sele_dist2d, theta_ns, pts = None)
    fig = plt.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_Poisson(model_fs, data_fs, vmin=1)
    fig.savefig('plots/' + figname + '.png')
def plot_gamma_dfe(popt, dfe_func, sele_dist1d, theta_ns, data_fs, figname):
    model_fs = dfe_func(popt, None, sele_dist1d, theta_ns, None)
    # Plot
    fig = plt.figure(1)
    fig.clear()
    dadi.Plotting.plot_2d_comp_Poisson(model_fs, data_fs, vmin=1)
    fig.savefig('plots/' + figname + '.png')
def main():
    #loading the cache, creating the extrapolation function, specifying the distribution, and loading our raw data frequency spectrum
    cache2d = pickle.load(open('mmd_HEL_mmd_GER_2d_cache.bpkl', 'rb'))
    dfe_func2d = cache2d.integrate
    dfe_func2d = dadi.Numerics.make_anc_state_misid_func(dfe_func2d)
    sele_dist2d = DFE.PDFs.biv_lognormal
    
    cache1d = pickle.load(open('mmd_HEL_mmd_GER_1d_cache.bpkl', 'rb'))
    dfe_func1d = cache1d.integrate
    dfe_func1d = dadi.Numerics.make_anc_state_misid_func(dfe_func1d)
    sele_dist1d = DFE.PDFs.lognormal
    
    data_fs = dadi.Spectrum.from_file('GER_HEL_nsyn_unfolded.fs')
    #synonymous theta from IM_pre_inbreeding demography fit
    theta0 = 21074.66312426641
    theta_ns = theta0 * 2.4
    
    popt_1d_cli, ll_1d_cli  = [5.702287980423937, 6.375569692079467, 0.11959079544441137], -950.3913048832605
    plot_lognormal_dfe(popt_1d_cli, dfe_func1d, sele_dist1d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_1D_lognormal_DFE_dadi_cli')
    print('\n')
    print('1d bivariate lognormal best fit params: ' + str(popt_1d_cli) + ' with log likelihood: ' + str(ll_1d_cli))
    print('\n')

    popt_2d_cli, ll_2d_cli  = [5.756526275212594, 6.3339827034054945, 0.9824880561428179, 0.11966465449191582], -947.93248741037
    plot_lognormal_dfe(popt_2d_cli, dfe_func2d, sele_dist2d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_2D_biv_lognormal_sym_DFE_dadi_cli')
    print('\n')
    print('2d bivariate symmetric lognormal best fit params: ' + str(popt_2d_cli) + ' with log likelihood: ' + str(ll_2d_cli))
    print('\n')
    
    popt_2d_cli_asym, ll_2d_cli_asym  = [8.91928220094154, 0.9959701906513141, 10.0, 3.6055970369426085, 0.5695847682533904, 0.11907311913341892], -756.4094618130488
    plot_lognormal_dfe(popt_2d_cli_asym, dfe_func2d, sele_dist2d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_2D_biv_lognormal_asym_DFE_dadi_cli')
    print('\n')
    print('2d bivariate asymmetric lognormal best fit params: ' + str(popt_2d_cli_asym) + ' with log likelihood: ' + str(ll_2d_cli_asym))
    print('\n')
    
    # popt_1d, ll_1d  = get_params('mmd_HEL_mmd_GER_1D_lognormal_DFE.txt')
    # plot_lognormal_dfe(popt_1d, dfe_func1d, sele_dist1d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_1D_lognormal_DFE')
    # print('\n')
    # print('1d lognormal best fit params: ' + str(popt_1d) + ' with log likelihood: ' + str(ll_1d))
    # print('\n')

    # popt_2d, ll_2d = get_params('mmd_HEL_mmd_GER_2D_biv_lognormal_DFE.txt')
    # plot_biv_lognormal_dfe(popt_2d, dfe_func2d, sele_dist2d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_2D_biv_lognormal_DFE')
    # print('\n')
    # print('2d bivariate lognormal best fit params: ' + str(popt_2d)+ ' with log likelihood: ' + str(ll_2d))
    # print('\n')

    # popt_mixture, ll_mix = get_params('mmd_HEL_mmd_GER_mixture_DFE.txt')
    # plot_mixture_simple_dfe(popt_mixture, cache1d, cache2d, sele_dist1d, sele_dist2d, theta_ns, data_fs, 'mmd_HEL_mmd_GER_simple_mixture_DFE')
    # print('\n')
    # print('mixture with shared selection best fit params: ' + str(popt_mixture) + ' with log likelihood: ' + str(ll_mix))
    # print('\n')
    
    # popt_gamma, ll_gamma = get_params('mmd_HEL_mmd_GER_gamma_DFE.txt')
    # gamma_dist = DFE.PDFs.gamma
    # plot_gamma_dfe(popt_gamma, dfe_func1d, gamma_dist, theta_ns, data_fs, 'mmd_HEL_mmd_GER_gamma_DFE')
    # print('\n')
    # print('gamma with shared selection best fit params: ' + str(popt_gamma) + ' with log likelihood: ' + str(ll_gamma))
    # print('\n')
if __name__ == '__main__':
    main()