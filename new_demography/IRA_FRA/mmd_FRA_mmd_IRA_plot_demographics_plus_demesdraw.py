import dadi
import matplotlib.pyplot as plt
import numpy as np
import csv
import matplotlib.pylab
import demes
import demesdraw

def get_popt(fname):
    popt = []
    fits_file = open(fname, 'r')
    fits = fits_file.readlines()
    popt = fits[0]
    popt = popt.split('\t')
    for f in fits:
        if f.split('\t')[0] < popt[0]:
            popt = f.split('\t')
    for i in range(len(popt)):
        popt[i] = float(popt[i])
    return popt

def plot_im_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (im): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.IM
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_IM_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_im_demography.png')
    out = dadi.Demes.output(deme_mapping = {'ancestral':['d1_1'], 'FRANCE': ['d1_2'], 'IRAN': ['d2_2']})
    demes.dump(out, 'demo_results/IRA_FRA_im_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_im_demography_demesdraw.png')
def plot_im_pre_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (im_pre): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.IM_pre
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_IM_pre_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_im_pre_demography.png')
    out = dadi.Demes.output(deme_mapping = {'ancestral':['d1_1'], 'FRANCE': ['d1_2'], 'IRAN': ['d2_2']})
    demes.dump(out, 'demo_results/IRA_FRA_im_pre_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_im_pre_demography_demesdraw.png')
def plot_bottlegrowth_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (bottlegrowth): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.bottlegrowth_2d
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_Bottlegrowth_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_bottlegrowth_demography.png')
    out = dadi.Demes.output()
    demes.dump(out, 'demo_results/IRA_FRA_bottlegrowth_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_bottlegrowth_demography_demesdraw.png')
def plot_bottlegrowth_split_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (bottlegrowth_split): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.bottlegrowth_split
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_Bottlegrowth_Split_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_bottlegrowth_split_demography.png')
    out = dadi.Demes.output(deme_mapping = {'ancestral':['d1_1'], 'FRANCE': ['d1_2'], 'IRAN': ['d2_2']})
    demes.dump(out, 'demo_results/IRA_FRA_bottlegrowth_split_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_bottlegrowth_split_demography_demesdraw.png')
def plot_bottlegrowth_split_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (bottlegrowth_split_mig): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.bottlegrowth_split_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_Bottlegrowth_Split_Migration_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_bottlegrowth_split_mig_demography.png')
    out = dadi.Demes.output(deme_mapping = {'ancestral':['d1_1'], 'FRANCE': ['d1_2'], 'IRAN': ['d2_2']})
    demes.dump(out, 'demo_results/IRA_FRA_bottlegrowth_split_mig_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_bottlegrowth_split_mig_demography_demesdraw.png')
def plot_split_asym_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (split_asym_mig): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.split_asym_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_Split_Asymmetric_Migration_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_split_asym_mig_demography.png')
    out = dadi.Demes.output(deme_mapping = {'ancestral':['d1_1'], 'FRANCE': ['d1_2'], 'IRAN': ['d2_2']})
    demes.dump(out, 'demo_results/IRA_FRA_split_asym_mig_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_split_asym_mig_demography_demesdraw.png')
def plot_split_delay_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (split_delay_mig): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.split_delay_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_Split_Delay_Migration_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_split_delay_mig_demography.png')
    out = dadi.Demes.output(deme_mapping = {'ancestral':['d1_1'], 'FRANCE': ['d1_2'], 'IRAN': ['d2_2']})
    demes.dump(out, 'demo_results/IRA_FRA_split_delay_mig_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_split_delay_mig_demography_demesdraw.png')
def plot_split_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (split_mig): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = dadi.Demographics2D.split_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_Split_Migration_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/IRA_FRA_split_mig_demography.png')
    out = dadi.Demes.output(deme_mapping = {'ancestral':['d1_1'], 'FRANCE': ['d1_2'], 'IRAN': ['d2_2']})
    demes.dump(out, 'demo_results/IRA_FRA_split_mig_demo_fits_demes.yaml')
    ax = demesdraw.tubes(out)
    ax.figure.savefig('plots/IRA_FRA_split_mig_demography_demesdraw.png')
def main():
    fs = dadi.Spectrum.from_file('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/SFS/IRA_FRA_syn_unfolded.fs')
    ns = fs.sample_sizes
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

    # #im_demography
    # im_popt = get_popt('demo_results/IRA_FRA_im_demo_fits_combined.txt')
    # plot_im_demography(fs, ns, pts_l, im_popt)

    # #im_pre_demography
    # im_pre_popt = get_popt('demo_results/IRA_FRA_im_pre_demo_fits_combined.txt')
    # plot_im_pre_demography(fs, ns, pts_l, im_pre_popt)

    #bottlegrowth_demography
    # bottlegrowth_popt = get_popt('demo_results/IRA_FRA_bottlegrowth_demo_fits_combined.txt')
    # plot_bottlegrowth_demography(fs, ns, pts_l, bottlegrowth_popt)

    #bottlegrowth_split_demography
    # bottlegrowth_split_popt = get_popt('demo_results/IRA_FRA_bottlegrowth_split_demo_fits_combined.txt')
    # plot_bottlegrowth_split_demography(fs, ns, pts_l, bottlegrowth_split_popt)

    #bottlegrowth_split_mig_demography
    # bottlegrowth_split_mig_popt = get_popt('demo_results/IRA_FRA_bottlegrowth_split_mig_demo_fits_combined.txt')
    # plot_bottlegrowth_split_mig_demography(fs, ns, pts_l, bottlegrowth_split_mig_popt)

    # #split_asym_mig_demography
    # split_asym_mig_popt = get_popt('demo_results/IRA_FRA_split_asym_mig_demo_fits_combined.txt')
    # plot_split_asym_mig_demography(fs, ns, pts_l, split_asym_mig_popt)

    #split_delay_mig_demography
    # split_delay_mig_popt = get_popt('demo_results/IRA_FRA_split_delay_mig_demo_fits_combined.txt')
    # plot_split_delay_mig_demography(fs, ns, pts_l, split_delay_mig_popt)

    #split_mig_demography
    split_mig_popt = get_popt('demo_results/IRA_FRA_split_mig_demo_fits_combined.txt')
    plot_split_mig_demography(fs, ns, pts_l, split_mig_popt)
if __name__ == '__main__':
    main()