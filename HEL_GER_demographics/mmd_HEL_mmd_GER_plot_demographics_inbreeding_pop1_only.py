import dadi
import matplotlib.pyplot as plt
import numpy as np
import csv
import inbreeding_models_pop1_only as mods

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
    print('Best Fit Parameters (im + inbreeding_pop1_only): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = mods.IM
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('GER_HEL_IM_Inbreeding_pop1_only_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/GER_HEL_im_inbreeding_pop1_only_demography.png')
def plot_im_pre_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (im_pre + inbreeding_pop1_only): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = mods.IM_pre
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('GER_HEL_IM_pre_Inbreeding_pop1_only_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/GER_HEL_im_pre_inbreeding_pop1_only_demography.png')
# def plot_bottlegrowth_demography(fs, ns, pts, popt):
#     print('Best Fit Parameters (bottlegrowth + inbreeding_pop1_only): ' + str(popt))
#     popt.pop(0)
#     popt.pop(-1)
#     demo_model = mods.bottlegrowth_2d
#     demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
#     demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
#     model = demo_model_ex(popt, ns, pts)
#     fig = plt.figure('GER_HEL_Bottlegrowth_Inbreeding_pop1_only_Demography')
#     fig.clear()
#     dadi.Plotting.plot_2d_comp_multinom(model, fs)
#     fig.savefig('plots/GER_HEL_bottlegrowth_inbreeding_pop1_only_demography.png')
def plot_bottlegrowth_split_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (bottlegrowth_split + inbreeding_pop1_only): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = mods.bottlegrowth_split
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('GER_HEL_Bottlegrowth_Split_Inbreeding_pop1_only_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/GER_HEL_bottlegrowth_split_inbreeding_pop1_only_demography.png')
def plot_bottlegrowth_split_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (bottlegrowth_split_mig + inbreeding_pop1_only): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = mods.bottlegrowth_split_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('GER_HEL_Bottlegrowth_Split_Migration_Inbreeding_pop1_only_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/GER_HEL_bottlegrowth_split_mig_inbreeding_pop1_only_demography.png')
def plot_split_asym_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (split_asym_mig + inbreeding_pop1_only): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = mods.split_asym_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('GER_HEL_Split_Asymmetric_Migration_Inbreeding_pop1_only_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/GER_HEL_split_asym_mig_inbreeding_pop1_only_demography.png')
def plot_split_delay_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (split_delay_mig + inbreeding_pop1_only): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = mods.split_delay_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('GER_HEL_Split_Delay_Migration_Inbreeding_pop1_only_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/GER_HEL_split_delay_mig_inbreeding_pop1_only_demography.png')
def plot_split_mig_demography(fs, ns, pts, popt):
    print('Best Fit Parameters (split_mig + inbreeding_pop1_only): ' + str(popt))
    popt.pop(0)
    popt.pop(-1)
    demo_model = mods.split_mig
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('GER_HEL_Split_Migration_Inbreeding_pop1_only_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('plots/GER_HEL_split_mig_inbreeding_pop1_only_demography.png')
def main():
    fs = dadi.Spectrum.from_file('GER_HEL_syn_unfolded.fs')
    ns = fs.sample_sizes
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]

    #im_demography
    im_popt = get_popt('results/GER_HEL_im_inbreeding_pop1_only_demo_fits_combined.txt')
    plot_im_demography(fs, ns, pts_l, im_popt)

    #im_pre_demography
    im_pre_popt = get_popt('results/GER_HEL_im_pre_inbreeding_pop1_only_demo_fits_combined.txt')
    plot_im_pre_demography(fs, ns, pts_l, im_pre_popt)

    # #bottlegrowth_demography
    # bottlegrowth_popt = get_popt('results/')
    # plot_bottlegrowth_demography(fs, ns, pts_l, bottlegrowth_popt)

    #bottlegrowth_split_demography
    bottlegrowth_split_popt = get_popt('results/GER_HEL_bottlegrowth_split_inbreeding_pop1_only_demo_fits_combined.txt')
    plot_bottlegrowth_split_demography(fs, ns, pts_l, bottlegrowth_split_popt)

    #bottlegrowth_split_mig_demography
    bottlegrowth_split_mig_popt = get_popt('results/GER_HEL_bottlegrowth_split_mig_inbreeding_pop1_only_demo_fits_combined.txt')
    plot_bottlegrowth_split_mig_demography(fs, ns, pts_l, bottlegrowth_split_mig_popt)

    #split_asym_mig_demography
    split_asym_mig_popt = get_popt('results/GER_HEL_split_asym_mig_inbreeding_pop1_only_demo_fits_combined.txt')
    plot_split_asym_mig_demography(fs, ns, pts_l, split_asym_mig_popt)

    #split_delay_mig_demography
    split_delay_mig_popt = get_popt('results/GER_HEL_split_delay_mig_inbreeding_pop1_only_demo_fits_combined.txt')
    plot_split_delay_mig_demography(fs, ns, pts_l, split_delay_mig_popt)

    #split_mig_demography
    split_mig_popt = get_popt('results/GER_HEL_split_mig_inbreeding_pop1_only_demo_fits_combined.txt')
    plot_split_mig_demography(fs, ns, pts_l, split_mig_popt)
if __name__ == '__main__':
    main()