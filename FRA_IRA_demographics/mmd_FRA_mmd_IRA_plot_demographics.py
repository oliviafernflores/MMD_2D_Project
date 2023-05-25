import dadi
import matplotlib.pyplot as plt
import numpy as np
import csv

def get_popt(fname):
    popt = []
    with open(fname, newline = '') as fits_file:
        fits = csv.reader(fits_file, delimiter = '\t')
        popt = fits[0]
        for fit in fits:
            if fit[0] > popt[0]:
                popt = fit
    return popt

def plot_im_demography(fs, ns, pts, popt):
    demo_model = dadi.Demographics2D.IM
    demo_model = dadi.Numerics.make_anc_state_misid_func(demo_model)
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
    model = demo_model_ex(popt, ns, pts)
    fig = plt.figure('IRA_FRA_IM_Demography')
    fig.clear()
    dadi.Plotting.plot_2d_comp_multinom(model, fs)
    fig.savefig('IRA_FRA_im_demography.png')
def main():
    fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')
    ns = fs.sample_sizes
    pts_l = [max(ns)+20, max(ns)+30, max(ns)+40]
    im_popt = get_popt('results/IRA_FRA_im_demo_fits_combined.txt')
    print(im_popt)
if __name__ == '__main__':
    main()