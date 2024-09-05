#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="FRA_GER_im_inbreeding_cache"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1

import dadi
import dadi.DFE as DFE
import pickle
# import numpy

def IM_pre_sel_inbreeding(params, ns, pts):
    """
    nuPre: Size after first size change
    TPre: Time before split of first size change.
    s: Fraction of nuPre that goes to pop1. (Pop 2 has size nuPre*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    gamma1: Scaled selection coefficient in pop 1 *and* ancestral pop.
    gamma2: Scaled selection coefficient in pop 2
    """
    nuPre,TPre,s,nu1,nu2,T,m12,m21,F1, F2, gamma1,gamma2 = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma1)
    phi = dadi.Integration.one_pop(phi, xx, TPre, nu=nuPre, gamma=gamma1)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    nu1_0 = nuPre*s
    nu2_0 = nuPre*(1-s)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1_func, nu2_func, m12=m12, m21=m21, gamma1=gamma1, gamma2=gamma2)

    fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2, 2))
    return fs
def IM_sel_inbreeding(params, ns, pts):
    s,nu1,nu2,T,m12,m21,F1,F2,gamma1,gamma2 = params
    return IM_pre_sel_inbreeding((1,0,s,nu1,nu2,T,m12,m21,F1,F2,gamma1,gamma2), ns, pts)
def IM_pre_sel_single_gamma_inbreeding(params, ns, pts):
    nuPre,TPre, s,nu1,nu2,T,m12,m21,F1,F2,gamma1 = params
    return IM_pre_sel_inbreeding((nuPre,TPre,s,nu1,nu2,T,m12,m21,F1,F2,gamma1,gamma1), ns, pts)
def IM_sel_single_gamma_inbreeding(params, ns, pts):
    s,nu1,nu2,T,m12,m21,F1,F2,gamma1 = params
    return IM_pre_sel_inbreeding((1,0,s,nu1,nu2,T,m12,m21,F1,F2,gamma1,gamma1), ns, pts)
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

def main():
    data_fs = dadi.Spectrum.from_file('FRA_GER_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns) + 200, max(ns) + 210, max(ns) + 220]

    # using the IM parameters with inbreeding

    demo_sel_model = IM_sel_inbreeding

    demo_params = get_popt('demo_results/FRA_GER_im_inbreeding_demo_fits_combined.txt')
    demo_params.pop(-1)
    demo_params.pop(-1)
    demo_params.pop(0)
    print(demo_params)

    spectra = DFE.Cache2D(demo_params, ns, demo_sel_model, pts = pts_l, gamma_bounds = (1e-5, 2000), gamma_pts = 50, additional_gammas = [1.2, 4.3])

    fid = open('mmd_GER_mmd_FRA_2d_cache.bpkl','wb')
    pickle.dump(spectra, fid, protocol = 2)
    
    demo_sel_1d = IM_sel_single_gamma_inbreeding
    spectra1d = DFE.Cache1D(demo_params, ns, demo_sel_1d, pts = pts_l, gamma_bounds = (1e-5, 2000), gamma_pts = 50, additional_gammas = [1.2, 4.3])
    fid = open('mmd_GER_mmd_FRA_1d_cache.bpkl','wb')
    pickle.dump(spectra1d, fid, protocol = 2)
    
    # numpy.sum(spectra.spectra < 0)
    # numpy.sum(spectra1d.spectra < 0)

if __name__ == '__main__':
    main()