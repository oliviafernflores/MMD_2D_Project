import pickle, random
import numpy as np
import nlopt
import dadi
import dadi.DFE as DFE
import matplotlib.pyplot as plt
from dadi.DFE import *

def gamma_dfe(theta_ns, spectra, data):
    
    
    sel_params = [0.2, 1000]
    lower_bound, upper_bound = [1e-3, 1e-2], [1, 50000]
    p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,
                                upper_bound=upper_bound)
    popt = dadi.Inference.opt(p0, data, spectra.integrate, pts=None,
                            func_args=[DFE.PDFs.gamma, theta_ns],
                            lower_bound=lower_bound, upper_bound=upper_bound, 
                            verbose=len(sel_params), maxiter=10, multinom=False)

    model_sfs = spectra.integrate(popt, None, DFE.PDFs.gamma, theta_ns, None)
    
    return model_sfs

def plot_gamma(model, data):
    print()
    
def main():
    theta_0 = 26681.214011656422
    theta_ns = theta_0 * 2.4
    spectra = pickle.load(open('mmd_FRA_mmd_IRA_2d_cache.bpkl', 'rb'))
    data = dadi.Spectrum.from_file('IRA_FRA_nsyn_unfolded.fs')
    
    gamma_dfe(theta_ns, spectra, data)
if __name__ == '__main__':
    main()