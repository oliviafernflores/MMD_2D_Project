#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="IRA_FRA_im_pre_inbreeding_cache"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1

import dadi
import dadi.DFE as DFE
import pickle


def main():
    data_fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns) + 200, max(ns) + 210, max(ns) + 220]

    # using the IM_pre parameters with inbreeding

    demo_sel_model = DFE.DemogSelModels.IM_pre

    demo_params = [2.3892200757761812, 0.09865800019618395, 0.059546476404449086, 0.29863880148384436, 2.9956042751467993, 0.1765401593467742, 2.1701288461162407, 0.19321473756551447]


    spectra = DFE.Cache2D(demo_params, ns, demo_sel_model, pts = pts_l, gamma_bounds = (1e-5, 2000), gamma_pts = 50, additional_gammas = [1.2, 4.3])

    fid = open('mmd_FRA_mmd_IRA_2d_cache.bpkl','wb')
    pickle.dump(spectra, fid, protocol = 2)
    
    demo_sel_1d = DFE.DemogSelModels.IM_pre_single_gamma
    spectra1d = DFE.Cache1D(demo_params, ns, demo_sel_1d, pts = pts_l, gamma_bounds = (1e-5, 2000), gamma_pts = 50, additional_gammas = [1.2, 4.3])
    fid = open('mmd_FRA_mmd_IRA_1d_cache.bpkl','wb')
    pickle.dump(spectra1d, fid, protocol = 2)
    
    # numpy.sum(cache.spectra < 0)

if __name__ == '__main__':
    main()