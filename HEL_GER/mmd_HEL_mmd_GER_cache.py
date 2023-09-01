#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="GER_HEL_im_pre_inbreeding_cache"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1

import dadi
import dadi.DFE as DFE
import pickle

if __name__ == '__main__':

    data_fs = dadi.Spectrum.from_file('GER_HEL_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns) + 110, max(ns) + 120, max(ns) + 130]

    # using the IM_pre parameters with inbreeding only in the HEL population

    demo_sel_model = DFE.DemogSelModels.IM_pre

    demo_params = [0.3060088614080329, 0.05367337350992898, 0.14744426867687283, 0.03573416869342573, 0.6340741046187689, 0.045314335100249334, 0.11760025775130044, 0.19431372592724705]

    spectra = DFE.Cache2D(demo_params, ns, demo_sel_model, pts = pts_l, gamma_bounds = (1e-5, 2000), gamma_pts = 50)

    fid = open('mmd_HEL_mmd_GER_2d_cache.bpkl','wb')
    pickle.dump(spectra, fid, protocol = 2)