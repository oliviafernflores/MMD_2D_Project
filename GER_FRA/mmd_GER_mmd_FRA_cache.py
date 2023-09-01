#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="FRA_GER_im_pre_inbreeding_cache"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1

import dadi
import dadi.DFE as DFE
import pickle

if __name__ == '__main__':

    data_fs = dadi.Spectrum.from_file('FRA_GER_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns) + 110, max(ns) + 120, max(ns) + 130]

    # using the IM_pre parameters with inbreeding only in the FRA population

    demo_sel_model = DFE.DemogSelModels.IM_pre

    demo_params = [0.4875358502474117, 0.4556158612301682, 0.4383815917563917, 0.16341077205662932, 0.3121766734685466, 0.12404523011097193, 5.197064094647251, 2.3002893927337307]

    spectra = DFE.Cache2D(demo_params, ns, demo_sel_model, pts = pts_l, gamma_bounds = (1e-5, 2000), gamma_pts = 50)

    fid = open('mmd_GER_mmd_FRA_2d_cache.bpkl','wb')
    pickle.dump(spectra, fid, protocol = 2)