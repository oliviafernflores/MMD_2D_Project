from inbreeding_models import IM, IM_pre
import dadi
import dadi.DFE as DFE
import pickle

if __name__ == '__main__':

    data_fs = dadi.Spectrum.from_file('IRA_FRA_syn_unfolded.fs')
    ns = data_fs.sample_sizes
    pts_l = [max(ns) + 10, max(ns) + 20, max(ns) + 30]
    demo_sel_model = DFE.DemogSelModels.IM_pre

    demo_params = [2.3892200757761812, 0.09865800019618395, 0.059546476404449086, 0.29863880148384436, 2.9956042751467993, 0.1765401593467742, 2.1701288461162407, 0.19321473756551447]


    spectra = DFE.Cache2D(demo_params, ns, demo_sel_model, pts = pts_l, gamma_bounds = (1e-5, 2000), gamma_pts = 50)

    fid = open('mmd_FRA_mmd_IRA_2d_cache.bpkl','wb')
    pickle.dump(spectra, fid, protocol = 2)