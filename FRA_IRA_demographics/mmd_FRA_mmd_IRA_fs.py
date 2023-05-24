import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt

dd = pickle.load(open('/groups/rgutenk/mouse_all_pops.dd_syn.bpkl', 'rb'))
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [14, 14]
fs = dadi.Spectrum_.from_data_dict(dd, pop_ids, ns, polarized = False)
plt.figure('IRA_FRA_fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('IRA_FRA_fs.png')