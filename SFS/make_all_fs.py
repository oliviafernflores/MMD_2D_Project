from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt


'''
Synonymous
'''
dd_syn = pickle.load(open('/Users/oliviafernflores/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_syn_with_ancestral.bpkl', 'rb'))
# France vs. Iran
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [16, 10]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('IRA_FRA_syn_unfolded.fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('IRA_FRA_syn_unfolded_fs.png')
plt.clf()
# Germany vs. France
pop_ids = ['Mmd_GER', 'Mmd_FRA']
ns = [16, 16]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('FRA_GER_syn_unfolded.fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('FRA_GER_syn_unfolded_fs.png')
plt.clf()
# Heligoland vs. Germany
pop_ids = ['Mmd_HEL', 'Mmd_GER']
ns = [6, 16]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('GER_HEL_syn_unfolded.fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('GER_HEL_syn_unfolded_fs.png')
plt.clf()
##########################################################
'''
Non Synonymous
'''
dd_nsyn = pickle.load(open('/Users/oliviafernflores/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_nsyn_with_ancestral.bpkl', 'rb'))
# France vs. Iran
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [16, 10]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('IRA_FRA_nsyn_unfolded.fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('IRA_FRA_nsyn_unfolded_fs.png')
plt.clf()
# Germany vs. France
pop_ids = ['Mmd_GER', 'Mmd_FRA']
ns = [16, 16]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('FRA_GER_nsyn_unfolded.fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('FRA_GER_nsyn_unfolded_fs.png')
plt.clf()
# Heligoland vs. Germany
pop_ids = ['Mmd_HEL', 'Mmd_GER']
ns = [6, 16]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('GER_HEL_nsyn_unfolded.fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('GER_HEL_nsyn_unfolded_fs.png')
plt.clf()