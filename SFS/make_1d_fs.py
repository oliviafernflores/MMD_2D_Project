from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt


'''
Synonymous
'''
dd_syn = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_syn_with_ancestral.bpkl', 'rb'))
# Iran
pop_ids = ['Mmd_IRA']
ns = [10]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('IRA_syn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('IRA_syn_unfolded_fs.png')
plt.clf()
# France
pop_ids = ['Mmd_FRA']
ns = [16]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('FRA_syn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('FRA_syn_unfolded_fs.png')
plt.clf()
# Germany
pop_ids = ['Mmd_GER']
ns = [16]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('GER_syn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('GER_syn_unfolded_fs.png')
plt.clf()
# Heligoland
pop_ids = ['Mmd_HEL']
ns = [6]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('HEL_syn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('HEL_syn_unfolded_fs.png')
plt.clf()
##########################################################
'''
Non Synonymous
'''
dd_nsyn = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_nsyn_with_ancestral.bpkl', 'rb'))
# Iran
pop_ids = ['Mmd_IRA']
ns = [10]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('IRA_nsyn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('IRA_nsyn_unfolded_fs.png')
plt.clf()
# France
pop_ids = ['Mmd_FRA']
ns = [16]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('FRA_nsyn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('FRA_nsyn_unfolded_fs.png')
plt.clf()
# Germany
pop_ids = ['Mmd_GER']
ns = [16]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('GER_nsyn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('GER_nsyn_unfolded_fs.png')
plt.clf()
# Heligoland
pop_ids = ['Mmd_HEL']
ns = [6]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('HEL_nsyn_unfolded.fs')
dadi.Plotting.plot_1d_fs(fs)
plt.savefig('HEL_nsyn_unfolded_fs.png')
plt.clf()