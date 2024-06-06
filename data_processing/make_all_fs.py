from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt


'''
Synonymous
'''
dd_syn = pickle.load(open('mouse_all_pops_no_relatives_yes_outgroup.dd_syn.bpkl', 'rb'))
# France vs. Iran
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [16, 10]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('IRA_FRA_syn_unfolded.fs')
# Germany vs. France
pop_ids = ['Mmd_GER', 'Mmd_FRA']
ns = [16, 16]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('FRA_GER_syn_unfolded.fs')
# Heligoland vs. Germany
pop_ids = ['Mmd_HEL', 'Mmd_GER']
ns = [6, 16]
fs = dadi.Spectrum.from_data_dict(dd_syn, pop_ids, ns)
fs.to_file('GER_HEL_syn_unfolded.fs')
##########################################################
'''
Non Synonymous
'''
dd_nsyn = pickle.load(open('mouse_all_pops_no_relatives_yes_outgroup.dd_nsyn.bpkl', 'rb'))
# France vs. Iran
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [16, 10]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('IRA_FRA_nsyn_unfolded.fs')
# Germany vs. France
pop_ids = ['Mmd_GER', 'Mmd_FRA']
ns = [16, 16]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('FRA_GER_nsyn_unfolded.fs')
# Heligoland vs. Germany
pop_ids = ['Mmd_HEL', 'Mmd_GER']
ns = [6, 16]
fs = dadi.Spectrum.from_data_dict(dd_nsyn, pop_ids, ns)
fs.to_file('GER_HEL_nsyn_unfolded.fs')