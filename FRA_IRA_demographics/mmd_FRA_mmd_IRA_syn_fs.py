#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="all_mouse_pops_syn_fs"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt

dd = pickle.load(open('/groups/rgutenk/mice/mouse_all_pops.dd_syn.bpkl', 'rb'))
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [14, 14]
fs = dadi.Spectrum.from_data_dict(dd, pop_ids, ns)
fs.to_file('IRA_FRA_syn_unfolded.fs')

nboot, chunk_size = 100, 2e6
chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, nboot, pop_ids, ns)
for i in range(len(boots)):
    boots[i].to_file('syn_bootstraps/IRA_FRA{0}.fs'.format(str(i)))

plt.figure('IRA_FRA_fs')
dadi.Plotting.plot_single_2d_sfs(fs)
plt.savefig('IRA_FRA_syn_fs.png')

plt.figure('IRA_FRA_fs_bootstraps')
dadi.Plotting.plot_single_2d_sfs(boots[0])
plt.savefig('IRA_FRA_syn_bootstraps_fs.png')