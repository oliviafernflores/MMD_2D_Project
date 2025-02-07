#!/usr/bin/env python
#!/usr/bin/env python
#SBATCH --job-name=syn_bootstraps
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1-5
from contextlib import AsyncExitStack
import dadi
import pickle
import nlopt
import matplotlib.pyplot as plt
import sys, os

dd_syn = pickle.load(open('/groups/rgutenk/oliviafernflores/mus_all_pops.dd_syn_with_ancestral.bpkl', 'rb'))
pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [16, 10]

chunks = dadi.Misc.fragment_data_dict(dd_syn, 1e7)

boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, 100, pop_ids, ns)

for i in range(len(boots)):
  boots[i].to_file('/groups/rgutenk/oliviafernflores/IRA_FRA/bootstraps/IRA_FRA_syn.boot_{0}.fs'.format(str(i)))
