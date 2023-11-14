#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="all_mouse_pops_dd"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

from contextlib import AsyncExitStack
import pickle
import dadi
import matplotlib.pyplot as plt

dd = pickle.load(open('all_pops_all_snps_mice_dd.bpkl', 'rb'))
dd_syn = {}
dd_nsyn = {}
ls_syn = []
ls_nsyn = []
ref_file = open('all_mouse_output.exonic_variant_function', 'r')
# full = ref_file.readline()
# split = full.split('\t')
# print(split)
for line in ref_file:
    single_line = line.split('\t')
    if single_line[1] == 'synonymous SNV' and 'Y' not in single_line[3] and 'X' not in single_line[3]:
        ls_syn.append(single_line[3] + '_' + single_line[4])
    if single_line[1] == 'nonsynonymous SNV' and 'Y' not in single_line[3] and 'X' not in single_line[3]:
        ls_nsyn.append(single_line[3] + '_' + single_line[4])
for i in range(len(ls_syn)):
    dd_syn[ls_syn[i]] = dd[ls_syn[i]]
for i in range(len(ls_nsyn)):
    dd_nsyn[ls_nsyn[i]] = dd[ls_nsyn[i]]
pick = open('mouse_all_pops.dd_syn.bpkl','wb')
pickle.dump(dd_syn, pick, 2)
pick = open('mouse_all_pops.dd_nsyn.bpkl','wb')
pickle.dump(dd_nsyn, pick, 2)