#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="mice_sort_snps"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

import pickle, dadi, pandas as pd, numpy as np

dd_all = pickle.load(open('mus_all_pops_all_snps_dd.bpkl','rb'))
print('data dictionary loaded')

dd_syn = {}
dd_nsyn = {}
ls_syn = []
ls_nsyn = []

ref_file = open('new_gene_annotation.exonic_variant_function', 'r')

for line in ref_file:
    single_line = line.split('\t')
    if single_line[1] == 'synonymous SNV' and 'Y' not in single_line[3] and 'X' not in single_line[3]:
        ls_syn.append(single_line[3] + '_' + single_line[4])
    if single_line[1] == 'nonsynonymous SNV' and 'Y' not in single_line[3] and 'X' not in single_line[3]:
        ls_nsyn.append(single_line[3] + '_' + single_line[4])
for i in range(len(ls_syn)):
    dd_syn[ls_syn[i]] = dd_all[ls_syn[i]]
for i in range(len(ls_nsyn)):
    dd_nsyn[ls_nsyn[i]] = dd_all[ls_nsyn[i]]

# Add SNPs to the respective dictionaries
for pos in ls_syn:
    dd_syn[pos] = dd_all.get(pos)

for pos in ls_nsyn:
    dd_nsyn[pos] = dd_all.get(pos)

print('Finished sorting into synonymous and nonsynonymous SNPs')

pick = open('mus_all_pops.dd_syn.bpkl','wb')
pickle.dump(dd_syn, pick, 2)
pick = open('mus_all_pops.dd_nsyn.bpkl','wb')
pickle.dump(dd_nsyn, pick, 2)