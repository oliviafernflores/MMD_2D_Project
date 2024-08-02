#!/usr/bin/env python

#SBATCH --job-name="mmd_ancestral"
#SBATCH --output=%x-%A_%a.out
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=1
#SBATCH --ntasks=90
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

import pickle, dadi

print('loading dd')
dd = pickle.load(open('mus_all_pops.dd_nsyn.bpkl','rb'))
# dd = pickle.load(open('../../mice/mouse_all_pops.dd_syn.bpkl','rb'))
print('dd loaded')

dd_new = {}

for chrom_pos in dd:
    if 16 in dd[chrom_pos]['calls']['Ms_SPRE']:
        if 20 in dd[chrom_pos]['calls']['Mmc_CAST']:
            for i in range(2):
                if dd[chrom_pos]['calls']['Ms_SPRE'][i] == 16 and dd[chrom_pos]['calls']['Mmc_CAST'][i] == 20:
                    dd_new[chrom_pos] = dd[chrom_pos]
                    dd_new[chrom_pos]['outgroup_allele'] = dd[chrom_pos]['segregating'][i]
                    dd_new[chrom_pos]['outgroup_context'] = '-' + dd[chrom_pos]['segregating'][i] + '-'
pick = open('mus_all_pops.dd_nsyn_with_ancestral.bpkl','wb')
pickle.dump(dd_new, pick, 2)