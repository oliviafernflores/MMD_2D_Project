#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="mmd_add_ancestral"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

import pickle, dadi, pandas as pd, numpy as np

mmd_dd = pickle.load(open('mus_domesticus_all_snps_all_pops_no_relatives_dd.bpkl','rb'))


anc_info_df = pd.read_table('Ms_SPRE_ancestral_annotation.txt', delimiter='\t')
anc_info_df.columns = ['chromosome', 'position', 'ref', 'alt', 'ancestral']
anc_info_df.set_index(['chromosome', 'position'], inplace=True)
print('Mus Spretus ancestral info loaded')
print('*' * 20)

mmd_dd_with_anc_info = {}
print('Beginning to add ancestral info for all chromosomes')
print('*' * 20)

for chrom_pos in mmd_dd:
    print(chrom_pos)
    chrom, pos = chrom_pos.split('_')
    ancestral_allele = anc_info_df.loc[(chrom, int(pos)), 'ancestral'] if (chrom, int(pos)) in anc_info_df.index else None

    if isinstance(ancestral_allele, str):
        mmd_dd_with_anc_info[chrom_pos] = {
            **mmd_dd[chrom_pos],
            'outgroup_allele': ancestral_allele,
            'outgroup_context': f'-{ancestral_allele}-'
        }

print("Finished adding ancestral info for all chromosomes")
print('*' * 20)

pick = open('mus_domesticus_all_snps_all_pops_no_relatives_with_ancestral_dd.bpkl','wb')
pickle.dump(mmd_dd_with_anc_info, pick, 2)

print("Finished saving ancestral info data dictionary")