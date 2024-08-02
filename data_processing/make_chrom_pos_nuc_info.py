#!/usr/bin/env python

#SBATCH --job-name=nuc_info
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=1
#SBATCH --ntasks=90
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

import pickle, dadi

print('loading dd')
dd = pickle.load(open('mus_all_pops_all_snps_dd.bpkl','rb'))
# dd = pickle.load(open('../../mice/mouse_all_pops.dd_syn.bpkl','rb'))
print('dd loaded')

fi = open('new_Ms_SPRE_ancestral_annotation.txt','w')
fi.write('#chromosome\tposition\tref\talt\tancestral\n')

fi_strict = open('new_Ms_SPRE_Mmc_CAST_ancestral_annotation.txt','w')
fi_strict.write('#chromosome\tposition\tref\talt\tancestral\n')

for chrom_pos in dd:
    fi.write('\t'.join(chrom_pos.split('_')+list(dd[chrom_pos]['segregating']))+'\t')
    fi_strict.write('\t'.join(chrom_pos.split('_')+list(dd[chrom_pos]['segregating']))+'\t')
    if 16 in dd[chrom_pos]['calls']['Ms_SPRE']:
        for i in range(2):
            if dd[chrom_pos]['calls']['Ms_SPRE'][i] == 16:
                fi.write(dd[chrom_pos]['segregating'][i]+'\n')
        if 20 in dd[chrom_pos]['calls']['Mmc_CAST']:
            for i in range(2):
                if dd[chrom_pos]['calls']['Ms_SPRE'][i] == 16 and dd[chrom_pos]['calls']['Mmc_CAST'][i] == 20:
                    fi_strict.write(dd[chrom_pos]['segregating'][i]+'\n')
