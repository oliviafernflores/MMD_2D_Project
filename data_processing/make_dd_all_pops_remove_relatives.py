#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="mice_dd"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=90
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

import dadi, pickle


dd = dadi.Misc.make_data_dict_vcf('/groups/rgutenk/mice/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz', '/groups/rgutenk/mice/mouse_pop_by_geography.remove_relatives.popfile.txt')
fid = open('all_pops_all_snps_mice_dd.bpkl','wb')
pickle.dump(dd, fid, protocol = 2)


print('loading dd')
dd = pickle.load(open('all_pops_all_snps_mice_dd.bpkl','rb'))

fi = open('AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.nucleotide_info_Ms_SPRE_ancestral_annotation.txt','w')
fi.write('#chromosome\tposition\tref\talt\tancestral\n')

fi_strict = open('AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.nucleotide_info_Ms_SPRE_and_Mmc_CAST_ancestral_annotation.txt','w')
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