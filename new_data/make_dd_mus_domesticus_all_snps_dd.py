#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="mmd_all_snps_dd"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

from contextlib import AsyncExitStack
import pickle
import dadi
import matplotlib.pyplot as plt

#this data dictionary was created June 26, 2024

#using the popfile to remove the relatives at this stage
#this dd will have all SNPs for the four populations of domesticus and will not include individuals who are relatives
dd = dadi.Misc.make_data_dict_vcf('mmd_annovar_merged.vcf.gz', 'mouse_pop_by_geography.remove_relatives.popfile.txt')
pick = open('mus_domesticus_all_snps_all_pops_no_relatives_dd.bpkl','wb')
pickle.dump(dd, pick, 2)