#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=standard
#SBATCH --job-name="mice_dd"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem-per-cpu=41gb
#SBATCH --constraint=hi_mem
#SBATCH --time=72:00:00

import dadi, pickle


dd = dadi.Misc.make_data_dict_vcf('/groups/rgutenk/mice/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz', '/groups/rgutenk/mice/mouse_pop_by_geography.remove_relatives.popfile.txt')
fid = open('test_mice_dd.bpkl','wb')
pickle.dump(dd, fid, protocol = 2)

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
pick = open('test_mouse_all_pops.dd_syn.bpkl','wb')
pickle.dump(dd_syn, pick, 2)
pick = open('test_mouse_all_pops.dd_nsyn.bpkl','wb')
pickle.dump(dd_nsyn, pick, 2)