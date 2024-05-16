#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=standard
#SBATCH --job-name="mice_dd"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem-per-cpu=41gb
#SBATCH --constraint=hi_mem
#SBATCH --time=72:00:00

import dadi, pickle


dd = dadi.Misc.make_data_dict_vcf('/groups/rgutenk/mice/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz', '/groups/rgutenk/mice/mouse_pop_by_geography.remove_relatives.popfile.txt')
fid = open('test_mice_dd.bpkl','wb')
pickle.dump(dd, fid, protocol = 2)