#!/usr/bin/env python
#SBATCH --account=rgutenk
#SBATCH --qos=user_qos_rgutenk
#SBATCH --partition=high_priority
#SBATCH --job-name="mmd_sort_snps"
#SBATCH --output=%x-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=500gb
#SBATCH --time=72:00:00

import pickle, dadi, pandas as pd, numpy as np

mmd_dd_with_anc_info = pickle.load(open('mus_domesticus_all_snps_all_pops_no_relatives_with_ancestral_dd.bpkl','rb'))


# Prepare to categorize SNPs
dd_syn = {}
dd_nsyn = {}

# Using sets for faster membership tests
ls_syn = set()
ls_nsyn = set()

# Load synonymous and nonsynonymous data
with open('mmd_table_gene_annotation.refGene.exonic_variant_function', 'r') as ref_file:
    print('Reference synonymous and nonsynonymous file loaded')
    print('*' * 20)
    print("Beginning sort into synonymous and nonsynonymous")
    print('*' * 20)

    for line in ref_file:
        single_line = line.split('\t')
        chrom = single_line[3]
        print(chrom)
        if chrom in {'X', 'Y'}:
            continue  # Skip X and Y chromosomes

        position = f"{chrom}_{single_line[4]}"
        if single_line[1] == 'synonymous SNV':
            ls_syn.add(position)
        elif single_line[1] == 'nonsynonymous SNV':
            ls_nsyn.add(position)

# Add SNPs to the respective dictionaries
for pos in ls_syn:
    dd_syn[pos] = mmd_dd_with_anc_info.get(pos)

for pos in ls_nsyn:
    dd_nsyn[pos] = mmd_dd_with_anc_info.get(pos)

print('Finished sorting into synonymous and nonsynonymous SNPs')

# Save the synonymous and nonsynonymous data dictionaries
with open('mus_musculus_domesticus.dd_syn.bpkl', 'wb') as pick:
    pickle.dump(dd_syn, pick, 2)

with open('mus_musculus_domesticus.dd_nsyn.bpkl', 'wb') as pick:
    pickle.dump(dd_nsyn, pick, 2)

print('Finished saving synonymous and nonsynonymous data dictionaries')