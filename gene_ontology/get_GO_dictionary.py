#!/usr/bin/env python
#SBATCH --job-name=get_GO_dictionary
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=20
#SBATCH --ntasks=100
#SBATCH --time=24:00:00

import os
import pickle
import csv
import pandas as pd

data_dict = pickle.load(open('/groups/rgutenk/oliviafernflores/mus_all_pops.dd_nsyn_with_ancestral.bpkl', 'rb'))
biomart = pd.read_table('/groups/rgutenk/oliviafernflores/mart_export.txt', sep = ',')

biomart_dd = {}

for row in biomart.index:
    chrom = biomart.loc[row, 'Chromosome/scaffold name']
    pos_start = biomart.loc[row, 'Gene start (bp)']
    pos_end = biomart.loc[row, 'Gene end (bp)']
    key = str(chrom) + '_' + str(pos_start) + '_' + str(pos_end)
    if key not in biomart_dd:
        if type(biomart.loc[row, 'GO term accession']) == type('hello'):
            biomart_dd[key] = [biomart.loc[row, 'GO term accession']]
    else:
        if type(biomart.loc[row, 'GO term accession']) == type('hello'):
            biomart_dd[key].append(biomart.loc[row, 'GO term accession'])
            
data_dict_GO = data_dict.copy()

biomart_chroms_only = {}

for key in biomart_dd.keys():
    if key[0] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
        biomart_chroms_only[key] = biomart_dd[key]

for key in data_dict_GO.keys():
    chr_num, position = key.split('_')
    position = int(position)
    data_dict_GO[key]['GO_terms'] = []
    for dict2_key, go_terms in biomart_chroms_only.items():
        chr_num_2, start_pos, end_pos = dict2_key.split('_')
        start_pos = int(start_pos)
        end_pos = int(end_pos)
        if chr_num[3:] == chr_num_2 and start_pos <= position <= end_pos:
            data_dict_GO[key]['GO_terms'].extend(go_terms)

with open('mus_all_pops_nsyn_with_ancestral_and_GO_terms.pkl', 'wb') as pickle_file:
    pickle.dump(data_dict_GO, pickle_file)

print("Dictionary saved as 'mus_all_pops_nsyn_with_ancestral_and_GO_terms.pkl'.")