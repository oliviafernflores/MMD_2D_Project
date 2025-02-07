#!/usr/bin/env python
#SBATCH --job-name=GO_bootstraps
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=10
#SBATCH --ntasks=50
#SBATCH --time=24:00:00

import os
import pickle
import numpy as np
import dadi
import matplotlib.pyplot as plt

data_folder = '/groups/rgutenk/oliviafernflores/gene_ontology/GO_term_dictionaries'
go_term_files = [f for f in os.listdir(data_folder) if f.endswith('_dict.pkl')]

pop_ids = ['Mmd_FRA', 'Mmd_IRA']
ns = [16, 10]

for go_term_file in go_term_files:
    go_term = os.path.splitext(go_term_file)[0]
    with open(os.path.join(data_folder, go_term_file), 'rb') as file:
        data_dict = pickle.load(file)

    if not data_dict:
        print(f"Skipping empty dictionary for GO term: {go_term}")
        continue 
    
    GO_index = go_term.split('_')[0]
    fname = f'GO_{GO_index}_bootstraps'
    go_term_dir = os.path.join('/groups/rgutenk/oliviafernflores/gene_ontology/GO_bootstraps', fname)
    os.makedirs(go_term_dir, exist_ok=True)

    chunks = dadi.Misc.fragment_data_dict(data_dict, 1e7)
    boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, 100, pop_ids, ns)
    
    boots_dir = os.path.join(go_term_dir, f'{GO_index}_IRA_FRA_boot')

    for i in range(len(boots)):
        boots[i].to_file(f'{boots_dir}_{str(i)}.fs')
