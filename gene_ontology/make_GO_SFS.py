#!/usr/bin/env python
#SBATCH --job-name=make_GO_SFS
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

population_pairs = [
    ('Mmd_FRA', 'Mmd_IRA'),
    ('Mmd_GER', 'Mmd_FRA'),
    ('Mmd_HEL', 'Mmd_GER')
]

data_folder = '/groups/rgutenk/oliviafernflores/gene_ontology/GO_term_dictionaries'
go_term_files = [f for f in os.listdir(data_folder) if f.endswith('_dict.pkl')]

for go_term_file in go_term_files:
    go_term = os.path.splitext(go_term_file)[0]
    with open(os.path.join(data_folder, go_term_file), 'rb') as file:
        data_dict = pickle.load(file)

    if not data_dict:
        print(f"Skipping empty dictionary for GO term: {go_term}")
        continue 
    
    GO_index = go_term.split('_')[0]
    fname = f'GO_{GO_index}_SFS'
    go_term_dir = os.path.join('/groups/rgutenk/oliviafernflores/gene_ontology/GO_term_SFS', fname)
    os.makedirs(go_term_dir, exist_ok=True)


    for pop1, pop2 in population_pairs:
        if pop2 == 'Mmd_IRA':
            sfs = dadi.Spectrum.from_data_dict(data_dict, [pop1, pop2], [16, 10])
        elif pop1 == 'Mmd_HEL':
            sfs = dadi.Spectrum.from_data_dict(data_dict, [pop1, pop2], [6, 16])
        else:
            sfs = dadi.Spectrum.from_data_dict(data_dict, [pop1, pop2], [16, 16])
        
        sfs_filepath = os.path.join(go_term_dir, f'{go_term}_{pop1}_{pop2}_SFS')
        sfs.to_file(f'{sfs_filepath}_nsyn_unfolded.fs')
        dadi.Plotting.plot_single_2d_sfs(sfs, vmin = 1)
        plt.savefig(f'{sfs_filepath}_nsyn_unfolded.png')
        plt.clf()

print("SFS calculations and plots are complete.")
