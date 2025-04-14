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

# import matplotlib
# matplotlib.use('Agg')  # Use a non-interactive backend

import matplotlib.pyplot as plt

import pandas as pd


population_pairs = [
    ('Mmd_FRA', 'Mmd_IRA'),
    ('Mmd_GER', 'Mmd_FRA'),
    ('Mmd_HEL', 'Mmd_GER')
]

folders_dir = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_term_dictionaries_descendants'
go_term_folders = [fol for fol in os.listdir(folders_dir)]

size_info_file = pd.read_csv('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_descendants_with_names_sorted.csv', index_col = 0)

sample_sizes = [50000, 25000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 500, 100]


for folder in go_term_folders:
    directory = os.path.join(folders_dir, folder)
    go_term_files = []
    if 'Store' not in folder and 'Store' not in folders_dir:
        for f in os.listdir(directory):
            if f.endswith('_dict.pkl'):
                go_term_files.append(os.path.join(directory, f))
                

    for go_term_file in go_term_files:
        go_term_file_name = go_term_file.split('/')[-1]
        if go_term_file_name.startswith('GO'):
            go_term = go_term_file_name.split(':')[1].split('_')[0]
        else:
            go_term = go_term_file.split('_')[0]
        with open(go_term_file, 'rb') as file:
            data_dict = pickle.load(file)

        if not data_dict:
            print(f"Skipping empty dictionary for GO term: {go_term}")
            continue 
        else:
            print(f'Made SFS for GO term: {go_term}')
        
        GO_index = go_term.split('_')[0]
        fname = f'GO_{GO_index}_SFS'
        go_term_dir = os.path.join(directory, fname)
        go_term_dir = go_term_dir.replace('GO_term_dictionaries_descendants', 'GO_SFS_descendants')
        go_term_dir_lst = go_term_dir.split('/')
        del go_term_dir_lst[len(go_term_dir_lst)-2]
        
        go_term_size = size_info_file.loc['GO:'+str(go_term), 'Number of Entries']
        go_dir_size_index = ''
        for i in range(len(sample_sizes)):
            if i == 0:
                if go_term_size >= sample_sizes[0]:
                    go_dir_size_index = str(sample_sizes[0]) + '_sites_or_more'
            elif i == len(sample_sizes)-1:
                if go_term_size < sample_sizes[-1]:
                    go_dir_size_index = 'less_than_' + str(sample_sizes[-1]) + '_sites'
            else:
                if go_term_size < sample_sizes[i] and go_term_size >= sample_sizes[i+1]:
                    go_dir_size_index = str(sample_sizes[i+1]) + '_sites_or_more'
        size_info_file.loc['GO:'+str(go_term), 'size directory'] = go_dir_size_index

        go_term_dir_lst.insert(-1, go_dir_size_index)
        
        go_term_dir = '/'.join(go_term_dir_lst)

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
            dadi.Plotting.plot_single_2d_sfs(sfs, vmin = 1, show = False)
            plt.savefig(f'{sfs_filepath}_nsyn_unfolded.png')
            plt.close()

size_info_file.to_csv('size_info_file_with_sfs_size_directories.csv')

print("SFS calculations and plots are complete.")
