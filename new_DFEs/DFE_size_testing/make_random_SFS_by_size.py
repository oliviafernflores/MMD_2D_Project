#!/usr/bin/env python
#SBATCH --job-name=make_GO_SFS_1D
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
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend

import matplotlib.pyplot as plt
import pandas as pd
import random

pops = ['Mmd_IRA', 'Mmd_FRA', 'Mmd_GER', 'Mmd_HEL', ['Mmd_FRA', 'Mmd_IRA'], ['Mmd_GER', 'Mmd_FRA'], ['Mmd_HEL', 'Mmd_GER']]

sample_sizes = [50000, 25000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1000, 500, 100]

data_dict = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_nsyn_with_ancestral.bpkl', 'rb'))
lst = list(data_dict.items())

for size in sample_sizes:
    for i in range(10):
        '''
        random sample from big dictionary
        need as many sites as there are in the df - this does not mean that the SFS will have the same sums
        '''
        data = dict(random.sample(lst, size))
        for p in pops:
            if p == ['Mmd_FRA', 'Mmd_IRA']:
                sfs = dadi.Spectrum.from_data_dict(data, p, [16, 10])
            elif p == ['Mmd_GER', 'Mmd_FRA']:
                sfs = dadi.Spectrum.from_data_dict(data, p, [16, 16])
            elif p == ['Mmd_HEL', 'Mmd_GER']:
                sfs = dadi.Spectrum.from_data_dict(data, p, [6, 16])
            elif p == 'Mmd_IRA':
                sfs = dadi.Spectrum.from_data_dict(data, [p], [10])
            elif p == 'Mmd_HEL':
                sfs = dadi.Spectrum.from_data_dict(data, [p], [6])
            else:
                sfs = dadi.Spectrum.from_data_dict(data, [p], [16])
            size_index = f'{size}_random_samples_SFS'
            size_dir = os.path.join('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/new_DFEs/DFE_size_testing', size_index)
            os.makedirs(size_dir, exist_ok=True)
            if p in ['Mmd_IRA', 'Mmd_FRA', 'Mmd_GER', 'Mmd_HEL']:
                sfs_filepath = os.path.join(size_dir, f'{p}_1d')
                os.makedirs(sfs_filepath, exist_ok=True)
                sfs_file_name = f'{str(size)}_random_samples_{p}_1D_random_SFS_{str(i)}'
                dadi.Plotting.plot_1d_fs(sfs)
            else:
                sfs_filepath = os.path.join(size_dir, f'{p[0]}_{p[1]}_2d')
                os.makedirs(sfs_filepath, exist_ok=True)
                sfs_file_name = f'{str(size)}_random_samples_{p[0]}_{p[1]}_2D_random_SFS_{str(i)}'
                dadi.Plotting.plot_single_2d_sfs(sfs, vmin = 1)
            sfs_save_string = f'{sfs_file_name}_nsyn_unfolded.fs'
            sfs_save_loc = os.path.join(sfs_filepath, sfs_save_string)
            sfs.to_file(sfs_save_loc)
            plt_save_string = f'{sfs_file_name}_nsyn_unfolded.png'
            plt_save_loc = os.path.join(sfs_filepath, plt_save_string)
            plt.savefig(plt_save_loc)
            plt.clf()
    print('Done making random SFS for GO ' + str(size))

print("SFS calculations and plots are complete.")
