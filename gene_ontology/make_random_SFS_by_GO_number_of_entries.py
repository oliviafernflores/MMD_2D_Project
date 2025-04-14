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
import matplotlib.pyplot as plt
import pandas as pd
import random

pops = ['Mmd_IRA', 'Mmd_FRA', 'Mmd_GER', 'Mmd_HEL']

info_file = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/go_term_summary.csv'

df = pd.read_csv(info_file, index_col = 0)
print(df)

data_dict = pickle.load(open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/data_dictionaries/mus_all_pops.dd_nsyn_with_ancestral.bpkl', 'rb'))


for row in df.index:
    if not np.isnan(row):
        if int(df.loc[row, 'Number of Entries']) >= 10000:
            number = int(df.loc[row, 'Number of Entries'])
            lst = list(data_dict.items())
            for i in range(10):
                '''
                random sample from big dictionary
                need as many sites as there are in the df - this does not mean that the SFS will have the same sums
                '''
                data = dict(random.sample(lst, number))
                for p in pops:
                    if p == 'Mmd_IRA':
                        sfs = dadi.Spectrum.from_data_dict(data, [p], [10])
                    elif p == 'Mmd_HEL':
                        sfs = dadi.Spectrum.from_data_dict(data, [p], [6])
                    else:
                        sfs = dadi.Spectrum.from_data_dict(data, [p], [16])
                    GO_index = f'GO_{row}_SFS'
                    go_term_dir = os.path.join('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_term_SFS_random_1D', GO_index)
                    os.makedirs(go_term_dir, exist_ok=True)
                    sfs_filepath = os.path.join(go_term_dir, f'{str(row)}_{p}_1D_random_SFS_{str(i)}')
                    sfs.to_file(f'{sfs_filepath}_nsyn_unfolded.fs')
                    dadi.Plotting.plot_1d_fs(sfs)
                    plt.savefig(f'{sfs_filepath}_nsyn_unfolded.png')
                    plt.clf()
            print('Done making random SFS for GO ' + str(row))
        else:
            print('Not enough entries for GO ' + str(row))

print("SFS calculations and plots are complete.")
