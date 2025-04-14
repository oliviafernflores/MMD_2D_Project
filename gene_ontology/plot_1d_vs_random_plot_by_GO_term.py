import os
import pickle
import numpy as np
import dadi
import matplotlib.pyplot as plt
import pandas as pd
import random

pops = ['Mmd_IRA', 'Mmd_FRA', 'Mmd_GER', 'Mmd_HEL']

info_file = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/go_term_summary.csv'

df = pd.read_csv(info_file, index_col=0)

over_cutoff = []

for row in df.index:
    if not np.isnan(row):
        if int(df.loc[row, 'Number of Entries']) >= 10000:
            over_cutoff.append(row)

random_dir = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_term_SFS_random_1D'
source_dir = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_term_SFS_1D'

source_dir_lst = [f for f in os.listdir(source_dir)]



for dir in source_dir_lst:
    GO_num = dir.split('_')[1]
    if GO_num != 'Store':
        if int(GO_num) in over_cutoff:
            # Create the figure and axes for all four populations in one plot
            # fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(10, 6), sharey = 'row', sharex = 'col')  # 3 rows, 4 columns
            fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(16, 7), sharex = 'col') 

            handles = []
            labels = []

            for p in range(len(pops)):
                source_sfs = dadi.Spectrum.from_file(os.path.join(source_dir, dir, f'{GO_num}_dict_{pops[p]}_SFS_nsyn_unfolded.fs'))
                random_sfs_lst = []
                print(pops[p])
                print('Source: ' + str(np.nansum(source_sfs)))
                for i in range(10):
                    fs_random = dadi.Spectrum.from_file(os.path.join(random_dir, f'GO_{int(GO_num)}.0_SFS', f'{int(GO_num)}.0_{pops[p]}_1D_random_SFS_{str(i)}_nsyn_unfolded.fs'))
                    random_sfs_lst.append(fs_random)
                    
                    if i < 1:
                        rand, = axes[0, p].plot(fs_random, marker='s', linestyle='-', color='blue', label='random')
                        print('Random: ' + str(np.nansum(fs_random)))
                        if p == 0:
                            handles.append(rand)
                            labels.append('random')
                    else:
                        axes[0, p].plot(fs_random, marker='s', linestyle='-', color='blue')

                # Plot the "truth" (source SFS)
                truth, = axes[0, p].plot(source_sfs, marker='o', linestyle='-', color='red', label='GO Data')
                if p == 0:
                    handles.append(truth)
                    labels.append('GO Data')
                # axes[0, p].legend()

                averages = []
                randoms = pd.DataFrame(random_sfs_lst)
                randoms[0] = np.nan
                randoms[len(randoms.columns)-1] = np.nan
                averages = list(randoms.mean())
                avgs, = axes[1, p].plot(averages, marker='s', linestyle='-', color='green', label='random average')
                if p == 0:
                    handles.append(avgs)
                    labels.append('random average')
                axes[1, p].plot(source_sfs, marker='o', linestyle='-', color='red', label='truth')
                # axes[1, p].legend()

                # Plot the residuals (difference between truth and average)
                resids, = axes[2, p].plot(source_sfs - averages, marker='*', linestyle='-', color='black', label='residuals')
                if p == 0:
                    handles.append(resids)
                    labels.append('residuals \n(GO data - average)')
                # axes[2, p].legend()

                # Set custom xticks for each population
                if pops[p] == 'Mmd_HEL':
                    for a in range(len(axes)):
                        axes[a, 3].set_xticks(range(0, 7))
                elif pops[p] == 'Mmd_IRA':
                    for a in range(len(axes)):
                        axes[a, 0].set_xticks(range(0, 11))
                else:
                    for a in range(len(axes)):
                        axes[a, 1].set_xticks(range(0, 17))
                        axes[a, 2].set_xticks(range(0, 17))

                axes[0, p].set_title(pops[p])

            # Set the title for the entire figure
            fig.suptitle(f'GO {GO_num}: {df.loc[int(GO_num), "GO Term Name"]}', fontsize=16)

            

            # Adjust layout to prevent overlap
            plt.tight_layout()
            plt.subplots_adjust(top=0.9, right=0.85)  # Adjust the top to make space for the title

            fig.legend(handles=handles, labels=labels, loc='center right', bbox_to_anchor=(1, 0.5), fontsize=12)
            
            plt.savefig(f'/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_SFS_1D_random_vs_truth/GO_{GO_num}_{df.loc[int(GO_num), "GO Term Name"]}.pdf')

            # Show the plot
            # plt.show()
