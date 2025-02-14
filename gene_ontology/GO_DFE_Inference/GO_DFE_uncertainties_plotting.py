import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import pandas as pd
import seaborn as sns

step_0_1_error_percents = []
step_0_01_error_percents = []
step_0_001_error_percents = []

def plot_1d_lognormal(root_dir):
    path_0_1 = root_dir + '0.1_step_size/1d lognormal_GO_results_plus_uncertainties_0.1_step.csv'
    df_0_1 = pd.read_csv(path_0_1, index_col='GO Term')
    
    path_0_01 = root_dir + '0.01_step_size/1d lognormal_GO_results_plus_uncertainties_0.01_step.csv'
    df_0_01 = pd.read_csv(path_0_01, index_col='GO Term')
    
    path_0_001 = root_dir + '0.001_step_size/1d lognormal_GO_results_plus_uncertainties_0.001_step.csv'
    df_0_001 = pd.read_csv(path_0_001, index_col='GO Term')
    
    params = [' log_mu', ' log_sigma', ' misid']

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    fig.suptitle('1D Lognormal DFE - IRA vs. FRA', fontsize=16)
    
    for i, p in enumerate(params):
        ax = axes[i]
        label = p.strip() + ' uncertainty'
        if 'mu' in label:
            ax.set_title('Mu')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=np.log(df_0_1['log_mu uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=np.log(df_0_01['log_mu uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=np.log(df_0_001['log_mu uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'sigma' in label:
            ax.set_title('Sigma')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=np.log(df_0_1['log_sigma uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=np.log(df_0_01['log_sigma uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=np.log(df_0_001['log_sigma uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'misid' in label:
            ax.set_title('Misid')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['misid uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['misid uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['misid uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        
        ax.scatter(df_0_1[p], df_0_1['GO Term Name'], color = 'black')
        ax.set_yticklabels(df_0_1['GO Term Name'])
        ax.set_xlabel(p)
        if i != 0:
            ax.get_yaxis().set_visible(False)

        legend_labels = ['0.1', '0.01', '0.001']
        legend_colors = ['blue', 'green', 'red']
        handles = [mlines.Line2D([0], [0], color=color, lw=2) for color in legend_colors]
    fig.legend(handles=handles, labels=legend_labels, title="Step Size", loc='upper left', ncol=3)
        


    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    # plt.show()
    plt.savefig(root_dir + '1d_lognormal_GO_DFE_uncertainties')

def plot_2d_biv_ind_gamma_asym(root_dir):
    path_0_1 = root_dir + '0.1_step_size/2d_biv_ind_gamma_asymmetric_GO_results_plus_uncertainties_0.1_step.csv'
    df_0_1 = pd.read_csv(path_0_1, index_col='GO Term')
    
    path_0_01 = root_dir + '0.01_step_size/2d_biv_ind_gamma_asymmetric_GO_results_plus_uncertainties_0.01_step.csv'
    df_0_01 = pd.read_csv(path_0_01, index_col='GO Term')
    
    path_0_001 = root_dir + '0.001_step_size/2d_biv_ind_gamma_asymmetric_GO_results_plus_uncertainties_0.001_step.csv'
    df_0_001 = pd.read_csv(path_0_001, index_col='GO Term')
    
    params = [' shape1', ' scale1', ' shape2', ' scale2', ' misid']

    fig, axes = plt.subplots(1, 5, figsize=(18, 6))

    fig.suptitle('2D Bivariate Independent Gamma Asymmetric DFE - IRA vs. FRA', fontsize=16)
    
    for i, p in enumerate(params):
        ax = axes[i]
        label = p.strip() + ' uncertainty'
        if 'shape1' in label:
            ax.set_title('Shape 1')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['shape1 uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['shape1 uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['shape1 uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'scale1' in label:
            ax.set_title('Scale 1')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['scale1 uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['scale1 uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['scale1 uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'shape2' in label:
            ax.set_title('Shape 2')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['shape2 uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['shape2 uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['shape2 uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'scale2' in label:
            ax.set_title('Scale 2')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['scale2 uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['scale2 uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['scale2 uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'misid' in label:
            ax.set_title('Misid')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['misid uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_01[p], df_0_01['GO Term Name'], xerr=(df_0_01['misid uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_001[p], df_0_001['GO Term Name'], xerr=(df_0_001['misid uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        
        ax.scatter(df_0_1[p], df_0_1['GO Term Name'], color = 'black')
        ax.set_xlabel(p)
        ax.set_yticklabels(df_0_1['GO Term Name'])
        
        if i != 0:
            ax.get_yaxis().set_visible(False)

        legend_labels = ['0.1', '0.01', '0.001']
        legend_colors = ['blue', 'green', 'red']
        handles = [mlines.Line2D([0], [0], color=color, lw=2) for color in legend_colors]
    fig.legend(handles=handles, labels=legend_labels, title="Step Size", loc='upper left', ncol=3)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    # plt.show()
    plt.savefig(root_dir + '2d_biv_ind_gamma_asymmetric_GO_DFE_uncertainties')

def plot_2d_biv_ind_gamma_sym(root_dir):
    path_0_1 = root_dir + '0.1_step_size/2d_biv_ind_gamma_symmetric_GO_results_plus_uncertainties_0.1_step.csv'
    df_0_1 = pd.read_csv(path_0_1, index_col='GO Term')
    
    path_0_01 = root_dir + '0.01_step_size/2d_biv_ind_gamma_symmetric_GO_results_plus_uncertainties_0.01_step.csv'
    df_0_01 = pd.read_csv(path_0_01, index_col='GO Term')
    
    path_0_001 = root_dir + '0.001_step_size/2d_biv_ind_gamma_symmetric_GO_results_plus_uncertainties_0.001_step.csv'
    df_0_001 = pd.read_csv(path_0_001, index_col='GO Term')
    
    params = [' shape', ' scale',  ' misid']

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    fig.suptitle('2D Bivariate Independent Gamma Symmetric DFE - IRA vs. FRA', fontsize=16)
    
    for i, p in enumerate(params):
        ax = axes[i]
        label = p.strip() + ' uncertainty'
        if 'shape' in label:
            ax.set_title('Shape')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['shape uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['shape uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['shape uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'scale' in label:
            ax.set_title('Scale')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['scale uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['scale uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['scale uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'misid' in label:
            ax.set_title('Misid')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['misid uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_01[p], df_0_01['GO Term Name'], xerr=(df_0_01['misid uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_001[p], df_0_001['GO Term Name'], xerr=(df_0_001['misid uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        
        ax.scatter(df_0_1[p], df_0_1['GO Term Name'], color = 'black')
        ax.set_xlabel(p)
        ax.set_yticklabels(df_0_1['GO Term Name'])
        
        if i != 0:
            ax.get_yaxis().set_visible(False)

        legend_labels = ['0.1', '0.01', '0.001']
        legend_colors = ['blue', 'green', 'red']
        handles = [mlines.Line2D([0], [0], color=color, lw=2) for color in legend_colors]
    fig.legend(handles=handles, labels=legend_labels, title="Step Size", loc='upper left', ncol=3)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    # plt.show()
    plt.savefig(root_dir + '2d_biv_ind_gamma_symmetric_GO_DFE_uncertainties')

def plot_2d_gamma(root_dir):
    path_0_1 = root_dir + '0.1_step_size/2d_gamma_GO_results_plus_uncertainties_0.1_step.csv'
    df_0_1 = pd.read_csv(path_0_1, index_col='GO Term')
    
    path_0_01 = root_dir + '0.01_step_size/2d_gamma_GO_results_plus_uncertainties_0.01_step.csv'
    df_0_01 = pd.read_csv(path_0_01, index_col='GO Term')
    
    path_0_001 = root_dir + '0.001_step_size/2d_gamma_GO_results_plus_uncertainties_0.001_step.csv'
    df_0_001 = pd.read_csv(path_0_001, index_col='GO Term')
    
    params = [' shape', ' scale',  ' misid']

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    fig.suptitle('2D Gamma DFE - IRA vs. FRA', fontsize=16)
    
    for i, p in enumerate(params):
        ax = axes[i]
        label = p.strip() + ' uncertainty'
        if 'shape' in label:
            ax.set_title('Shape')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['shape uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['shape uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['shape uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'scale' in label:
            ax.set_title('Scale')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['scale uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_01['GO Term Name'], xerr=(df_0_01['scale uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_1[p], df_0_001['GO Term Name'], xerr=(df_0_001['scale uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        elif 'misid' in label:
            ax.set_title('Misid')
            ax.errorbar(df_0_1[p], df_0_1['GO Term Name'], xerr=(df_0_1['misid uncertainty']), color = 'blue', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_01[p], df_0_01['GO Term Name'], xerr=(df_0_01['misid uncertainty']), color = 'green', fmt='o', capsize=5, markersize=0)
            ax.errorbar(df_0_001[p], df_0_001['GO Term Name'], xerr=(df_0_001['misid uncertainty']), color = 'red', fmt='o', capsize=5, markersize=0)
        
        ax.scatter(df_0_1[p], df_0_1['GO Term Name'], color = 'black')
        ax.set_xlabel(p)
        ax.set_yticklabels(df_0_1['GO Term Name'])
        
        if i != 0:
            ax.get_yaxis().set_visible(False)

        legend_labels = ['0.1', '0.01', '0.001']
        legend_colors = ['blue', 'green', 'red']
        handles = [mlines.Line2D([0], [0], color=color, lw=2) for color in legend_colors]
    fig.legend(handles=handles, labels=legend_labels, title="Step Size", loc='upper left', ncol=3)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    # plt.show()
    plt.savefig(root_dir + '2d_gamma_GO_DFE_uncertainties')


if __name__ == "__main__":
    root_dir = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/'
    plot_1d_lognormal(root_dir)
    plot_2d_biv_ind_gamma_asym(root_dir)
    plot_2d_biv_ind_gamma_sym(root_dir)
    plot_2d_gamma(root_dir)