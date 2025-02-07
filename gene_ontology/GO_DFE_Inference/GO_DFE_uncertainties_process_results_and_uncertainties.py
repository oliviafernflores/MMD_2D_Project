import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def extract_info_from_filename(filename):
    '''
    Get the model name and GO number from each filename. 
    Important because everything is in nested in directories sorted
    by GO number. This function is used in process_dictionaries().
    '''
    base_name = os.path.basename(filename)
    parts = base_name.split('_')
    model_name = "_".join(parts[4:-1])
    go_number = float(parts[-1].split('.')[0])
    return model_name, go_number

def process_file(file_path):
    '''
    For each results file, need to get the log likelihood, 
    parameter values, list of parameter names, and know
    if the optimizations converged or not. This function
    is used in process_dictionaries().
    '''
    with open(file_path, 'r') as f:
        lines = f.readlines()
    start_index = None
    for i, line in enumerate(lines):
        if line.strip() == "# Top 100 results":
            start_index = i
            break
    if start_index is None:
        raise ValueError(f"Could not find '# Top 100 results' in {file_path}")
    converged = "no"
    for i, line in enumerate(lines):
        if line.strip() == "# Converged results":
            converged = "yes"
            break
    param_names_line = lines[start_index + 1].strip().split('\t')
    result_line = lines[start_index + 2].strip().split('\t')
    params_df = pd.DataFrame([result_line], columns=param_names_line, index=None)
    for col in params_df.columns:
        params_df.loc[0, col] = round(float(params_df.loc[0, col]), 3)
    log_likelihood = params_df.iloc[0, 0]
    params_dict = params_df.to_dict('index')
    return log_likelihood, params_dict[0], param_names_line, converged


def process_directories(root_folder, go_csv_path, output_csv_path):
    '''
    For each of the files in the root folder:
    1. If it's a best fits file, get all the paramter info using process_file()
    and the model info from extract_info_from_filename()
    
    Turn this all into a csv file and save (has all info for all files)
    '''
    go_df = pd.read_csv(go_csv_path, index_col='GO Term')
    results = {}
    for subdir, _, files in os.walk(root_folder):
        for file in files:
            if file.startswith("mmd_FRA_mmd_IRA") and file.endswith(".InferDFE.bestfits"):
                file_path = os.path.join(subdir, file)
                try:
                    log_likelihood, params_dict, param_names, converged = process_file(file_path)
                    if log_likelihood is None:
                        continue
                    model_name, go_number = extract_info_from_filename(file)
                    if go_number not in results:
                        results[go_number] = {}
                    model_data = f"LL: {log_likelihood}, "
                    filtered_params = {param: value for param, value in params_dict.items() if param not in ['# Log(likelihood)', 'theta']}
                    model_data += ", ".join([f"{param}: {value}" for param, value in filtered_params.items()])
                    model_data += f", Converged: {converged}"
                    results[go_number][model_name] = model_data
                except ValueError as e:
                    print(f"Error processing {file}: {e}")
                    continue
    final_data = []
    for go_number, models in results.items():
        row = {'GO Term': go_number}
        for model_name, model_data in models.items():
            row[model_name] = model_data
        final_data.append(row)
    final_df = pd.DataFrame(final_data)
    final_df = go_df.join(final_df.set_index('GO Term'), how='left', on='GO Term')
    final_df.to_csv(output_csv_path, index = True)
    
def add_uncerts(root_folder, go_plus_results_dir, output_csv_path):
    '''
    Using the csv file from process_dictionaries(), add all the uncertainties
    info. This is not formatted, raw text from the uncerts file is added.
    '''
    df = pd.read_csv(go_plus_results_dir,index_col='GO Term')
    results = {}
    for subdir, _, files in os.walk(root_folder):
        for file in files:
            if file.startswith("mmd_FRA_mmd_IRA") and file.endswith(".InferDFE.godambe"):
                file_path = os.path.join(subdir, file)
                model_name, go_number = extract_info_from_filename(file)
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                if go_number not in results:
                        results[go_number] = {}
                results[go_number][model_name] = lines
    for go_num in df.index:
        for model in df.columns:
            if go_num in results:
                if model in results[go_num]:
                    try:
                        for i in range(len(results[go_num][model])):
                            if str(results[go_num][model][i]).startswith('Estimated'):
                                df.loc[go_num, model] += '\n' + str(results[go_num][model][i]).strip('\n')
                    except:
                        continue
    df.to_csv(output_csv_path, index = True)
    
def get_data_by_model(uncerts_csv_path, root_folder):
    '''
    Using the csv file from add_uncerts() that has all models and all uncertainties,
    split that data into columns based on step size of the uncerts. 
    '''
    df = pd.read_csv(uncerts_csv_path, index_col="GO Term")
    model_names = list(df.columns[2:])
    for model in model_names:
        model_results = df[model].dropna()
        go_names = df["GO Term Name"].dropna()
        data = pd.DataFrame([model_results, go_names]).T
        
        data_dict = {}
        for i in range(len(data.index)):
            if data.index[i] not in data_dict:
                data_dict[data.index[i]] = []
                vals = data.iloc[i, 0].split('\n')
                if len(vals) > 1:
                    data_dict[data.index[i]].append(vals[0].split(','))
                    data_dict[data.index[i]].append(vals[1].split(':')[1])
                    data_dict[data.index[i]].append(vals[2].split(':')[1])
                    data_dict[data.index[i]].append(vals[3].split(':')[1])
                else:
                    data_dict[data.index[i]].append(vals[0].split(','))
        col_names = [ele.split(":")[0] for ele in (data_dict[9987.0][0])]
        col_names.append('step size 0.1')
        col_names.append('step size 0.01')
        col_names.append('step size 0.001')
        clean_df = pd.DataFrame(index = data_dict.keys(), columns = col_names)
        for r in range(len(clean_df.index)):
            step_0_1 = None
            step_0_01 = None
            step_0_001 = None
            row = clean_df.index[r]
            if len(data_dict[row]) > 1:
                param_vals = [ele.split(":")[1] for ele in (data_dict[row][0])]
                step_0_1 = data_dict[row][1]
                step_0_01 = data_dict[row][2]
                step_0_001 = data_dict[row][3]
            else:
                param_vals = [ele.split(":")[1] for ele in (data_dict[row][0])]
            for i in range(len(param_vals)):
                clean_df.iloc[r, i] = param_vals[i]
            try:
                clean_df.loc[row, 'step size 0.1'] = step_0_1
                clean_df.loc[row, 'step size 0.01'] = step_0_01
                clean_df.loc[row, 'step size 0.001'] = step_0_001
            except:
                clean_df.loc[row, 'step size 0.1'] = []
                clean_df.loc[row, 'step size 0.01'] = []
                clean_df.loc[row, 'step size 0.001'] = []
        clean_df['GO Term Name'] = go_names

        output_location = root_folder + "/" + model + '_GO_results_plus_uncertainties.csv'
        
        clean_df.to_csv(output_location, index_label='GO Term')
        

def separate_df_by_step_sizes(big_df_path, model_name):
    '''
    Given the path to the big dataframe, last modified in get_data_by_model(),
    and the model name, create new dataframes for each model with each different 
    step size. Ex. 1d_lognormal will have a file for 0.1, 0.01, and 0.001 step sizes.
    '''
    big_df = pd.read_csv(big_df_path, index_col='GO Term')
    df_0_1 = big_df.drop('step size 0.01', axis = 1).drop('step size 0.001', axis = 1)
    df_0_01 = big_df.drop('step size 0.1', axis = 1).drop('step size 0.001', axis = 1)
    df_0_001 = big_df.drop('step size 0.1', axis = 1).drop('step size 0.01', axis = 1)
    
    
    
    param_names = []
    for val in big_df.columns:
        if val not in ['LL', ' Converged','step size 0.1', 'step size 0.01', 'step size 0.001', 'GO Term Name', 'Unnamed: 0']:
            param_names.append(val.strip())
    for row in df_0_1.index:
        for param in param_names:
            try:
                df_0_1.loc[row, param + ' uncertainty'] = float(big_df.loc[row, 'step size 0.1'].split()[param_names.index(param)].strip('[').strip(']'))
            except:
                df_0_1.loc[row, param + ' uncertainty'] = None
    for row in df_0_01.index:
        for param in param_names:
            try:
                df_0_01.loc[row, param + ' uncertainty'] = float(big_df.loc[row, 'step size 0.01'].split()[param_names.index(param)].strip('[').strip(']'))
            except:
                df_0_01.loc[row, param + ' uncertainty'] = None
    for row in df_0_001.index:
        for param in param_names:
            try:
                df_0_001.loc[row, param + ' uncertainty'] = float(big_df.loc[row, 'step size 0.001'].split()[param_names.index(param)].strip('[').strip(']'))
            except:
                df_0_001.loc[row, param + ' uncertainty'] = None
    df_0_1 = df_0_1.drop('step size 0.1', axis = 1)
    df_0_01 = df_0_01.drop('step size 0.01', axis = 1)
    df_0_001 = df_0_001.drop('step size 0.001', axis = 1)
    
    output_location_0_1 = "/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/0.1_step_size/" + model_name + '_GO_results_plus_uncertainties_0.1_step.csv'
    output_location_0_01 = "/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/0.01_step_size/" + model_name + '_GO_results_plus_uncertainties_0.01_step.csv'
    output_location_0_001 = "/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/0.001_step_size/" + model_name + '_GO_results_plus_uncertainties_0.001_step.csv'
    
    df_0_1.to_csv(output_location_0_1)
    df_0_01.to_csv(output_location_0_01)
    df_0_001.to_csv(output_location_0_001)


if __name__ == "__main__":
    root_folder = "/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference"
    go_csv_path = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/go_term_DFE_inference_results.csv'
    output_csv_path = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/go_term_DFE_inference_results_output_all_models.csv'

    process_directories(root_folder, go_csv_path, output_csv_path)
    
    go_results_path = output_csv_path
    
    uncerts_output_csv_path = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/go_term_DFE_inference_results_output_all_models_and_uncertainties.csv'
    
    add_uncerts(root_folder, go_results_path, uncerts_output_csv_path)
    
    get_data_by_model(uncerts_output_csv_path, root_folder)
    
    path_1d_lognormal = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/1d_lognormal_DFE_GO_results_plus_uncertainties.csv'
    path_2d_biv_ind_gamma_asymmetric = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/2d_biv_ind_gamma_asymmetric_DFE_GO_results_plus_uncertainties.csv'
    path_2d_biv_ind_gamma_symmetric = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/2d_biv_ind_gamma_symmetric_DFE_GO_results_plus_uncertainties.csv'
    path_2d_biv_lognormal_asymmetric = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/2d_biv_lognormal_asymmetric_DFE_GO_results_plus_uncertainties.csv'
    path_2d_biv_lognormal = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/2d_biv_lognormal_DFE_GO_results_plus_uncertainties.csv'
    path_2d_gamma = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/uncertainties_by_model/2d_gamma_DFE_GO_results_plus_uncertainties.csv'

    separate_df_by_step_sizes(path_1d_lognormal, '1d lognormal')
    separate_df_by_step_sizes(path_2d_biv_ind_gamma_asymmetric, '2d_biv_ind_gamma_asymmetric')
    separate_df_by_step_sizes(path_2d_biv_ind_gamma_symmetric, '2d_biv_ind_gamma_symmetric')
    separate_df_by_step_sizes(path_2d_biv_lognormal_asymmetric, '2d_biv_lognormal_asymmetric')
    separate_df_by_step_sizes(path_2d_biv_lognormal, '2d_biv_lognormal')
    separate_df_by_step_sizes(path_2d_gamma, '2d_gamma')
    