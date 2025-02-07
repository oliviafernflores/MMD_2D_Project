#!/usr/bin/env python
#SBATCH --job-name=GO_DFE_process_results
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --time=24:00:00
#SBATCH --array=1

import os
import pandas as pd



# Function to extract information from the file name
def extract_info_from_filename(filename):
    # Extract the model name and the row locator (GO number) from the file name
    base_name = os.path.basename(filename)
    parts = base_name.split('_')
    model_name = "_".join(parts[4:-1])  # Everything between "mmd_FRA_mmd_IRA" and ".InferDFE.bestfits"
    go_number = float(parts[-1].split('.')[0])  # The number right before ".InferDFE.bestfits"
    
    return model_name, go_number

def process_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Look for the '# Top 100 results' marker
    start_index = None
    for i, line in enumerate(lines):
        if line.strip() == "# Top 100 results":
            start_index = i
            break

    if start_index is None:
        raise ValueError(f"Could not find '# Top 100 results' in {file_path}")
    
    # Check if '# Converged results' is present in the file
    converged = "no"
    for i, line in enumerate(lines):
        if line.strip() == "# Converged results":
            converged = "yes"
            break

    # The next line after the '# Top 100 results' line contains the column headers
    param_names_line = lines[start_index + 1].strip().split('\t')
    
    # The line after that contains the actual results
    result_line = lines[start_index + 2].strip().split('\t')

    
    params_df = pd.DataFrame([result_line], columns = param_names_line, index = None)

    for col in params_df.columns:
        params_df.loc[0, col] = round(float(params_df.loc[0, col]), 3)
    
    log_likelihood = params_df.iloc[0, 0]
    params_dict = params_df.to_dict('index')
    
    return log_likelihood, params_dict[0], param_names_line, converged


def process_directories(root_folder, go_csv_path, output_csv_path):
    # Read the GO terms file
    go_df = pd.read_csv(go_csv_path, index_col='GO Term')

    # Initialize a dictionary to store the best model information for each row locator
    best_fit_data = {}

    # Iterate through each directory in the root folder
    for subdir, _, files in os.walk(root_folder):
        for file in files:
            if file.startswith("mmd_FRA_mmd_IRA") and file.endswith(".InferDFE.bestfits"):
                # Get the full path to the file
                file_path = os.path.join(subdir, file)
                
                try:
                    # Process the file to extract log likelihood, parameters, and converged status
                    log_likelihood, params_dict, param_names, converged = process_file(file_path)
                    
                    # Skip files where processing failed (i.e., where we returned None)
                    if log_likelihood is None:
                        continue
                    
                    # Extract model name and GO number (row locator)
                    model_name, go_number = extract_info_from_filename(file)
                    
                    # print(f"Best fit so far for GO {go_number}: {best_fit_data.get(go_number)}")
                    
                    
                    # Check if this is the best fit for the given GO number
                    if go_number not in best_fit_data:
                        params_dict.pop("# Log(likelihood)")
                        params_dict.pop('theta')
                        best_fit_data[go_number] = {
                            'log_likelihood': round(float(log_likelihood), 3),
                            'model_name': model_name,
                            'params_dict': params_dict,
                            'converged': converged
                        }
                    else:
                        if float(log_likelihood) > (best_fit_data[go_number]['log_likelihood']):
                            params_dict.pop("# Log(likelihood)")
                            params_dict.pop('theta')
                            best_fit_data[go_number] = {
                                'log_likelihood': round(float(log_likelihood), 3),
                                'model_name': model_name,
                                'params_dict': params_dict,
                                'converged': converged
                            }
                except ValueError as e:
                    print(f"Error processing {file}: {e}")
                    continue

    # Add the best fit data to the GO terms DataFrame
    best_fit_df = pd.DataFrame.from_dict(best_fit_data, orient='index')
    best_fit_df.index.name = 'GO Term'
    # print(best_fit_df)
    
    # Merge with the existing GO DataFrame
    final_df = go_df.join(best_fit_df, how='left')

    # Save the final DataFrame to CSV
    final_df.to_csv(output_csv_path)





# Example usage
if __name__ == "__main__":
    # root_folder = "/groups/rgutenk/oliviafernflores/gene_ontology/GO_DFE_Inference"  # Replace with the path to your folder containing subdirectories
    # go_csv_path = '/groups/rgutenk/oliviafernflores/gene_ontology/GO_DFE_Inference/go_term_DFE_inference_results.csv'  # Replace with the path to your GO numbers CSV file
    # output_csv_path = '/groups/rgutenk/oliviafernflores/gene_ontology/GO_DFE_Inference/go_term_DFE_inference_results_output.csv'  # Path where the final CSV will be saved

    root_folder = "/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference"  # Replace with the path to your folder containing subdirectories
    go_csv_path = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/go_term_DFE_inference_results.csv'  # Replace with the path to your GO numbers CSV file
    output_csv_path = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_DFE_Inference/go_term_DFE_inference_results_output.csv'  # Path where the final CSV will be saved



    process_directories(root_folder, go_csv_path, output_csv_path)
    print("Process completed, and the data is saved to", output_csv_path)
