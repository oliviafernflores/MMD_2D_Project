#!/usr/bin/env python
#SBATCH --job-name=sort_GO_dictionaries
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32000

import os
import csv
import pickle
import argparse

# Function to read GO terms from the specified file
def read_go_terms(filename):
    go_terms = set()
    with open(filename) as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                go_terms.add(parts[0])
    return go_terms

# Set up argument parser to specify which GO file to process
def parse_args():
    parser = argparse.ArgumentParser(description="Process a specific GO file and its descendants.")
    parser.add_argument('go_file', help="The GO file to process (e.g., GO_00000122.txt)")
    return parser.parse_args()

# Main processing function
def process_go_file(go_file):
    # Read the GO terms from the specified file
    go_terms_from_file = read_go_terms(go_file)
    base_go_term = os.path.basename(go_file).split('_')[-1].split('.')[0]
    go_terms_from_file.add(base_go_term)

    # Output CSV for the summary of this specific GO file
    summary_file = f'go_term_summary_{base_go_term}.csv'
    
    # Open the summary CSV file to write results incrementally
    with open(summary_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['GO Term', 'Number of Entries'])

        # Open the pickle file incrementally
        with open('/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/mus_all_pops_nsyn_with_ancestral_and_GO_terms.pkl', 'rb') as pickle_file:
            big_dict = pickle.load(pickle_file)

            # Process each GO term one by one to minimize memory usage
            for term in go_terms_from_file:
                # Prepare the dictionary for this term
                term_entries = {}

                # Iterate over `big_dict` items, process and write entries as we go
                for key, value in big_dict.items():
                    # Check if the specific term is in the 'GO_terms' of the current big_dict entry
                    if term in value['GO_terms']:
                        term_entries[key] = value


                # Save the dictionary for each term
                term_dict_path = f'/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_term_dictionaries_descendants/{term}/{term}_dict.pkl'
                os.makedirs(os.path.dirname(term_dict_path), exist_ok=True)
                
                with open(term_dict_path, 'wb') as term_pickle_file:
                    pickle.dump(term_entries, term_pickle_file)

                # Write the summary incrementally
                writer.writerow([term, len(term_entries)])

                print(f"Processed GO term: {term} - {len(term_entries)} entries. Dictionary saved.")

    print(f"Summary CSV file '{summary_file}' created.")

# Entry point of the script
if __name__ == '__main__':
    # Parse command-line arguments
    args = parse_args()

    # Process the specific GO file provided as an argument
    process_go_file(args.go_file)
