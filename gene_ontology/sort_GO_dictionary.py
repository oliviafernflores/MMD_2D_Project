#!/usr/bin/env python
#SBATCH --job-name=sort_GO_dictionaries
#SBATCH --output=hpc_outfiles/%x-%j.out
#SBATCH --error=hpc_outfiles/%x-%j.err
#SBATCH --account=rgutenk
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_rgutenk
#SBATCH --nodes=10
#SBATCH --ntasks=50
#SBATCH --time=24:00:00

import os
import glob
import csv
import pickle


big_dict = pickle.load(open('/groups/rgutenk/oliviafernflores/gene_ontology/mus_all_pops_nsyn_with_ancestral_and_GO_terms.pkl', 'rb'))


def read_go_terms(filename):
    go_terms = set()
    with open(filename) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                go_terms.add(parts[0])
    return go_terms

directory_path = '/groups/rgutenk/oliviafernflores/gene_ontology/GO_term_descendants' 
go_files = glob.glob(os.path.join(directory_path, "all_descendants_of_GO_*.txt"))

go_dicts = {}


for go_file in go_files:
    go_terms_from_file = read_go_terms(go_file)
    base_go_term = os.path.basename(go_file).split('_')[-1].split('.')[0]
    go_terms_from_file.add(base_go_term)

    go_dicts[base_go_term] = {}

    for key, value in big_dict.items():
        if any(go_term in go_terms_from_file for go_term in value['GO_terms']):
            go_dicts[base_go_term][key] = value

summary_file = 'go_term_summary.csv'
with open(summary_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['GO Term', 'Number of Entries'])
    
    for go_term, entries in go_dicts.items():
        writer.writerow([go_term, len(entries)]) 
        
        with open(f'/groups/rgutenk/oliviafernflores/gene_ontology/GO_term_dictionaries/{go_term}_dict.pkl', 'wb') as pickle_file:
            pickle.dump(entries, pickle_file)

print(f"Summary CSV file '{summary_file}' created.")
for go_term in go_dicts.keys():
    print(f"Dictionary for {go_term} saved as '{go_term}_dict.pkl'.")
