import os
import pickle
import numpy as np
import dadi
import matplotlib.pyplot as plt
import pandas as pd

data_folder = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_term_descendants'
summary_files = [f for f in os.listdir(data_folder) if f.startswith('go_term_summary')]

df = pd.DataFrame(columns = ['GO Term', 'Number of Entries'])

df_lst = []

for file in summary_files:
    print(file)
    path = '/Users/olivia/Documents/2D_demographics_DFE/MMD_2D_Project/gene_ontology/GO_term_descendants/' + file
    df_temp = pd.read_csv(path)
    df_lst.append(df_temp)


df = pd.concat(df_lst)
df = df.set_index('GO Term')

# Define the path to the ontology file
ontology_file_path = '/Users/olivia/Desktop/go-basic.obo'

# Initialize an empty dictionary to hold GO terms and their names
go_terms_dict = {}

# Read the ontology file
with open(ontology_file_path, 'r') as f:
    current_go_id = None
    for line in f:
        line = line.strip()
        if line.startswith("id:"):
            current_go_id = line.split(" ")[1]
        elif line.startswith("name:") and current_go_id:
            go_name = line.split(":")[1]
            go_terms_dict[current_go_id] = go_name

# Convert the dictionary to a DataFrame
go_terms_df = pd.DataFrame(list(go_terms_dict.items()), columns=['GO Term', 'Name']).set_index('GO Term')

# Now merge this DataFrame with your existing DataFrame
merged_df = df.join(go_terms_df, how='left')
merged_df = merged_df[~merged_df.index.duplicated(keep='first')]

# Sort the merged DataFrame by 'Number of Entries' in descending order
sorted_df = merged_df.sort_values(by='Number of Entries', ascending=False)

# Display the sorted DataFrame
print(sorted_df)

# Save the sorted DataFrame to a CSV file
sorted_df.to_csv('GO_descendants_with_names_sorted.csv')

