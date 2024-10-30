import pandas as pd
from goatools import obo_parser

# Load the GO ontology
go_obo_file = "/Users/olivia/Desktop/go-basic.obo"  # Replace with the path to your GO OBO file
go = obo_parser.GODag(go_obo_file)


# Read the input text file
go_terms = pd.read_csv("direct_biological_process_terms.txt", sep="\t")

# Function to get all descendants of a given GO term
def get_all_descendants(go_id):
    descendants = set()
    queue = [go_id]  # Initialize the queue with the specified GO term

    while queue:
        current_term = queue.pop(0)

        if current_term in go:
            current_children = go[current_term].children
            for child in current_children:
                descendants.add(child.id)
                queue.append(child.id)

    return descendants

# Loop through each GO term in the input file
for index, row in go_terms.iterrows():
    goid = row['GOID']
    all_descendants = get_all_descendants(goid)

    # Create a DataFrame for the descendants
    descendant_terms = [(descendant_id, go[descendant_id].name) for descendant_id in all_descendants if descendant_id in go]
    descendant_df = pd.DataFrame(descendant_terms, columns=['GOID', 'TERM'])

    # Save to a new text file named after the GO term
    output_file = f"all_descendants_of_{goid.replace(':', '_')}.txt"
    descendant_df.to_csv(output_file, sep='\t', index=False)

    print(f"File created successfully for {goid}")

print("All files created.")
