from goatools import obo_parser

# Load the GO ontology
go_obo_file = "/Users/olivia/Desktop/go-basic.obo"  # Replace with the path to your GO OBO file
go = obo_parser.GODag(go_obo_file)

# Find the GO term for "biological process"
biological_process_term = "GO:0008150"  # The GO ID for biological process
biological_process = go[biological_process_term]

# Get direct children
direct_children = biological_process.children


print(go['GO:0048518'].parents)

# Prepare the output
results = []
missing_terms = []

for child in direct_children:
    child_id = child.id
    if child_id in go:
        child_term = go[child_id]
        results.append(f"{child_term.id}\t{child_term.name}")
    else:
        missing_terms.append(child_id)

# Save results to a text file
output_file = "direct_children_biological_process.txt"
with open(output_file, 'w') as f:
    f.write("GO Term\tDescription\n")  # Header
    f.write("\n".join(results))

print(f"Direct children of biological process saved to {output_file}.")

# Output missing terms for debugging
if missing_terms:
    print("Missing terms:", missing_terms)
