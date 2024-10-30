import matplotlib.pyplot as plt
import networkx as nx
from goatools import obo_parser

# Load the GO ontology
go_obo_file = "/Users/olivia/Desktop/go-basic.obo"  # Replace with your path
go = obo_parser.GODag(go_obo_file)

# Find the GO term for "biological process"
biological_process_term = "GO:0008150"  # The GO ID for biological process
biological_process = go.get(biological_process_term)

# Create a directed graph
G = nx.DiGraph()

# Function to recursively add edges to the graph with a depth limit
def add_edges(term, current_depth, max_depth=1):
    if current_depth > max_depth:
        return
    for child_id in term.children:
        child_id_str = str(child_id).split()[0]  # Get only the GO ID
        if child_id_str in go:
            child_term = go[child_id_str]
            G.add_edge(term.id, child_term.id)
            add_edges(child_term, current_depth + 1, max_depth)

# Start adding edges from the biological process term
add_edges(biological_process, current_depth=0)

# Debugging: Print number of nodes and edges
print(f"Number of nodes: {len(G.nodes)}")
print(f"Number of edges: {len(G.edges)}")

# Draw the graph only if it has nodes
if len(G.nodes) > 0:
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, k=0.5, iterations=50)  # You can reduce iterations for faster layout
    nx.draw_networkx_nodes(G, pos, node_size=100, node_color='lightblue', alpha=0.6)
    nx.draw_networkx_edges(G, pos, arrowstyle='-|>', arrowsize=10, alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8, font_family="sans-serif")
    
    # Set title and show the plot
    plt.title("GO Hierarchy for Biological Process (Limited Depth)")
    plt.axis('off')  # Turn off the axis
    plt.tight_layout()  # Adjust layout to fit everything
    plt.show()
else:
    print("No nodes to display in the graph.")

# plt.savefig('biological_process_GO.png')
