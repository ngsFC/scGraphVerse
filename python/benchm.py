import pandas as pd
import numpy as np
from GENIE3 import GENIE3, get_link_list
import networkx as nx
from sklearn.metrics import precision_recall_curve, f1_score

# Step 1: Load the original GRN from interaction data
def load_grn_from_file(interaction_file):
    # Load interaction file
    data = pd.read_csv(interaction_file, sep='\t')  # Adjust separator if needed (e.g., ',' for CSV)
    
    # Create a graph and add edges based on the interaction data
    G = nx.Graph()
    
    for index, row in data.iterrows():
        gene1 = row['Gene 1'].strip().lower()  # Normalize gene names to lowercase
        gene2 = row['Gene 2'].strip().lower()  # Normalize gene names to lowercase
        weight = row['Weight'] if 'Weight' in row else np.random.uniform(-1, 1)  # Assign weights if available
        G.add_edge(gene1, gene2, weight=weight)
    
    return G

# Step 2: Apply GENIE3 to infer a GRN using 8 threads
def infer_grn_with_genie3(expression_data, gene_names):
    # Debugging: Check the shape of the input data
    print("Shape of expression data:", expression_data.shape)
    print("Sample of expression data:")
    print(expression_data[:5, :5])  # Check a sample of the first few rows and columns
    
    # Run GENIE3 to infer the regulatory network from gene expression data using 8 threads
    inferred_network = GENIE3(expression_data, nthreads=8)
    
    # Debugging: Check if the inferred network is valid
    if inferred_network is None:
        print("GENIE3 did not return a valid inferred network.")
        return None
    
    print("Inferred network matrix shape:", inferred_network.shape)
    
    # Convert the list of edges with their weights into a pandas DataFrame
    links = get_link_list(inferred_network)
    
    # Debugging: Check if links are valid
    if links is None:
        print("No links were inferred from the network.")
        return None
    
    print("Sample inferred edges:", links[:5])
    
    # Replace the generic identifiers with the original gene names using the gene_names list
    for i in range(len(links)):
        links[i][0] = gene_names[int(links[i][0][1:]) - 1]  # Regulator (convert Gx to gene name)
        links[i][1] = gene_names[int(links[i][1][1:]) - 1]  # Target (convert Gx to gene name)
    
    # Convert the inferred network into an edge list for easier comparison
    inferred_edges = pd.DataFrame(links, columns=["Regulator", "Target", "Weight"])
    
    # Normalize gene names to lowercase for consistent comparison
    inferred_edges['Regulator'] = inferred_edges['Regulator'].str.strip().str.lower()
    inferred_edges['Target'] = inferred_edges['Target'].str.strip().str.lower()
    
    return inferred_edges

def run_benchmark(expression_file, original_grn_file, output_file="inferred_interactions.csv"):
    # Load the expression data (cells x genes matrix) and normalize gene names
    expression_data = pd.read_csv(expression_file, index_col=0)
    
    # Debugging: Check the first 5 column (gene) names and row (cell) names
    print("First 5 gene names (columns):", expression_data.columns[:5])
    print("First 5 sample names (rows):", expression_data.index[:5])
    
    # Get the gene names from the columns and normalize them
    gene_names = expression_data.columns.str.strip().str.lower().tolist()
    
    # Convert the dataframe to a NumPy array for GENIE3
    expression_data = expression_data.values
    
    # Check for NaNs or zero rows/columns
    if np.isnan(expression_data).any():
        print("Warning: Expression data contains NaNs. This may cause issues.")
    if not np.any(expression_data):
        print("Warning: Expression data contains all zeros or is empty.")
    
    # Load the original GRN
    original_grn = load_grn_from_file(original_grn_file)
    
    # Infer GRN with GENIE3 using 8 threads and map back to original gene names
    inferred_edges = infer_grn_with_genie3(expression_data, gene_names)
    
    if inferred_edges is None:
        print("No valid edges inferred.")
        return
    
    # Save inferred edges to a CSV file
    save_inferred_edges(inferred_edges, output_file)
    
    # Evaluate the inferred GRN compared to the original one
    true_labels, predicted_scores = evaluate_inferred_grn(inferred_edges, original_grn)
    
    if true_labels and predicted_scores:
        # Calculate benchmark metrics
        precision, recall, f1 = calculate_benchmark_metrics(true_labels, predicted_scores)
        
        # Print the results
        print("Precision-Recall Curve: Precision =", precision, "Recall =", recall)
        print("F1-Score:", f1)
    else:
        print("No valid true labels or predicted scores to compute metrics.")

# Example usage:
run_benchmark('simulated_scrna_seq.csv', './../data/genemania-interactions.txt', 'inferred_interactions_with_names.csv')

