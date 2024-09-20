import pandas as pd
from GENIE3 import GENIE3, get_link_list  # Assuming you have GENIE3.py in your folder

# Step 1: Load the simulated scRNA-seq count matrix
def load_simulated_matrix(file):
    # Load the matrix where genes are columns and cells are rows
    return pd.read_csv(file, index_col=0)

# Step 2: Use GENIE3 to reconstruct the gene regulatory network
def reconstruct_network_genie3(simulated_matrix, nthreads=8):
    # Transpose the simulated matrix to fit GENIE3 input (genes as rows, cells as columns)
    print("Simulated matrix shape (before transpose):", simulated_matrix.shape)
    transposed_matrix = simulated_matrix.T
    print("Simulated matrix shape (after transpose):", transposed_matrix.shape)

    # Run GENIE3 with the transposed matrix (genes as rows, cells as columns)
    inferred_network = GENIE3(transposed_matrix.values, nthreads=nthreads)  # Use 8 threads
    
    if inferred_network is None:
        print("GENIE3 returned None, check your input matrix.")
        return None
    
    print("GENIE3 inference completed.")
    
    # Get the list of gene interactions from the inferred network
    link_list = get_link_list(inferred_network)
    
    if link_list is None:
        print("get_link_list returned None, check your GENIE3 results.")
        return None
    
    return link_list

# Step 3: Save the inferred gene regulatory network to a file
def save_inferred_network(link_list, output_file='inferred_network.txt'):
    # Convert the link list to a DataFrame and save it
    inferred_df = pd.DataFrame(link_list, columns=['Gene 1', 'Gene 2', 'Weight'])
    inferred_df.to_csv(output_file, sep='\t', index=False)
    print(f'Inferred network saved to {output_file}')

# Main function to run the pipeline
def main(simulated_file, output_file, nthreads=8):
    # Step 1: Load the simulated scRNA-seq matrix
    simulated_matrix = load_simulated_matrix(simulated_file)
    
    # Step 2: Reconstruct the gene regulatory network using GENIE3 with multiple threads
    reconstructed_network = reconstruct_network_genie3(simulated_matrix, nthreads=nthreads)
    
    if reconstructed_network is None:
        print("Error: Reconstructed network is None. Exiting.")
        return
    
    # Step 3: Save the inferred network to a file
    save_inferred_network(reconstructed_network, output_file)

# Example usage
if __name__ == "__main__":
    # Replace with your actual file paths
    simulated_matrix_file = 'simulated_scRNAseq_counts.csv'
    output_file = 'inferred_network.txt'
    
    # Run the network reconstruction with GENIE3
    main(simulated_matrix_file, output_file, nthreads=8)

