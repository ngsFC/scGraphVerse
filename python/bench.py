import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from arboreto.algo import grnboost2
import sys
sys.path.append("GENIE3.py")
from GENIE3 import GENIE3
from scipy.stats import spearmanr

def run_grnboost(expression_matrix):
    """Run GRNBoost2 to infer gene regulatory network."""
    return grnboost2(expression_matrix)

def run_genie3(expression_matrix):
    """Run GENIE3 to infer gene regulatory network."""
    inferred_network = GENIE3(expression_matrix)
    return inferred_network

def compare_networks(grnboost_df, genie3_df):
    """Compare GRNBoost2 and GENIE3 results."""
    merged_df = pd.merge(grnboost_df, genie3_df, on=['regulator', 'target'], suffixes=('_grnboost', '_genie3'))
    
    # Calculate Spearman correlation between the interaction strengths
    spearman_corr, p_value = spearmanr(merged_df['importance_grnboost'], merged_df['importance_genie3'])
    
    print(f"Spearman correlation between GRNBoost2 and GENIE3 networks: {spearman_corr:.4f}")
    
    # Plot a scatter plot comparing the two networks
    plt.figure(figsize=(8, 6))
    plt.scatter(merged_df['importance_grnboost'], merged_df['importance_genie3'], alpha=0.5)
    plt.xlabel('GRNBoost2 Importance')
    plt.ylabel('GENIE3 Importance')
    plt.title('Comparison of GRNBoost2 and GENIE3 Networks')
    plt.show()

if __name__ == '__main__':
    # Step 1: Load the simulated gene expression data
    expression_noisy_df = pd.read_csv("simulated_expression_tcell_with_dropout.csv", index_col=0)  # Read the file with cell names as row indices

    # Prepare the expression data (transpose to fit GRNBoost2 and GENIE3 format)
    expression_matrix = expression_noisy_df.T  # Genes x cells format
    genes = list(expression_matrix.index)

    # Step 2: Run GRNBoost2
    print("Running GRNBoost2...")
    grnboost_network = run_grnboost(expression_matrix)

    # Save GRNBoost2 results
    grnboost_network_df = pd.DataFrame(grnboost_network, columns=['regulator', 'target', 'importance'])
    grnboost_network_df.to_csv("inferred_network_grnboost2.csv", index=False)

    # Step 3: Run GENIE3
    print("Running GENIE3...")
    genie3_network = run_genie3(expression_matrix.values)

    # Convert GENIE3 results into a dataframe
    genie3_network_df = pd.DataFrame(genie3_network, columns=genes, index=genes)
    genie3_network_df = genie3_network_df.stack().reset_index()
    genie3_network_df.columns = ['regulator', 'target', 'importance']
    genie3_network_df = genie3_network_df[genie3_network_df['importance'] > 0]  # Filter out zero interactions

    # Save GENIE3 results
    genie3_network_df.to_csv("inferred_network_genie3.csv", index=False)

    # Step 4: Compare the results
    compare_networks(grnboost_network_df, genie3_network_df)
