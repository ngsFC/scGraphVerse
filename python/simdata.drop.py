# Apply dropout to simulate scRNA-seq-like data
def apply_dropout(expression_matrix, dropout_rate=0.1):
    dropout_mask = np.random.rand(*expression_matrix.shape) > dropout_rate
    return expression_matrix * dropout_mask

# Apply noise to the simulated gene expression
def add_noise(expression_matrix, noise_level=0.1):
    noise = np.random.normal(0, noise_level, expression_matrix.shape)
    return expression_matrix + noise

# Apply dropout and noise
expression_with_dropout = apply_dropout(simulated_expression, dropout_rate=0.3)
expression_noisy = add_noise(expression_with_dropout, noise_level=0.1)

# Convert to pandas dataframe for easier manipulation
expression_noisy_df = pd.DataFrame(expression_noisy, columns=nodes)

# Save the noisy expression data
expression_noisy_df.to_csv("simulated_expression_tcell_with_dropout.csv", index=False)

