#!/bin/bash

# List of R scripts
scripts=(
  "simp100n500.R"
  "simp100.R"
  "simp200n500.R"
  "simp200.R"
  "simp500n500.R"
  "simp500.R"
)

# Loop through and run each script
for script in "${scripts[@]}"; do
  echo "Running $script..."
  Rscript "$script"
  echo "$script completed."
  echo "-----------------------------"
done

echo "All scripts have been executed."

