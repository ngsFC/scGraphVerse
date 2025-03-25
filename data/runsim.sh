#!/bin/bash



# List of Rmd files to process

FILES=(

  "simulations.Rmd"

  "simulations.p100.Rmd"

  "simulations.p100.n500.Rmd"

  "simulations.p200.Rmd"

  "simulations.p200.n500.Rmd"

  "simulations.p500.n500.Rmd"

)



for file in "${FILES[@]}"; do

  echo "Processing $file ..."

  

  # Get the R script filename

  Rfile="${file%.Rmd}.R"



  # Extract R code from the Rmd file

  Rscript -e "knitr::purl('$file', output = '$Rfile')"



  # Run the extracted R script

  echo "Running $Rfile ..."

  Rscript "$Rfile"



  # Remove the generated R script

  rm -f "$Rfile"

  

  echo "Cleaned up $Rfile"

  echo "------------------------"

done


