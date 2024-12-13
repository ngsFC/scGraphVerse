download_Atlas <- function(file_url) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("The 'httr' package is required but not installed.")
  }
  
  response <- httr::GET(file_url)
  
  if (httr::status_code(response) != 200) {
    stop("Failed to download the data. HTTP status code: ", httr::status_code(response))
  }
  
  raw_data <- httr::content(response, as = "raw")
  readRDS(gzcon(rawConnection(raw_data)))
}
