#' Download and Load an RDS File from a URL
#'
#' This function downloads an RDS file from a given URL and reads it into R.
#'
#' @param file_url A character string specifying the URL of the RDS file.
#' @return The contents of the RDS file as an R object.
#' @examples
#' \dontrun{
#' data <- download_Atlas("https://example.com/data.rds")
#' }
#' @export
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
