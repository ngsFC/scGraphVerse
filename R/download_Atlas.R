#' Download and Load an RDS File from a URL
#'
#' Downloads an RDS file from a specified URL and reads its contents into R.
#'
#' @param file_url A character string specifying the URL of the RDS file to download.
#'
#' @return
#' An R object loaded from the downloaded RDS file.
#'
#' @details
#' This function uses \pkg{httr} to perform the file download. The RDS file is read directly
#' from the raw connection without saving it to disk. An internet connection is required.
#'
#' If the download fails (e.g., invalid URL, server error), an informative error message is returned.
#'
#' @importFrom httr GET status_code content
#' @export
#'
#' @examples
#' \dontrun{
#' # Example: Load an RDS file from a public URL
#' atlas_data <- download_Atlas("https://example.com/mydata.rds")
#' }
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
