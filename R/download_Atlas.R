#' Download and Load an RDS File from a URL
#'
#' Downloads an RDS file from a specified URL (or loads from a local file URL) and reads its contents into R.
#'
#' @param file_url A character string specifying the URL of the RDS file to download.
#'
#' @return
#' An R object loaded from the downloaded RDS file.
#'
#' @details
#' Supports both \code{http(s)} URLs (download via \pkg{httr}) and \code{file} URLs (read from local file).
#' An internet connection is required for remote URLs.
#'
#' @importFrom httr GET status_code content timeout
#' @export
#'
#' @examples
#' # Create a temporary RDS file for demonstration
#' temp_file <- tempfile(fileext = ".rds")
#' saveRDS(mtcars, temp_file)
#' 
#' # Convert local path to a file URL
#' temp_url <- paste0("file://", temp_file)
#' 
#' # Load the RDS file using download_Atlas
#' loaded_data <- download_Atlas(temp_url)
#' 
#' # Check the loaded object
#' head(loaded_data)
download_Atlas <- function(file_url) {
  if (missing(file_url) || !is.character(file_url) || length(file_url) != 1 || file_url == "") {
    stop("'file_url' must be a non-empty character string.")
  }
  
  if (grepl("^file://", file_url)) {
    # Local file case
    local_path <- sub("^file://", "", file_url)
    if (!file.exists(local_path)) {
      stop("The specified local file does not exist: ", local_path)
    }
    return(readRDS(local_path))
  }
  
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("The 'httr' package is required but not installed.")
  }
  
  response <- httr::GET(file_url, httr::timeout(60))
  
  if (httr::status_code(response) != 200) {
    stop("Failed to download the data. HTTP status code: ", httr::status_code(response))
  }
  
  raw_data <- httr::content(response, as = "raw")
  readRDS(gzcon(rawConnection(raw_data)))
}
