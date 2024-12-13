download_Atlas <- function(file_url, dest_file) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("The 'httr' package is required but not installed.")
  }
  
  httr::GET(file_url, httr::write_disk(dest_file, overwrite = TRUE))
  
  if (file.exists(dest_file)) {
    message("File downloaded successfully: ", dest_file)
  } else {
    stop("File download failed.")
  }
}
