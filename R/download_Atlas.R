#' Download and Load an RDS File from a URL
#'
#' Downloads an RDS file from a specified URL and reads its contents into R.
#' We used it for https://www.singlecellatlas.org
#'
#' @param file_url Character; URL of the RDS file to download.
#'
#' @return
#' An R object loaded from the downloaded RDS file.
#'
#' @details
#' This function uses \pkg{httr} to perform the download. The RDS file is read
#' directly from a raw connection without saving to disk. An internet connection
#' is required.
#'
#' If the download fails (e.g., invalid URL, server error), an informative error
#' message is returned.
#'
#' @importFrom httr GET status_code content
#' @export
#'
#' @examples
#' url <- paste0(
#'     "https://www.dropbox.com/s/r8qwsng79rhp9gf/",
#'     "SCA_scRNASEQ_TISSUE_WHOLE_BLOOD.RDS?dl=1"
#' )
#' atlas_data <- download_Atlas(url)
download_Atlas <- function(file_url) {
    if (!requireNamespace("httr", quietly = TRUE)) {
        stop("'httr' is required but not installed.")
    }

    response <- httr::GET(file_url)
    if (httr::status_code(response) != 200) {
        stop(
            "Failed to download data. HTTP status code: ",
            httr::status_code(response)
        )
    }

    raw_data <- httr::content(response, as = "raw")
    readRDS(gzcon(rawConnection(raw_data)))
}
