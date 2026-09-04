# Download and Load an RDS File from a URL

Downloads an RDS file from a specified URL and reads its contents into
R. We used it for https://www.singlecellatlas.org

## Usage

``` r
download_Atlas(file_url)
```

## Arguments

- file_url:

  Character; URL of the RDS file to download.

## Value

An R object loaded from the downloaded RDS file.

## Details

This function uses httr to perform the download. The RDS file is read
directly from a raw connection without saving to disk. An internet
connection is required.

If the download fails (e.g., invalid URL, server error), an informative
error message is returned.

## Examples

``` r
url <- "https://zenodo.org/records/15511027/files/sce_obj.rds?download=1"
atlas_data <- download_Atlas(url)
```
