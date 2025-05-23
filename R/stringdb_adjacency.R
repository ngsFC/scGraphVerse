#' Build Adjacency Matrices for Physical Interactions from STRING (POST API)
#'
#' Constructs weighted and binary adj matrices for physical protein-protein
#' interactions using a POST request to the STRING database API.
#'
#' @param genes A character vector of gene symbols or identifiers, e.g.,
#'   \code{c("TP53", "BRCA1", ...)}.
#' @param species Integer. NCBI taxonomy ID of the species. Default is
#'   \code{9606} (human).
#' @param required_score Integer in [0,1000]. Minimum confidence score for
#'   interactions. Default is \code{400}.
#' @param keep_all_genes Logical. If \code{TRUE} (default), includes all input
#'   genes in the final matrix even if unmapped.
#' @param verbose Logical. If \code{TRUE}, displays progress messages. Default
#'   is \code{TRUE}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{weighted}: A square numeric adjacency matrix with scores as
#'       weights.
#'     \item \code{binary}: A corresponding binary (0/1) adjacency matrix.
#'   }
#'
#' @details This function:
#'   \enumerate{
#'     \item Maps input genes to STRING internal IDs.
#'     \item Uses a POST request to retrieve physical protein-protein
#'       interactions from STRING.
#'     \item Builds a weighted adjacency matrix using the STRING combined score.
#'     \item Builds a binary adjacency matrix indicating presence/absence.
#'   }
#'
#' Genes not mapped to STRING are optionally retained as zero rows/columns if
#' \code{keep_all_genes = TRUE}.
#'
#' @note Requires packages: \pkg{STRINGdb}, \pkg{httr}, \pkg{jsonlite}.
#'
#' @importFrom STRINGdb STRINGdb
#' @importFrom httr POST content
#' @importFrom jsonlite fromJSON
#' @export
#'
#' @examples
#' url <- paste0(
#'     "https://www.dropbox.com/s/r8qwsng79rhp9gf/",
#'     "SCA_scRNASEQ_TISSUE_WHOLE_BLOOD.RDS?dl=1"
#' )
#' seu <- download_Atlas(file_url = url)
#' genes <- selgene(
#'     object = seu,
#'     top_n = 50,
#'     cell_type = "T_cells",
#'     cell_type_col = "CELL_TYPE",
#'     remove_rib = TRUE,
#'     remove_mt = TRUE
#' )
#' str_res <- stringdb_adjacency(
#'     genes = genes,
#'     species = 9606,
#'     required_score = 900,
#'     keep_all_genes = FALSE
#' )
#' wadj_truth <- str_res$weighted
#' adj_truth <- str_res$binary
stringdb_adjacency <- function(
    genes,
    species = 9606,
    required_score = 400,
    keep_all_genes = TRUE,
    verbose = TRUE) {
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
        stop("Package 'STRINGdb' is required. Install via Bioconductor.")
    }
    if (!requireNamespace("httr", quietly = TRUE)) {
        stop("Package 'httr' is required. Please install it.")
    }
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required. Please install it.")
    }
    if (length(genes) == 0) {
        stop("Please provide at least one gene in 'genes'.")
    }
    if (verbose) {
        message("Initializing STRINGdb...")
    }

    string_db <- STRINGdb$new(
        version         = "11.5",
        species         = species,
        score_threshold = required_score,
        input_directory = ""
    )

    if (verbose) {
        message("Mapping genes to STRING IDs...")
    }
    mapping <- .map_genes_to_string(string_db, genes)
    mapped_genes <- mapping$mapped
    unmapped_genes <- mapping$unmapped

    if (verbose) {
        message(
            "Mapped ", nrow(mapped_genes),
            " genes to STRING IDs."
        )
        if (length(unmapped_genes) > 0 && keep_all_genes) {
            message(
                length(unmapped_genes),
                " genes not found in STRING; included as zero rows/cols."
            )
        }
    }
    if (nrow(mapped_genes) == 0) {
        stop("No valid STRING IDs found for the provided genes.")
    }
    if (verbose) {
        message("Retrieving physical interactions from STRING API...")
    }
    interactions <- .query_string_api(
        mapped_genes$STRING_id,
        species,
        required_score
    )
    if (!is.data.frame(interactions) || nrow(interactions) == 0) {
        if (verbose) {
            message("No STRING physical interactions found.")
        }
        return(.zero_matrix_result(genes))
    }
    if (verbose) {
        message(
            "Found ", nrow(interactions),
            " STRING physical interactions."
        )
    }

    matrices <- .build_adjacency_matrices(
        interactions,
        mapped_genes,
        genes,
        keep_all_genes
    )
    if (verbose) {
        message("Adjacency matrices constructed successfully.")
    }

    return(matrices)
}
