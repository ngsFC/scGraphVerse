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
#' data(count_matrices)
#' genes <- selgene(
#'     object = count_matrices[[1]],
#'     top_n = 20,
#'     cell_type = "T_cells",
#'     cell_type_col = "CELL_TYPE",
#'     remove_rib = TRUE,
#'     remove_mt = TRUE,
#'     assay = "counts"
#' )
#' 
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

    # Alternative approach: Try to use a direct API method without STRINGdb initialization
    if (verbose) {
        message("Attempting direct STRING API access to bypass SSL initialization issues...")
    }
    
    # Try direct API call first
    direct_result <- tryCatch({
        .direct_string_api_call(genes, species, required_score, verbose, keep_all_genes)
    }, error = function(e) {
        if (verbose) message("Direct API approach failed, trying STRINGdb initialization...")
        NULL
    })
    
    if (!is.null(direct_result)) {
        return(direct_result)
    }
    
    # Fallback to STRINGdb initialization with SSL workarounds
    old_method <- getOption("download.file.method")
    old_extra <- getOption("download.file.extra")
    old_curl_bundle <- Sys.getenv("CURL_CA_BUNDLE")
    old_curl_opts <- getOption("RCurlOptions")
    
    on.exit({
        options(download.file.method = old_method)
        options(download.file.extra = old_extra)
        Sys.setenv(CURL_CA_BUNDLE = old_curl_bundle)
        options(RCurlOptions = old_curl_opts)
    })
    
    # Configure SSL bypass methods
    options(download.file.method = "libcurl")
    options(download.file.extra = c("-k", "--insecure"))
    Sys.setenv(CURL_CA_BUNDLE = "")
    Sys.setenv(CURL_INSECURE = "1")
    options(RCurlOptions = list(ssl.verifypeer = FALSE, ssl.verifyhost = FALSE))
    
    string_db <- tryCatch({
        STRINGdb$new(
            version         = "11.5",
            species         = species,
            score_threshold = required_score,
            input_directory = ""
        )
    }, error = function(e1) {
        if (verbose) message("LibCurl failed, trying alternative methods...")
        
        # Try with internal method override
        Sys.setenv(DOWNLOAD_FILE_METHOD = "internal")
        options(download.file.method = "internal")
        
        tryCatch({
            STRINGdb$new(
                version         = "11.5",
                species         = species,
                score_threshold = required_score,
                input_directory = ""
            )
        }, error = function(e2) {
            if (verbose) {
                message("STRINGdb initialization failed due to SSL certificate issues.")
                message("Error: ", e1$message)
                message("This is a known issue with some network configurations.")
                message("Workaround: Use pre-downloaded STRING data or alternative networks.")
            }
            stop("Unable to initialize STRINGdb. SSL certificate verification failed. ",
                 "Consider using cached data or alternative protein interaction databases.")
        })
    })

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

# Direct STRING API call that bypasses STRINGdb R package initialization
.direct_string_api_call <- function(genes, species, required_score, verbose, keep_all_genes = TRUE) {
    if (verbose) {
        message("Using direct STRING API approach...")
    }
    
    # First, try to get STRING IDs using direct API
    mapping_url <- "https://string-db.org/api/json/get_string_ids"
    
    # Prepare gene list for POST request
    gene_list <- paste(genes, collapse = "\n")
    
    # Try to get STRING IDs
    string_ids <- tryCatch({
        if (verbose) message("Mapping genes to STRING IDs via direct API...")
        
        response <- httr::POST(
            url = mapping_url,
            body = list(
                identifiers = gene_list,
                species = species,
                limit = 1,
                echo_query = 1,
                caller_identity = "scGraphVerse"
            ),
            encode = "form",
            httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE)
        )
        
        if (httr::status_code(response) != 200) {
            stop("STRING API request failed with status: ", httr::status_code(response))
        }
        
        result <- httr::content(response, as = "text", encoding = "UTF-8")
        jsonlite::fromJSON(result)
    }, error = function(e) {
        if (verbose) message("Direct STRING ID mapping failed: ", e$message)
        return(NULL)
    })
    
    if (is.null(string_ids) || nrow(string_ids) == 0) {
        if (verbose) message("No STRING IDs obtained via direct API")
        return(NULL)
    }
    
    # Get interactions using direct API
    interactions <- tryCatch({
        if (verbose) message("Retrieving interactions via direct API...")
        
        string_protein_ids <- unique(string_ids$stringId)
        protein_list <- paste(string_protein_ids, collapse = "%0d")
        
        interaction_url <- "https://string-db.org/api/json/network"
        
        response <- httr::POST(
            url = interaction_url,
            body = list(
                identifiers = protein_list,
                species = species,
                required_score = required_score,
                add_white_nodes = 0,
                caller_identity = "scGraphVerse"
            ),
            encode = "form",
            httr::config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE)
        )
        
        if (httr::status_code(response) != 200) {
            stop("STRING interaction API request failed with status: ", httr::status_code(response))
        }
        
        result <- httr::content(response, as = "text", encoding = "UTF-8")
        jsonlite::fromJSON(result)
    }, error = function(e) {
        if (verbose) message("Direct STRING interaction retrieval failed: ", e$message)
        return(NULL)
    })
    
    if (is.null(interactions) || nrow(interactions) == 0) {
        if (verbose) message("No interactions obtained via direct API")
        return(.zero_matrix_result(genes))
    }
    
    # Convert STRING IDs back to gene names and build matrices
    # Create mapping from STRING ID to gene name
    id_to_gene <- setNames(string_ids$queryItem, string_ids$stringId)
    
    # Map interaction partners to gene names
    interactions$gene1 <- id_to_gene[interactions$stringId_A]
    interactions$gene2 <- id_to_gene[interactions$stringId_B]
    
    # Remove rows where genes couldn't be mapped back
    interactions <- interactions[!is.na(interactions$gene1) & !is.na(interactions$gene2), ]
    
    if (nrow(interactions) == 0) {
        if (verbose) message("No valid gene mappings for interactions")
        return(.zero_matrix_result(genes))
    }
    
    # Build adjacency matrices directly
    all_genes <- unique(c(interactions$gene1, interactions$gene2))
    if (keep_all_genes) {
        missing_genes <- setdiff(genes, all_genes)
        all_genes <- c(all_genes, missing_genes)
    }
    all_genes <- sort(all_genes)
    
    n_genes <- length(all_genes)
    weighted_adj <- matrix(0, nrow = n_genes, ncol = n_genes)
    rownames(weighted_adj) <- colnames(weighted_adj) <- all_genes
    
    # Fill in the weights (use combined_score)
    for (i in seq_len(nrow(interactions))) {
        g1 <- interactions$gene1[i]
        g2 <- interactions$gene2[i]
        score <- as.numeric(interactions$score[i])
        
        if (g1 %in% all_genes && g2 %in% all_genes) {
            weighted_adj[g1, g2] <- score
            weighted_adj[g2, g1] <- score  # Make symmetric
        }
    }
    
    # Create binary matrix
    binary_adj <- (weighted_adj > 0) * 1
    
    if (verbose) {
        message("Direct API approach successful: ", nrow(interactions), " interactions found")
    }
    
    return(list(
        weighted = weighted_adj,
        binary = binary_adj
    ))
}
