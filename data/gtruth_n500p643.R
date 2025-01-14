# ---- Load Library ----
setwd("/home/francescoc/Desktop/GRN_project/data")

library(tidyverse)
library(biomaRt)
library(igraph)
library(learn2count)
library(ggraph)

source("./../R/download_Atlas.R")
source("./../R/exploreCells.R")
source("./../R/pathg.R")
source("./../R/BioGRID_Adj.R")

# ---- Load BioGRID ----

biogrid_data <- read.delim("/home/francescoc/Downloads/BIOGRID-ALL-4.4.241.tab3.txt", header = TRUE, stringsAsFactors = FALSE)

biogrid_data <- biogrid_data[, c("Official.Symbol.Interactor.A", 
                                 "Official.Symbol.Interactor.B", 
                                 "Score",
                                 "Experimental.System.Type",
                                 "Organism.Name.Interactor.A",
                                 "Organism.Name.Interactor.B")]
names(biogrid_data) <- c("Interactor_A", "Interactor_B", "Score", "type", "org1", "org2")

biogrid_data <- biogrid_data %>%
  mutate(Score = as.numeric(Score)) %>%
  filter(type == "physical") %>%       
  filter(!is.na(Score)) %>%            
  filter(abs(Score) > 300) %>%
  filter(org1 == "Homo sapiens" & org2 == "Homo sapiens")

write.table(biogrid_data, "biogrid_physical_s300.txt", sep = "\t", quote = F, col.names = T, row.names = T)

# ---- Load Atlas data ----

options(timeout = 600)
seurat_url <- "https://datasets.cellxgene.cziscience.com/8e64f5c1-e56c-4b0a-bc83-1447fed2e7a4.rds"

seurat_object <- download_Atlas(seurat_url)
celltype <- exploreCells(seurat_object)

mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

pathgenes <- pathg(seurat_object = seurat_object, cell_type = "mature NK T cell", mart = mart, top_n = 1000)

#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

top_genes <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = pathgenes,
  mart = mart
)

top_genes <- top_genes$external_gene_name[top_genes$external_gene_name != ""]

result <- BioGRID_Adj(top_genes, biogrid_data)
wadjm <- result$weighted
adjm <- result$binary

common_names <- intersect(rownames(adjm), colnames(adjm))
adjm <- adjm[common_names, common_names, drop = FALSE]

print(dim(wadjm))
print(dim(adjm))

write.table(adjm, "./../analysis/adjm_p500n643.txt", sep = "\t", quote = F, col.names = T, row.names = T)
write.table(wadjm, "./../analysis/wadjm_p500n643.txt", sep = "\t", quote = F, col.names = T, row.names = T)

# ---- Create Ground Truth ----

gtruth <- igraph::graph_from_adjacency_matrix(adjm, mode = "undirected", diag = F)

num_nodes <- vcount(gtruth)
num_edges <- ecount(gtruth)

set.seed(1234)

p1 <- ggraph(gtruth, layout = "fr") + 
  geom_edge_link(color = "gray", width = 0.5) +  # Set edge color and width
  geom_node_point(color = "steelblue", size = 0.7) +  # Set node color and size
  labs(title = paste("Ground Truth\nNodes:", vcount(gtruth), "Edges:", ecount(gtruth))) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave("./../analysis/plots/gtruth_top1000_s300.png", p1, width = 8, height = 6, dpi = 300, bg = "white")

# ---- Simulate Single-cell expression data ----

ncell <- 500
nodes <- nrow(adjm)

set.seed(1130)
mu_values <- c(3, 6, 9)
theta_values <- c(1, 0.7, 0.5)

count_matrices <- lapply(1:3, function(i) {
  set.seed(1130 + i)
  mu_i <- mu_values[i]
  theta_i <- theta_values[i]
  
  count_matrix_i <- simdata(n = ncell, p = nodes, B = adjm, family = "ZINB", 
                            mu = mu_i, mu_noise = 1, theta = theta_i, pi = 0.2)
  
  count_matrix_df <- as.data.frame(count_matrix_i)
  colnames(count_matrix_df) <- colnames(adjm)
  rownames(count_matrix_df) <- paste("cell", 1:nrow(count_matrix_df), sep = "")
  
  return(count_matrix_df)
})

saveRDS(count_matrices, "./../analysis/count_matrices_n500xp643.RDS")
