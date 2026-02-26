# Load required libraries
library(SpiecEasi)
library(phyloseq)
library(igraph)

# Load OTU, taxonomic, and sample data
otu <- read.csv("all_rarefy_sparcc_counts18.csv", header = T, row.names = 1)
tax <- read.csv("all_rarefy_sparcc_counts18_info.csv", header = T, row.names = 1)
sample <- read.csv("sample.csv", header = T, row.names = 1)

# Create phyloseq object
data <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(as.matrix(tax)), sample_data(sample))

# Construct co-occurrence network using SPIEC-EASI
t1 <- spiec.easi(data, method = 'glasso', 
                 sel.criterion = "bstars", 
                 lambda.log = TRUE, 
                 lambda.min.ratio = 1e-2, 
                 nlambda = 20, 
                 pulsar.select = TRUE, 
                 pulsar.params = list(thresh = 0.05, 
                                      subsample.ratio = 0.8, 
                                      rep.num = 20, 
                                      seed = 123))

# Convert network to igraph object
t2 <- adj2igraph(getRefit(t1), rmEmptyNodes = TRUE, diag = FALSE, 
                 vertex.attr = list(name = taxa_names(data)))

# Get the precision matrix (sparse inverse covariance matrix)
prec_matrix <- as.matrix(getOptCov(t1))
colnames(prec_matrix) <- taxa_names(data)
rownames(prec_matrix) <- taxa_names(data)

# Standardize to correlation matrix
weighted_matrix <- cov2cor(prec_matrix)
weighted_matrix <- -weighted_matrix
diag(weighted_matrix) <- 1

# Convert to adjacency matrix with weighted edges
adj_structure <- as.matrix(getRefit(t1))
final_network_matrix <- weighted_matrix * adj_structure
diag(final_network_matrix) <- 0
final_network_matrix <- (final_network_matrix + t(final_network_matrix)) / 2
colnames(final_network_matrix) <- taxa_names(data)
rownames(final_network_matrix) <- taxa_names(data)

# Save the network matrix
write.csv(final_network_matrix, "final_network_matrix_SpiecEasi.csv", row.names = TRUE)

# Create binary network for further analysis
mat_binary <- ifelse(final_network_matrix != 0, 1, 0)
g <- graph_from_adjacency_matrix(mat_binary, mode = "undirected", weighted = NULL, diag = FALSE)
g <- delete_vertices(g, which(degree(g) == 0))

# Save adjacency matrix of the binary network
adjacency_otu <- as.matrix(get.adjacency(g, sparse = FALSE))
write.csv(adjacency_otu, "adjacency_otu.csv", row.names = TRUE)

# Prepare OTU presence/absence matrix (original counts > 0)
otu_raw <- read.csv("all_rarefy_sparcc_counts18.csv", row.names = 1, check.names = FALSE)
otu_present <- ifelse(otu_raw > 0, 1, 0)
write.csv(otu_present, "otu_present.csv")

# Filter OTU data based on network and abundance table
adjacency_otu_df <- read.csv("adjacency_otu.csv", row.names = 1, check.names = FALSE)
otu_present_df <- read.csv("otu_present.csv", row.names = 1, check.names = FALSE)

common_otus <- intersect(rownames(adjacency_otu_df), rownames(otu_present_df))

# Filter OTU presence/absence matrix based on common OTUs
otu_present_filtered <- otu_present_df[common_otus, ]
write.csv(otu_present_filtered, "otu_present_filtered.csv", row.names = TRUE)

# Extract subnetwork for each sample
g_filtered <- induced_subgraph(g, common_otus)
sub_graph <- list()
sample_names <- colnames(otu_present_filtered)

for (i in sample_names) {
  sample_otus <- rownames(otu_present_filtered)[otu_present_filtered[, i] == 1]
  if (length(sample_otus) >= 2) {
    sub_graph[[i]] <- induced_subgraph(g_filtered, sample_otus)
  }
}

# Calculate topological properties for each subnetwork
sub_graph_stat <- data.frame(
  nodes_num = sapply(sub_graph, function(g) vcount(g)),
  edges_num = sapply(sub_graph, function(g) ecount(g)),
  average_degree = sapply(sub_graph, function(g) mean(degree(g))),
  average_path_length = sapply(sub_graph, function(g) mean_distance(g, directed = FALSE)),
  betweenness_centralization = sapply(sub_graph, function(g) centr_betw(g, directed = FALSE)$centralization),
  graph_density = sapply(sub_graph, function(g) edge_density(g)),
  clustering_coefficient = sapply(sub_graph, function(g) transitivity(g, type = "global")),
  row.names = sample_names
)

# Print summary of subnetwork statistics
print(head(sub_graph_stat))

# Save final results
write.csv(sub_graph_stat, "sub_graph_statistics.csv", row.names = TRUE)