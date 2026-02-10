Network_analysis
#########################################
# Microbial Co-occurrence Network Analysis 
#########################################

library(Hmisc)       
library(igraph)
library(dplyr)
library(reshape2)
library(psych)

rm(list=ls())

# ----------------------------
# 1. Import ASV table
# ----------------------------
asv_table <- read.csv("E:/aaa2019zydna/net1/net_adata.csv", row.names = 1, check.names = FALSE)
head(asv_table)

# ----------------------------
# 2. Spearman correlations
# ----------------------------
cor_res <- corr.test(t(asv_table), method = "spearman", adjust = "none") 

r <- cor_res$r      
p <- cor_res$p      

# Adjust p-values (BH)
p_adj <- p.adjust(p, method = "BH")

# Threshold: significant and strong correlations
r[p_adj > 0.05 | abs(r) < 0.8] <- 0

# ----------------------------
# 3. Construct network
# ----------------------------
net_igraph <- graph_from_adjacency_matrix(r,
                                          mode = "undirected",
                                          weighted = TRUE,
                                          diag = FALSE)

net_igraph <- simplify(net_igraph)
net_igraph <- delete.vertices(net_igraph, V(net_igraph)[degree(net_igraph) == 0])
E(net_igraph)$correlation <- E(net_igraph)$weight

# Export adjacency
adjacency_asv <- as.data.frame(get.adjacency(net_igraph, sparse = FALSE))
adjacency_cor <- as.data.frame(get.adjacency(net_igraph, sparse = FALSE, attr = "correlation"))

write.csv(adjacency_asv, "adjacency_asv.csv")
write.csv(adjacency_cor, "adjacency_cor_asv.csv")

# ----------------------------
# 4. Network summary
# ----------------------------
cat("Number of ASVs (nodes):", vcount(net_igraph), "\n")
cat("Number of edges:", ecount(net_igraph), "\n")

# ----------------------------
# 5. Identify representative ASVs
# ----------------------------
asv_represent <- ifelse(asv_table >= 0.001, 1, 0)
asv_represent <- asv_represent[rowSums(asv_represent) > 0, ]
write.csv(asv_represent, "asv_represent.csv")

# ----------------------------
# 6. Filter ASVs present in network
# ----------------------------
adjacency_asv <- read.csv("adjacency_asv.csv")
asv_represent1 <- read.csv("asv_represent.csv")

asv_represent_result <- asv_represent1[asv_represent1$X %in% adjacency_asv$X, ]
write.csv(asv_represent_result, "asv_represent_result.csv", row.names = FALSE)

# ----------------------------
# 7. Subgraph extraction per sample 
# ----------------------------
asv <- read.csv("asv_represent_result.csv", row.names = 1, check.names = FALSE)

sub_graph <- list()
for (i in names(asv)) {
  sample_i <- asv[i]
  select_node <- rownames(sample_i)[which(sample_i != 0)]
  sub_graph[[i]] <- subgraph(net_igraph, select_node) 
}

# ----------------------------
# 8. Compute topological properties
# ----------------------------
sub_graph_stat <- data.frame(
  nodes_num = numeric(),
  edges_num = numeric(),
  average_degree = numeric(),
  average_path_length = numeric(),
  betweenness_centralization = numeric(),
  graph_diameter = numeric(),
  graph_density = numeric(),
  clustering_coefficient = numeric(),
  degree_centralization = numeric(),
  average_degree_centralization = numeric(),
  modularity = numeric()
)

for(i in 1:length(sub_graph)){
  g <- sub_graph[[i]]
  sample_name <- names(sub_graph[i])
  
  sub_graph_stat <- rbind(sub_graph_stat, data.frame(
    sample_name = sample_name,
    nodes_num = length(V(g)),
    edges_num = length(E(g)),
    average_degree = mean(degree(g)),
    average_path_length = average.path.length(g, directed = FALSE),
    betweenness_centralization = centralization.betweenness(g)$centralization,
    graph_diameter = diameter(g, directed = FALSE),
    graph_density = graph.density(g),
    clustering_coefficient = transitivity(g),
    degree_centralization = centralization.degree(g)$centralization,
    average_degree_centralization = mean(centralization.degree(g)$centralization),
    modularity = modularity(g, membership(cluster_fast_greedy(g)))
  ))
}

write.csv(sub_graph_stat, "sub_graph_stat_asv.csv", row.names = FALSE)
