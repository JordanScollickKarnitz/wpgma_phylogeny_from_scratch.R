############################################################
# WPGMA Phylogenetic Tree Construction (From Scratch)
# Author: Jordan Scollick-Karnitz
# Language: R
#
# Description:
# Implements the WPGMA (Weighted Pair Group Method with
# Arithmetic Mean) clustering algorithm using Jukes–Cantor
# distances to construct a phylogenetic tree.
############################################################

############################
# Input sequences
############################

seq1 <- "ACTGGGCT"
seq2 <- "CGTGAGCT"

############################
# Jukes–Cantor distance
############################

jukes_cantor_distance <- function(seq1, seq2) {
  
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  
  # Proportion of mismatches
  p <- sum(s1 != s2) / length(s1)
  
  # Jukes–Cantor correction
  d <- (-3 / 4) * log(1 - (4 * p / 3))
  return(d)
}

############################
# Build initial distance matrix
############################

all_sequences <- list(seq1 = seq1, seq2 = seq2)
n <- length(all_sequences)
seq_names <- names(all_sequences)

distance_matrix <- matrix(0, nrow = n, ncol = n,
                          dimnames = list(seq_names, seq_names))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      distance_matrix[i, j] <-
        jukes_cantor_distance(all_sequences[[i]],
                              all_sequences[[j]])
    }
  }
}

cat("Initial Distance Matrix:\n")
print(distance_matrix)

############################
# WPGMA clustering setup
############################

active_clusters <- seq_names

tree_edges <- data.frame(
  parent = character(),
  child1 = character(),
  child2 = character(),
  height = numeric(),
  stringsAsFactors = FALSE
)

next_node_id <- n + 1
current_matrix <- distance_matrix

############################
# Find minimum distance
############################

find_min_distance <- function(mat) {
  
  diag(mat) <- Inf
  mat[lower.tri(mat)] <- Inf
  
  min_val <- min(mat)
  min_pos <- which(mat == min_val, arr.ind = TRUE)[1, ]
  
  return(list(
    value = min_val,
    row = min_pos[1],
    col = min_pos[2]
  ))
}

############################
# Update distance matrix (WPGMA)
############################

update_matrix_wpgma <- function(mat, cluster1, cluster2, new_cluster) {
  
  idx1 <- match(cluster1, colnames(mat))
  idx2 <- match(cluster2, colnames(mat))
  
  other_clusters <- colnames(mat)[-c(idx1, idx2)]
  new_distances <- numeric(length(other_clusters))
  
  # WPGMA: equal weighting of merged clusters
  for (i in seq_along(other_clusters)) {
    oc <- other_clusters[i]
    new_distances[i] <- (mat[cluster1, oc] + mat[cluster2, oc]) / 2
  }
  
  updated <- mat[other_clusters, other_clusters, drop = FALSE]
  
  updated <- rbind(updated, new_distances)
  rownames(updated)[nrow(updated)] <- new_cluster
  
  updated <- cbind(updated, c(new_distances, 0))
  colnames(updated)[ncol(updated)] <- new_cluster
  
  return(updated)
}

############################
# Main WPGMA loop
############################

while (length(active_clusters) > 1) {
  
  min_res <- find_min_distance(current_matrix)
  c1 <- rownames(current_matrix)[min_res$row]
  c2 <- colnames(current_matrix)[min_res$col]
  min_dist <- min_res$value
  
  new_height <- min_dist / 2
  new_node <- paste0("Node", next_node_id)
  next_node_id <- next_node_id + 1
  
  tree_edges <- rbind(
    tree_edges,
    data.frame(
      parent = new_node,
      child1 = c1,
      child2 = c2,
      height = new_height,
      stringsAsFactors = FALSE
    )
  )
  
  active_clusters <- setdiff(active_clusters, c(c1, c2))
  active_clusters <- c(active_clusters, new_node)
  
  current_matrix <-
    update_matrix_wpgma(current_matrix, c1, c2, new_node)
}

############################
# Tree plotting function
############################

plot_wpgma_tree <- function(tree_edges, tips) {
  
  nodes <- unique(c(tree_edges$parent,
                    tree_edges$child1,
                    tree_edges$child2))
  
  x_pos <- setNames(rep(NA, length(nodes)), nodes)
  y_pos <- setNames(rep(0, length(nodes)), nodes)
  
  for (i in seq_along(tips)) {
    x_pos[tips[i]] <- i
  }
  
  for (i in 1:nrow(tree_edges)) {
    y_pos[tree_edges$parent[i]] <- tree_edges$height[i]
  }
  
  for (i in 1:nrow(tree_edges)) {
    p <- tree_edges$parent[i]
    x_pos[p] <- mean(c(x_pos[tree_edges$child1[i]],
                       x_pos[tree_edges$child2[i]]))
  }
  
  plot(NULL,
       xlim = c(0.5, length(tips) + 0.5),
       ylim = c(0, max(tree_edges$height) * 1.1),
       xlab = "",
       ylab = "Distance",
       main = "WPGMA Phylogenetic Tree",
       xaxt = "n")
  
  axis(1, at = seq_along(tips), labels = tips, las = 2)
  
  for (i in 1:nrow(tree_edges)) {
    p <- tree_edges$parent[i]
    c1 <- tree_edges$child1[i]
    c2 <- tree_edges$child2[i]
    
    segments(x_pos[c1], y_pos[p], x_pos[c2], y_pos[p], lwd = 2)
    segments(x_pos[c1], y_pos[c1], x_pos[c1], y_pos[p], lwd = 2)
    segments(x_pos[c2], y_pos[c2], x_pos[c2], y_pos[p], lwd = 2)
  }
  
  points(x_pos[tips], y_pos[tips], pch = 19, col = "blue")
  internal <- setdiff(names(x_pos), tips)
  points(x_pos[internal], y_pos[internal], pch = 19, col = "red")
}

############################
# Run and visualize
############################

plot_wpgma_tree(tree_edges, seq_names)

cat("Tree Structure:\n")
print(tree_edges)
