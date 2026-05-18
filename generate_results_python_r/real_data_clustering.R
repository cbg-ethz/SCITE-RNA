# Script used to run DENDRO and SClineager on cancer datasets and cluster SCITE-RNA's genotype matrix.
library(DENDRO)
library(SClineager)

base_path <- file.path("data")

sample <- "BT_S2"
clones <- c(2,3,4,5,6,7)

read.matrix <- function(path){
  mat <- as.matrix(read.table(path, header=FALSE, sep=" "))
  dimnames(mat) <- NULL
  return(mat)
}

merge.to.parent <- function(merge.mat){
  n.leaves <- nrow(merge.mat) + 1
  result <- rep(-1, 2*n.leaves - 1)

  for (i in 1:nrow(merge.mat)){
    for (child in merge.mat[i,]){
      parent.id <- n.leaves + i - 1
      if (child < 0) result[-child] <- parent.id
      else result[n.leaves + child] <- parent.id
    }
  }

  return(result)
}

dir.create(file.path(base_path, "results", sample, "dendro", "dendro_clones"), recursive = TRUE)
dir.create(file.path(base_path, "results", sample, "sclineager", "sclineager_clones"), recursive = TRUE)
dir.create(file.path(base_path, "results", sample, "sciterna", "sciterna_clones"), recursive = TRUE)
dir.create(file.path(base_path, "results", sample, "sciterna_bootstrap", "sciterna_bootstrap_clones"), recursive = TRUE)
dir.create(file.path(base_path, "results", sample, "sciterna_bootstrap", "sciterna_clustering_parent_vec"), recursive = TRUE)

ref <- t(as.matrix(read.csv(file.path(base_path, "input_data", sample, "ref.csv"), row.names = 1)))
ref[is.na(ref)] <- 0
dimnames(ref) <- NULL
alt <- t(as.matrix(read.csv(file.path(base_path, "input_data", sample, "alt.csv"), row.names = 1)))
alt[is.na(alt)] <- 0
dimnames(alt) <- NULL

coverage <- ref + alt
mutations_mat <- alt / coverage

n_round <- 1
genotype_matrix <- as.matrix(read.table(file.path(base_path, "results", sample, "sciterna", "sciterna_genotype",
                                                  paste0("sciterna_genotype_", n_round, "r", 0, ".txt")),
                                        stringsAsFactors = FALSE))

mapping_dict <- c("A" = 1.0, "H" = 0.5, "R" = 0)

# Apply the mapping
genotype_sciterna <- matrix(mapping_dict[genotype_matrix],
                      nrow = nrow(genotype_matrix),
                      ncol = ncol(genotype_matrix))
dist_scite <- dist(genotype_sciterna) # use same method for other models
hc_scite <- hclust(dist_scite, method='ward.D')

genotype_matrix_consensus <- as.matrix(read.matrix(file.path(base_path, "results", sample, "sciterna_bootstrap", "sciterna_consensus_genotype",
                                                  paste0("median_consensus_genotype_", n_round, ".txt"))))

dist_consensus <- dist(genotype_matrix_consensus)
hc_consensus <- hclust(dist_consensus, method = "ward.D")
cluster_consensus <- DENDRO.cluster(dist_consensus, plot = FALSE, type = "phylogram")
parent_vec_consensus <- merge.to.parent(cluster_consensus$merge)
write.table(parent_vec_consensus, file.path(base_path, "results", sample, "sciterna_bootstrap", "sciterna_clustering_parent_vec", paste0("sciterna_clustering_parent_vec.txt")), row.names=FALSE, col.names=FALSE)



genotype_sclineager <- as.matrix(read.table(file.path(base_path, "results", sample, "sclineager", "sclineager_vaf",
                                                  paste0("sclineager_vaf.txt"))))
dist_scl <- dist(genotype_sclineager)
hc_scl <- hclust(dist_scl, method='ward.D')

mut_indicator <- matrix(ifelse(round(mutations_mat, 1) > 0.3, 1, 0), nrow = nrow(mutations_mat),
                        ncol = ncol(mutations_mat))
filtered <- FilterCellMutation(alt, coverage, mut_indicator, cut.off.VAF = 0.02, cut.off.sd = 10, plot=FALSE)
dist_dendro <- DENDRO.dist(filtered$X,filtered$N,filtered$Z,show.progress=TRUE)
hc_dendro <- hclust(dist_dendro, method='ward.D')

for (clone in clones) {
    memb_pred_scite <- cutree(hc_scite, k = clone)
    write.table(memb_pred_scite, file.path(base_path, "results", sample, "sciterna", "sciterna_clones", paste0("sciterna_clones_", clone, ".txt")), row.names=FALSE, col.names=FALSE)

    memb_pred_consensus <- cutree(hc_consensus, k = clone)
    write.table(memb_pred_consensus, file.path(base_path, "results", sample, "sciterna_bootstrap", "sciterna_bootstrap_clones", paste0("sciterna_bootstrap_clones_", clone, ".txt")), row.names=FALSE, col.names=FALSE)

    memb_pred_scl <- cutree(hc_scl, k = clone)
    write.table(memb_pred_scl, file.path(base_path, "results", sample, "sclineager", "sclineager_clones", paste0("sclineager_clones_", clone, ".txt")), row.names=FALSE, col.names=FALSE)

    memb_pred <- cutree(hc_dendro, k = clone)
    write.table(memb_pred, file.path(base_path, "results", sample, "dendro", "dendro_clones", paste0("dendro_clones_", clone, ".txt")), row.names=FALSE, col.names=FALSE)
}