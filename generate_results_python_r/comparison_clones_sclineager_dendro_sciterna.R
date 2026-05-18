# Script used to run DENDRO and SClineager on simulated datasets and cluster
# cells based on SCITE-RNA genotype predictions.

library(DENDRO)
library(SClineager)

n_tests <- 100
n_cells_list <- c(200)
n_mut_list <- c(5000)
clones_list <- c("")

base_dir <- file.path("data", "simulated_data")

read.matrix <- function(path) {
  mat <- as.matrix(read.table(path, header = FALSE, sep = " "))
  dimnames(mat) <- NULL
  return(mat)
}

read.str.matrix <- function(path) {
  mat <- as.matrix(read.table(path, header = FALSE, sep = " "))
  mat[mat == "A"] <- 1
  mat[mat == "H"] <- 0.5
  mat[mat == "R"] <- 0
  dimnames(mat) <- NULL
  return(mat)
}

merge.to.parent <- function(merge.mat) {
  n.leaves <- nrow(merge.mat) + 1
  result <- rep(-1, 2 * n.leaves - 1)

  for (i in seq_len(nrow(merge.mat))){
    for (child in merge.mat[i, ]){
      parent.id <- n.leaves + i - 1
      if (child < 0) result[-child] <- parent.id
      else result[n.leaves + child] <- parent.id
    }
  }

  return(result)
}
generate.parent.vec <- function(base_path, n_tests = 10, clones = 5,
run_sclineager = TRUE, run_sciterna_clustering = TRUE, run_dendro = TRUE,
test_index = NULL) {
  if (run_sclineager) {
    dir.create(file.path(base_path, "sclineager", "sclineager_vaf"), recursive = TRUE)
    dir.create(file.path(base_path, "sclineager", "sclineager_parent_vec"), recursive = TRUE)
    dir.create(file.path(base_path, "sclineager", "sclineager_clones"), recursive = TRUE)
    dir.create(file.path(base_path, "sclineager", "sclineager_selected"), recursive = TRUE)
  }

  if (run_sciterna_clustering) {
    dir.create(file.path(base_path, "sciterna", "sciterna_parent_vec_clustering"), recursive = TRUE)
    dir.create(file.path(base_path, "sciterna", "sciterna_clones"), recursive = TRUE)
  }

  if (run_dendro) {
    dir.create(file.path(base_path, "dendro", "dendro_parent_vec"), recursive = TRUE)
    dir.create(file.path(base_path, "dendro", "dendro_clones"), recursive = TRUE)
  }

  cat("Theoretical number of clones: ", clones, "\n")
  sclineager_runtimes <- c()
  dendro_runtimes <- c()

  # If a specific test index was provided, only run that test. Otherwise run all tests 0..(n_tests-1)
  if (!is.null(test_index)) {
    indices <- as.integer(test_index)
    if (any(is.na(indices))) stop("Provided test_index is not an integer")
    # validate range
    if (any(indices < 0) || any(indices > (n_tests - 1))) {
      stop(sprintf("Provided test_index out of range. Must be between 0 and %d", n_tests - 1))
    }
  } else {
    indices <- 0:(n_tests - 1)
  }

  for (i in indices) {
    n_round <- 0
    ref <- t(read.matrix(file.path(base_path, sprintf("ref/ref_%d.txt", i))))
    alt <- t(read.matrix(file.path(base_path, sprintf("alt/alt_%d.txt", i))))

    if (run_sciterna_clustering || run_sclineager) {
      genotype_matrix <- as.matrix(read.table(
        file.path(
          base_path, "sciterna", "sciterna_genotype",
          paste0("sciterna_genotype_", n_round, "r", i, ".txt")
        ),
        stringsAsFactors = FALSE
      ))

      mapping_dict <- c("A" = 1.0, "H" = 0.5, "R" = 0)
      genotype_sciterna <- matrix(mapping_dict[genotype_matrix],
                                  nrow = nrow(genotype_matrix),
                                  ncol = ncol(genotype_matrix))
    }

    coverage <- ref + alt
    mut_indicator <- read.matrix(file.path(base_path, sprintf("mut_indicator/mut_indicator_%d.txt", i)))
    mutations_mat <- alt / coverage

    unique_rows <- unique(t(mut_indicator))
    n_clones <- nrow(unique_rows)
    cat(" Actual number of clones: ", n_clones, "\n")

    if (run_sciterna_clustering) {
      dist_scite <- dist(genotype_sciterna)
      hc_scite <- hclust(dist_scite, method = "ward.D")
      memb_pred_scite <- cutree(hc_scite, k = n_clones)
      cluster_scite <- DENDRO.cluster(dist_scite, plot = FALSE, type = "phylogram")
      parent_vec_scite <- merge.to.parent(cluster_scite$merge)
      write.table(
        parent_vec_scite,
        file.path(
          base_path,
          sprintf(
            "sciterna/sciterna_parent_vec_clustering/sciterna_parent_vec_clustering_%d.txt",
            i
          )
        ),
        row.names = FALSE,
        col.names = FALSE
      )
      write.table(
        memb_pred_scite,
        file.path(
          base_path,
          sprintf("sciterna/sciterna_clones/sciterna_clones_%d.txt", i)
        ),
        row.names = FALSE,
        col.names = FALSE
      )
    }

    if (run_sclineager) {
      start_time_sclineager <- Sys.time()

      keep <- apply(mutations_mat, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) > 0.01)
      keep_numeric <- as.numeric(keep)
      write.table(
        keep_numeric,
        file = file.path(
          base_path,
          sprintf(
            "sclineager/sclineager_selected/sclineager_selected_%d.txt",
            i
          )
        ),
        row.names = FALSE,
        col.names = FALSE
      )
      mutations_mat_sclineager <- mutations_mat[keep, ]
      coverage_sclineager <- coverage[keep, ]

      res_scl <- sclineager_internal(
        mutations_mat = mutations_mat_sclineager,
        coverage_mat = coverage_sclineager,
        max_iter = 2000,
        vaf_offset = 0.01,
        dfreedom = ncol(mutations_mat_sclineager),
        psi = diag(10, ncol(mutations_mat_sclineager)),
        save = FALSE
      )

      dist_scl <- dist(t(res_scl[["genotype_mat"]]))
      hc_scl <- hclust(dist_scl, method = "ward.D")
      memb_pred_scl <- cutree(hc_scl, k = n_clones)
      cluster_scl <- DENDRO.cluster(dist_scl, plot = FALSE, type = "phylogram")
      parent_vec_scl <- merge.to.parent(cluster_scl$merge)

      end_time_sclineager <- Sys.time()
      runtime_sclineager <- as.numeric(difftime(end_time_sclineager, start_time_sclineager, units = "secs"))
      sclineager_runtimes <- c(sclineager_runtimes, runtime_sclineager)

      write.table(
        t(res_scl[["genotype_mat"]]),
        file = file.path(
          base_path,
          sprintf("sclineager/sclineager_vaf/sclineager_vaf_%d.txt", i)
        ),
        row.names = FALSE,
        col.names = FALSE
      )
      write.table(
        parent_vec_scl,
        file.path(
          base_path,
          sprintf(
            "sclineager/sclineager_parent_vec/sclineager_parent_vec_%d.txt",
            i
          )
        ),
        row.names = FALSE,
        col.names = FALSE
      )
      write.table(
        memb_pred_scl,
        file.path(
          base_path,
          sprintf("sclineager/sclineager_clones/sclineager_clones_%d.txt", i)
        ),
        row.names = FALSE,
        col.names = FALSE
      )
    }

    if (run_dendro) {
      start_time_dendro <- Sys.time()

      filtered <- FilterCellMutation(
        alt,
        coverage,
        mut_indicator,
        cut.off.VAF = 0.0,
        cut.off.sd = 10,
        plot = FALSE
      )
      dist <- DENDRO.dist(
        filtered$X,
        filtered$N,
        filtered$Z,
        show.progress = FALSE
      )

      hc <- hclust(dist, method = "ward.D")
      memb_pred <- cutree(hc, k = n_clones)
      cluster <- DENDRO.cluster(dist, plot = FALSE, type = "phylogram")
      dendro_parent_vec <- merge.to.parent(cluster$merge)

      end_time_dendro <- Sys.time()
      runtime_dendro <- as.numeric(difftime(end_time_dendro, start_time_dendro, units = "secs"))
      dendro_runtimes <- c(dendro_runtimes, runtime_dendro)

      write.table(
        dendro_parent_vec,
        file.path(
          base_path,
          sprintf("dendro/dendro_parent_vec/dendro_parent_vec_%d.txt", i)
        ),
        row.names = FALSE,
        col.names = FALSE
      )
      write.table(
        memb_pred,
        file.path(
          base_path,
          sprintf("dendro/dendro_clones/dendro_clones_%d.txt", i)
        ),
        row.names = FALSE,
        col.names = FALSE
      )
    }
  }

  if (run_sclineager) {
    # If a single test index was requested, include it in the runtime filename
    idx_suffix <- ""
    if (!is.null(test_index) && length(indices) == 1) {
      idx_suffix <- sprintf("_%d", indices)
    }
    write.table(
      sclineager_runtimes,
      file.path(base_path, sprintf("sclineager/sclineager_runtimes%s.txt", idx_suffix)),
      row.names = FALSE,
      col.names = FALSE
    )
  }

  if (run_dendro) {
    # If a single test index was requested, include it in the runtime filename
    idx_suffix <- ""
    if (!is.null(test_index) && length(indices) == 1) {
      idx_suffix <- sprintf("_%d", indices)
    }
    write.table(
      dendro_runtimes,
      file.path(base_path, sprintf("dendro/dendro_runtimes%s.txt", idx_suffix)),
      row.names = FALSE,
      col.names = FALSE
    )
  }
}


paths <- NULL
clones <- NULL

for (i in seq_along(n_cells_list)) {
  n_cells <- n_cells_list[i]
  n_mut <- n_mut_list[i]

  for (clone in clones_list) {
    path <- paste(base_dir, "/", n_cells, "c", n_mut, "m", clone, sep = "")
    paths <- c(paths, path)
    clones <- c(clones, clone)

  }
}

print(paths)

# Parse command line arguments to allow running a specific test index.
# Usage examples:
# Rscript comparison_clones_sclineager_dendro_sciterna.R --test=5
# Rscript comparison_clones_sclineager_dendro_sciterna.R --test=2,4,7
args <- commandArgs(trailingOnly = TRUE)
test_index <- NULL
if (length(args) > 0) {
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      cat("Usage: Rscript comparison_clones_sclineager_dendro_sciterna.R [--test=IDX[,IDX...]]\n")
      quit(status = 0)
    }

    # accept --test=NUM, --test-index=NUM, or -t=NUM
    if (grepl('^--test(?:-index)?=', arg) || grepl('^-t=', arg)) {
      val <- sub('^[^=]*=', '', arg)
      parts <- strsplit(val, ",")[[1]]
      parts <- trimws(parts)
      idx <- as.integer(parts)
      if (any(is.na(idx))) stop("Provided test index contains non-integer values: ", val)
      test_index <- idx
    }
  }
}

for (i in seq_along(paths)) {
  n_clones <- clones[i]
  path <- paths[i]
  generate.parent.vec(path, n_tests, n_clones, run_sclineager = TRUE,
                      run_sciterna_clustering = FALSE, run_dendro = TRUE,
                      test_index = test_index)
}