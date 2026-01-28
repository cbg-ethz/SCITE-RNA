REF.FILE <- snakemake@input[[1]]
ALT.FILE <- snakemake@input[[2]]
SCORE.FILE <- snakemake@output[[1]]
K <- snakemake@params$K
NUM.SAMPLES <- snakemake@params$num_samples
RNG.SEED <- snakemake@params$rng_seed
START <- snakemake@params$start
SKIP <- snakemake@params$skip
NUM.CORES <- snakemake@params$num_cores


library(parallel)
source("phylinsic_scripts/beastlib.R")

my.write <- function(X, filename, row.names=FALSE, col.names=FALSE) {
  data.out <- as.matrix(X)
  if(is.logical(row.names)) {
    if(row.names)
      row.names <- rownames(X)
    else
      row.names <- c()
  }
  if(is.logical(col.names)) {
    if(col.names)
      col.names <- colnames(X)
    else
      col.names <- c()
  }
  if(length(col.names))
    data.out <- rbind(col.names, data.out)
  if(length(row.names)) {
    if(length(col.names))
      row.names <- c("", row.names)
    data.out <- cbind(row.names, data.out)
  }
  write.table(data.out, filename, quote=FALSE, sep="\t",
    row.names=FALSE, col.names=FALSE)
}


# Adapted to read in the SCITERNA simulated data
data.ref <- read.table(REF.FILE, header = FALSE, sep = "", colClasses = "integer", comment.char = "")
data.alt <- read.table(ALT.FILE, header = FALSE, sep = "", colClasses = "integer", comment.char = "")

data.ref <- as.data.frame(t(data.ref))
data.alt <- as.data.frame(t(data.alt))

data.ref <- cbind(SNV = paste0("SNV", seq_len(nrow(data.ref))), data.ref)
data.alt <- cbind(SNV = paste0("SNV", seq_len(nrow(data.alt))), data.alt)

colnames(data.ref)[-1] <- paste0("Cell", seq_len(ncol(data.ref) - 1))
colnames(data.alt)[-1] <- paste0("Cell", seq_len(ncol(data.alt) - 1))


if(nrow(data.ref) != nrow(data.alt)) stop("bad 1")
if(ncol(data.ref) != ncol(data.alt)) stop("bad 2")
if(any(data.ref[,1] != data.alt[,1])) stop("bad 3")
if(any(colnames(data.ref) != colnames(data.alt))) stop("bad 4")
  
M.ref <- as.matrix(data.ref[,2:ncol(data.ref)])
M.alt <- as.matrix(data.alt[,2:ncol(data.alt)])
M.ref[is.na(M.ref)] <- 0
M.alt[is.na(M.alt)] <- 0


# Calculate the probability distributions from cell j.
# m-length list of n x 3 matrices of probabilities.
p.j <- mclapply(1:ncol(M.ref), function(j) {
  calc.p.dist(M.ref[,j], M.alt[,j])
  }, mc.cores=NUM.CORES)

# This loop is really slow.
# Make a ncol x ncol matrix of cell similarities.  Optimization:
# calculate the upper diagonal only.

I <- seq(START, ncol(M.ref), SKIP)
x <- mclapply(I, function(j1) {
  set.seed(RNG.SEED)  # make sure this is reproducible.
  S.j <- rep(0, ncol(M.ref))
  # Will score lower diagonal of the matrix.
  for(j2 in j1:ncol(M.ref)) {
    S.j[j2] <- calc.cell.similarity(p.j[[j1]], p.j[[j2]], NUM.SAMPLES)
  }
  return(S.j)
}, mc.cores=NUM.CORES)
# Columns are the cells in this batch.
S <- matrix(unlist(x), nrow=ncol(M.ref), ncol=length(I))

# Write out the score file.
cell.names <- colnames(M.ref)
data.out <- cbind(cell.names, S)
colnames(data.out) <- c("Score", cell.names[I])
my.write(data.out, SCORE.FILE, col.names=TRUE)

