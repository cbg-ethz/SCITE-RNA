# Snakemake workflow for running PhylinSic on our processed real cancer single-cell RNA-seq data.
# See https://github.com/U54Bioinformatics/PhylinSic_Project for the full PhylinSic pipeline

import os
from os.path import join as opj


APPENDIX = str(config.get("APPENDIX", ""))
NUM_CELLS = int(config.get("NUM_CELLS", 1169))


DATA_DIR = rf"data/results/{APPENDIX}"
GENOTYPE_DIR = opj(DATA_DIR, "phylinsic", "output/genotype")
PHYLO_DIR = opj(DATA_DIR, "phylinsic", "output/phylogeny")
PHYLO_LOG_DIR = opj(DATA_DIR, "phylinsic", "logs/phylogeny")
DEMUX_DIR = opj(DATA_DIR, "phylinsic", "output/demux")
BEAST_OUTPUT = opj(DATA_DIR, "phylinsic", "output/beast2")
BEAST2_DIR = r"/beast2/beast" #TODO replace

PHYLINSIC_GENOTYPE_DIR = opj(DATA_DIR, "phylinsic", "phylinsic_genotype")
PHYLINSIC_PARENT_VEC_DIR = opj(DATA_DIR, "phylinsic", "phylinsic_parent_vec")

os.makedirs(GENOTYPE_DIR, exist_ok=True)
os.makedirs(PHYLO_DIR, exist_ok=True)
os.makedirs(PHYLO_LOG_DIR, exist_ok=True)
os.makedirs(DEMUX_DIR, exist_ok=True)
os.makedirs(BEAST_OUTPUT, exist_ok=True)
os.makedirs(PHYLINSIC_GENOTYPE_DIR, exist_ok=True)
os.makedirs(PHYLINSIC_PARENT_VEC_DIR, exist_ok=True)

cells = [f"Cell{i}\tno\tA" for i in range(1, NUM_CELLS+1)]
with open(opj(DEMUX_DIR, "cells.txt"), "w") as f:
    f.write("Cell\tOutgroup\tCategory\n" + "\n".join(cells))

LOGCOMBINER = r"/beast2/beast/bin/logcombiner" #TODO replace
RSCRIPT = "Rscript"
JAVA = "java"

# PARAMETERS

# This is the number of neighbors to use to smooth the genotype
# calls.
TERNARY_K = 10

# We compare two cells by comparing sequences drawn from
# distributions determined by their read counts and scoring the
# similarity of their sequences.  This is the number of samples
# to generate. Higher numbers will lead to more accurate scores,
# but takes longer.
TERNARY_SAMPLES = 100

# We smooth the genotypes for each cell by taking a weighted
# average of the probability distribution of the genotypes of a
# cell and its neighbors.  This controls how much weight to give
# to its neighbors.  The default value will give each cell the
# same weight.  Set to 0 to turn off smoothing.
TERNARY_DELTA = 0

# If a cells has no reads at a site, this controls whether to
# impute the genotype at that site, or to leave it as a missing
# value.
TERNARY_IMPUTE = False

# Set the seed of the random number generator.  Can use to
# make sure results are reproducible.
PHYLO_RNG_SEED = 1




# PARAMETERS - PHYLOGENETIC INFERENCE
# -----------------------------------

BEAST2_SITE_MODEL = "gtr"
BEAST2_CLOCK_MODEL = "rln"
BEAST2_TREE_PRIOR = "yule"
BEAST2_ITERATIONS = 10000000


# How many iterations to discard for burn-in.
BEAST2_BURNIN = int(BEAST2_ITERATIONS*0.50)

BEAST2_SAMPLE_INTERVAL = 5000
BEAST2_RNG_SEED = 1


KNN_BATCHES = ["00"]
batch = "00"


rule all:
    input:
        opj(PHYLINSIC_GENOTYPE_DIR, f"phylinsic_genotype.txt"),
        opj(PHYLINSIC_PARENT_VEC_DIR, f"phylinsic_parent_vec.txt")

rule calc_neighbor_scores2:
    input:
        opj(DATA_DIR, "ref.txt"),
        opj(DATA_DIR, "alt.txt"),
    output:
        opj(GENOTYPE_DIR, "knn.neighbors.0.txt"),
    params:
        K=TERNARY_K,
        num_samples=TERNARY_SAMPLES,
        rng_seed=PHYLO_RNG_SEED,
        start=1,
        skip=1,
        num_cores=1
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    script:
        "phylinsic_scripts/calc_neighbor_scores1.R"


rule select_neighbors2:
    input:
        opj(GENOTYPE_DIR, "knn.neighbors.0.txt"),
    output:
        opj(GENOTYPE_DIR, "knn.neighbors.txt"),
    params:
        K=TERNARY_K,
    script:
        "phylinsic_scripts/select_neighbors1.py"


rule call_genotypes2:
    input:
        opj(DATA_DIR, "ref.txt"),
        opj(DATA_DIR, "alt.txt"),
        opj(GENOTYPE_DIR, "knn.neighbors.txt"),
    output:
        opj(GENOTYPE_DIR, "genotypes.txt"),
        opj(GENOTYPE_DIR, "probabilities.txt"),
    params:
        delta=TERNARY_DELTA,
        impute=TERNARY_IMPUTE,
        num_cores=1 #workflow.cores,
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    script:
        "phylinsic_scripts/call_genotypes1.R"


rule make_fasta_file:
    input:
        opj(GENOTYPE_DIR, "genotypes.txt"),
    output:
        opj(PHYLO_DIR, "mutations.fa")
    script:
        "phylinsic_scripts/make_fasta_file.py"


rule make_phylinsic_genotype_file:
    input:
        opj(GENOTYPE_DIR, "genotypes.txt"),
    output:
        opj(PHYLINSIC_GENOTYPE_DIR, "phylinsic_genotype.txt")
    script:
        "phylinsic_scripts/fasta_to_phylinsic_genotype.py"


rule run_beast2:
    input:
        opj(PHYLO_DIR, "mutations.fa")
    output:
        beast_output_dir = directory(opj(BEAST_OUTPUT,"beast2")),
        beast_model=opj(BEAST_OUTPUT, "beast2", "beast2.model.RDS"),
        tree_log=opj(BEAST_OUTPUT, "beast2", "tree.log"),
    log:
        opj(PHYLO_LOG_DIR,"beast2.log"),
    params:
        RSCRIPT=RSCRIPT,
        beast2_dir=BEAST2_DIR,
        site_model=BEAST2_SITE_MODEL,
        clock_model=BEAST2_CLOCK_MODEL,
        tree_prior=BEAST2_TREE_PRIOR,
        iterations=BEAST2_ITERATIONS,
        sample_interval=BEAST2_SAMPLE_INTERVAL,
        rng_seed=BEAST2_RNG_SEED,
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    shell:
        """{RSCRIPT} phylinsic_scripts/run_beast2.R -i {input} -o {output[0]} \
            --beast2_dir {params.beast2_dir} \
            --site_model {params.site_model} \
            --clock_model {params.clock_model} \
            --tree_prior {params.tree_prior} \
            --iterations {params.iterations} \
            --sample_interval {params.sample_interval} \
            --rng_seed {params.rng_seed} >& {log}
        """

rule summarize_beast2:
    input:
        opj(BEAST_OUTPUT, "beast2", "beast2.model.RDS"),
    output:
        opj(PHYLO_DIR, "summary.txt"),
        opj(PHYLO_DIR, "summary.ess.txt"),
    log:
        opj(PHYLO_LOG_DIR, "summary.log")
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    params:
        RSCRIPT=RSCRIPT,
        sample_interval=BEAST2_SAMPLE_INTERVAL,
        burnin=BEAST2_BURNIN,
    shell:
        """{params.RSCRIPT} phylinsic_scripts/summarize_beast2.R \
            {input[0]} {params.sample_interval} {params.burnin} \
            {output[0]} {output[1]} >& {log}
         """


perc_burnin = int(round(float(BEAST2_BURNIN) / BEAST2_ITERATIONS * 100))
assert perc_burnin > 0 and perc_burnin < 100

rule combine_trees:
    input:
        opj(BEAST_OUTPUT, "beast2", "tree.log")
    output:
        opj(PHYLO_DIR, "beast2.trees.nexus.txt")
    log:
        opj(PHYLO_LOG_DIR, "beast2.trees.nexus.log")
    params:
        LOGCOMBINER=LOGCOMBINER,
        perc_burnin=perc_burnin,
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    shell:
        """
        {params.LOGCOMBINER} -b {params.perc_burnin} \
            -log {input[0]} -o {output} >& {log}
        """

rule make_mcc_tree:
    input:
        opj(PHYLO_DIR, "beast2.trees.nexus.txt")
    output:
        opj(PHYLO_DIR, "max_clade_cred.nexus.txt")
    log:
        opj(PHYLO_LOG_DIR, "max_clade_cred.nexus.log")
    params:
        JAVA=JAVA,
        BEAST2_DIR=BEAST2_DIR,
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    shell:
        # -burnin 0 because burning already removed by tree combiner
        """
        {params.JAVA} -Xms1g -Xmx8g \
            -Dlauncher.wait.for.exit=true \
            -Duser.language=en \
            -Djava.library.path={params.BEAST2_DIR}/lib \
            -cp {params.BEAST2_DIR}/lib/launcher.jar \
            beast.app.treeannotator.TreeAnnotatorLauncher \
            -burnin 0 \
            -heights mean \
            -lowMem \
            {input} {output} >& {log}
        """


rule analyze_mcc_tree:
    input:
        opj(PHYLO_DIR, "max_clade_cred.nexus.txt"),
        opj(DEMUX_DIR, "cells.txt"),
    output:
        # analyze_mcc_tree.R will generate the rerooted trees only if
        # a root node can be determined (e.g. there is an outgroup).
        # If not, then the files will not be generated.
        #
        # Snakemake considers missing output files to be errors.  So
        # make the rerooted trees parameters rather than output files.
        opj(PHYLO_DIR, "max_clade_cred.newick.txt"),
        #opj(PHYLO_DIR, "max_clade_cred.rerooted.newick.txt"),
        opj(PHYLO_DIR, "max_clade_cred.dist.txt"),
        opj(PHYLO_DIR, "max_clade_cred.metadata.txt"),
        #opj(PHYLO_DIR, "max_clade_cred.rerooted.metadata.txt"),
    log:
        opj(PHYLO_LOG_DIR, "max_clade_cred.analysis.log")
    params:
        RSCRIPT=RSCRIPT,
        REROOT_NEWICK=opj(PHYLO_DIR, "max_clade_cred.rerooted.newick.txt").replace("\\", "/"),
        REROOT_META=opj(PHYLO_DIR, "max_clade_cred.rerooted.metadata.txt").replace("\\", "/"),
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    shell:
        """{params.RSCRIPT} phylinsic_scripts/analyze_mcc_tree.R \
            {input[0]} \
            {input[1]} \
            {output[0]} \
            {params.REROOT_NEWICK} \
            {output[1]} \
            {output[2]} \
            {params.REROOT_META} >& {log}
        """

rule newick_to_parent_vec:
    input:
        opj(PHYLO_DIR, "max_clade_cred.newick.txt"),
        opj(DEMUX_DIR, "cells.txt"),
    output:
        opj(PHYLINSIC_PARENT_VEC_DIR, "phylinsic_parent_vec.txt"),
    conda:
        "phylinsic_scripts/R_scripts.yaml"
    shell:
        """python phylinsic_scripts/newick_to_parent_vec.py \
            --newick_file {input[0]} \
            --cell_file {input[1]} \
            --output_file {output[0]}
        """
