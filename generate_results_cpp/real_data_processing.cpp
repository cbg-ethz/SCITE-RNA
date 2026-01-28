/*
Script used to run SCITE-RNA on real cancer datasets.
*/

#include <vector>
#include <string>
#include <mutation_filter.h>

#include "utils.h"
#include <mutation_filter.h>
#include "generate_results.h"


int main() {
    int bootstrap_samples = 100; // how many bootstrap samples to use for tree inference, only used if use_bootstrap is true
    bool use_bootstrap = false; // if true, bootstrap samples of the SNVs are used, otherwise the results are computed from the original data
    int n_snps = 3000; // how many mutations to select for tree inference, only used if method is "first_k"
    double posterior_threshold = 0.05; // threshold for filtering mutations, only used if method is "threshold"
    std::string method = "threshold"; // criterion for filtering mutations, can be "threshold" (loci with posterior probability of being mutated > threshold are selected) , "first_k" (k loci with highest posterior probability of being mutated are selected) or "highest_post" (loci where highest probability is a mutation and not no mutation are selected)
    int n_rounds = 2; // how many rounds of tree inference and parameter optimization to perform
    std::string sample = "BT_S2"; // specify the sample name, e.g. "mm16", "mm34"
    bool flipped_mutation_direction = true; // if true, we allow the model to switch the root genotype/mutation direction during tree inference
    bool only_preprocessing = false; // if true, only the mutation filtering step is performed, no tree inference or parameter optimization
    std::vector<std::string> tree_space = {"c", "m"}; // which tree spaces to use during optimization, "c" for cell lineage tree, "m" for mutation tree, the order determines in which space the optimization starts
    bool reshuffle_nodes = false; // false makes optimization faster for large numbers of mutations as individual nodes are not pruned and reinserted, instead only pruning and reattaching of subtrees is performed
    bool load_from_file = false; // if true, the selected mutations and genotypes are loaded from file, otherwise they are computed from the data
    bool reduced_output = false; // if true the number of output files is reduced, i.e. the genotype and mutation indicator files are not saved

    std::string input_path = "../data/input_data/" + sample;
    std::string output_path = "../data/results/" + sample + "/sciterna";

    std::vector<std::vector<int>> ref = read_csv(input_path + "/ref.csv");
    std::vector<std::vector<int>> alt = read_csv(input_path + "/alt.csv");

//    for (auto& row : ref) {
//        for (auto& entry : row) {
//            entry += 1;
//        }
//    }

    generate_sciterna_results(ref, alt, input_path, output_path,
                            bootstrap_samples, use_bootstrap, tree_space,
                            flipped_mutation_direction, n_snps, posterior_threshold,
                            n_rounds, only_preprocessing, method, reshuffle_nodes, load_from_file, reduced_output);
    return 0;
}
