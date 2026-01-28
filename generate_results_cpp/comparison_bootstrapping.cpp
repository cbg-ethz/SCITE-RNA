/*
Run bootstrapping comparison for simulated data.
*/

#include <vector>
#include <string>
#include <sstream>
#include <mutation_filter.h>

#include "utils.h"
#include <mutation_filter.h>
#include "generate_results.h"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_dir> <index>" << std::endl;
        return 1;
    }

    std::string sample_dir = argv[1]; // e.g., "50c500m" (only for simulated data), otherwise use real_data_processing.cpp
    int idx = std::stoi(argv[2]); // simulated sample id

    std::string input_path = "../data/simulated_data/" + sample_dir + "/";
    std::string output_path = "../data/results/" + sample_dir + "/sciterna_" + std::to_string(idx);

    std::string alt_file = input_path + "alt/alt_" + std::to_string(idx) + ".txt";
    std::string ref_file = input_path + "ref/ref_" + std::to_string(idx) + ".txt";

    std::vector<std::vector<int>> alt = load_txt(alt_file);
    std::vector<std::vector<int>> ref = load_txt(ref_file);

    int bootstrap_samples = 1000; // number of bootstrap samples
    bool use_bootstrap = true; // use bootstrap samples or not
    int n_snvs = std::stoi(sample_dir.substr(sample_dir.find('c') + 1, sample_dir.find('m') - sample_dir.find('c') - 1)); // number of mutations
    double posterior_threshold = 0.05; // posterior threshold for filtering mutations, only used if method is "threshold"
    std::string method = "first_k"; // "threshold", "first_k", or "highest_post"
    int n_rounds = 2; // number of rounds for tree inference and parameter optimization
//    bool flipped_mutation_direction = true; // whether to allow to change the root genotype during tree inference
    bool only_preprocessing = false; // if true, only preprocess the data and do not run the inference
    std::vector<std::string> tree_space = {"c", "m"}; // tree space to use, "c" for cell space, "m" for mutation space, or both
    bool reshuffle_nodes = true; // whether to reshuffle the nodes in the mutation tree space by pruning and re-inserting individual nodes
    bool load_from_file = false; // if true, load the selected mutations and genotypes from file, otherwise preprocess the data
    bool reduced_output = true; // if true doesn't save the genotype and mutation indicator files

    generate_sciterna_results(ref, alt, input_path, output_path,
                              bootstrap_samples, use_bootstrap, tree_space,
                              flipped_mutation_direction, n_snvs, posterior_threshold,
                              n_rounds, only_preprocessing, method, reshuffle_nodes, load_from_file, reduced_output);

    return 0;
}