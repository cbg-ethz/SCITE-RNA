#ifndef SCITE_RNA_GENERATE_RESULTS_H
#define SCITE_RNA_GENERATE_RESULTS_H


#include "mutation_filter.h"
#include "swap_optimizer.h"

void create_directories(const std::string& pathout, bool reduced_ouput);
void process_rounds(MutationFilter &mf, SwapOptimizer &optimizer, const std::vector<std::vector<int>> &ref,
                    const std::vector<std::vector<int>> &alt, int n_snvs, int n_rounds, const std::string &pathout,
                    int i, std::vector<int> selected, const std::vector<char> &gt1, const std::vector<char> &gt2,
                    const std::vector<char> &not_selected_genotypes, int max_loops, bool insert_nodes,
                    bool reduced_output);
void generate_sciterna_simulation_results(const std::string& path,
                                          const std::string& pathout,
                                          int n_tests,
                                          const std::vector<std::string>& tree_space,
                                          bool flipped_mutation_direction,
                                          int n_keep,
                                          int n_rounds,
                                          bool insert_nodes);
void generate_sciterna_results(
        const std::vector<std::vector<int>>& ref,
        const std::vector<std::vector<int>>& alt,
        const std::string& path,
        std::string pathout,
        int n_bootstrap,
        bool use_bootstrap,
        const std::vector<std::string>& tree_space,
        bool flipped_mutation_direction,
        int n_keep,
        double posterior_threshold,
        int n_rounds,
        bool only_preprocessing,
        const std::string& method,
        bool insert_nodes,
        bool load_from_file,
        bool reduced_output);

#endif //SCITE_RNA_GENERATE_RESULTS_H
