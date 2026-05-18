/*
Script used to run SCITE-RNA on simulated data with a variable number of cells, SNVs and clones.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>

#include "utils.h"
#include "generate_results.h"


int main() {
    int n_tests = 100; //number of runs
    int n_rounds = 2; //number of optimization rounds
    std::vector<int> n_cells_list = {50}; // number of cells in the simulation
    std::vector<int> n_mut_list = {500}; // number of mutations in the simulation
    std::vector<std::string> tree_space = {"c", "m"}; // {"c"} {"m"} for cell lineage and mutation tree, respectively, starting optimization in the first tree space
    std::vector<std::string> clones = {""}; //{"", "5", "10", "20"}; // list of clones to simulate, empty string for randomly placed mutations
    bool flipped_mutation_direction = true; // if true, we allow the model to switch the root genotype/mutation direction during tree inference
    bool insert_nodes = true; // if true, we allow SCITE-RNA to reinsert single pruned mutations at their optimal location during optimization

    load_config("../config/config.yaml");
    std::cout << "Random seed: " << config_variables["random_seed"] << std::endl;

//    std::vector<std::string> param_names = {
//            "dropout", "overdispersion_Het", "overdispersion_Hom", "error_rate",
//            "coverage_mean", "coverage_zero_inflation", "coverage_dispersion"
//    };
    std::vector<std::string> param_names = {
            "CNV_fraction"
    };

    std::vector<std::vector<std::string>> param_values = {
            {"0", "0_2", "0_5", "0_8"}
    };

//    std::vector<std::vector<std::string>> param_values = {
//            {"0", "0_2", "0_4", "0_6"},
//            {"3", "6", "10", "100"},
//            {"3", "6", "10", "100"},
//            {"0_001", "0_01", "0_05", "0_1"},
//            {"10", "30", "60", "100"},
//            {"0", "0_2", "0_4", "0_6"},
//            {"1", "2", "5", "10"}
//    };


    for (const auto& clone : clones) {
        for (int i = 0; i < n_cells_list.size(); ++i) {
            int n_cells = n_cells_list[i];
            int n_mut = n_mut_list[i];

            if (n_mut >= 500){
                insert_nodes = false; // disable node insertion for large datasets to speed up inference
            }

            for (size_t p = 0; p < param_names.size(); ++p) {
                for (const auto& val : param_values[p]) {
                    std::string param_str = param_names[p] + "_" + val;
                    std::string base_path = "../data/simulated_data/" + std::to_string(n_cells) + "c" +
                                            std::to_string(n_mut) + "m" + clone + "_param_testing/";
                    std::string path = base_path + param_str + "/";
                    std::string pathout = path + "sciterna/";

                    generate_sciterna_simulation_results(
                            path, pathout, n_tests, tree_space,
                            flipped_mutation_direction, n_mut, n_rounds,
                            insert_nodes
                    );
                }
            }
        }
    }
    return 0;
}