"""
Script used to generate simulated datasets. The number of cells and SNVs and the number of clones simulated can be set.
Additionally, you have the option to set various other simulation parameters. Tree inference can be run as well.
However, you might want to use generate_results_cpp/comparison_num_clones.cpp for tree inference as it is faster.
"""

import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from src_python.data_generator import DataGenerator
from src_python.mutation_tree import MutationTree
from src_python.utils import load_config_and_set_random_seed
from src_python.generate_results import generate_sciterna_simulation_results

from tqdm import tqdm
import numpy as np

config = load_config_and_set_random_seed()


def generate_comparison_data(n_cells: int, n_mut: int, size=100, path="./comparison_data/", random_seed=None,
                             n_clones=None, coverage_distribution="zinb", **kwargs):
    """
    Generates simulated data.
    """
    if random_seed is not None:
        np.random.seed(random_seed)

    if os.path.exists(path):
        while True:
            ans = input(f"Directory {path} already exists. Existing files will be overwritten. Continue? [Y/N] ")
            match ans:
                case "Y" | "y" | "Yes" | "yes":
                    break
                case "N" | "n" | "No" | "no":
                    return
    else:
        os.makedirs(path)

    os.makedirs(os.path.join(path, "ref"), exist_ok=True)
    os.makedirs(os.path.join(path, "alt"), exist_ok=True)
    os.makedirs(os.path.join(path, "parent_vec"), exist_ok=True)
    os.makedirs(os.path.join(path, "mut_indicator"), exist_ok=True)
    os.makedirs(os.path.join(path, "genotype"), exist_ok=True)
    os.makedirs(os.path.join(path, "dropout_probs"), exist_ok=True)
    os.makedirs(os.path.join(path, "overdispersions_H"), exist_ok=True)
    os.makedirs(os.path.join(path, "mutation_location"), exist_ok=True)

    generator = DataGenerator(n_cells, n_mut, coverage_method=coverage_distribution, **kwargs)

    for i in tqdm(range(size)):
        ref, alt, dropout_probs, overdispersions_H = generator.generate_reads(new_tree=True, new_mut_type=True,
                                                                                                  num_clones=n_clones)

        mut_indicator = np.zeros((n_mut, n_cells), dtype=bool)
        for j in range(generator.n_mut):
            if generator.gt1[j] == generator.gt2[j]:
                continue  # not mutated
            for leaf in generator.ct.leaves(generator.ct.mut_loc[j]):  # cells below the mutation
                mut_indicator[j, leaf] = True

        mt = MutationTree(n_mut=n_mut, n_cells=n_cells)
        mt.fit_cell_tree(generator.ct)

        np.savetxt(os.path.join(path, f"ref/ref_{i}.txt"), ref, fmt="%i")
        np.savetxt(os.path.join(path, f"alt/alt_{i}.txt"), alt, fmt="%i")
        np.savetxt(os.path.join(path, f"parent_vec/parent_vec_{i}.txt"), generator.ct.parent_vec, fmt="%i")
        np.savetxt(os.path.join(path, f"mut_indicator/mut_indicator_{i}.txt"), mut_indicator, fmt="%i")
        np.savetxt(os.path.join(path, f"genotype/genotype_{i}.txt"), generator.genotype, fmt="%s")
        # np.savetxt(os.path.join(path, f"genotype/alt_alleles_{i}.txt"), generator.alt_alleles, fmt="%i")
        # np.savetxt(os.path.join(path, f"genotype/ref_alleles_{i}.txt"), generator.ref_alleles, fmt="%i")
        np.savetxt(os.path.join(path, f"dropout_probs/dropout_probs_{i}.txt"), dropout_probs)
        np.savetxt(os.path.join(path, f"overdispersions_H/overd"
                                      f"ispersions_H_{i}.txt"), overdispersions_H)
        np.savetxt(os.path.join(path, f"mutation_location/mutation_location_{i}.txt"), generator.ct.mut_loc, fmt="%i")


num_tests = 1  # Number of simulated samples
n_rounds = 2  # Number of rounds of SCITE-RNA to optimize the SNV specific parameters like dropout probabilities
n_cells_list = [50]  # Number of cells in the simulated dataset
n_mut_list = [50]  # Number of SNVs in the simulated dataset
clones = [""]  # Number of clones in the simulated dataset, empty string means random mutation placement
flipped_mutation_direction = True  # Whether to allow to flip the mutation direction (change root genotype)
tree_space = ["c", "m"]  # Tree spaces to use during tree inference and which space to start. ["c", "m] means we start optimizing a cell tree then switch to optimizing a mutation tree
coverage_method = "zinb"  # Determines from which distribution the coverage is sampled. Can be "zinb", "poisson", "geometric" or from a real data sample
run_tree_inference = True  # Whether to run tree inference after generating the simulated data

default_params = {
    'dropout': config["dropout_alpha"]/(config["dropout_alpha"] + config["dropout_beta"]),
    'overdispersion_Het': config["overdispersion_h"],
    'overdispersion_Hom': config["overdispersion"],
    'error_rate': config["error_rate"],
    'coverage_mean': 60,
    'coverage_zero_inflation': 0.39,
    'coverage_dispersion': 5.88,
    "CNV_fraction": 0.0,
    "homoplasy_fraction": 0.0
}

param_sets = {
    'dropout': [0, 0.2, 0.4, 0.6],
    'overdispersion_Het': [3, 6, 10, 100],
    'overdispersion_Hom': [3, 6, 10, 100],
    'error_rate': [0.001, 0.01, 0.05, 0.1],
    'coverage_mean': [10, 30, 60, 100],
    'coverage_zero_inflation': [0, 0.2, 0.4, 0.6],
    'coverage_dispersion': [1, 2, 5, 10],
    'CNV_fraction': [0, 0.2, 0.5, 0.8],
    "homoplasy_fraction": [0, 0.1, 0.2, 0.5] # Fraction of loci that are affected twice by independent mutations
}

# uncomment to run data simulations for different parameter settings
# for clone in clones:
#     for num_cells, num_mut in zip(n_cells_list, n_mut_list):
#         for param_name, param_values in param_sets.items():
#             for param_value in param_values:
#                 params = default_params.copy()
#                 params[param_name] = param_value
#
#                 param_str = f"{param_name}_{param_value}".replace('.', '_')
#                 data_path = f"../data/simulated_data/{num_cells}c{num_mut}m{clone}_param_testing/{param_str}"
#
#                 generate_comparison_data(
#                     num_cells, num_mut, num_tests, path=data_path, n_clones=clone,
#                     coverage_distribution=coverage_method, **params
#                 )
#
#                 if run_tree_inference:
#                     path_results = os.path.join(data_path, "sciterna")
#                     generate_sciterna_simulation_results(
#                         path=data_path, pathout=path_results, n_tests=num_tests,
#                         tree_space=tree_space,
#                         flipped_mutation_direction=flipped_mutation_direction,
#                         n_keep=num_mut, n_rounds=n_rounds
#                     )
for clone in clones:
    for num_cells, num_mut in zip(n_cells_list, n_mut_list):
        data_path = f"../data/simulated_data/{num_cells}c{num_mut}m{clone}"
        generate_comparison_data(num_cells, num_mut, num_tests, path=data_path, n_clones=clone, coverage_distribution=coverage_method)
        if run_tree_inference:
            path_results = os.path.join(data_path, "sciterna")
            generate_sciterna_simulation_results(path=data_path, pathout=path_results, n_tests=num_tests,
                                                 tree_space=tree_space,
                                                 flipped_mutation_direction=flipped_mutation_direction,
                                                 n_keep=num_mut, n_rounds=n_rounds)
