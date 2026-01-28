# SCITE-RNA

This repository contains the code and data for the paper **Phylogenetic tree inference from single-cell RNA sequencing data**. 
The code and datasets provided here enable users to replicate the experiments and figures presented in the paper, as well as to run SCITE-RNA on new data.

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
  - [Requirements](#requirements)
  - [Cloning the Repository](#cloning-the-repository)
- [Running the Model](#running-the-model)
  - [Set Model Parameters](#set-model-parameters)
  - [Simulated Data](#simulated-data)
  - [Cancer Data](#multiple-myeloma-data)
  - [Run on New Data](#run-on-new-data)
- [Generating Figures](#generating-figures)

---

## Overview

We implement a new method for reconstructing phylogenetic trees from single-cell RNA sequencing data. 
SCITE-RNA selects single-nucleotide variants (SNVs), and reconstructs a phylogenetic tree of the sequenced cells. 
We maximize the likelihood of the inferred tree by alternating between the cell lineage and mutation tree spaces until convergence is achieved in both.
This repository provides:
1. Scripts to execute SCITE-RNA. The model is split into C++ `src_cpp` and Python files `src_python`. Especially for large numbers of cells and SNVs it is recommended to use the C++ code, as it is significantly faster. The inferred trees should be comparable between the C++ and Python implementations, but as the method is stochastic likely won't produce the exact same tree. 
2. Data used in the paper are available in the `data_summary/`and `data/` directories, which contain all necessary files to reproduce the figures.
3. Visualization scripts/notebooks to generate plots as presented in the paper.


## Repository Structure

**SCITE-RNA**<br>
├── data_summary/               # Summary data files, as the raw output is quite large <br>
├── data/                       # Input data files and results <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── input_data               # Alternative and reference read counts among other files. <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;├── results                  # Inferred trees of the cancer datasets and consensus tree results <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;└── simulated_data           # Simulated data and inference results <br>
├── generate_results_cpp/       # C++ scripts to run SCITE-RNA on various datasets <br>
├── generate_results_python_r/  # Python and R scripts for simulating data, inferring trees and visualization <br>
├── phylinsic_scripts/          # Slightly adapted code to run PhylinSic https://github.com/U54Bioinformatics/PhylinSic_Project on our simulated data <br>
├── src_cpp/                    # C++ source files for SCITE-RNA <br>
├── src_python/                 # Python source files for SCITE-RNA <br>
├── config/                     # Model parameters <br>
├── CMakeLists.txt              # Primary configuration file for CMake <br>
├── README.md                   # Project overview and setup instructions <br>
├── environment.yml             # Conda environment file with required Python packages <br>
├── Snakefile                   # Snakemake file to run the PhylinSic pipeline on simulated data <br>
└── Snakefile_real_data.smk     # Snakemake file to run the PhylinSic pipeline on cancer data

## Installation

### Requirements

#### Python Libraries:
- numpy
- pandas
- scikit-learn
- matplotlib
- seaborn
- scipy
- numba
- jupyter
- pyyaml
- python-graphviz
- pygraphviz
- dendropy 

Using conda:

    conda env create -f environment.yml


#### C++ Requirements

- CMake (>= 3.27)
- LBFGS https://github.com/yixuan/LBFGSpp
- Eigen https://gitlab.com/libeigen/eigen

### Cloning the Repository

To set up the SCITE-RNA project locally:

    git clone https://github.com/cbg-ethz/SCITE-RNA.git
    cd SCITE-RNA


## Running the Model

### Set Model Parameters

If desired adjust model parameters in `config/config.yaml`.

### Simulated Data

To **generate new simulated data** with various parameter settings and different numbers of clones execute:

    generate_results_python_r/comparison_data_generation.py

It offers the option to set the number of cells and SNVs and the number of clones simulated.
The same file can also be used for tree inference. Alternatively, run the C++ version:

    generate_results_cpp/comparison_num_clones.cpp 

for tree inference (not data generation), which is a lot faster.

To **compare different optimization strategies** of tree space switching run:

    generate_results_cpp/comparison_tree_spaces_switching.cpp

Otherwise, the model will by default alternate between cell lineage and
mutation tree optimization, starting from a random cell lineage tree. 

All simulated results will be saved in `data/simulated_data/`.

### Cancer Data

To run SCITE-RNA on the cancer datasets:

Run either 

       generate_results_python_r/real_data_processing.py

or (recommended) run the faster C++ version:

       generate_results_cpp/real_data_processing.cpp

Results will be saved in `data/results/`.

### Run on New Data
To use SCITE-RNA on new data:

1. Prepare reference and alternative allele count files in `.csv` format. 
Use the format provided in `data/input_data/` as a reference, 
where columns represent cells and rows represent SNVs.

2. Set the number of bootstrap samples (optional) and run SCITE-RNA tree inference with the following script:

       generate_results_cpp/real_data_processing.cpp

3. The results are saved in `data/results`.


## Generating Figures

### Data Preparation

To reproduce the figures quickly you can use the files provided in `data` and `data_summary`. 
As the size and the number of raw data files was quite large, we produced summary statistics. To reproduce the plots presented in the paper, follow the instructions below:
        
- **Figure 3: Comparison of tree optimization strategies**

      generate_results_python_r/comparison_tree_spaces_switching.py

  If you want to rerun the full analysis first run with the desired number of cells and SNVs
  
        generate_results_python_r/comparison_data_generation.py
        generate_results_cpp/comparison_tree_spaces_switching.cpp
        generate_results_cpp/space_switching_results_postprocessing.cpp
<br>

- **Figure 4: Comparison to SClineager, PhylinSic and DENDRO including runtimes**
        
        generate_results_python_r/comparison_num_clones.ipynb
   
  If you want to rerun the full analysis first run with the desired number of clones, SNVs, cells 
 
          generate_results_python_r/comparison_data_generation.py
          generate_results_cpp/comparison_num_clones.cpp
          generate_results_python_r/comparison_clones_sclineager_dendro_sciterna.R       
      
- **Figure 5/6/7: Cancer datasets**
          
          generate_results_python_r/results_real_data.ipynb

  Rerun the full analysis with and without bootstrapping

          generate_results_cpp/real_data_processing.cpp
          generate_results_python_r/generate_consensus_parent_vector.py
 
Figures will be saved in `data/results/figures/`.