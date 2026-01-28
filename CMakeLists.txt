cmake_minimum_required(VERSION 3.27)
project(SCITE-RNA)

set(CMAKE_TOOLCHAIN_FILE "../data/")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(src_cpp)
include_directories("./data/yixuan-LBFGSpp-9db7a6d/include")
include_directories("./data/eigen3")

add_executable(SCITE-RNA
        src_cpp/cell_tree.cpp
        src_cpp/cell_tree.h
        src_cpp/mutation_tree.cpp
        src_cpp/mutation_tree.h
        src_cpp/swap_optimizer.cpp
        src_cpp/swap_optimizer.h
        src_cpp/mutation_filter.cpp
        src_cpp/mutation_filter.h
        src_cpp/generate_results.cpp
        src_cpp/generate_results.h
        src_cpp/utils.cpp
        src_cpp/utils.h
#        generate_results_cpp/comparison_parameters.cpp
       generate_results_cpp/comparison_num_clones.cpp
#        generate_results_cpp/real_data_processing.cpp
#        generate_results_cpp/comparison_bootstrapping.cpp
#        generate_results_cpp/space_switching_results_postprocessing.cpp
#        generate_results_cpp/comparison_tree_spaces_switching.cpp
)
