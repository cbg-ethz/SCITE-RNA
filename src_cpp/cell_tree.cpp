/*
Defines cell lineage trees and how they are optimized.
*/

#include <chrono>
#include <random>
#include <vector>
#include <algorithm>
#include <stack>
#include <stdexcept>
#include <cassert>
#include <limits>

#include "cell_tree.h"
#include "utils.h"

// CELL TREE FUNCTIONS

// initialize the cell tree
CellTree::CellTree(int n_cells, int n_mut, bool reversible_mut)
        : n_cells(n_cells), n_mut(n_mut), reversible(reversible_mut) {
    if (n_cells < 3) {
        std::cerr << "Warning: Cell tree too small, nothing to explore" << std::endl;
    }
    n_vtx = 2 * n_cells - 1;
    parent_vector_ct = std::vector<int>(n_vtx, -1);
    flipped = std::vector<bool>(n_mut, false);
    children_list_ct.resize(n_vtx + 1);
    for (size_t i = n_cells; i <= n_vtx; ++i) {
        children_list_ct[i].reserve(2); // Reserve space for exactly two children per hidden node
    }
    mut_loc.resize(n_mut);
    bool initialTree = true;
    randSubtree(initialTree); // Initialize with a completely random structure
//    std::vector<int> new_parent_vector = {53,56,88,73,63,53,55,57,64,79,69,58,57,68,75,67,62,51,64,72,74,65,66,58,60,54,72,74,61,61,
//    77,56,81,73,52,54,50,69,52,51,76,62,59,65,63,86,55,97,50,78,60,70,67,77,59,71,70,80,88,82,
//    66,91,71,68,86,93,80,85,84,79,83,85,83,90,75,76,78,92,81,89,82,84,87,96,87,94,90,89,95,91,
//    95,92,93,94,96,97,98,98,-1};
//    useParentVec(new_parent_vector);
}

// use parent vector to initialize cell tree
void CellTree::useParentVec(const std::vector<int>& new_parent_vec, int main_r) {
    if (new_parent_vec.size() != static_cast<size_t>(n_vtx)) {
        throw std::invalid_argument("Parent vector must have the same length as number of vertices.");
    }

    parent_vector_ct = new_parent_vec;
    children_list_ct = std::vector<std::vector<int>>(n_vtx + 1); // Reinitialize the children_list

    for (int vtx = 0; vtx < n_vtx; ++vtx) {
        addChild(new_parent_vec[vtx],vtx);
    }

    if (main_r == -1) {
        main_root_ct = roots()[0];
    } else if (parent_vector_ct[main_root_ct] != -1) {
        throw std::invalid_argument("Provided main root is not a root.");
    } else {
        main_root_ct = main_r;
    }
}


// create random tree
void CellTree::randSubtree(bool initialTree, std::vector<int>* leaves, std::vector<int>* internals) {
    // Determine the leaf and internal vertices
    std::vector<int> local_leaves;
    std::vector<int> local_internals;

    if (leaves == nullptr || internals == nullptr) {
        local_leaves.reserve(n_cells);
        local_internals.reserve(n_cells - 1);
        for (int i = 0; i < n_cells; ++i) {
            local_leaves.push_back(i);
        }
        for (int i = n_cells; i < n_vtx; ++i) {
            local_internals.push_back(i);
        }
    } else {
        if (leaves->size() != internals->size() + 1) {
            throw std::invalid_argument("There must be exactly one more leaf than internals.");
        }
        local_leaves = *leaves;
        local_internals = *internals;
    }

//    std::sort(local_leaves.begin(), local_leaves.end());
//    std::sort(local_internals.begin(), local_internals.end());

    // Randomly assign two children to an internal vertex
    std::mt19937 rng = create_rng();
    for (int parent : local_internals) {
        std::shuffle(local_leaves.begin(), local_leaves.end(), rng);
        int child1 = local_leaves.front();
        local_leaves.erase(local_leaves.begin());
        int child2 = local_leaves.front();
        local_leaves.erase(local_leaves.begin());

        assignParent(child1, parent);
        assignParent(child2, parent);
        local_leaves.push_back(parent);
    }

    reroot(local_internals.back());
    if (initialTree) {
        assignParent(local_leaves[0], -1);
    }
}

// convert mutation tree to cell tree
void CellTree::fitMutationTree(MutationTree mt) {
    assert(this->n_cells == mt.n_cells);

    // reinitialize empty children list
    for (size_t i = n_cells; i <= n_vtx; ++i) {
        children_list_ct[i].clear();
        children_list_ct[i].reserve(2);
    }

    std::vector<int> mrca(mt.n_vtx, -1); // Most recent common ancestor of cells below a mutation node
    int next_internal = n_cells;

    for (int mtvx : mt.rdfs(mt.main_root_mt)) {
        std::vector<int> leaves;
        for (int child : mt.children_list_mt[mtvx]) {
            if (mrca[child] != -1) {
                leaves.push_back(mrca[child]);
            }
        }

        for (int i = 0; i < mt.cell_loc.size(); ++i) {
            if (mt.cell_loc[i] == mtvx) {
                leaves.push_back(i);
            }
        }

        if (leaves.empty()) { // No cell below, nothing to do
            continue;
        } else if (leaves.size() == 1) { // One cell below, no internal node added
            mrca[mtvx] = leaves[0];
        } else if (leaves.size() > 1) { // More than one cell below, add new internal node(s)
            std::vector<int> internals(leaves.size() - 1);
            std::iota(internals.begin(), internals.end(), next_internal);
            randSubtree(false, &leaves, &internals);
            mrca[mtvx] = internals.back();
            next_internal += static_cast<int>(internals.size());
        }
    }
    children_list_ct[n_vtx] = std::vector<int>{main_root_ct};
}

// initialize LLR, joint LL, LLR internal nodes, mutation locations
void CellTree::fitLlh(const std::vector<std::vector<double>>& llh_1,
                       const std::vector<std::vector<double>>& llh_2) {
    if (llh_1.size() != llh_2.size() || llh_1[0].size() != llh_2[0].size()) {
        throw std::invalid_argument("llh_1 and llh_2 must have the same shape.");
    }

    if (llh_1.size() != n_cells) {
        std::cerr << "Number of cells does not match the row count of the llh matrix." << std::endl;
    } else if (llh_1[0].size() != n_mut) {
        std::cerr << "Number of mutations does not match the column count of the llh matrix." << std::endl;
    }

    llr = std::vector<std::vector<double>>(n_vtx, std::vector<double>(n_mut));

    // Calculate llr of genotype 1 and 2
    for (int i = 0; i < n_cells; ++i) {
        for (int j = 0; j < n_mut; ++j) {
            llr[i][j] = llh_2[i][j] - llh_1[i][j];
        }
    }

    loc_joint_1 = std::vector<double>(n_mut);
    loc_joint_2 = std::vector<double>(n_mut);

    // Calculate joint likelihood of each locus when all cells have genotype 1 or 2
    for (int j = 0; j < n_mut; ++j) {
        loc_joint_1[j] = 0;
        loc_joint_2[j] = 0;
        for (int i = 0; i < n_cells; ++i) {
            loc_joint_1[j] += llh_1[i][j];
            loc_joint_2[j] += llh_2[i][j];
        }
    }

    // pre compute joint 1 genotype for tree joint llh calculation
    joint_1_sum = std::accumulate(loc_joint_1.begin(), loc_joint_1.end(), 0.0);

    // Assign mutations to optimal locations
    updateAll();
}

// update LLR and Mutation location
void CellTree::updateAll() {
    updateLlr();
    updateMutLoc();
}

// update the LLR
void CellTree::updateLlr() {
    for (int rt : roots()) {
        for (int vtx : rdfs(rt)) {
            if (isLeaf(vtx)) {
                continue;
            }
            // LLR at internal vertex is the sum of LLR of both children
            std::vector<int> children = children_list_ct[vtx];
            for (int j = 0; j < llr[0].size(); ++j) {
                llr[vtx][j] = 0.0;
                for (int child : children) {
                    llr[vtx][j] += llr[child][j];
                }
            }
        }
    }
}

// update the optimal mutation locations
void CellTree::updateMutLoc() {
    if (reversible) {
        std::vector<int> locs_neg(n_mut);
        std::vector<int> locs_pos(n_mut);
        std::vector<double> llhs_neg(n_mut);
        std::vector<double> llhs_pos(n_mut);
        std::vector<double> loc_joint(n_mut);

        for (int j = 0; j < n_mut; ++j) {
            double minValue = std::numeric_limits<double>::infinity();
            double maxValue = -std::numeric_limits<double>::infinity();
            for (int i = 0; i < n_vtx; ++i) {
                if (llr[i][j] < minValue) {
                    minValue = llr[i][j];
                    locs_neg[j]  = i;
                }
                if (llr[i][j] > maxValue) {
                    maxValue = llr[i][j];
                    locs_pos[j] = i;
                }
            }
            llhs_neg[j] = loc_joint_2[j] - minValue;
            llhs_pos[j] = loc_joint_1[j] + maxValue;
        }

        std::vector<bool> neg_larger(n_mut);
        for (int j = 0; j < n_mut; ++j) {
            neg_larger[j] = llhs_neg[j] > llhs_pos[j];
            mut_loc[j] = neg_larger[j] ? locs_neg[j] : locs_pos[j];
            flipped[j] = neg_larger[j];
            loc_joint[j] = neg_larger[j] ? llhs_neg[j] : llhs_pos[j];
        }

        joint = std::accumulate(loc_joint.begin(), loc_joint.end(), 0.0);
    }
    else {
        joint = 0.0;
        for (int j = 0; j < n_mut; ++j) {
            // Find the index of the maximum element in the column
            auto max_it = std::max_element(llr.begin(), llr.end(), [j](const std::vector<double>& a, const std::vector<double>& b) {
                return a[j] < b[j];
            });

            // Compute the index of the maximum element in the column
            int max_idx = static_cast<int>(std::distance(llr.begin(), max_it));
            mut_loc[j] = max_idx;

            joint += llr[max_idx][j] + loc_joint_1[j];
        }
    }
}

// exhaustive greedy pruning and reattaching of subtrees
void CellTree::exhaustiveOptimize(bool leaf_only) {
    updateAll();

    std::vector<int> sr_candidates;
    sr_candidates.reserve(n_vtx);
    if (leaf_only) {
        for (int i = 0; i < n_cells; ++i) {
            sr_candidates.push_back(i);
        }
    } else {
        for (int i = 0; i < n_vtx; ++i) {
            sr_candidates.push_back(i);
        }
    }

    // randomize order of candidate roots

    std::mt19937 rng = create_rng(); // Persistent RNG for this function

    std::shuffle(sr_candidates.begin(), sr_candidates.end(), rng);

    for (int sr : sr_candidates) {
        if (sr == main_root_ct) {
            continue;
        }

//        std::cout << sr << " " << parent_vector_ct[sr] << " " << sibling(sr) << std::endl;
        binaryPrune(sr);
        updateLlr();
        greedyInsert(parent_vector_ct[sr], sr);
    }
    updateAll();
}

// Translation of greedy_insert_experimental function
void CellTree::greedyInsert(int anchor, int subroot) {
    auto [bestTarget, bestJoint] = searchInsertionLoc(main_root_ct, anchor, subroot);
    insert(anchor, bestTarget);
}

// search optimal insertion location of pruned subroot with anchor
std::pair<int, double> CellTree::searchInsertionLoc(int target, int anchor, int subroot) {
    std::stack<std::tuple<int, bool, std::vector<double>, std::vector<double>, std::vector<double>>> stack;
    stack.push({target, false, {}, {}, {}});

    std::vector<int> best_targets;
    double best_joint = -std::numeric_limits<double>::infinity();
    int check = 0;
    bool was_leaf = true; // only anchor row of llr changes for attaching above the initial target
    std::vector<double> current_llr_max_without_anchor;
    current_llr_max_without_anchor.resize(llr[0].size());
    std::vector<double> current_llr_min_without_anchor;
    current_llr_min_without_anchor.resize(llr[0].size());


    // Calculate the maximum and minimum value for each column, skipping the anchor row
    for (size_t j = 0; j < llr[0].size(); ++j) {
        for (size_t i = 0; i < llr.size(); ++i) {
            if (static_cast<int>(i) != anchor) {
                if (llr[i][j] > current_llr_max_without_anchor[j]) {
                    current_llr_max_without_anchor[j] = llr[i][j];
                }
                if (llr[i][j] < current_llr_min_without_anchor[j]) {
                    current_llr_min_without_anchor[j] = llr[i][j];
                }
            }
        }
    }

    while (!stack.empty()) {
        check += 1;
        if (check > (n_vtx + n_cells)){
            throw std::runtime_error("greedy insert function encountered an issue");
        }
        auto [currentTarget, visited, originalLlr, originalMaxWithoutAnchor, originalMinWithoutAnchor] = stack.top();
        stack.pop();

        if (!visited) {
            llr[anchor] = addVectors(llr[subroot], llr[currentTarget]);

            double currentJoint;
            if (reversible) {
                // much faster, but maybe less intuitive joint tree llh calculation than updateMutLoc
                auto [updated_currentJoint, updated_was_leaf, updated_current_llr_max_without_anchor, updated_current_llr_min_without_anchor] =
                        jointCalculationFlipped(was_leaf, anchor, current_llr_max_without_anchor, current_llr_min_without_anchor,
                                                   currentTarget, originalLlr);
                currentJoint = updated_currentJoint;
                was_leaf = updated_was_leaf;
                current_llr_max_without_anchor = updated_current_llr_max_without_anchor;
                current_llr_min_without_anchor = updated_current_llr_min_without_anchor;
            } else {
                auto [updated_currentJoint, updated_was_leaf, updated_current_llr_max_without_anchor] =
                        jointCalculationNotFlipped(was_leaf, anchor, current_llr_max_without_anchor,
                                                   currentTarget, originalLlr);
                currentJoint = updated_currentJoint;
                was_leaf = updated_was_leaf;
                current_llr_max_without_anchor = updated_current_llr_max_without_anchor;
            }

//            updateMutLoc();
//            double currentJoint = joint;

            if (currentJoint == best_joint) {
                best_targets.push_back(currentTarget);
            }
            if (currentJoint > best_joint) {
                best_targets.clear();
                best_targets.push_back(currentTarget);
                best_joint = currentJoint;
            }

            if (isLeaf(currentTarget)){
                was_leaf = true;
            } else {
                std::vector<double> currentLlr = llr[currentTarget];
                llr[currentTarget] = addVectors(llr[currentTarget], llr[subroot]);

                stack.emplace(currentTarget, true, currentLlr, current_llr_max_without_anchor, current_llr_min_without_anchor);
                std::vector<int> childrenNodes = children_list_ct[currentTarget];

                for (int child : childrenNodes) {
                    stack.push({child, false, currentLlr, {}, {}});
                }
            }
        }

        if (visited) {
            llr[currentTarget] = originalLlr;
            current_llr_max_without_anchor = originalMaxWithoutAnchor;
            current_llr_min_without_anchor = originalMinWithoutAnchor;
        }
    }

    // choose random best target

//    std::random_device rd;   // Seed generator
//    std::mt19937 gen(rd());  // Mersenne Twister RNG
    std::mt19937 rng = create_rng();
    std::uniform_int_distribution<std::vector<int>::size_type> dis(0, best_targets.size() - 1); // Uniform distribution over the vector indices
//    for (int element : best_targets) {
//        std::cout << element << " ";
//    }
//    std::cout << std::endl;
    int random_best_target = best_targets[dis(rng)];
//    auto min_iter = std::min_element(best_targets.begin(), best_targets.end());
//    int random_best_target = *min_iter;
    return {random_best_target, best_joint};
}

// calculate the tree joint llh if flipped mutation directions are allowed
std::tuple <double, bool, std::vector<double>, std::vector<double>> CellTree::jointCalculationFlipped(bool was_leaf, int anchor,
                                                                                    std::vector<double> current_llr_max_without_anchor,
                                                                                    std::vector<double> current_llr_min_without_anchor,
                                                                                    int currentTarget,
                                                                                    std::vector<double> originalLlr) {
    double Joint = 0.0;
    if (was_leaf) { // only the anchor llr has changed in that case
        for (int i = 0; i < n_mut; ++i) {
            double max_val = loc_joint_1[i] + std::max(current_llr_max_without_anchor[i], llr[anchor][i]);
            double min_val = loc_joint_2[i] - std::min(current_llr_min_without_anchor[i], llr[anchor][i]);

            if (max_val > min_val){
                Joint += max_val;
            } else{ //flipped mutation
                Joint += min_val;
            }
        }
        was_leaf = false;
    } else { // the llr of anchor and previous current_target changed. Update current_llr_max_without_anchor as well
        std::vector<double> new_llr_node = llr[parent_vector_ct[currentTarget]];
        for (int i = 0; i < n_mut; ++i) {
            double pos_value = loc_joint_1[i] + std::max({current_llr_max_without_anchor[i], new_llr_node[i], llr[anchor][i]});
            double neg_value = loc_joint_2[i] - std::min({current_llr_min_without_anchor[i], new_llr_node[i], llr[anchor][i]});
            bool neg_larger = (neg_value) > (pos_value);
            if ((originalLlr[i] < current_llr_max_without_anchor[i]) and not neg_larger) { // if the maximum was not in the node row that changed and meximum is larger than the minimum->not flipped
                Joint += pos_value;
            }
            else if ((originalLlr[i] > current_llr_min_without_anchor[i]) and neg_larger){// if the minimum was not in the node row that changed and flipped genotypes
                Joint += neg_value;
            }
            else { // if it was, we have to recalculate the maximum of the row with and without anchor
                double max_in_column = -std::numeric_limits<double>::infinity();
                double min_in_column = std::numeric_limits<double>::infinity();

                for (size_t j = 0; j < n_vtx; ++j) {
                    double value = llr[j][i];
                    max_in_column = std::max(max_in_column, value);
                    min_in_column = std::min(min_in_column, value);
                }
                double max_val = loc_joint_1[i] + max_in_column;
                double min_val = loc_joint_2[i] - min_in_column;

                if (max_val > min_val){
                    Joint += max_val;
                } else{ //flipped mutation
                    Joint += min_val;
                }

                // if the value at anchor is smaller than the maximum, the maximum without anchor stays the same
                if (llr[anchor][i] < max_in_column){
                    current_llr_max_without_anchor[i] = max_in_column;
                }
                //  otherwise we have to recalculate without anchor
                else{
                    double max_in_column_without_anchor = -std::numeric_limits<double>::infinity();
                    for (size_t j = 0; j < n_vtx; ++j) {
                        if (j != anchor){
                            double value = llr[j][i];
                            max_in_column_without_anchor = std::max(max_in_column_without_anchor, value);
                        }
                    }
                    current_llr_max_without_anchor[i] = max_in_column_without_anchor;
                }
                // if the value at anchor is larger than the minimum, the minimum without anchor stays the same
                if (llr[anchor][i] > min_in_column){
                    current_llr_min_without_anchor[i] = min_in_column;
                }
                //  otherwise we have to recalculate without anchor
                else{
                    double min_in_column_without_anchor = std::numeric_limits<double>::infinity();
                    for (size_t j = 0; j < n_vtx; ++j) {
                        if (j != anchor){
                            double value = llr[j][i];
                            min_in_column_without_anchor = std::min(min_in_column_without_anchor, value);
                        }
                    }
                    current_llr_min_without_anchor[i] = min_in_column_without_anchor;
                }
            }
            if (new_llr_node[i] > current_llr_max_without_anchor[i]) {
                current_llr_max_without_anchor[i] = new_llr_node[i];
            }
            if (new_llr_node[i] < current_llr_min_without_anchor[i]) {
                current_llr_min_without_anchor[i] = new_llr_node[i];
            }
        }
    }
    double currentJoint = Joint;
    return {currentJoint, was_leaf, current_llr_max_without_anchor, current_llr_min_without_anchor};
}


// calculate the joint likelihood of the tree and update variables
std::tuple <double, bool, std::vector<double>> CellTree::jointCalculationNotFlipped(bool was_leaf, int anchor,
                                            std::vector<double> current_llr_max_without_anchor, int currentTarget,
                                            std::vector<double> originalLlr) {
    double Joint = 0.0;
    if (was_leaf) { // only the anchor llr has changed in that case
        for (size_t i = 0; i < llr[anchor].size(); ++i) {
            double max_val = std::max(current_llr_max_without_anchor[i], llr[anchor][i]);
            Joint += max_val;
        }
        was_leaf = false;
    } else { // the llr of anchor and previous current_target changed. Update current_llr_max_without_anchor as well
        std::vector<double> new_llr_node = llr[parent_vector_ct[currentTarget]];
        for (size_t i = 0; i < n_mut; ++i) {
            if ((originalLlr)[i] < current_llr_max_without_anchor[i]) { // if the maximum was not in the node row that changed
                Joint += std::max({current_llr_max_without_anchor[i], new_llr_node[i], llr[anchor][i]});
            }
            else { // if it was, we have to recalculate the maximum of the row with and without anchor
                double max_in_column = -std::numeric_limits<double>::infinity();

                for (size_t j = 0; j < n_vtx; ++j) {
                    double value = llr[j][i];
                    max_in_column = std::max(max_in_column, value);
                }
                Joint += max_in_column;

                // if the value at anchor is smaller than the maximum, the maximum without anchor stays the same
                if (llr[anchor][i] < max_in_column){
                    current_llr_max_without_anchor[i] = max_in_column;
                }
                //  otherwise we have to recalculate without anchor
                else{
                    max_in_column = -std::numeric_limits<double>::infinity();
                    for (size_t j = 0; j < n_vtx; ++j) {
                        if (j != anchor){
                            double value = llr[j][i];
                            max_in_column = std::max(max_in_column, value);
                        }
                    }
                    current_llr_max_without_anchor[i] = max_in_column;
                }

            }
            if (new_llr_node[i] > current_llr_max_without_anchor[i]) {
                current_llr_max_without_anchor[i] = new_llr_node[i];
            }
        }
    }
    double currentJoint = Joint + joint_1_sum;
    return {currentJoint, was_leaf, current_llr_max_without_anchor};
}

// BASIC TREE FUNCTIONS

// assign child to new parent
void CellTree::assignParent(int child, int new_parent) {
    removeChild(parent_vector_ct[child], child);
    parent_vector_ct[child] = new_parent;
    addChild(new_parent, child);
}

// add new main root
void CellTree::reroot(int new_main_root) {
    if (!isRoot(new_main_root)) {
        prune(new_main_root);
    }
    main_root_ct = new_main_root;
}

// remove child from children list of binary tree
void CellTree::removeChild(int parent, int child) {
    if (parent == -1) {
        parent = n_vtx;
    }

    auto& childList = children_list_ct[parent];

    // Check if the first element is the child to be removed
    if (!childList.empty() && childList[0] == child) {
        childList[0] = childList.back(); // Replace it with the other child
        childList.pop_back();
    }
    // Check if the second element is the child to be removed
    else if (childList.size() == 2 && childList[1] == child) {
        childList.pop_back(); // Remove the last child
    }
}

// add child to children_list of parent node
void CellTree::addChild(int parent, int child) {
    if (parent == -1) {
        parent = n_vtx;
    }
    children_list_ct[parent].push_back(child);
}

// check if node is a leaf
bool CellTree::isLeaf(int node) {
    return children_list_ct[node].empty();
}

// roots have parent -1
bool CellTree::isRoot(int vtx) {
    return parent_vector_ct[vtx] == -1;
}

// prune subroot
void CellTree::prune(int subroot){
    assignParent(subroot, -1);
}

// prune the parent of subroot and reattach the sibling to the grandparent
void CellTree::binaryPrune(int subroot) {
    int sib = sibling(subroot);
    splice(sib);
}

// Splice a subtree
void CellTree::splice(int subroot) {
    if (isRoot(subroot)) {
        throw std::invalid_argument("Cannot splice a root");
    }

    int parent = parent_vector_ct[subroot];
    int grandparent = parent_vector_ct[parent];

    prune(parent);
    assignParent(subroot, grandparent);

    if (grandparent == -1) {
        main_root_ct = subroot;
    }
}

// Insert anchor into the tree above target and assign target to anchor
void CellTree::insert(int anchor, int target) {
    int parent = parent_vector_ct[target];
    assignParent(target, anchor);          // target becomes the child of subroot
    assignParent(anchor, parent); // subroot becomes the parent of target's current parent

    if (target == main_root_ct) {
        main_root_ct = anchor;                // Update main_root if necessary
    }
}

// Get the siblings of a node
int CellTree::sibling(int node) {
    int sib = 0;
    int parent = parent_vector_ct[node];
    if (parent != -1) { // If the node has a parent, that is not the root
        for (int s : children_list_ct[parent]) {
            if (s != node) {
                sib = s;
            }
        }
    }
    return sib;
}

// return all roots
std::vector<int> CellTree::roots() {
    return children_list_ct.back();
}

// depth first search
std::vector<int> CellTree::dfs(int subroot) {
    std::vector<int> stack = {subroot};
    std::vector<int> output;
    output.reserve(n_vtx);
    int check = 0;

    while (!stack.empty()) {
        check += 1;
        if (check > n_vtx){
            throw std::runtime_error("Problem with dfs");
        }

        int node = stack.back();
        stack.pop_back();
        output.push_back(node);

        // Add children of the current node to the stack
        for (int child : children_list_ct[node]) {
            stack.push_back(child);
        }
    }
    return output;
}

// reverse depth first search
std::vector<int> CellTree::rdfs(int subroot) {
    std::vector<int> stack = {subroot};
    std::vector<int> output;
    output.reserve(n_vtx);
    int check = 0;

    while (!stack.empty()) {
        check += 1;
        if (check > n_vtx){
            std::cerr << "Problem with rdfs." << std::endl;
            throw std::runtime_error("rdfs function encountered an issue");
        }

        int node = stack.back();
        stack.pop_back();
        output.push_back(node);

        // Add children of the current node to the stack
        for (int child : children_list_ct[node]) {
            stack.push_back(child);
        }
    }

    // Reverse the output to get nodes in reversed DFS order
    std::reverse(output.begin(), output.end());

    return output;
}

// get all leaves of the tree starting from subroot
std::vector<int> CellTree::leaves(int subroot) {
    std::vector<int> leaf_vertices;
    for (int vtx : dfs(subroot)) {
        if (isLeaf(vtx)) {
            leaf_vertices.push_back(vtx);
        }
    }
    return leaf_vertices;
}


// HELPER FUNCTIONS


std::mt19937& CellTree::create_rng() {
    static std::mt19937 rng([] {
        load_config("../config/config.yaml");
        return std::mt19937(std::stoi(config_variables["random_seed"]));
    }());
    return rng;
}
