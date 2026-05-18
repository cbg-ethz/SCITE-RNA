"""
This script processes bootstrap samples of trees to generate a consensus parent vector from simulated data.
"""

def normalize_split(split):
    A, B = map(set, split)

    # Remove trivial splits (less than 2 leaves on either side)
    if len(A) < 2 or len(B) < 2:
        return None

    # Create frozensets for each side
    side1 = frozenset(A)
    side2 = frozenset(B)

    # Return a canonical bipartition as frozenset of frozensets (order insensitive)
    return frozenset([side1, side2])

def get_splits(ct, node, labels=None):
    splits = []
    for child in ct.dfs(node):
        if ct.isleaf(child):
            continue

        child_leaves = [leaf for leaf in ct.leaves(child)]
        other_partition = [leaf for leaf in labels if leaf not in child_leaves]

        if 0 < len(child_leaves) < len(labels):
            splits.append((child_leaves, other_partition))

    return splits

def annotate_clade_frequencies(tree, split_counter, total_trees):
    """
    Annotates each internal node in the consensus tree with its clade frequency
    based on provided split_counter.
    Frequencies are stored in node.label
    """
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            continue

        # Get leaf labels on each side of the split
        node_leaves = set(int(leaf.taxon.label) for leaf in node.leaf_iter())
        rest_leaves = set(int(leaf.taxon.label) for leaf in tree.leaf_node_iter()) - node_leaves

        if len(node_leaves) < 2 or len(rest_leaves) < 2:
            continue  # Not a valid split

        split = frozenset([frozenset(node_leaves), frozenset(rest_leaves)])
        freq = split_counter.get(split, 0) / total_trees
        node.label = f"{freq:.3f}"

def tree_to_parent_vector(tree):
    """
    Converts a DendroPy Tree object to a parent vector.
    Leaf nodes retain their labels (assumed to be string representations of integers).
    Internal nodes are assigned new labels (as strings of integers) starting from the max leaf label + 1.

    Returns:
        parent_vector: list where index i contains the parent of node i (both as integers)
        node_index_map: dict mapping each node object to its string label
    """
    node_index_map = {}

    leaf_labels = [int(node.taxon.label) for node in tree.leaf_node_iter()]
    max_leaf_label = max(leaf_labels)

    for node in tree.leaf_node_iter():
        node_index_map[node] = node.taxon.label  # Keep leaf labels (as str)

    next_internal_label = max_leaf_label + 1
    for node in tree.postorder_node_iter():
        if node not in node_index_map:
            node_index_map[node] = str(next_internal_label)
            next_internal_label += 1

    parent_dict = {}
    for parent in tree.preorder_node_iter():
        parent_id = int(node_index_map[parent])
        for child in parent.child_node_iter():
            child_id = int(node_index_map[child])
            parent_dict[child_id] = parent_id

    total_nodes = next_internal_label
    parent_vector = [-1] * total_nodes  # Initialize with -1 for all nodes

    for child_id, parent_id in parent_dict.items():
        parent_vector[child_id] = parent_id

    return parent_vector, node_index_map

def to_newick(ct, node):
    if ct.isleaf(node):
        return f"{node}"
    children = ct.children(node)
    return "(" + ",".join(to_newick(ct, child) for child in children) + f"){node}"

def main():
    import dendropy
    import os
    import json
    import argparse
    import numpy as np
    from collections import Counter
    from dendropy import Tree, TreeList, TaxonNamespace
    from src_python.cell_tree import CellTree

    parser = argparse.ArgumentParser(description="Convert bootstrap samples of trees to consensus parent vectors")
    parser.add_argument("--input_folder", type=str, help="Path to the input folder containing bootstrap tree files.", default="BT_S2") # "50c500m"
    parser.add_argument("--base_path", type=str, help="Base path for the files", default=r"../data/results")
    parser.add_argument("--model", type=str, help="Model used for the bootstrap samples.", default="sciterna")
    parser.add_argument("--simulated", type=bool, help="Run on simulated data.", default=False)
    parser.add_argument("--n_samples", type=int, help="Number of simulated samples to process.", default=100)
    parser.add_argument("--round", type=int, help="Which round to use. Each round updates optimized SNV specific and global parameters like dropout probabilities", default=1)
    parser.add_argument("--n_bootstrap", type=int, help="Number of bootstrap samples to process.", default=200)
    args = parser.parse_args()

    model = args.model
    n_samples = args.n_samples
    round = args.round
    n_bootstrap = args.n_bootstrap
    input_folder = args.input_folder
    base_path = args.base_path
    simulated = args.simulated

    for s in range(n_samples):

        if simulated:
            path = os.path.join(base_path, rf"{input_folder}/{model}_{s}_bootstrap")
        else:
            path = os.path.join(base_path, rf"{input_folder}/{model}_bootstrap")

        taxa = TaxonNamespace()
        trees = TreeList(taxon_namespace=taxa)
        split_counter = Counter()
        for test in range(n_bootstrap):
            path_parent = os.path.join(path, f"{model}_parent_vec", f"{model}_parent_vec_{round}r{test}.txt")
            path_selected = os.path.join(path, f"{model}_selected_loci", f"{model}_selected_loci_{round}r{test}.txt")

            if not os.path.exists(path_parent):
                print(f"Missing files for sample {s}, test {test}. Skipping.")
                continue

            parent_vec = np.loadtxt(path_parent, dtype=int)

            if not os.path.exists(path_selected):
                selected_mutations = []
            else:
                selected_mutations = np.loadtxt(path_selected, dtype=int)

            n_cells = int(((len(parent_vec) + 1) / 2))

            ct = CellTree(n_cells=n_cells, n_mut=len(selected_mutations))
            ct.use_parent_vec(parent_vec)

            newick_str = to_newick(ct, ct.main_root) + ";"
            tree = Tree.get(data=newick_str, schema="newick", taxon_namespace=taxa)
            trees.append(tree)

            # Count normalized splits
            labels = list(ct.leaves(ct.main_root))
            raw_splits = get_splits(ct, ct.main_root, labels)
            for split in raw_splits:
                norm = normalize_split(split)
                if norm is not None:
                    split_counter[norm] += 1

        consensus_tree = trees.consensus(min_freq=0.001, resolve_polytomies=True, suppress_unifurcations=True)
        consensus_tree.resolve_polytomies(update_bipartitions=False)
        annotate_clade_frequencies(consensus_tree, split_counter, total_trees=n_bootstrap)
        print(consensus_tree.as_ascii_plot(show_internal_node_labels=True))

        consensus_parent_vec, node_index_map = tree_to_parent_vector(consensus_tree)

        support_values = {}
        for node in consensus_tree.postorder_node_iter():
            if not node.is_leaf():
                node_id = int(node_index_map[node])
                try:
                    support = float(node.label)
                    support_values[str(node_id)] = support
                except (ValueError, TypeError):
                    continue

        if simulated:
            support_values_file = os.path.join(path, "..", f"{model}_consensus_parent_vec",
                                               f"{model}_support_values_{round}r{s}.json")
        else:
            support_values_file = os.path.join(path, f"{model}_consensus_parent_vec",
                                               f"{model}_support_values_{round}r.json")

        os.makedirs(os.path.dirname(support_values_file), exist_ok=True)

        with open(support_values_file, "w") as f:
            json.dump(support_values, f, indent=4)

        if simulated:
            np.savetxt(
                os.path.join(path, "..", f"{model}_consensus_parent_vec", f"{model}_parent_vec_{round}r{s}.txt"),
                consensus_parent_vec, fmt='%d')
        else:
            np.savetxt(
                os.path.join(path, f"{model}_consensus_parent_vec", f"{model}_parent_vec_{round}r.txt"),
                consensus_parent_vec, fmt='%d')
            break

if __name__ == "__main__":
    main()