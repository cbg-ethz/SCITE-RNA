import argparse
import numpy as np
from Bio import Phylo

def read_cell_order(cell_file):
    cells = np.loadtxt(cell_file, dtype=str)
    if cells[0, 0] != "Cell":
        raise SystemExit(f"Expected first line of cell file {cell_file} to start with 'Cell'")
    cells = list(cells[1:, 0])
    return cells

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--newick_file", required=True)
    p.add_argument("--cell_file", required=True)
    p.add_argument("--output_file", required=True)
    args = p.parse_args()

    cells = read_cell_order(args.cell_file)
    n_leaves = len(cells)
    leaf_name_to_index = {name: i for i, name in enumerate(cells)}

    tree = Phylo.read(args.newick_file, "newick")

    terminals = tree.get_terminals()
    for t in terminals:
        if t.name is None:
            raise SystemExit("Found terminal node without a name in the tree.")
        if t.name not in leaf_name_to_index:
            raise SystemExit(f"Leaf name '{t.name}' not found in cell file {args.cell_file}.")

    internal_nodes = list(tree.get_nonterminals(order="postorder"))
    internal_map = {node: n_leaves + i for i, node in enumerate(internal_nodes)}

    # build clade -> index map (terminals + internals)
    clade_to_index = {}
    for term in terminals:
        clade_to_index[term] = leaf_name_to_index[term.name]
    for node, idx in internal_map.items():
        clade_to_index[node] = idx

    total_nodes = n_leaves + len(internal_nodes)
    parent = [None] * total_nodes

    def walk(par):
        for ch in par.clades:
            child_idx = clade_to_index[ch]
            parent_idx = clade_to_index[par]
            parent[child_idx] = parent_idx
            if not ch.is_terminal():
                walk(ch)

    root = tree.root

    walk(root)

    # mark root parent as -1
    root_idx = clade_to_index[root]
    parent[root_idx] = -1

    with open(args.output_file, "w", encoding="utf-8") as out:
        for v in parent:
            out.write(str(v) + "\n")

    print(f"Wrote parent vector for {total_nodes} nodes ({n_leaves} leaves) to {args.output_file}")

if __name__ == "__main__":
    main()