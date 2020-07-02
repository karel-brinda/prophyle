#! /usr/bin/env python3

import sys
import os

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


class TreeIndex:
    def __init__(self, tree_newick_fn):
        self.tree_newick_fn = tree_newick_fn
        self.tree = pro.load_nhx_tree(tree_newick_fn)

    def process_node(self, node):
        if node.is_leaf():
            path = False
            if hasattr(node, "path"):
                fastas_fn = node.path.split("@")
                path = True
            elif hasattr(node, "fastapath"):
                fastas_fn = node.fastapath.split("@")
                path = True
            if path:
                for fasta_fn in fastas_fn:
                    print(fasta_fn)

        else:
            children = node.get_children()

            for child in children:
                self.process_node(child)


if __name__ == "__main__":
    assert (len(sys.argv) == 2)
    newick_fn = sys.argv[1]

    ti = TreeIndex(tree_newick_fn=newick_fn, )
    ti.process_node(ti.tree.get_tree_root())
