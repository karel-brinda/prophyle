#! /usr/bin/env python3
"""K-mer propagation pre-processing.

Extract subtrees, merge ProPhyle Newick/NHX trees, possibly run autocomplete
(add path and create names for internal nodes) and subsample the resulting tree.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Examples:
    $ prophyle_propagation_preprocessing.py ~/prophyle/bacteria.nw ~/prophyle/viruses.nw bv.nw
    $ prophyle_propagation_preprocessing.py ~/prophyle/bacteria.nw@562 ecoli.nw
"""

import argparse
import ete3
import os
import random
import re
import sys

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro

random.seed(42)


def add_prefix(tree, prefix):
    for node in tree.traverse("postorder"):
        node.name = "{}-{}".format(prefix, node.name)
    return tree


def autocomplete_fasta_path(tree):
    print("Autocompleting FASTA paths", file=sys.stderr)
    for n in tree.traverse():
        if len(n.children) == 0:
            n.add_features(path="{}.fa".format(n.name))
    return tree


def autocomplete_internal_node_names(tree):
    print("Autocompleting internal node names", file=sys.stderr)

    re_inferred = re.compile(r'^(.*)-up(\d+)$')

    for n in tree.traverse("postorder"):
        if len(n.children) == 0:
            assert hasattr(n, "name")
        else:
            for x in n.children:
                assert hasattr(x, "name")

            if not hasattr(n, "name") or n.name == "" or n.name is None:
                names = [x.name for x in n.children]
                lmin_name = sorted(names)[0]

                m = re_inferred.match(lmin_name)
                if m is not None:
                    left, right = m.groups()
                    right = int(right) + 1
                    n.name = "{}-up{}".format(left, right)
                else:
                    n.name = lmin_name + "-up1"

    return tree


def merge_trees(input_trees_fn, output_tree_fn, verbose, add_prefixes, sampling_rate, autocomplete):
    assert sampling_rate is None or 0.0 <= float(sampling_rate) <= 1.0

    t = ete3.Tree(name="merge_root", )

    if len(input_trees_fn) == 1:
        if verbose:
            print("Only one tree, don't add any prefix", file=sys.stderr)
        add_prefixes = False

    for i, x in enumerate(input_trees_fn, 1):
        if verbose:
            print("Loading '{}'".format(x), file=sys.stderr)

        tree_fn, _, root_name = x.partition("@")
        tree_to_add = pro.load_nhx_tree(tree_fn, validate=False)

        # subtree extraction required
        if root_name != '':
            tree_to_add = tree_to_add & root_name

        # prepend prefixes to node names
        if add_prefixes:
            tree_to_add = add_prefix(tree_to_add, i)

        t.add_child(tree_to_add)

    if autocomplete:
        if not (pro.has_attribute(t, "path") or pro.has_attribute(t, "fastapath")):
            t = autocomplete_fasta_path(t)
        t = autocomplete_internal_node_names(t)

    if sampling_rate is not None:
        sampling_rate = float(sampling_rate)

        leaves_1 = []

        for node in t.traverse("postorder"):
            if len(node.children) == 0:
                leaves_1.append(node)

        leaves_1.sort(key=lambda x: x.name)

        leaves_2 = random.sample(leaves_1, max(round(sampling_rate * len(leaves_1)), 1))
        leaves_2.sort(key=lambda x: x.name)

        leaves_to_remove = list(set(leaves_1) - set(leaves_2))
        leaves_to_remove.sort(key=lambda x: x.name)

        if verbose:
            print(
                "Removing the following leaves: {}".format(", ".join(map(lambda x: x.name, leaves_to_remove))),
                file=sys.stderr
            )

        for node in leaves_to_remove:
            while len(node.up.children) == 1:
                node = node.up
            node.detach()

        print(
            "Subsampling the tree with rate {:.4f}, {} leaves were kept (out of {})".format(
                sampling_rate, len(leaves_2), len(leaves_1)
            ), file=sys.stderr
        )

    for node in t.traverse("postorder"):
        if hasattr(node, "fastapath"):
            node.add_feature("path", node.fastapath)
            node.del_feature("fastapath")

    if verbose:
        print("Writing to '{}'".format(output_tree_fn), file=sys.stderr)

    pro.save_nhx_tree(t, output_tree_fn)


def parse_args():
    parser = argparse.ArgumentParser(
        description="\n".join(
            [
                'Merge multiple ProPhyle trees. Specific subtrees might be extracted before merging. Examples:',
                '\t$ prophyle_merge_trees.py ~/prophyle/bacteria.nw ~/prophyle/viruses.nw bv.nw',
                '\t$ prophyle_merge_trees.py ~/prophyle/bacteria.nw@562 ecoli.nw'
            ]
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        'in_tree',
        metavar='<in_tree.nw{@node_name}>',
        type=str,
        help='input tree',
        nargs='+',
    )

    parser.add_argument(
        'out_tree',
        metavar='<out_tree.nw>',
        type=str,
        help='output tree',
    )

    parser.add_argument(
        '-s',
        help='rate of sampling the tree [no sampling]',
        dest='sampling_rate',
        metavar='FLOAT',
        type=str,
        default=None,
    )

    parser.add_argument(
        '-A',
        help='autocomplete tree (names of internal nodes and FASTA paths)',
        dest='autocomplete',
        action='store_true',
    )

    parser.add_argument(
        '-V',
        help='verbose',
        dest='verbose',
        action='store_true',
    )

    parser.add_argument(
        '-P',
        help='do not add prefixes to node names',
        dest='add_prefixes',
        action='store_false',
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    merge_trees(
        input_trees_fn=args.in_tree,
        output_tree_fn=args.out_tree,
        verbose=args.verbose,
        add_prefixes=args.add_prefixes,
        sampling_rate=args.sampling_rate,
        autocomplete=args.autocomplete,
    )


if __name__ == "__main__":
    main()
