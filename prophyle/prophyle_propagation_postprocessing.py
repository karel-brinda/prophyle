#! /usr/bin/env python3
"""K-mer propagation post-processing.

Create the main index FASTA file and a tree with k-mer annotations (number of
k-mers during propagation - size of the  full k-mer set and reduced k-mer set).

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT
"""

import argparse
import glob
import os
import re
import sys

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


def load_kmer_stats(tsv_fn):
    re_fa = re.compile(r'(.*)\.([a-z]*)\.fa')

    counts = {"reduced": {}, "full": {}}
    with open(tsv_fn) as f:
        for x in f:
            x = x.strip()
            if len(x) == 0:
                continue
            if x[0] == "#":
                continue
            else:
                (fa, count) = x.split("\t")
                count = int(count)
                fa_short = os.path.basename(fa)
                m = re_fa.match(fa_short)

                cat = m.group(2)
                nname = m.group(1)

                try:
                    assert counts[cat][
                        nname
                    ] == count, "Different k-mer sizes reported for the same node '{}', probably a bug of prophyle_assembler".format(
                        nname
                    )
                except KeyError:
                    counts[cat][nname] = count

    return counts


def enrich_tree(tree, count_tb):
    for node in tree.traverse("postorder"):
        node.del_feature('kmers_full')
        node.del_feature('kmers_reduced')

        nname = node.name

        assert nname != "", "There is a node without any name ('')"

        singleton = len(node.children) == 1

        try:
            node.add_features(kmers_full=count_tb["full"][nname])
        except KeyError:
            if singleton:
                # kmer_full the same as in the child
                node.add_features(kmers_full=node.children[0].kmers_full)
            else:
                print("Warning: full-{} is missing".format(nname), file=sys.stderr)
                node.add_features(kmers_full=0)

        try:
            node.add_features(kmers_reduced=count_tb["reduced"][nname])
        except KeyError:
            if not singleton:
                print("Warning: reduced-{} is missing".format(nname), file=sys.stderr)
            node.add_features(kmers_reduced=0)

        assert 0 <= node.kmers_reduced
        assert node.kmers_reduced <= node.kmers_full


def create_fasta(propagation_dir, index_fasta_fn, suffix, verbose=False):
    re_name = re.compile(r'(.*)\.' + suffix.replace(r'.', '\.') + r'$')

    fa_fns = glob.glob("{}/*".format(propagation_dir))
    fa_fns.sort()

    with open(index_fasta_fn, "w+") as fa_fo:
        for fn in fa_fns:
            if verbose:
                print("Processing '{}'".format(fn), file=sys.stderr)

            fn2 = fn.split("/")[-1]

            m = re_name.match(fn2)

            if not m:
                continue

            node_name = m.group(1)

            with open(fn) as f:
                for x in f:
                    if len(x) == 0:
                        continue
                    if x[0] == ">":
                        contig_name = x[1:]
                        print(">{}@{}".format(node_name, contig_name), end="", file=fa_fo)
                    else:
                        print(x, end="", file=fa_fo)


def main():

    parser = argparse.ArgumentParser(
        description='K-mer propagation postprocessing: merging FASTA files and k-mer annotation.'
    )

    parser.add_argument('dir', metavar='<propagation.dir>', type=str, help='directory with FASTA files')

    parser.add_argument(
        'index_fasta_fn',
        type=str,
        metavar='<index.fa>',
        help='output fast file',
    )

    parser.add_argument(
        'in_tree_fn',
        type=str,
        metavar='<in.tree.nw>',
        help='input phylogenetic tree',
    )

    parser.add_argument(
        'counts_fn',
        type=str,
        metavar='<counts.tsv>',
        help='input phylogenetic tree',
    )

    parser.add_argument(
        'out_tree_fn',
        type=str,
        metavar='<out.tree.nw>',
        help='output phylogenetic tree',
    )

    #parser.add_argument (
    #   '-D',
    #   dest='nondel',
    #   action='store_true',
    #   help='Non-deleting propagation',
    #)

    args = parser.parse_args()

    dir_fn = args.dir
    index_fasta_fn = args.index_fasta_fn
    in_tree_fn = args.in_tree_fn
    out_tree_fn = args.out_tree_fn
    tsv_fn = args.counts_fn

    suffix = "reduced.fa"

    #if args.nondel:
    #   suffix = "full.fa"
    #else:
    #   suffix = "reduced.fa"

    create_fasta(dir_fn, index_fasta_fn, suffix)

    tree = pro.load_nhx_tree(in_tree_fn)
    stats = load_kmer_stats(tsv_fn)
    enrich_tree(tree, stats)
    pro.save_nhx_tree(tree, out_tree_fn)


if __name__ == "__main__":
    main()
