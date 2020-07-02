#! /usr/bin/env python3
"""Build taxonomic tree for genomic libraries using etetoolkit.

Author: Simone Pignotti <pignottisimone@gmail.com>

Licence: MIT

TODO(Simone):
    * use NCBI eutils to translate accesion numbers to TaxID
    * find solution for internal nodes with fasta file associated
TODO(Karel):
    * add a unit test
    * unify cli interface with the main prophyle program (parameters, etc.)
"""

# ! /usr/bin/env python3

import sys
import os
import argparse

from ete3 import NCBITaxa

desc = """\
    Program: prophyle_ncbi_tree

    Build a taxonomic tree in the New Hampshire newick format #1 for NCBI sequences
"""


def acquire_sequences(library, library_dir, log):
    seqs = {}
    acquired = 0
    skipped = 0
    for dirpath, _, filenames in os.walk(os.path.join(library_dir, library)):
        for filename in (f for f in filenames if f.endswith('.fai')):
            fn = os.path.join(dirpath, filename)
            rel_fn = fn[(len(library_dir) + 1):-4]
            with open(fn, 'r') as faidx:
                for seq in faidx:
                    try:
                        seqname, seqlen, offset, _, _ = seq.split('\t')
                        # RefSeq filenames start with accession numbers
                        f = (filename.split('.')[0]).split('_')
                        acc = f[0].strip() + '_' + f[1].strip()
                        try:
                            seqs[acc]['fn'] += '@' + rel_fn
                            seqs[acc]['seqname'] += '@' + seqname
                            seqs[acc]['seqlen'] += '@' + seqlen
                            seqs[acc]['offset'] += '@' + offset
                        except KeyError:
                            seqs[acc] = {'fn': rel_fn, 'seqname': seqname, 'seqlen': seqlen, 'offset': offset}
                            acquired += 1
                            pass
                    except:
                        if log:
                            print(
                                '[prophyle_ncbi_tree] ERROR: Sequence \"' + seqname + '\" in file \"' + rel_fn +
                                '\" not acquired', file=log
                            )
                        skipped += 1
                        pass
    return seqs, acquired, skipped


def assign_sequences(taxid_map_f, seqs, log):
    taxa2acc = {}
    acc2taxa = {}
    assigned = 0
    skipped = 0

    with open(taxid_map_f, 'r') as taxid_map:
        for line in taxid_map:
            acc, taxid = line.split('\t')
            taxid = int(taxid)
            try:
                taxa2acc[taxid].append(acc)
            except KeyError:
                taxa2acc[taxid] = [acc]
                pass
            acc2taxa[acc] = taxid

    for acc, dic in seqs.items():
        try:
            taxid = acc2taxa[acc]
            dic['taxid'] = taxid
            assigned += 1
        except KeyError:
            if log:
                print('[prophyle_ncbi_tree] ERROR: TaxID not found for sequence ' + acc, file=log)
            skipped += 1
            pass

    return taxa2acc, assigned, skipped


def build_tree(seqs, taxa2acc, red_factor, root, log):
    # Important: you should update ETE DB before running this script.
    # This is done automatically only if it has not been downloaded yet.
    ncbi = NCBITaxa()
    taxa = []
    for s in seqs.values():
        try:
            taxa.append(s['taxid'])
        except KeyError:
            continue
    built = False
    while not built:
        try:
            t = ncbi.get_topology(taxa, intermediate_nodes=True)
            built = True
        except KeyError as e:
            taxid_not_found = int(e.args[0])
            taxa.remove(taxid_not_found)
            if log:
                print(
                    '[prophyle_ncbi_tree] ERROR: TaxID ' + str(taxid_not_found) +
                    ' not found in ETE DB (try updating it)', file=log
                )
            pass

    # [Issue #153] Ignore internal nodes with fasta associated till we find a solution for it
    if log:
        internal_with_fasta = 0
        for node in t.traverse('postorder'):
            if not node.is_leaf() and node.taxid in taxa:
                internal_with_fasta += len([acc for acc in taxa2acc[node.taxid] if acc in seqs.keys()])
        print(
            '[prophyle_ncbi_tree] ' + str(internal_with_fasta) + ' sequences' +
            ' ignored because associated to internal node (see issue #153)', file=log
        )
    leaves_taxa = [leaf.taxid for leaf in t if leaf.taxid in taxa2acc]
    t = ncbi.get_topology(leaves_taxa, intermediate_nodes=True)

    if red_factor:
        i = 0
        red_taxa = []
        for leaf in t:
            if i % red_factor == 0:
                red_taxa.append(leaf.taxid)
            i += 1
        t = ncbi.get_topology(red_taxa, intermediate_nodes=True)

    if root:
        taxa_to_keep = []
        for leaf in t:
            if root in leaf.named_lineage:
                taxa_to_keep.append(leaf.taxid)
        t = ncbi.get_topology(taxa_to_keep, intermediate_nodes=True)

    node_count = len(t.get_descendants()) + 1
    seq_count = 0

    for node in t.traverse('postorder'):
        node.name = node.taxid
        if node.is_leaf():
            first = True
            for acc in taxa2acc[node.taxid]:
                try:
                    s = seqs[acc]
                    if first:
                        accession = '@'.join([acc] * (s['offset'].count('@') + 1))
                        path = s['fn']
                        base_len = s['seqlen']
                        infasta_offset = s['offset']
                        first = False
                    else:
                        accession += ('@' + acc) * (s['offset'].count('@') + 1)
                        path += '@' + s['fn']
                        base_len += '@' + s['seqlen']
                        infasta_offset += '@' + s['offset']
                    seq_count += 1
                except KeyError:
                    pass
            node.add_features(path=path, base_len=base_len, infasta_offset=infasta_offset, accession=accession)

    if not hasattr(t, 'taxid'):
        t.add_features(taxid=0)
    t.name = t.taxid

    return t, seq_count, node_count


def main_fun(library, library_dir, output_f, taxid_map, red_factor, root, log_file):
    library_dir = os.path.abspath(library_dir)

    seqs, acquired, skipped = acquire_sequences(library, library_dir, log_file)

    print(
        '[prophyle_ncbi_tree] Acquired ' + str(acquired) + ' sequences (' + str(skipped) + ' skipped)', file=sys.stderr
    )
    if log_file:
        print(
            '[prophyle_ncbi_tree] Acquired ' + str(acquired) + ' sequences (' + str(skipped) + ' skipped)',
            file=log_file
        )

    taxa2seqid, assigned, skipped = assign_sequences(taxid_map, seqs, log_file)

    print(
        '[prophyle_ncbi_tree] TaxID assigned to ' + str(assigned) + ' sequences (' + str(skipped) + ' skipped)',
        file=sys.stderr
    )
    if log_file:
        print(
            '[prophyle_ncbi_tree] TaxID assigned to ' + str(assigned) + ' sequences (' + str(skipped) + ' skipped)',
            file=log_file
        )

    tax_tree, seq_count, node_count = build_tree(seqs, taxa2seqid, red_factor, root, log_file)

    print(
        '[prophyle_ncbi_tree] Built taxonomic tree for ' + str(seq_count) + ' sequences (' + str(node_count) +
        ' nodes)', file=sys.stderr
    )
    if log_file:
        print(
            '[prophyle_ncbi_tree] Built taxonomic tree for ' + str(seq_count) + ' sequences (' + str(node_count) +
            ' nodes)', file=log_file
        )

    tax_tree.write(features=['accession', 'path', 'rank'], format=1, format_root_node=True, outfile=output_f)


def main():
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument(
        'library',
        type=str,
        metavar='<library>',
        help='directory with the library sequences (e.g. bacteria, viruses etc.)',
    )
    parser.add_argument(
        'library_dir', type=str, metavar='<library_dir>',
        help='library path (parent of library, e.g. main ProPhyle directory)'
    )
    parser.add_argument('output_file', type=str, metavar='<output_file>', help='output file')
    parser.add_argument(
        'taxid_map_file', type=str, metavar='<taxid_map>', help='tab separated accession number to taxid map'
    )
    parser.add_argument(
        '-l', type=argparse.FileType('a+'), metavar='log_file', dest='log_file', help='log file [stderr]'
    )
    parser.add_argument(
        '-r', type=int, metavar='red_factor', dest='red_factor', help='build reduced tree (one sequence every n)'
    )
    parser.add_argument(
        '-u', type=str, metavar='root', dest='root',
        help='root of the tree (e.g. Bacteria); will exclude sequences which are not its descendants'
    )

    args = parser.parse_args()

    main_fun(
        args.library, args.library_dir, args.output_file, args.taxid_map_file, args.red_factor, args.root, args.log_file
    )


if __name__ == '__main__':
    main()
