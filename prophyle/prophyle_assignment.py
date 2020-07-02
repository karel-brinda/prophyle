#! /usr/bin/env python3
"""ProPhyle assignment algorithm (reference implementation).

Example: ./prophyle/prophyle_index/prophyle_index query -k 5 -u -b _index_test/index.fa tests/simulation_bacteria.1000.fq |./prophyle/prophyle_assignment.py -m h1 -f sam _index_test/tree.nw 5 -

Author:  Karel Brinda <kbrinda@hsph.harvard.edu>

License: MIT

Todo:
    - test with CRAM
"""

import argparse
#import bitarray
import collections
import functools
import itertools
import os
import sys

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

###############################################################################################
###############################################################################################

CONFIG = {
    # this should be longer than any possible read, because of CRAM (non-tested yet)
    'FAKE_CONTIG_LENGTH': 42424242,
    # print diagnostics messages
    'DIAGNOSTICS': False,
    # sort nodes alphabetically when reporting assignments
    'SORT_NODES': False,
}

###############################################################################################
###############################################################################################

from bitarray import bitarray as _bitarray


class bitarray(_bitarray):
    def __hash__(self):
        return self.tobytes().__hash__()


class Assignment:
    """Class for handling a single read.

    Args:
        output_fo (file): Output file object.
        tree_index (TreeIndex): Tree index.
        kmer_lca (bool): Simulate LCA on the k-mer level.
        tie_lca (bool): If a tie (multiple winning nodes), compute LCA.
        annotate (bool): If taxonomic info present in the tree, annotate the assignments using SAM tags.
        mask_unmatched_bases (bool): Mask unmatched bases.

    Attributes:
        output_fo (file): Output file object.
        tree_index (TreeIndex): Tree index.
        k (int): k-mer size.
        kmer_lca (bool): Simulate LCA on the k-mer level.
        tie_lca (bool): If a tie (multiple winning nodes), compute LCA.
        mask_unmatched_bases (bool): Mask unmatched bases.
        annotate (bool): If taxonomic info present in the tree, annotate the assignments using SAM tags.
        krakline_parser (KraklineParser): Parser of kraklines.
        hitmasks_dict (dict): Hit masks
        covmasks_dict (dict): Cov masks.
        ass_dict (dict): Assignment dictionary.
        max_nodenames (list of str): List of nodenames of winners.
        max_val (int/float): Maximal value of the measure.
    """

    def __init__(
        self, output_fo, tree_index, kmer_lca=False, tie_lca=False, annotate=False, mask_unmatched_bases=False
    ):
        self.output_fo = output_fo
        self.tree_index = tree_index
        self.k = self.tree_index.k
        self.kmer_lca = kmer_lca
        self.annotate = annotate
        self.tie_lca = tie_lca
        self.mask_unmatched_bases = mask_unmatched_bases

        self.krakline_parser = KraklineParser(k=self.k)

        self.hitmasks_dict = {}
        self.covmasks_dict = {}
        self.ass_dict = {}

        self.max_nodenames = []
        self.max_val = 0

    def process_read(self, krakline, form, measure):
        """Process one Kraken-like line.

        Args:
            krakline (str): Kraken-like line.
            form (str): Expected output format (sam/kraken).
            measure (str): Measure (h1/c1/h2/c2).
        """

        self.krakline_parser.parse_krakline(krakline)
        if CONFIG['DIAGNOSTICS']:
            self.krakline_parser.diagnostics()

        self.blocks_to_masks(self.krakline_parser.kmer_blocks, self.kmer_lca)
        if CONFIG['DIAGNOSTICS']:
            self.diagnostics()

        self.compute_assignments()
        if CONFIG['DIAGNOSTICS']:
            self.diagnostics()

        self.select_best_assignments(measure, self.tie_lca)
        if CONFIG['DIAGNOSTICS']:
            self.diagnostics()

        self.print_selected_assignments(form)

    def blocks_to_masks(self, kmer_blocks, kmer_lca):
        """Extract hit and coverage masks from krakline blocks (without propagation) and store them in self.{hit,cov}masks_dict.

        Args:
            kmer_blocks (list): List of assigned k-mers, i.e., list of (list of node_names, count).
            kmer_lca (bool): Simulate k-mer LCA on the k-mer level.
        """

        readlen = sum([x[1] for x in kmer_blocks])

        hitmask_len = readlen
        covmask_len = readlen + self.k - 1

        # zero masks
        hitmask_empty = self.bitarray_block(hitmask_len, 0, 0)
        covmask_empty = self.bitarray_block(covmask_len, 0, 0)

        # fast copying (bitarray trick)
        hitmasks_dict = collections.defaultdict(lambda: bitarray(hitmask_empty))
        covmasks_dict = collections.defaultdict(lambda: bitarray(covmask_empty))

        pos = 0

        for (node_names, count) in kmer_blocks:

            # Kraken output format: 0 and A have special meanings, no blocks
            if node_names != ["0"] and node_names != ["A"]:

                hitmask_1block = self.bitarray_block(hitmask_len, count, pos)
                covmask_1block = self.bitarray_block(covmask_len, count + self.k - 1, pos)

                if kmer_lca:
                    node_names = [self.tree_index.lca(node_names)]

                for node_name in node_names:
                    hitmasks_dict[node_name] |= hitmask_1block
                    covmasks_dict[node_name] |= covmask_1block

            pos += count

        self.readlen = readlen
        self.hitmasks_dict = hitmasks_dict
        self.covmasks_dict = covmasks_dict

    def compute_assignments(self):
        """Compute assignments & characteristics.

        Compute and their characteristics from hitmasks and store
        them in self.ass_dict.
        """

        nodenames = self.hitmasks_dict.keys()
        self.ass_dict = {nodename: self.evaluate_single_assignment(nodename) for nodename in nodenames}

    def evaluate_single_assignment(self, nodename):
        """Evaluate a single assignment.

        Args:
            nodename (str): Name of the node for which we will compute characteristics.

        Returns:
            assignment (dict): Assignment dictionary.
        """

        #################################
        # A) Start with the current masks
        #################################
        hitmask = bitarray(self.hitmasks_dict[nodename])
        covmask = bitarray(self.covmasks_dict[nodename])

        ##########################
        # B) Update from ancestors
        ##########################
        ancestors = self.tree_index.nodename_to_upnodenames[nodename] & self.hitmasks_dict.keys()
        for anc_nodename in ancestors:
            hitmask |= self.hitmasks_dict[anc_nodename]
            covmask |= self.covmasks_dict[anc_nodename]

        ##############################
        # C) Calculate characteristics
        ##############################
        hit = hitmask.count()
        cov = covmask.count()
        readlen = self.krakline_parser.readlen
        kcount = self.tree_index.nodename_to_kmercount[nodename]

        assignment = {
            # 1. hit count
            'hitmask': hitmask,
            #'hitcigar': self.cigar_from_bitmask(hitmask),
            'h1': [hit],
            'hf': [hit / (readlen - self.k + 1)],
            'h2': [hit / kcount if kcount > 0 else 0],

            # 2. coverage
            'covmask': covmask,
            #'covcigar': self.cigar_from_bitmask(covmask),
            'c1': [cov],
            'cf': [cov / readlen],
            'c2': [cov / kcount if kcount > 0 else 0],

            # 3. other values
            'ln': readlen,
        }

        return assignment

    def select_best_assignments(self, measure, tie_lca=False):
        """Find the best assignments and save it to self.max_nodenames (max value: self.max_val).

        Args:
            measure (str): Measure (h1/c1/h2/c2).
            tie_lca (bool): Compute LCA in case of tie.
        """

        self.max_val = 0
        self.max_nodenames = []

        for nodename, ass in self.ass_dict.items():

            if ass[measure][0] > self.max_val:
                self.max_val = ass[measure][0]
                self.max_nodenames = [nodename]

            elif ass[measure][0] == self.max_val:
                self.max_nodenames.append(nodename)

        if tie_lca:
            self.make_lca_from_winners()

        if CONFIG['SORT_NODES']:
            self.max_nodenames.sort()

    def make_lca_from_winners(self):
        """Create LCA from winners.

        Assemble the characteristics of the LCA from the characteristics
        of the winning node. Output mutliple tag measures (e.g., h1 or c1).
        """

        # all characteristics are already computed
        if len(self.max_nodenames) <= 1:
            return

        lca = self.tree_index.lca(self.max_nodenames)

        ass = {
            'covmask': None,
            'covcigar': None,
            'hitmask': None,
            'hitcigar': None,
            'ln': self.krakline_parser.readlen,
        }

        for tag in ['h1', 'hf', 'h2', 'c1', 'cf', 'c2']:
            # ass[tag] = [self.ass_dict[nodename][tag] for nodename in self.max_nodenames]
            ass[tag] = [0]

        self.max_nodenames = [lca]
        self.ass_dict[lca] = ass

    def print_selected_assignments(self, form):
        """Print the best assignments.

        Args:
            form (str): Expected output format (sam/kraken).
        """

        if form == "sam":
            if len(self.max_nodenames):
                tag_is = len(self.max_nodenames)
                for tag_ii, nodename in enumerate(self.max_nodenames, 1):
                    ass = self.ass_dict[nodename]
                    # compute cigar
                    if ass['covmask'] is None:
                        ass['covcigar'] = None
                    else:
                        ass['covcigar'] = self.cigar_from_bitmask(ass['covmask'])

                    if ass['hitmask'] is None:
                        ass['hitcigar'] = None
                    else:
                        ass['hitcigar'] = self.cigar_from_bitmask(ass['hitmask'])

                    suffix_parts = ["ii:i:{}".format(tag_ii), "is:i:{}".format(tag_is)]
                    if self.annotate:
                        suffix_parts.append(self.tree_index.nodename_to_samannot[nodename])
                    self.print_sam_line(nodename, "\t".join([""] + suffix_parts))
            else:
                self.print_sam_line(None)
        elif form == "kraken":
            self.print_kraken_line(*self.max_nodenames)

    @staticmethod
    @functools.lru_cache(maxsize=5)
    def cigar_from_bitmask(bitmask):
        """Create a CIGAR from a binary mask.

        Args:
            mask (list): Bitmask.

        Return:
            cigar (str): Cigar string (X/= ops).
        """

        c = []
        bitmask = map(bool, bitmask)
        runs = itertools.groupby(bitmask)
        for run in runs:
            c.append(str(len(list(run[1]))))
            c.append('=' if run[0] else 'X')
        return "".join(c)

    def print_sam_line(self, node_name, suffix=""):
        """Print a single SAM record.

        Args:
            node_name (str): Node name. None if unassigned.
            suffix (str): Suffix with additional tags.
        """

        qname = self.krakline_parser.readname

        if self.mask_unmatched_bases and self.krakline_parser.seq and node_name is not None:
            covmask = self.ass_dict[node_name]["covmask"]
            seq = "".join(("N" if covmask[i] == 0 else ch for i, ch in enumerate(self.krakline_parser.seq)))
        else:
            seq = self.krakline_parser.seq

        if node_name is not None:
            flag = 0
            mapq = "60"
            pos = "1"
            cigar = self.ass_dict[node_name]['covcigar']
        else:
            flag = 4
            node_name = None
            pos = None
            mapq = None
            cigar = None

        columns = [
            qname,  # QNAME
            str(flag),  # FLAG
            node_name if node_name else "*",  # RNAME
            pos if pos else "0",  # POS
            mapq if mapq else "0",  # MAPQ
            cigar if cigar else "*",  # CIGAR
            "*",  # RNEXT
            "0",  # PNEXT
            "0",  # TLEN
            seq if self.krakline_parser.seq else "*",  # SEQ
            self.krakline_parser.qual if self.krakline_parser.qual else "*",  # QUAL
        ]

        columns.append("{}:{}:{}".format("ln", "i", self.readlen))

        if node_name is not None:
            asg = self.ass_dict[node_name]

            for tag, datatype in [
                ('h1', 'i'),
                ('h2', 'f'),
                ('hf', 'f'),
                ('c1', 'i'),
                ('c2', 'f'),
                ('cf', 'f'),
            ]:
                for val in asg[tag]:
                    columns.append("{}:{}:{}".format(tag, datatype, val))

            if asg['hitcigar']:
                columns.append("hc:Z:{}".format(asg['hitcigar']))

        print("\t".join(columns), suffix, file=self.output_fo, sep="")

    def print_sam_header(self):
        """Print SAM headers.
        """
        print("@HD", "VN:1.5", "SO:unsorted", sep="\t", file=self.output_fo)
        print("@PG", "PN:prophyle", "ID:prophyle", "VN:{}".format(version.VERSION), sep="\t", file=self.output_fo)
        for node in self.tree_index.tree.traverse("postorder"):

            try:
                ur = "\tUR:{}".format(node.path)
            except:
                ur = ""

            try:
                sp = "\tSP:{}".format(node.sci_name)
            except:
                sp = ""

            try:
                as_ = "\tAS:{}".format(node.gi)
            except:
                as_ = ""

            if node.name != '':
                print(
                    "@SQ\tSN:{rname}\tLN:{rlen}{as_}{ur}{sp}".format(
                        rname=node.name,
                        rlen=CONFIG['FAKE_CONTIG_LENGTH'],
                        as_=as_,
                        ur=ur,
                        sp=sp,
                    ), file=self.output_fo
                )

    def print_kraken_line(self, *nodenames):
        """Print a single record in the Kraken-like format.

        Args:
            *nodenames (list of str): Node names of assignments to report.
        """

        if len(nodenames) == 0:
            stat = "U"
            krak_ass = "0"
        else:
            stat = "C"
            krak_ass = ",".join(nodenames)

        # recompute krakenmers
        if self.kmer_lca:
            nodenames_lca_seq = []
            for [nodenames, count] in self.krakline_parser.kmer_blocks:
                if len(nodenames) == 1:
                    nodename = nodenames[0]
                    #if nodename == "A" or nodename == "0":
                    #    pass
                    #else:
                    #    pass
                else:
                    nodename = self.tree_index.lca(nodenames)
                nodenames_lca_seq.extend(count * [nodename])
            c = []
            runs = itertools.groupby(nodenames_lca_seq)
            for run in runs:
                c.append("{}:{}".format(str(run[0]), len(list(run[1]))))
            krakmers = " ".join(c)
        else:
            krakmers = self.krakline_parser.krakmers

        columns = [stat, self.krakline_parser.readname, krak_ass, str(self.krakline_parser.readlen), krakmers]
        print("\t".join(columns), file=self.output_fo)

    @staticmethod
    def bitarray_block(alen, blen, pos):
        """Create a bitarray containing a block of one's.

        Args:
            alen (int): Array length.
            blen (int): Block length (within the array).
            pos (int): Position of the block (0-based).

        Return:
            bitarray (bitarray)
        """
        return bitarray(pos * "0" + blen * "1" + (alen - pos - blen) * "0")

    def diagnostics(self):
        """Print debug messages.
        """
        print("---------------------", file=sys.stderr)
        print("Alignment diagnostics", file=sys.stderr)
        print("---------------------", file=sys.stderr)
        print("Alignment.hitmasks_dict: ", dict(self.hitmasks_dict), file=sys.stderr)
        print("Alignment.covmasks_dict: ", dict(self.covmasks_dict), file=sys.stderr)
        print("Alignment.ass_dict:      ", dict(self.ass_dict), file=sys.stderr)
        print("Alignment.max_nodenames: ", self.max_nodenames, file=sys.stderr)
        print("Alignment.max_val:       ", self.max_val, file=sys.stderr)


###############################################################################################
###############################################################################################


class TreeIndex:
    """Class for an indexed Phylogenetic tree.

    Args:
        tree_newick_fn (str): Filename of the phylogenetic tree.
        k (int): K-mer size.

    Attributes:
        tree (ete3.Tree): Minimal subtree of the original phylogenetic tree.
        k (int): K-mer size.
        nodename_to_node (dict): node name => node.
        nodename_to_samannot (dict): node name => string to append in SAM.
        nodename_to_upnodenames (dict): node name => set of node names of ancestors.
        nodename_to_kmercount (dict): nname => number of k-mers (full set).
    """

    def __init__(self, tree_newick_fn, k):
        tree = pro.load_nhx_tree(tree_newick_fn)
        self.tree = pro.minimal_subtree(tree)

        self.k = k

        self.nodename_to_node = {}
        self.nodename_to_kmercount = {}

        self.nodename_to_samannot = {}

        self.nodename_to_upnodenames = collections.defaultdict(lambda: set())

        for node in self.tree.traverse("postorder"):
            nodename = node.name
            self.nodename_to_node[nodename] = node
            self.nodename_to_kmercount[nodename] = int(node.kmers_full)

            # annotations
            tags_parts = []
            try:
                tags_parts.append("gi:Z:{}".format(node.gi))
            except AttributeError:
                pass

            try:
                tags_parts.append("sn:Z:{}".format(node.sci_name))
            except AttributeError:
                pass

            try:
                tags_parts.append("ra:Z:{}".format(node.rank))
            except AttributeError:
                pass

            self.nodename_to_samannot[nodename] = "\t".join(tags_parts)

            # set of upper nodes
            while node.up:
                node = node.up
                self.nodename_to_upnodenames[nodename].add(node.name)

    def lca(self, node_names):
        """Return LCA for a given list of nodes.

        node_names (list of str): List of node names.

        Returns:
            str: Name of the LCA.
        """
        assert len(node_names) > 0
        if len(node_names) == 1:
            return node_names[0]
        nodes = [self.nodename_to_node[n] for n in node_names]
        lca = nodes[0].get_common_ancestor(nodes)

        if lca.is_root() and len(lca.children) == 1:
            lca = lca.children[0]
        assert lca.name != ""  # , [x.name for x in lca.children]

        return lca.name

    def diagnostics(self):
        """Print debug messages.
        """
        print("---------------------", file=sys.stderr)
        print("TreeIndex diagnostics", file=sys.stderr)
        print("---------------------", file=sys.stderr)
        print("TreeIndex.k:                       ", self.k, file=sys.stderr)
        print("TreeIndex.nodename_to_node:        ", self.nodename_to_node, file=sys.stderr)
        print("TreeIndex.nodename_to_samannot:    ", self.nodename_to_samannot, file=sys.stderr)
        print("TreeIndex.nodename_to_upnodenames: ", self.nodename_to_upnodenames, file=sys.stderr)
        print("TreeIndex.nodename_to_kmercount:   ", self.nodename_to_kmercount, file=sys.stderr)


###############################################################################################
###############################################################################################


class KraklineParser():
    """Class for parsing Kraken-like input into a structure.

    Args:
        k (int): K-mer size.

    Attributes:
        krakline (str): Original krakline.
        readname (str): Name of the read.
        raedlen (str): Length of the read.
        seq (str): Sequence of nucleotides. None if unknown.
        qual (str): Sequence of qualities. None if unknown.
        kmer_blocks (list of (list of str, int)): Assigned k-mer blocks, list of (nodenames, count).
    """

    def __init__(self, k):
        self.k = k
        self.krakline = None
        self.readname = None
        self.readlen = None
        self.seq = None
        self.qual = None
        self.kmer_blocks = []

    def parse_krakline(self, krakline):
        """Load a krakline to the current object.

        Args:
            krakline (str): Kraken-like line.
        """

        self.krakline = krakline
        parts = krakline.strip().split("\t")
        self.readname, _, readlen, self.krakmers = parts[1:5]
        self.readlen = int(readlen)
        if len(parts) == 7:
            self.seq = parts[5]
            self.qual = parts[6]
        else:
            self.seq = None
            self.qual = None

        # list of (count,list of nodes)
        self.kmer_blocks = []
        kmer_countdown = self.readlen - self.k + 1

        for block in self.krakmers.split(" "):
            try:
                (ids, count) = block.split(":")
                count = int(count)
                kmer_countdown -= count
                nodenames = ids.split(",")
                self.kmer_blocks.append((nodenames, count))
            except ValueError:
                pro.message("Warning: prophex output for read '{}' has been truncated.".format(self.readname))
        self.kmer_blocks.append((['0'], kmer_countdown))
        #pro.message(str(self.kmer_blocks))

    def check_consistency(self, k):
        """Check consistency of the fields loaded from the krakline.

        Args:
            k (int): k-mer length.

        Returns:
            bool: Consistent.
        """
        if self.qlen < k:
            return True

        block_len_sum = sum([x[1] for x in self.kmer_blocks])

        if not self.readlen == block_len_sum + k - 1:
            return False

        if self.seq is not None and len(self.seq) != self.readlen():
            return False

        if self.qual is not None and len(self.qual) != self.readlen():
            return False

        return True

    def diagnostics(self):
        """Print debug messages.
        """
        print("--------------------------", file=sys.stderr)
        print("KraklineParser diagnostics", file=sys.stderr)
        print("--------------------------", file=sys.stderr)
        #print("KraklineParser.krakline:    ", self.krakline.strip(), file=sys.stderr)
        print("KraklineParser.readname:    ", self.readname, file=sys.stderr)
        print("KraklineParser.krakmers:    ", self.krakmers, file=sys.stderr)
        print("KraklineParser.readlen:     ", self.readlen, file=sys.stderr)
        print("KraklineParser.seq:         ", self.seq, file=sys.stderr)
        print("KraklineParser.qual:        ", self.qual, file=sys.stderr)
        print("KraklineParser.kmer_blocks: ", self.kmer_blocks, file=sys.stderr)


###############################################################################################
###############################################################################################


def assign_all_reads(
    tree_fn,
    inp_fo,
    kmer_lca,
    tie_lca,
    form,
    k,
    measure,
    annotate,
):
    assert form in ["kraken", "sam"]
    assert k > 1

    tree_index = TreeIndex(
        tree_newick_fn=tree_fn,
        k=k,
    )

    if CONFIG['DIAGNOSTICS']:
        tree_index.diagnostics()

    assignment = Assignment(
        output_fo=sys.stdout,
        tree_index=tree_index,
        kmer_lca=kmer_lca,
        tie_lca=tie_lca,
        annotate=annotate,
    )

    if form == "sam":
        assignment.print_sam_header()

    for krakline in inp_fo:
        assignment.process_read(krakline, form=form, measure=measure)


def parse_args():
    parser = argparse.ArgumentParser(description='Implementation of assignment algorithm')

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nhx>',
        help='phylogenetic tree (Newick/NHX)',
    )

    parser.add_argument(
        'k',
        type=int,
        metavar='<k>',
        help='k-mer length',
    )

    parser.add_argument(
        'input_file',
        type=argparse.FileType('r'),
        metavar='<assignments.txt>',
        help='assignments in generalized Kraken format',
    )

    parser.add_argument(
        '-f',
        choices=['kraken', 'sam'],
        default='sam',
        dest='format',
        help='format of output [sam]',
    )

    parser.add_argument(
        '-m',
        choices=['h1', 'c1', 'c2', 'h2'],
        default='h1',
        dest='measure',
        help='measure: h1=hit count, c1=coverage, h2=norm.hit count, c2=norm.coverage [h1]',
    )

    parser.add_argument(
        '-A',
        action='store_true',
        dest='annotate',
        help='annotate assignments',
    )

    parser.add_argument(
        '-L',
        action='store_true',
        dest='tie_lca',
        help='use LCA when tie (multiple assignments with the same score)',
    )

    parser.add_argument(
        '-X',
        action='store_true',
        dest='kmer_lca',
        help='use LCA for k-mers (multiple hits of a k-mer)',
    )

    parser.add_argument(
        '-c',
        dest='config',
        metavar='STR',
        nargs='*',
        type=str,
        default=[],
        help='configuration (a JSON dictionary)',
    )

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    global CONFIG
    prophyle_conf_string = pro.load_prophyle_conf(CONFIG, args.config)

    try:
        assign_all_reads(
            tree_fn=args.tree_fn,
            inp_fo=args.input_file,
            form=args.format,
            k=args.k,
            measure=args.measure,
            annotate=args.annotate,
            tie_lca=args.tie_lca,
            kmer_lca=args.kmer_lca,
        )

    # Karel: I don't remember why I was considering also IOError here
    # except (BrokenPipeError, IOError):
    except BrokenPipeError:
        # pipe error (e.g., when head is used)
        sys.stderr.close()
        sys.stdout.close()
        exit(0)

    except KeyboardInterrupt:
        pro.message("Error: Keyboard interrupt")
        pro.close_log()
        exit(1)

    finally:
        try:
            sys.stdout.flush()
        except BrokenPipeError:
            pass
        finally:
            try:
                sys.stderr.flush()
            except:
                pass


if __name__ == "__main__":
    main()
