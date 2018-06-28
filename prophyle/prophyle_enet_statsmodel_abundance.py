#! /usr/bin/env python3
"""Estimate abundances from ProPhyle's assignment using a regularized Poisson GLM

Author:  Simone Pignotti <pignottisimone@gmail.com>

License: MIT

"""

###############################################################################################
###############################################################################################

CONFIG = {
    # print diagnostics messages
    'DIAGNOSTICS': False,
}

###############################################################################################
###############################################################################################

import os
import sys
import pysam
import argparse
import numpy as np
import statsmodels.api as sm
from ete3 import Tree

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version


def analyse_assignments(ass_fn, nodes2leaves, vec_pos):

    assignments = np.zeros(len(vec_pos))

    ass_f, in_format = pro.open_asg(ass_fn)
    prev_read_name = ""
    cur_ref = []

    if in_format == 'kraken':
        for read in ass_f:
            fields = read.split('\t')
            if fields[0] == 'C':
                read_name = fields[2]
                read_ref = fields[3]
                if read_name == prev_read_name:
                    cur_ref.append(read_ref)
                else:
                    for ref in cur_ref:
                        try:
                            for leaf in nodes2leaves[ref]:
                                assignments[vec_pos[leaf]] += 1
                        except KeyError:
                            print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)
                    cur_ref = [read_ref]
                    prev_read_name = read_name
    else:
        for read in ass_f.fetch(until_eof=True):
            if not read.is_unmapped:
                read_name = read.qname
                read_ref = read.reference_name
                if read_name == prev_read_name:
                    cur_ref.append(read_ref)
                else:
                    for ref in cur_ref:
                        try:
                            for leaf in nodes2leaves[ref]:
                                assignments[vec_pos[leaf]] += 1
                        except KeyError:
                            print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)
                    cur_ref = [read_ref]
                    prev_read_name = read_name

    # last assignment
    if len(cur_ref) > 0:
        try:
            for ref in cur_ref:
                for leaf in nodes2leaves[ref]:
                    assignments[vec_pos[leaf]] += 1
        except KeyError:
            print('[prophyle_sim_matrix] Warning: assignments to {} ignored because not in the tree. Are you using the right tree/index?'.format(ref), file=sys.stderr)

    ass_f.close()

    return assignments


def estimate_abundances(tree_fn, asg_fn, sim_mat_fn, out_fn, alpha=0.1, l1_ratio=0.99):

    tree = Tree(tree_fn, format=1)
    leaves = [leaf.name for leaf in tree]
    nodes2leaves = {node.name: {leaf.name for leaf in node} for node in tree.traverse("postorder")}
    vec_pos = {leaf.name: i for i, leaf in enumerate(tree)}

    count_fn = '.'.join(asg_fn.split('.')[:-1]+['npy'])
    if os.path.isfile(count_fn):
        print("Loading counts from existing npy file", file=sys.stderr)
        map_counts = np.load(count_fn)
    else:
        map_counts = analyse_assignments(asg_fn, nodes2leaves, vec_pos)
        np.save(count_fn, map_counts)

    assert len(leaves) == len(map_counts), "Length of mappings different from #leaves...try to remove <asg.npy> and analyse assignments again using the right tree!"

    sim_mat = np.load(sim_mat_fn)
    assert len(leaves) == len(sim_mat), "Size of similarity matrix different from #leaves...have you used the right index/tree?"

    nonzero_pos = map_counts > 100
    map_counts = map_counts[nonzero_pos]
    print(len(map_counts))

    temp_sim_mat = np.empty((len(map_counts), len(map_counts)))
    j = 0
    for i, is_nonzero in enumerate(nonzero_pos):
        if is_nonzero:
            temp_sim_mat[j,] = sim_mat[i][nonzero_pos]
            j += 1
    sim_mat = temp_sim_mat

    sim_mat[sim_mat==0] = 1e-10
    print(sim_mat)
    print(len(sim_mat))
    print(len(sim_mat[0]))

    # glm = sm.GLM(map_counts, sim_mat, family=sm.families.Poisson(link=sm.families.links.identity))
    glm = sm.GLM(map_counts, sim_mat, family=sm.families.Gaussian(link=sm.families.links.identity))

    glm_results = glm.fit_regularized(alpha=alpha, L1_wt=l1_ratio, refit=False)
    print(glm_results.summary())

    with open(out_fn, 'w') as out_f:
        for leaf, ab in zip(leaves, glm_results.params):
            print(leaf, ab, sep='\t', file=out_f)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Estimate abundances from ProPhyle's assignment using a regularized Poisson GLM"
    )

    parser.add_argument(
        'tree_fn',
        type=str,
        metavar='<tree.nw>',
        help='taxonomic tree (tree.preliminary.nw)'
    )

    parser.add_argument(
        'asg_fn',
        type=str,
        metavar='<pseudo_aln.bam>',
        help='assignments (output of prophyle classify)'
    )

    parser.add_argument(
        'sim_mat_fn',
        type=str,
        metavar='<sim_matrix.npy>',
        help='similarity matrix'
    )

    parser.add_argument(
        'out_fn',
        type=str,
        metavar='<output_fn>',
        help='output file'
    )

    parser.add_argument(
        '-a',
        dest='alpha',
        type=float,
        default=0.1,
        metavar='FLOAT',
        help='regularization weight'
    )

    parser.add_argument(
        '-l',
        dest='l1_ratio',
        type=float,
        default=0.99,
        metavar='FLOAT',
        help='l1 ratio'
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
        estimate_abundances(
            tree_fn=args.tree_fn,
            asg_fn=args.asg_fn,
            sim_mat_fn=args.sim_mat_fn,
            out_fn=args.out_fn,
            alpha=args.alpha,
            l1_ratio=args.l1_ratio,
        )

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


if __name__ == '__main__':
    main()
