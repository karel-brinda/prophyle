#! /usr/bin/env python3

"""Create a Makefile for ProPhyle k-mer propagation.

This script first loads a phylogenetic tree, computes its minimal
subtree (i.e., removes all nodes with one child), traverses the
subtree and prints the corresponding k-mer propagation rules.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Example:

	prophyle_propagation_makefile


Propagation parameters (in the Makefile, can be changed through CL):
	* NONPROP: no k-mer propagation (sequences for leaves only)
	* REASM: re-assemble sequences in leaves
	* NONDEL: non-deletative propagation, implies REASM
	* MASKREP: mask repeats in leaves
"""

import os
import shutil
import datetime
import sys
import argparse
import textwrap

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro


def _compl(fn):
	"""Get complete marker file name.

	Args:
		fn (str): Original file name.
	"""
	return fn + ".complete"


def _compl_l(fns):
	"""Get complete marker file names.

	Args:
		fns (list of str): Original file names.
	"""
	return [_compl(x) for x in fns]


def merge_fasta_files(input_files_fn, output_file_fn, is_leaf, makefile_fo, nhx_file_fn=None):
	"""Print Makefile lines for merging FASTA files and removing empty lines.

	Args:
		input_files_fn (list of str): List of input files.
		output_file_fn (str): Output file.
		is_leaf (str): Is a leaf (i.e., copying must be done).
		makefile_fo (file): Output file.
		nhx_file_fn (str): File with the tree (for including in rule dependencies).
	"""

	if is_leaf:
		cmd = textwrap.dedent("""\

				{ocompl}: {i}
					cat $^ $(CMD_MASKING) $(CMD_REASM) > {o}
					@touch $@

			""".format(
			i=' '.join(input_files_fn),
			o=output_file_fn,
			ocompl=_compl(output_file_fn),
		))
	else:

		cmd = textwrap.dedent("""\

				{ocompl}: {icompl} {nhx}
					cat {i} > {o}
					@touch $@

			""".format(
			i=' '.join(input_files_fn),
			icomp=' '.join(_compl_l(input_files_fn)),
			o=output_file_fn,
			ocomp=_compl(output_file_fn),
			nhx=nhx_file_fn if nhx_file_fn is not None else "",
		))

	print(
		textwrap.dedent("""\
			#
			# Merging FASTA files: {output_file}
			#
			""".format(output_file=output_file_fn)
		)
		+ cmd,
		file=makefile_fo
	)


def assembly(input_files_fn, output_files_fn, intersection_file_fn, makefile_fo, counts_fn="/dev/null",
		nhx_file_fn=None):
	"""Print Makefile lines for running prophyle_assembler.

	Args:
		input_files_fn (list of str): List of input files.
		output_files_fn (list of str): List of output files.
		intersection_file_fn (str): File with intersection.
		makefile_fo (file): Output file.
		counts_fn (str): File with count statistics.
		nhx_file_fn (str): File with the tree (for including in rule dependencies).
	"""

	assert len(input_files_fn) == len(output_files_fn)
	# assert intersection_file not in input_files
	# print(intersection_file, input_files,file=sys.stderr)
	cmd = textwrap.dedent("""\
			ifdef NONDEL
			   CMD_ASM_OUT_{nid} =
			else
			   CMD_ASM_OUT_{nid} = -o {oo}
			endif

			ifdef NONPROP
			   CMD_ASM_{nid} = @touch {x} {o}
			else
			   CMD_ASM_{nid} = $(PRG_ASM) -S -k $(K) -x {x} -i {ii} $(CMD_ASM_OUT_{nid}) -s {c}
			endif

			{xcompl}: {icompl} {nhx}
				@echo starting propagation for $@
				$(CMD_ASM_{nid})
				@touch $@
			""".format(
		icompl=' '.join(_compl_l(input_files_fn)),
		o=' '.join(output_files_fn),
		ii=' -i '.join(input_files_fn),
		oo=' -o '.join(output_files_fn),
		x=intersection_file_fn,
		xcompl=_compl(intersection_file_fn),
		c=counts_fn,
		nid=intersection_file_fn,
		nhx=nhx_file_fn if nhx_file_fn is not None else "",
	)
	)

	print(
		textwrap.dedent("""\
			#
			# Assemblying FASTA files: {intersection_file}
			#
			""".format(intersection_file=intersection_file_fn)
		) + cmd,
		file=makefile_fo,
	)


class TreeIndex:
	"""Main class for k-mer propagation.
	"""

	def __init__(self, tree_newick_fn, index_dir, library_dir, makefile_fn):
		"""Init the class.

		Args:
			tree_newick_fn (str): Tree file name.
			index_dir (str): Directory of the index.
			library_dir (str): Directory with FASTA files.
			makefile_fn (str): Output Makefile.
		"""

		self.tree_newick_fn = tree_newick_fn
		tree = pro.load_nhx_tree(tree_newick_fn)
		self.tree = pro.minimal_subtree(tree)
		self.newick_dir = os.path.dirname(tree_newick_fn)
		self.index_dir = index_dir
		self.library_dir = library_dir
		self.makefile_fn = makefile_fn
		os.makedirs(self.index_dir, exist_ok=True)

	@staticmethod
	def _node_debug(node):
		if hasattr(node, "common_name") and node.common_name != "":
			return "{}_{}".format(node.name, node.common_name)
		elif hasattr(node, "sci_name") and node.sci_name != "":
			return "{}_{}".format(node.name, node.sci_name)
		else:
			return "{}".format(node.name)

	def nonreduced_fasta_fn(self, node):
		"""Get name of the full FASTA file (k-mer propagation).

		Args:
			node: Node of the tree.
		"""
		return os.path.join(self.index_dir, node.name + ".full.fa")

	def reduced_fasta_fn(self, node):
		"""Get name of the reduced FASTA file (k-mer propagation).

		Args:
			node: Node of the tree.
		"""
		return os.path.join(self.index_dir, node.name + ".reduced.fa")

	def count_fn(self, node):
		"""Get FASTA name of the file with k-mer counts.

		Args:
			node: Node of the tree.
		"""
		return os.path.join(self.index_dir, node.name + ".count.tsv")

	def process_node(self, node, makefile_fo):
		"""Recursive function for treating an individual node of the tree.

		Args:
			node: Node of the tree.
			makefile_fo: Output file.
		"""

		if node.is_leaf():

			if hasattr(node, "fastapath"):
				fastas_fn = node.fastapath.split("@")
				for i in range(len(fastas_fn)):
					fastas_fn[i] = os.path.join(self.library_dir, fastas_fn[i])
				merge_fasta_files(fastas_fn, self.nonreduced_fasta_fn(node), is_leaf=True, makefile_fo=makefile_fo)

		else:
			children = node.get_children()

			# 1) process children
			for child in children:
				self.process_node(child, makefile_fo=makefile_fo)
			# print(child.name, "processed",file=sys.stderr)

			# 2) k-mer propagation & assembly
			input_files = [self.nonreduced_fasta_fn(x) for x in children]
			output_files = [self.reduced_fasta_fn(x) for x in children]
			intersection_file = self.nonreduced_fasta_fn(node)
			count_file = self.count_fn(node)
			assembly(input_files, output_files, intersection_file, counts_fn=count_file, makefile_fo=makefile_fo)

	def build_index(self, k):
		"""Print Makefile for the tree.

		Args:
			k (int): K-mer size.
		"""

		with open(self.makefile_fn, 'w+') as f:
			print(textwrap.dedent("""\
					include params.mk\n

					.PHONY: all clean

					SHELL=/usr/bin/env bash -euc -o pipefail

					PRG_ASM?=prophyle_assembler
					PRG_DUST?=dustmasker

					$(info )
					$(info /------------------------------------------------------------------)

					ifdef K
					   $(info | K-mer length:           $(K))
					else
					   $(error | K-mer length is not specified)
					endif

					$(info | Assembler:              $(PRG_ASM))

					$(info | DustMasker:             $(PRG_DUST))

					ifdef MASKREP
					   $(info | Masking repeats:        On)
					   CMD_MASKING= | $(PRG_DUST) -infmt fasta -outfmt fasta | sed '/^>/! s/[^AGCT]/N/g'
					else
					   $(info | Masking repeats:        Off)
					   CMD_MASKING=
					endif

					ifdef NONPROP
					   $(info | K-mer propagation:      Off)
					else
					   $(info | K-mer propagation:      On)
					endif

					ifdef NONDEL
					   $(info | K-mer propagation mode: Non-deletative)
					   REASM=1
					else
					   $(info | K-mer propagation mode: Deletative)
					endif

					ifdef REASM
					   $(info | Re-assembling leaves:   On)
					   CMD_REASM= | $(PRG_ASM) -k $(K) -S -i - -o -
					else
					   $(info | Re-assembling leaves:   Off)
					   CMD_REASM=
					endif
					$(info \------------------------------------------------------------------)
					$(info )

					all: {root_red_compl}

					clean:
						rm -f *.complete
						rm -f *.fa
						rm -f *.tsv

					{root_red_compl}: {root_nonred_compl}
						ln -s {root_nonred} {root_red}
						@touch $@

					""".format(
				root_nonred=self.nonreduced_fasta_fn(self.tree.get_tree_root()),
				root_nonred_compl=_compl(self.nonreduced_fasta_fn(self.tree.get_tree_root())),
				root_red=self.reduced_fasta_fn(self.tree.get_tree_root()),
				root_red_compl=_compl(self.reduced_fasta_fn(self.tree.get_tree_root())),
			)
			), file=f)

			self.process_node(self.tree.get_tree_root(), makefile_fo=f)


def main():
	parser = argparse.ArgumentParser(description='Create Makefile for parallelized ProPhyle k-mer propagation.')
	parser.add_argument(
		'newick_fn',
		type=str,
		metavar='<tree.nw>',
		help='phylogenetic tree (in Newick/NHX).',
	)
	parser.add_argument(
		'-k',
		type=int,
		metavar='int',
		dest='k',
		required=True,
		help='k-mer length',
	)
	parser.add_argument(
		'library_dir_fn',
		metavar='<library.dir>',
		help='directory with the library',
	)
	parser.add_argument(
		'output_dir_fn',
		type=str,
		metavar='<output.dir>',
		help='output directory for the index',
	)
	parser.add_argument(
		'makefile_fn',
		type=str,
		metavar='<Makefile>',
		help='output Makefile',
	)

	args = parser.parse_args()

	k = args.k
	assert k > 0
	newick_fn = args.newick_fn
	output_dir_fn = args.output_dir_fn
	library_dir_fn = args.library_dir_fn
	makefile_fn = args.makefile_fn

	ti = TreeIndex(
		tree_newick_fn=newick_fn,
		library_dir=library_dir_fn,
		index_dir=output_dir_fn,
		makefile_fn=makefile_fn,
	)
	ti.build_index(
		k=k,
	)


if __name__ == "__main__":
	main()
