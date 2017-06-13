#! /usr/bin/env python3

"""Main ProPhyle file.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Example:

	Download sequences:

		$ prophyle download bacteria

	Create an index for k=10 and the small testing bacterial tree:

		$ prophyle index -k 10 ~/prophyle/test_bacteria.nw ~/prophyle/test_viruses.nw test_idx

	Classify some reads:

		$ prophyle classify test_idx reads.fq > result.sam

TODO:
	* save configuration (trees, k, etc.) into a json; if anything changed from the last time, remove all marks
	* _is_complete should be combined with a test of files: is_missing => remove mark
	* index: kmer annotation to the tree
	* classification: support for c2, h2
"""

import argparse
import hashlib
import multiprocessing
import os
import sys
import textwrap

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

GITDIR = os.path.basename(sys.argv[0])[-3:] == ".py"
if GITDIR:
	C_D = os.path.abspath(os.path.dirname(sys.argv[0]))
else:
	C_D = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

TREE_D = os.path.join(C_D, "trees")

BWA = os.path.join(C_D, "prophyle_index", "bwa", "bwa")
IND = os.path.join(C_D, "prophyle_index", "prophyle_index")
ASM = os.path.join(C_D, "prophyle_assembler", "prophyle_assembler")

# git
if GITDIR:
	ASSIGN = os.path.join(C_D, "prophyle_assignment.py")
	PROPAGATION_POSTPROCESSING = os.path.join(C_D, "prophyle_propagation_postprocessing.py")
	PROPAGATION_PREPROCESSING = os.path.join(C_D, "prophyle_propagation_preprocessing.py")
	NEWICK2MAKEFILE = os.path.join(C_D, "prophyle_propagation_makefile.py")
	READ = os.path.join(C_D, "prophyle_paired_end.py")
	TEST_TREE = os.path.join(C_D, "prophyle_validate_tree.py")

# package
else:
	ASSIGN = "prophyle_assignment.py"
	PROPAGATION_POSTPROCESSING = "prophyle_propagation_postprocessing.py"
	PROPAGATION_PREPROCESSING = "prophyle_propagation_preprocessing.py"
	NEWICK2MAKEFILE = "prophyle_propagation_makefile.py"
	READ = "prophyle_paired_end.py"
	TEST_TREE = "prophyle_validate_tree.py"

DEFAULT_K = 31
DEFAULT_THREADS = multiprocessing.cpu_count()
# DEFAULT_THREADS=1
DEFAULT_MEASURE = 'h1'
DEFAULT_OUTPUT_FORMAT = 'sam'
DEFAULT_HOME_DIR = os.path.join(os.path.expanduser('~'), 'prophyle')

LIBRARIES = ['bacteria', 'viruses', 'plasmids', 'hmp']

FTP_NCBI = 'https://ftp.ncbi.nlm.nih.gov'


def _file_md5(fn, block_size=2 ** 20):
	md5 = hashlib.md5()
	with open(fn, 'rb') as f:
		while True:
			data = f.read(block_size)
			if not data:
				break
			md5.update(data)
	return md5.hexdigest()


def _log_file_md5(fn, remark=None):
	md5 = _file_md5(fn)
	size = pro.file_sizes(fn)[0]
	m = "File {}{} has md5 checksum {} and size {} B".format(
		os.path.basename(fn),
		" ({})".format(remark) if remark is not None else "",
		md5,
		size,
	)
	pro.message(m, only_log=True)


def _test_tree(fn):
	"""Test if given tree is valid for ProPhyle.

	Args:
		fn (str): Newick/NHX tree.

	Raises:
		AssertionError: The tree is not valid.
	"""
	tree = pro.load_nhx_tree(fn, validate=False)
	assert pro.validate_prophyle_nhx_tree(tree, verbose=False, throw_exceptions=False, output_fo=sys.stderr)


def _compile_prophyle_bin():
	"""Compile ProPhyle binaries if they don't exist yet. Recompile if not up-to-date.
	"""

	try:
		command = ["make", "-C", C_D]
		pro.run_safe(command, output_fo=sys.stderr)
	except RuntimeError:
		if not os.path.isfile(IND) or not os.path.isfile(ASM):
			print("Error: ProPhyle executables could not be compiled. Please, the command '{}' manually.".format(
				" ".join(command)), file=sys.stderr)
			sys.exit(1)
		else:
			print("Warning: ProPhyle executables could not be recompiled. Going to use the old ones.", file=sys.stderr)


#####################
# PROPHYLE DOWNLOAD #
#####################

def __mark_fn(d, i, name):
	"""Create a mark name.

	Args:
		d (str): Directory.
		i (int): Number of the step.
		name (str): Name of the mark.
	"""
	if name is None:
		return os.path.join(d, ".complete.{}".format(i))
	else:
		return os.path.join(d, ".complete.{}.{}".format(name, i))


def _mark_complete(d, i=1, name=None):
	"""Create a mark file (an empty file to mark a finished step nb i).

	Args:
		d (str): Directory.
		i (int): Number of the step.
		name (str): Name of the mark.
	"""

	assert i > 0

	pro.touch(__mark_fn(d, i, name))


def _is_complete(d, i=1, name=None):
	"""Check if a mark file i exists AND is newer than the mark file (i-1).

	Args:
		d (str): Directory.
		i (int): Number of the step.
		name (str): Name of the mark.
	"""

	assert i > 0
	fn = __mark_fn(d, i, name)
	fn0 = __mark_fn(d, i - 1, name)

	if i == 1:
		return os.path.isfile(fn)
	else:
		return pro.existing_and_newer(fn0, fn)


def _missing_library(d):
	"""Check if library has been already downloaded.

	Args:
		d (str): Directory.
	"""

	l = os.path.dirname(d)
	pro.makedirs(d)
	if _is_complete(d, 1):
		pro.message("Skipping downloading library '{}' (already exists)".format(l))
		return False
	else:
		pro.message("Downloading library '{}'".format(l))
		return True


def _pseudo_fai(d):
	"""Generate a psedudofai file for given directory (directory/*.fa => directory.fai).

	Pseudofai format = TSV with 2 two columns: filename, sequence header (text after > in FASTA).

	Args:
		d (str): Directory.
	"""
	l = os.path.dirname(d)
	pseudofai_fn = d + ".pseudofai"
	pro.makedirs(d)
	if _is_complete(d, 2) and os.path.isfile(pseudofai_fn):
		pro.message("Skipping generating pseudofai for library '{}' (already exists)".format(l))
	else:
		pro.message("Generating pseudofai for library '{}'".format(l))
		assert d[-1] != "/"
		# cmd=['grep -r --include=\\*.{fa,ffn,fna}', '">"', d, '| sed "s/:>/\t/"']
		cmd = [
			'find', d, '-name', "'*.fa'", "-o", "-name", "'*.ffn'", "-o", "-name", "'*.fna'", "-exec", "grep", "-H",
			'">"', "{}", "\\;",
			"|", 'sed', '"s/\:>/\t/"']

		pro.run_safe(cmd, output_fn=pseudofai_fn)
		_mark_complete(d, 2)


def prophyle_download(library, library_dir, force=False):
	"""Create a library Download genomic library and copy the corresponding tree.

	Args:
		library (str): Library to download (bacteria / viruses / ...)
		library_dir (str): Directory where download files will be downloaded.

	TODO:
		* Add support for alternative URLs (http / ftp, backup refseq sites, etc.).
			* http://downloads.hmpdacc.org/data/HMREFG/all_seqs.fa.bz2
			* ftp://public-ftp.hmpdacc.org/HMREFG/all_seqs.fa.bz2
	"""

	if library == "all":
		for l in LIBRARIES:
			prophyle_download(l, library_dir, force)
		return
	else:
		assert library in LIBRARIES

	if library_dir is None:
		d = os.path.join(os.path.expanduser("~/prophyle"), library)
	else:
		d = os.path.join(library_dir, library)
	# print('making',d, file=sys.stderr)
	# os.makedirs(d, exist_ok=True)
	pro.makedirs(d)

	pro.message("Checking library '{}' in '{}'".format(library, d))

	lib_missing = _missing_library(d)
	if lib_missing or force:
		for test_prefix in ["", "test_"]:
			fn = "{}{}.nw".format(test_prefix, library, )
			nhx = os.path.join(TREE_D, fn)
			new_nhx = os.path.join(d, "..", fn)
			pro.test_files(nhx)
			pro.message("Copying Newick/NHX tree '{}' to '{}'".format(nhx, new_nhx))
			pro.cp_to_file(nhx, new_nhx)

	if library == 'bacteria':
		if lib_missing or force:
			cmd = ['cd', d, '&&', 'curl', FTP_NCBI + '/genomes/archive/old_refseq/Bacteria/all.fna.tar.gz', '|', 'tar',
				'xz']
			pro.run_safe(cmd)
			_mark_complete(d, 1)
		# _pseudo_fai(d)

	elif library == 'viruses':
		if lib_missing or force:
			# cmd=['cd', d, '&&', 'curl', FTP_NCBI+'/genomes/Viruses/all.ffn.tar.gz', '|', 'tar', 'xz']
			# pro.run_safe(cmd)
			cmd = ['cd', d, '&&', 'curl', FTP_NCBI + '/genomes/Viruses/all.fna.tar.gz', '|', 'tar', 'xz']
			pro.run_safe(cmd)
			_mark_complete(d, 1)
		# _pseudo_fai(d)

	elif library == 'plasmids':
		if lib_missing or force:
			cmd = ['cd', d, '&&', 'curl', FTP_NCBI + '/genomes/archive/old_refseq/Plasmids/plasmids.all.fna.tar.gz',
				'|', 'tar', 'xz', '--strip', '5']
			pro.run_safe(cmd)
			_mark_complete(d, 1)
		# _pseudo_fai(d)

	elif library == 'hmp':
		if lib_missing or force:
			# fix when error appears
			cmd = ['cd', d, '&&', 'curl', 'http://downloads.hmpdacc.org/data/HMREFG/all_seqs.fa.bz2', '|', 'bzip2',
				'-d']
			pro.run_safe(cmd, os.path.join(d, "all_seqs.fa"))
			_mark_complete(d, 1)
		# _pseudo_fai(d)

	else:
		raise ValueError('Unknown library "{}"'.format(library))


##################
# PROPHYLE INDEX #
##################

def _create_makefile(index_dir, k, library_dir, mask_repeats=False):
	"""Create a Makefile for k-mer propagation.

	Args:
		index_dir (str): Index directory.
		k (int): K-mer size.
		library_dir (library_dir): Library directory.
		mask_repeats (bool): Mask repeats using DustMasker.

	TODO:
		* Add checking of params.mk
	"""
	pro.message('Creating Makefile for k-mer propagation')
	propagation_dir = os.path.join(index_dir, 'propagation')
	pro.makedirs(propagation_dir)

	makefile = os.path.join(propagation_dir, 'Makefile')
	tree_fn = os.path.join(index_dir, 'tree.preliminary.nw')
	_test_tree(tree_fn)
	# pro.test_files(NEWICK2MAKEFILE, tree_fn)
	command = [NEWICK2MAKEFILE, '-k', k, tree_fn, os.path.abspath(library_dir), './', makefile]

	with open(os.path.join(propagation_dir, "params.mk"), "w+") as f:
		f.write('PRG_ASM="{}"\n'.format(ASM))
		f.write("K={}\n".format(k))
		if mask_repeats:
			f.write("MASKREP=1\n")
	pro.run_safe(command)
	_log_file_md5(makefile)


def _propagate(index_dir, threads):
	"""Run k-mer propagation.

	Args:
		index_dir (str): Index directory.
		threads (int): Number of threads for Makefile.
	"""
	pro.message('Running k-mer propagation')
	propagation_dir = os.path.join(index_dir, 'propagation')
	pro.test_files(os.path.join(propagation_dir, 'Makefile'), test_nonzero=True)

	# test if input files for propagation exist
	command = ['make', '-C', propagation_dir, '-n', '-s', '>', '/dev/null']
	pro.run_safe(
		command,
		err_msg="Some FASTA files needed for k-mer propagation are probably missing, see the messages above.",
		thr_exc=False,
	)

	# run propagation
	command = ['make', '-j', threads, '-C', propagation_dir, 'V=1']
	pro.run_safe(
		command,
		err_msg="K-mer propagation has not been finished because of an error. See messages above.",
		thr_exc=False,
	)


def _kmer_stats(index_dir):
	"""Create a file with k-mer statistics.

	Args:
		index_dir (str): Index directory.
	"""
	propagation_dir = os.path.join(index_dir, 'propagation')
	command = ["cat", propagation_dir + "/*.count.tsv", "|", "grep", "-v", "^#", "|", "sort", "|", "uniq", ">",
		os.path.join(index_dir, "index.fa.kmers.tsv")]
	pro.run_safe(
		command,
		err_msg="A file with k-mer statistics could not be created.",
		thr_exc=False,
	)


def _propagation_preprocessing(in_trees, out_tree, no_prefixes, sampling_rate):
	"""Merge input trees into a single tree.

	Args:
		in_trees (list of str): Input NHX trees (possibly with a root specifier).
		out_tree (str): Output NHX tree.
		no_prefixes (bool): Don't prepend prefixes to node names during tree merging.
		sampling rate (float): Sampling rate for subsampling the tree or None for no subsampling.
	"""

	pro.message('Generating index tree')
	# existence already checked
	# pro.test_files(*in_trees)
	command = [PROPAGATION_PREPROCESSING]
	if sampling_rate is not None:
		command += ['-s', sampling_rate]
	command += in_trees + [out_tree]
	if no_prefixes:
		command += ['-P']
	pro.run_safe(
		command,
		err_msg="The main tree could not be generated.",
		thr_exc=False,
	)
	_log_file_md5(out_tree)


def _remove_tmp_propagation_files(index_dir):
	"""Run k-mer propagation.

	Args:
		index_dir (str): Index directory.
	"""
	pro.message('Removing temporary files')
	propagation_dir = os.path.join(index_dir, 'propagation')

	command = ['make', '-C', propagation_dir, 'clean', '>', '/dev/null']
	pro.run_safe(command)


def _propagation_postprocessing(index_dir, in_tree_fn, out_tree_fn):
	"""Merge reduced FASTA files after k-mer propagation and create index.fa.

	Args:
		index_dir (str): Index directory.
		in_tree_fn (str): Input tree in Newick/NHX.
		out_tree_fn (str): Output tree in Newick/NHX.
	"""

	pro.message('Propagation post-processing')

	propagation_dir = os.path.join(index_dir, 'propagation')
	tsv_fn = os.path.join(index_dir, "index.fa.kmers.tsv")
	index_fa = os.path.join(index_dir, "index.fa")

	command = ["cat", os.path.join(propagation_dir, "*.tsv"), '>', tsv_fn]
	pro.run_safe(
		command,
		err_msg="K-mer statistics could not be created.",
		thr_exc=True,
	)

	command = [PROPAGATION_POSTPROCESSING, propagation_dir, index_fa, in_tree_fn, tsv_fn, out_tree_fn]
	pro.run_safe(
		command,
		err_msg="Main ProPhyle FASTA file could not be generated",
		thr_exc=True,
	)
	pro.touch(index_fa + ".complete")
	_log_file_md5(index_fa)
	_log_file_md5(in_tree_fn)
	_log_file_md5(out_tree_fn)


def _fa2pac(fa_fn):
	"""Run `bwa fa2pac` (FA => 2bit).

	Args:
		fa_fn (str): FASTA file.
	"""

	pro.message('Generating packed FASTA file')
	pro.test_files(BWA, fa_fn)
	command = [BWA, 'fa2pac', fa_fn, fa_fn]
	pro.run_safe(
		command,
		err_msg="Packaged file could not be created.",
		thr_exc=True,
	)
	_log_file_md5(fa_fn + ".pac")


def _pac2bwt(fa_fn):
	"""Run `bwa pac2bwtgen` (2bit => BWT).

	Args:
		fa_fn (str): FASTA file.
	"""

	pro.message('Generating BWT')
	pro.test_files(BWA, fa_fn + ".pac")
	command = [BWA, 'pac2bwtgen', fa_fn + ".pac", fa_fn + ".bwt"]
	pro.run_safe(
		command,
		err_msg="Burrows-Wheeler Transform could not be computed.",
		thr_exc=True,
	)
	_log_file_md5(fa_fn + ".bwt", remark="without OCC")


def _bwt2bwtocc(fa_fn):
	"""Run `bwa bwtupdate` (BWT => BWT+OCC).

	Args:
		fa_fn (str): FASTA file.
	"""

	pro.message('Generating sampled OCC array')
	pro.test_files(BWA, fa_fn + ".bwt")
	command = [BWA, 'bwtupdate', fa_fn + ".bwt"]
	pro.run_safe(
		command,
		err_msg="OCC array could not be computed.",
		thr_exc=True,
	)
	_log_file_md5(fa_fn + ".bwt", remark="with OCC")


def _bwtocc2sa(fa_fn):
	"""Run `bwa bwt2sa` (BWT+, remark="with OCC"OCC => SSA).

	Args:
		fa_fn (str): FASTA file.
	"""

	pro.message('Generating sampled SA')
	pro.test_files(BWA, fa_fn + ".bwt", remark="with OCC")
	command = [BWA, 'bwt2sa', fa_fn + ".bwt", fa_fn + ".sa"]
	pro.run_safe(
		command,
		err_msg="Sampled Suffix Array computation failed.",
		thr_exc=True,
	)
	_log_file_md5(fa_fn + ".sa")


def _bwtocc2klcp(fa_fn, k):
	"""Create k-LCP `` (BWT => k-LCP).

	Args:
		fa_fn (str): FASTA file.
		k (int): K-mer size.
	"""

	pro.message('Generating k-LCP array')
	pro.test_files(IND, fa_fn + ".bwt")
	command = [IND, 'build', '-k', k, fa_fn]
	pro.run_safe(
		command,
		err_msg="k-Longest Common Prefix array construction failed.",
		thr_exc=True,
	)
	_log_file_md5("{}.{}.klcp".format(fa_fn, k))


def _bwtocc2sa_klcp(fa_fn, k):
	"""Create k-LCP `` (BWT => k-LCP).

	Args:
		fa_fn (str): FASTA file.
		k (int): K-mer size.
	"""

	pro.message('Generating k-LCP array and SA in parallel')
	pro.test_files(IND, fa_fn + ".bwt")
	command = [IND, 'build', '-s', '-k', k, fa_fn]
	pro.run_safe(
		command,
		err_msg="Parallel construction of k-Longest Common Prefix array and Sampled Suffix Array failed.",
		thr_exc=True,
	)
	_log_file_md5(fa_fn + ".sa")
	_log_file_md5("{}.{}.klcp".format(fa_fn, k))


def prophyle_index(index_dir, threads, k, trees_fn, library_dir, construct_klcp, force, no_prefixes, mask_repeats,
		keep_tmp_files, sampling_rate):
	"""Build a ProPhyle index.

	Args:
		index_dir (str): Index directory.
		threads (int): Number of threads in k-mer propagation.
		k (int): K-mer size.
		trees_fn (list of str): Newick/NHX tree, possibly with a root spec (@root).
		library_dir (str): Library directory.
		klcp (bool): Generate klcp.
		force (bool): Rewrite files if they already exist.
		no_prefixes (bool): Don't prepend prefixes to node names during tree merging.
		mask_repeats (bool): Mask repeats using DustMasker.
		keep_tmp_files (bool): Keep temporary files from k-mer propagation.
		sampling rate (float): Sampling rate for subsampling the tree or None for no subsampling.
	"""

	assert isinstance(k, int)
	assert isinstance(threads, int)
	assert k > 1
	assert threads > 0
	assert sampling_rate is None or 0.0 <= float(sampling_rate) <= 1.0

	_compile_prophyle_bin()

	index_fa = os.path.join(index_dir, 'index.fa')
	index_tree_1 = os.path.join(index_dir, 'tree.preliminary.nw')
	index_tree_2 = os.path.join(index_dir, 'tree.nw')

	# recompute = recompute everything from now on
	# force==True => start to recompute everything from beginning
	recompute = force

	# make index dir
	pro.makedirs(index_dir)

	#
	# 1) Newick
	#

	if not _is_complete(index_dir, 1) or not pro.existing_and_newer_list(trees_fn, index_tree_1):
		recompute = True

	if recompute:
		pro.message('[1/5] Copying/merging trees', upper=True)
		for tree_fn in trees_fn:
			tree_fn, _, root = tree_fn.partition("@")
			tree = pro.load_nhx_tree(tree_fn)
			pro.validate_prophyle_nhx_tree(tree)
			if root != "":
				assert len(tree.search_nodes(name=root)) != 0, "Node '{}' does not exist in '{}'.".format(root, tree_fn)
		if len(trees_fn) != 1:
			pro.message('Merging {} trees{}'.format(len(trees_fn)))
		_propagation_preprocessing(trees_fn, index_tree_1, no_prefixes=no_prefixes, sampling_rate=sampling_rate)
		_mark_complete(index_dir, 1)
	else:
		pro.message('[1/5] Tree already exists, skipping copying', upper=True)

	#
	# 2) Create and run Makefile for propagation, and merge FASTA files
	#

	if not _is_complete(index_dir, 2):
		recompute = True

	if recompute:
		pro.message('[2/5] Running k-mer propagation', upper=True)
		_create_makefile(index_dir, k, library_dir, mask_repeats=mask_repeats)
		_propagate(index_dir, threads=threads)
		_propagation_postprocessing(index_dir, index_tree_1, index_tree_2)
		_kmer_stats(index_dir)
		if not keep_tmp_files:
			_remove_tmp_propagation_files(index_dir)
		else:
			pro.message('Keeping temporary files')
		_mark_complete(index_dir, 2)
	else:
		pro.message('[2/5] K-mers have already been propagated, skipping propagation', upper=True)

	#
	# 3) BWT + OCC
	#

	if not _is_complete(index_dir, 3):
		recompute = True

	# if ccontinue and os.path.isfile(index_fa+'.bwt') and os.path.isfile(index_fa+'.bwt.complete'):

	if recompute:
		pro.message('[3/5] Constructing BWT+OCC', upper=True)
		pro.rm(index_fa + '.bwt', index_fa + '.bwt.complete')
		_fa2pac(index_fa)
		_pac2bwt(index_fa)
		_bwt2bwtocc(index_fa)
		_mark_complete(index_dir, 3)
	else:
		pro.message('[3/5] BWT and OCC already exist, skipping their construction', upper=True)

	#
	# 4) SA + 5) KLCP (compute SA + KLCP in parallel)
	#

	klcp_fn = "{}.{}.klcp".format(index_fa, k)

	if construct_klcp:

		if not _is_complete(index_dir, 4):
			# SA not computed yet => compute it in parallel with KLCP
			recompute = True

		if recompute:
			pro.message('[4/5],[5/5] Constructing SA + KLCP in parallel ', upper=True)
			_bwtocc2sa_klcp(index_fa, k)
			_mark_complete(index_dir, 4)
			_mark_complete(index_dir, 5)
			return

	#
	# 4) SA (compute only SA)
	#

	if not _is_complete(index_dir, 4):
		recompute = True

	if recompute:
		pro.message('[4/5] Constructing SA', upper=True)
		_bwtocc2sa(index_fa)
	else:
		pro.message('[4/5] SA already exists, skipping its construction', upper=True)

	#
	# 5) KLCP (compute only KLCP)
	#

	if construct_klcp:
		if not _is_complete(index_dir, 5):
			recompute = True

		if recompute:
			pro.message('[5/5] Constructing k-LCP', upper=True)
			_bwtocc2klcp(index_fa, k)
			_mark_complete(index_dir, 5)
		else:
			pro.message('[5/5] k-LCP already exists, skipping its construction', upper=True)


#####################
# PROPHYLE CLASSIFY #
#####################

def prophyle_classify(index_dir, fq_fn, fq_pe_fn, k, use_rolling_window, out_format, mimic_kraken, measure, annotate,
		tie_lca,
		print_seq):
	"""Run ProPhyle classification.

	Args:
		index_dir (str): Index directory.
		fq_fn (str): Input reads (single-end or first of paired-end).
		fq_pe_fn (str): Input reads (second paired-end, None if single-end)
		k (int): K-mer size (None => detect automatically).
		use_rolling_window (bool): Use rolling window.
		out_format (str): Output format: sam / kraken.
		mimic_kraken (bool): Mimic Kraken algorithm (compute LCA for each k-mer).
		measure (str): Measure used for classification (h1 / h2 / c1 / c2).
		annotate (bool): Annotate assignments (insert annotations from Newick to SAM).
		tie_lca (bool): If multiple equally good assignments found, compute their LCA.
		print_seq (bool): Print sequencing in SAM.
	"""

	_compile_prophyle_bin()
	index_fa = os.path.join(index_dir, 'index.fa')
	index_tree = os.path.join(index_dir, 'tree.nw')

	if k is None:
		k = pro.detect_k_from_index(index_dir)
		pro.message("Automatic detection of k-mer length: k={}".format(k))

	_test_tree(index_tree)

	if fq_pe_fn:
		pro.test_files(fq_fn, fq_pe_fn)
	elif fq_fn != '-':
		pro.test_files(fq_fn)

	pro.test_files(index_fa, IND)

	pro.test_files(
		index_fa + '.bwt',
		index_fa + '.pac',
		index_fa + '.sa',
		index_fa + '.ann',
		index_fa + '.amb',
	)

	(bwt_s, sa_s, pac_s) = pro.file_sizes(index_fa + '.bwt', index_fa + '.sa', index_fa + '.pac')
	assert abs(bwt_s - 2 * sa_s) < 1000, 'Inconsistent index (SA vs. BWT)'
	assert abs(bwt_s - 2 * pac_s) < 1000, 'Inconsistent index (PAC vs. BWT)'

	if use_rolling_window:
		klcp_fn = "{}.{}.klcp".format(index_fa, k)
		pro.test_files(klcp_fn)
		(klcp_s,) = pro.file_sizes(klcp_fn)
		assert abs(bwt_s - 4 * klcp_s) < 1000, 'Inconsistent index (KLCP vs. BWT)'

	if mimic_kraken:
		cmd_assign = [ASSIGN, '-i', '-', '-k', k, '-n', index_tree, '-m', 'h1', '-f', 'kraken', '-l', '-t']
	else:
		cmd_assign = [ASSIGN, '-i', '-', '-k', k, '-n', index_tree, '-m', measure, '-f', out_format]
		if annotate:
			cmd_assign += ['--annotate']
		if tie_lca:
			cmd_assign += ['--tie-lca']

	if fq_pe_fn:
		cmd_read = [READ, fq_fn, fq_pe_fn, '|']
		in_read = '-'
	else:
		cmd_read = []
		# fq_fn can be '-' as well
		in_read = fq_fn

	cmd_query = [IND, 'query', '-k', k, '-u' if use_rolling_window else '', '-b' if print_seq else '', index_fa, in_read,
		'|']

	command = cmd_read + cmd_query + cmd_assign
	pro.run_safe(command)


########
# MAIN #
########

def parser():
	class MyParser(argparse.ArgumentParser):
		def error(self, message):
			if len(sys.argv) == 2:
				self.print_help()
			else:
				print('error: {}'.format(message), file=sys.stderr)
			sys.exit(2)

	desc = """\
		Program: prophyle (phylogeny-based metagenomic classification)
		Version: {V}
		Authors: Karel Brinda <kbrinda@hsph.harvard.edu>, Kamil Salikhov <kamil.salikhov@univ-mlv.fr>,
		         Simone Pignotti <pignottisimone@gmail.com>, Gregory Kucherov <gregory.kucherov@univ-mlv.fr>

		Usage:   prophyle <command> [options]
		""".format(V=version.VERSION)
	parser = MyParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(desc))

	parser.add_argument('-v', '--version',
		action='version',
		version='%(prog)s {}'.format(version.VERSION),
	)

	subparsers = parser.add_subparsers(help="", description=argparse.SUPPRESS, dest='subcommand', metavar="")
	fc = lambda prog: argparse.HelpFormatter(prog, max_help_position=27)

	##########

	parser_download = subparsers.add_parser(
		'download',
		help='download a genomic database',
		# description='Download RefSeq and HMP databases.',
		formatter_class=fc,
	)

	parser_download.add_argument(
		'library',
		metavar='<library>',
		nargs='+',
		choices=LIBRARIES + ['all'],
		help='genomic library {}'.format(LIBRARIES + ['all']),
	)

	parser_download.add_argument(
		'-d',
		metavar='DIR',
		dest='home_dir',
		type=str,
		default=None,
		help='directory for the tree and the sequences [~/prophyle]',
	)

	parser_download.add_argument(
		'-l',
		dest='log_fn',
		metavar='STR',
		type=str,
		help='log file',
		default=None,
	)

	parser_download.add_argument(
		'-F',
		dest='force',
		action='store_true',
		help='rewrite library files if they already exist',
	)

	##########

	parser_index = subparsers.add_parser(
		'index',
		help='build index',
		formatter_class=fc,
	)

	parser_index.add_argument(
		'tree',
		metavar='<tree.nw>',
		type=str,
		nargs='+',
		help='phylogenetic tree (in Newick/NHX)',
	)

	parser_index.add_argument(
		'index_dir',
		metavar='<index.dir>',
		type=str,
		help='index directory (will be created)',
	)

	parser_index.add_argument(
		'-g',
		metavar='DIR',
		dest='library_dir',
		type=str,
		help='directory with the library sequences [dir. of the first tree]',
		default=None,
		# required=True,
	)

	parser_index.add_argument(
		'-j',
		metavar='INT',
		dest='threads',
		type=int,
		help='number of threads [auto ({})]'.format(DEFAULT_THREADS),
		default=DEFAULT_THREADS,
	)

	parser_index.add_argument(
		'-k',
		dest='k',
		metavar='INT',
		type=int,
		help='k-mer length [{}]'.format(DEFAULT_K),
		default=DEFAULT_K,
	)

	parser_index.add_argument(
		'-l',
		dest='log_fn',
		metavar='STR',
		type=str,
		help='log file [<index.dir>/log.txt]',
		default=None,
	)

	parser_index.add_argument(
		'-s',
		metavar='FLOAT',
		help='rate of sampling of the tree [no sampling]',
		dest='sampling_rate',
		type=str,
		default=None,
	)

	parser_index.add_argument(
		'-F',
		dest='force',
		action='store_true',
		help='rewrite index files if they already exist',
	)

	parser_index.add_argument(
		'-M',
		action='store_true',
		dest='mask_repeats',
		help='mask repeats/low complexity regions (using DustMasker)',
	)

	parser_index.add_argument(
		'-P',
		dest='no_prefixes',
		action='store_true',
		help='do not add prefixes to node names when multiple trees are used',
	)

	parser_index.add_argument(
		'-K',
		dest='klcp',
		action='store_false',
		help='skip k-LCP construction',
	)

	parser_index.add_argument(
		'-T',
		dest='keep_tmp_files',
		action='store_true',
		help='keep temporary files from k-mer propagation',
	)

	##########

	parser_classify = subparsers.add_parser(
		'classify',
		help='classify reads',
		# description='Classify reads.',
		formatter_class=fc,
	)

	parser_classify.add_argument(
		'index_dir',
		metavar='<index.dir>',
		type=str,
		help='index directory',
	)

	parser_classify.add_argument(
		'reads',
		metavar='<reads1.fq>',
		type=str,
		help='first file with reads in FASTA or FASTQ (use - for standard input)',
	)

	parser_classify.add_argument(
		'reads_pe',
		metavar='<reads2.fq>',
		type=str,
		help='second file with reads in FASTA or FASTQ',
		nargs='?',
		default=None,
	)

	parser_classify.add_argument(
		'-k',
		dest='k',
		metavar='INT',
		type=int,
		help='k-mer length [detect automatically from the index]',
		default=None,
	)

	parser_classify.add_argument(
		'-R',
		dest='rolling_window',
		action='store_false',
		help='use restarted search for matching rather than rolling window (slower, but k-LCP is not needed)',
	)

	parser_classify.add_argument(
		'-m',
		dest='measure',
		choices=['h1', 'c1'],
		help='measure: h1=hit count, c1=coverage [{}]'.format(DEFAULT_MEASURE),
		default=DEFAULT_MEASURE,
	)

	parser_classify.add_argument(
		'-f',
		dest='oform',
		choices=['kraken', 'sam'],
		default=DEFAULT_OUTPUT_FORMAT,
		help='output format [{}]'.format(DEFAULT_OUTPUT_FORMAT),
	)

	parser_classify.add_argument(
		'-l',
		dest='log_fn',
		metavar='STR',
		type=str,
		help='log file',
		default=None,
	)

	parser_classify.add_argument(
		'-A',
		dest='annotate',
		action='store_true',
		help='annotate assignments',
	)

	parser_classify.add_argument(
		'-L',
		dest='tie',
		action='store_true',
		help='use LCA when tie (multiple hits with the same score)',
	)

	parser_classify.add_argument(
		'-M',
		dest='mimic',
		action='store_true',
		# help='mimic Kraken algorithm and output (for debugging purposes)',
		help=argparse.SUPPRESS,
	)

	parser_classify.add_argument(
		'-P',
		dest='print_seq',
		action='store_true',
		help='print sequences and qualities in SAM (otherwise \'*\' is used)',
	)

	##########

	return parser


def main():
	try:
		par = parser()
		args = par.parse_args()
		subcommand = args.subcommand

		if subcommand == "download":
			pro.open_log(args.log_fn)
			for single_lib in args.library:
				pro.message('Downloading "{}" started'.format(single_lib))
				prophyle_download(
					library=single_lib,
					library_dir=args.home_dir,
					force=args.force,
				)
				pro.message('Downloading "{}" finished'.format(single_lib))
			pro.close_log()

		elif subcommand == "index":
			if args.library_dir is None:
				library_dir = os.path.dirname(args.tree[0])
			else:
				library_dir = args.library_dir

			if args.log_fn is None:
				args.log_fn = os.path.join(args.index_dir, "log.txt")

			pro.open_log(args.log_fn)
			pro.message('Index construction started')
			prophyle_index(
				index_dir=args.index_dir,
				threads=args.threads,
				k=args.k,
				trees_fn=args.tree,
				library_dir=library_dir,
				force=args.force,
				construct_klcp=args.klcp,
				no_prefixes=args.no_prefixes,
				mask_repeats=args.mask_repeats,
				keep_tmp_files=args.keep_tmp_files,
				sampling_rate=args.sampling_rate,
			)
			pro.message('Index construction finished')
			pro.close_log()

		elif subcommand == "classify":
			# if args.log_fn is None:
			#	args.log_fn = os.path.join(args.index_dir, "log.txt")

			pro.open_log(args.log_fn)
			pro.message('Classification started')
			prophyle_classify(
				index_dir=args.index_dir,
				fq_fn=args.reads,
				fq_pe_fn=args.reads_pe,
				k=args.k,
				use_rolling_window=args.rolling_window,
				out_format=args.oform,
				mimic_kraken=args.mimic,
				measure=args.measure,
				tie_lca=args.tie,
				annotate=args.annotate,
				print_seq=args.print_seq,
			)
			pro.message('Classification finished')
			pro.close_log()

		else:
			msg_lns = par.format_help().split("\n")[2:]
			msg_lns = [x for x in msg_lns if x.find("optional arguments") == -1 and x.find("--") == -1]
			msg = "\n".join(msg_lns)
			msg = msg.replace("\n\n", '\n').replace("subcommands:\n", "Command:").replace("Usage", "\nUsage")
			print(file=sys.stderr)
			print(msg, file=sys.stderr)
			sys.exit(1)

	except BrokenPipeError:
		# pipe error (e.g., when head is used)
		sys.stderr.close()
		exit(0)

	except KeyboardInterrupt:
		pro.message("Error: Keyboard interrupt")
		pro.close_log()
		exit(1)


if __name__ == "__main__":
	main()
