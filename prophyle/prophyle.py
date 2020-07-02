#! /usr/bin/env python3
"""Main ProPhyle file.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT

Example:

    Download sequences:

        $ prophyle download bacteria

    Create an index for k=10 and the small testing bacterial tree:

        $ prophyle index -k 10 -s 0.1 ~/prophyle/bacteria.nw ~/prophyle/viruses.nw test_idx

    Classify some reads:

        $ prophyle classify test_idx reads.fq > result.sam
"""

import argparse
import collections
import hashlib
import multiprocessing
import os
import sys
import tarfile
import tempfile
import textwrap

sys.path.append(os.path.dirname(__file__))
import prophylelib as pro
import version

CONFIG = {}

GITDIR = os.path.basename(sys.argv[0])[-3:] == ".py"
if GITDIR:
    C_D = os.path.abspath(os.path.dirname(sys.argv[0]))
else:
    C_D = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

TREE_D = os.path.join(C_D, "trees")

BWA = os.path.join(C_D, "prophyle_index", "bwa", "bwa")
IND = os.path.join(C_D, "prophyle_index", "prophyle_index")
ASM = os.path.join(C_D, "prophyle_assembler", "prophyle_assembler")
C_ASSIGN = os.path.join(C_D, "prophyle_assignment", "prophyle_assignment")

# executed from the git repo
if GITDIR:
    PROPHYLE = os.path.join(C_D, "prophyle.py")
    PY_ASSIGN = os.path.join(C_D, "prophyle_assignment.py")
    ANALYZE = os.path.join(C_D, "prophyle_analyze.py")
    PROPAGATION_POSTPROCESSING = os.path.join(C_D, "prophyle_propagation_postprocessing.py")
    PROPAGATION_PREPROCESSING = os.path.join(C_D, "prophyle_propagation_preprocessing.py")
    NEWICK2MAKEFILE = os.path.join(C_D, "prophyle_propagation_makefile.py")
    READ = os.path.join(C_D, "prophyle_paired_end.py")
    TEST_TREE = os.path.join(C_D, "prophyle_validate_tree.py")
    SPLIT_FA = os.path.join(C_D, "prophyle_split_allseq.py")
# executed from a Python package
else:
    PROPHYLE = "prophyle"
    PY_ASSIGN = "prophyle_assignment.py"
    ANALYZE = "prophyle_analyze.py"
    PROPAGATION_POSTPROCESSING = "prophyle_propagation_postprocessing.py"
    PROPAGATION_PREPROCESSING = "prophyle_propagation_preprocessing.py"
    NEWICK2MAKEFILE = "prophyle_propagation_makefile.py"
    READ = "prophyle_paired_end.py"
    TEST_TREE = "prophyle_validate_tree.py"
    SPLIT_FA = "prophyle_split_allseq.py"

DEFAULT_K = 31
DEFAULT_THREADS = multiprocessing.cpu_count()
# DEFAULT_THREADS=1
DEFAULT_MEASURE = 'h1'
DEFAULT_OUTPUT_FORMAT = 'sam'
DEFAULT_HOME_DIR = os.path.join(os.path.expanduser('~'), 'prophyle')

LIBRARIES = ['bacteria', 'viruses', 'plasmids', 'hmp']

ZENODO_URL = 'https://zenodo.org/record/1054426'

ANALYZE_IN_FMTS = ['sam', 'bam', 'cram', 'uncompressed_bam', 'kraken', 'histo']
ANALYZE_STATS = ['w', 'u', 'wl', 'ul']

FILES_TO_ARCHIVE = [
    ".complete.1",
    ".complete.2",
    ".complete.3",
    "tree.nw",
    "tree.preliminary.nw",
    "index.json",
    "index.fa.bwt",
    "index.fa.ann",
    "index.fa.amb",  # but will be empty
    'index.fa.kmers.tsv'
]


def _file_md5(fn, block_size=2**20):
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
    """
    tree = pro.load_nhx_tree(fn, validate=False)
    if not pro.validate_prophyle_nhx_tree(tree, verbose=True, throw_exceptions=False, output_fo=sys.stderr):
        error("The tree '{}' could not be properly parsed.".format(fn))


def _compile_prophyle_bin(clean=False, parallel=False, silent=True, force=False):
    """Compile ProPhyle binaries if they don't exist yet. Recompile if not up-to-date.

    Args:
        clean (bool): Run make clean instead of make.
        parallel (bool): Run make in parallel.
        silent (bool): Run make silently.
        force (bool): Force recompile (make -B).
    """

    try:
        command = ["make"]

        if parallel:
            command += ['-j']

        if silent:
            command += ['-s']

        if force:
            command += ['-B']

        command += ["-C", C_D]

        if clean:
            command += ['clean']

        pro.run_safe(command, output_fo=sys.stderr)

    except RuntimeError:
        if not os.path.isfile(IND) or not os.path.isfile(ASM):
            pro.error(
                "Error: ProPhyle executables could not be compiled. Please, the command '{}' manually.".format(
                    " ".join(command)
                )
            )
        else:
            print("Warning: ProPhyle executables could not be recompiled. Going to use the old ones.", file=sys.stderr)


def _add_configuration_parameter(parser, visible=True):
    parser.add_argument(
        '-c',
        dest='config',
        metavar='STR',
        nargs='*',
        type=str,
        default=[],
        help='advanced configuration (a JSON dictionary)' if visible else argparse.SUPPRESS,
    )


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


def _is_complete(d, i=1, name=None, dont_check_previous=False):
    """Check if a mark file i exists AND is newer than the mark file (i-1).

    Args:
        d (str): Directory.
        i (int): Number of the step.
        name (str): Name of the mark.
    """

    assert i > 0
    fn = __mark_fn(d, i, name)
    fn0 = __mark_fn(d, i - 1, name)

    if i == 1 or dont_check_previous:
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
            '">"', "{}", "\\;", "|", 'sed', '"s/\:>/\t/"'
        ]

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

    #pro.message("Checking library '{}' in '{}'".format(library, d))
    lib_missing = _missing_library(d)

    if library == 'bacteria':
        if lib_missing or force:
            cmd = [
                'cd', d + "/..", '&&', 'curl', '-O', ZENODO_URL + '/files/bacteria.nw', '&&', 'curl',
                ZENODO_URL + '/files/bacteria.tar.gz', '|', 'tar', 'xz'
            ]
            pro.run_safe(cmd)
            _mark_complete(d, 1)
        # _pseudo_fai(d)

    elif library == 'viruses':
        if lib_missing or force:
            cmd = [
                'cd', d + "/..", '&&', 'curl', '-O', ZENODO_URL + '/files/viruses.nw', '&&', 'curl',
                ZENODO_URL + '/files/viruses.tar.gz', '|', 'tar', 'xz'
            ]
            pro.run_safe(cmd)
            _mark_complete(d, 1)
        # _pseudo_fai(d)

    elif library == 'plasmids':
        if lib_missing or force:
            cmd = [
                'cd', d + "/..", '&&', 'curl', '-O', ZENODO_URL + '/files/plasmids.nw', '&&', 'curl',
                ZENODO_URL + '/files/plasmids.tar.gz', '|', 'tar', 'xz'
            ]
            pro.run_safe(cmd)
            _mark_complete(d, 1)
        # _pseudo_fai(d)

    elif library == 'hmp':
        if lib_missing or force:
            # fix when error appears
            cmd = [
                'cd', d, '&&', 'curl', 'http://downloads.hmpdacc.org/data/HMREFG/all_seqs.fa.bz2', '|', 'bzip2', '-d',
                '|', SPLIT_FA,
                os.path.abspath(d)
            ]
            pro.run_safe(cmd)
            _mark_complete(d, 1)
        # _pseudo_fai(d)

    else:
        pro.error('Unknown library "{}"'.format(library))


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

    config = collections.OrderedDict()
    config['prophyle-version'] = version.VERSION
    config['prophyle-revision'] = version.REVCOUNT
    config['prophyle-commit'] = version.SHORTHASH
    config['k'] = k

    pro.save_index_config(index_dir, config)

    with open(os.path.join(propagation_dir, "params.mk"), "w+") as f:
        f.write('PRG_ASM="{}"\n'.format(ASM))
        f.write("K={}\n".format(k))
        if mask_repeats:
            f.write("MASKREP=1\n")
    pro.run_safe(command)
    _log_file_md5(makefile)


def _propagate(index_dir, threads, nonprop=0):
    """Run k-mer propagation.

    Args:
        index_dir (str): Index directory.
        threads (int): Number of threads for Makefile.
        nonprop (bool): Switch propagation off.
    """
    pro.message('Running k-mer propagation')
    propagation_dir = os.path.join(index_dir, 'propagation')
    pro.test_files(os.path.join(propagation_dir, 'Makefile'), test_nonzero=True)

    if nonprop:
        nonprop_cmd_str = "NONPROP=1"
    else:
        nonprop_cmd_str = ""

    # test if input files for propagation exist
    command = ['make', '-j', threads, '-C', propagation_dir, '-n', '-s', nonprop_cmd_str, '>', '/dev/null']
    pro.run_safe(
        command,
        err_msg="Some FASTA files needed for k-mer propagation are probably missing, see the messages above.",
        thr_exc=False,
        silent=True,
    )

    # run propagation
    command = ['make', '-j', threads, '-C', propagation_dir, nonprop_cmd_str, 'V=1']
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
    command = [
        "cat", propagation_dir + "/*.count.tsv", "|", "grep", "-v", "^#", "|", "sort", "|", "uniq", ">",
        os.path.join(index_dir, "index.fa.kmers.tsv")
    ]
    pro.run_safe(
        command,
        err_msg="A file with k-mer statistics could not be created.",
        thr_exc=False,
    )


def _propagation_preprocessing(in_trees, out_tree, no_prefixes, sampling_rate, autocomplete):
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
    if autocomplete:
        command += ['-A']
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

    #pro.message('Generating packed FASTA file')
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

    #pro.message('Generating BWT')
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

    #pro.message('Generating sampled OCC array')
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

    #pro.message('Generating sampled SA')
    pro.test_files(BWA, fa_fn + ".bwt")
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

    #pro.message('Generating k-LCP array')
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


def prophyle_index(
    index_dir,
    threads,
    k,
    trees_fn,
    library_dir,
    construct_klcp,
    force,
    no_prefixes,
    mask_repeats,
    keep_tmp_files,
    sampling_rate,
    autocomplete,
    nonprop,
):
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
        autocomplete (bool): Autocomplete names of internal nodes and fasta paths.
        nonprop (bool): Switch propagation off.
    """

    assert isinstance(k, int)
    assert isinstance(threads, int)
    assert k > 1
    assert threads > 0
    assert sampling_rate is None or 0.0 <= float(sampling_rate) <= 1.0

    _compile_prophyle_bin(parallel=True)

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

    #if not _is_complete(index_dir, 1) or not pro.existing_and_newer_list(trees_fn, index_tree_1):
    if not _is_complete(index_dir, 1):
        recompute = True

    if recompute:
        pro.message('[1/6] Copying/merging trees', upper=True)
        for tree_fn in trees_fn:
            tree_fn, _, root = tree_fn.partition("@")
            tree = pro.load_nhx_tree(tree_fn, validate=False)
            # postpone for autocomplete
            if not autocomplete:
                pro.validate_prophyle_nhx_tree(tree)
            if root != "":
                if len(tree.search_nodes(name=root)) == 0:
                    pro.error("Node '{}' does not exist in '{}'.".format(root, tree_fn))
        if len(trees_fn) != 1:
            pro.message('Merging {} trees'.format(len(trees_fn)))
        _propagation_preprocessing(
            trees_fn, index_tree_1, no_prefixes=no_prefixes, sampling_rate=sampling_rate, autocomplete=autocomplete
        )
        _test_tree(index_tree_1)
        _mark_complete(index_dir, 1)
    else:
        pro.message('[1/6] Tree already exists, skipping its creation', upper=True)

    #
    # 2) Create and run Makefile for propagation, and merge FASTA files
    #

    if not _is_complete(index_dir, 2):
        recompute = True

    if recompute:
        pro.message('[2/6] Running k-mer propagation', upper=True)
        _create_makefile(index_dir, k, library_dir, mask_repeats=mask_repeats)
        _propagate(index_dir, threads=threads, nonprop=nonprop)
        _propagation_postprocessing(index_dir, index_tree_1, index_tree_2)
        _test_tree(index_tree_2)
        _kmer_stats(index_dir)
        if not keep_tmp_files:
            _remove_tmp_propagation_files(index_dir)
        else:
            pro.message('Keeping temporary files')
        _mark_complete(index_dir, 2)
    else:
        pro.message('[2/6] K-mers have already been propagated, skipping propagation', upper=True)

    #
    # 3) BWT
    #

    if not _is_complete(index_dir, 3) and not _is_complete(index_dir, 4, dont_check_previous=True):
        recompute = True

    if recompute:
        pro.message('[3/6] Constructing BWT', upper=True)
        pro.rm(index_fa + '.bwt', index_fa + '.bwt.complete')
        _fa2pac(index_fa)
        _pac2bwt(index_fa)
        _mark_complete(index_dir, 3)
    else:
        pro.message('[3/6] BWT already exists, skipping its construction', upper=True)

    #
    # 3) OCC
    #

    if not _is_complete(index_dir, 4):
        recompute = True

    if recompute:
        pro.message('[4/6] Constructing OCC', upper=True)
        _bwt2bwtocc(index_fa)
        _mark_complete(index_dir, 4)
    else:
        pro.message('[4/6] OCC already exists, skipping their construction', upper=True)

    #
    # 4) SA + 5) KLCP (compute SA + KLCP in parallel)
    #

    klcp_fn = "{}.{}.klcp".format(index_fa, k)

    if construct_klcp:

        if not _is_complete(index_dir, 5):
            # SA not computed yet => compute it in parallel with KLCP
            recompute = True

        if recompute:
            pro.message('[5/6],[6/6] Constructing SA + KLCP in parallel ', upper=True)
            _bwtocc2sa_klcp(index_fa, k)
            _mark_complete(index_dir, 5)
            _mark_complete(index_dir, 6)
            return

    #
    # 5) SA (compute only SA)
    #

    if not _is_complete(index_dir, 5):
        recompute = True

    if recompute:
        pro.message('[5/6] Constructing SA', upper=True)
        _bwtocc2sa(index_fa)
    else:
        pro.message('[5/6] SA already exists, skipping its construction', upper=True)

    #
    # 6) KLCP (compute only KLCP)
    #

    if construct_klcp:
        if not _is_complete(index_dir, 6):
            recompute = True

        if recompute:
            pro.message('[6/6] Constructing k-LCP', upper=True)
            _bwtocc2klcp(index_fa, k)
            _mark_complete(index_dir, 6)
        else:
            pro.message('[6/6] k-LCP already exists, skipping its construction', upper=True)


#####################
# PROPHYLE CLASSIFY #
#####################


def prophyle_classify(
    index_dir, fq_fn, fq_pe_fn, k, out_format, mimic_kraken, measure, annotate, tie_lca, kmer_lca, print_seq, cimpl,
    force_restarted_search, prophyle_conf_string
):
    """Run ProPhyle classification.

    Args:
        index_dir (str): Index directory.
        fq_fn (str): Input reads (single-end or first of paired-end).
        fq_pe_fn (str): Input reads (second paired-end, None if single-end)
        k (int): K-mer size (None => detect automatically).
        out_format (str): Output format: sam / kraken.
        mimic_kraken (bool): Mimic Kraken algorithm (compute LCA for each k-mer).
        measure (str): Measure used for classification (h1 / h2 / c1 / c2).
        annotate (bool): Annotate assignments (insert annotations from Newick to SAM).
        tie_lca (bool): If multiple equally good assignments found, compute their LCA.
        kmer_lca (bool): Replace k-mer matches by their LCA.
        print_seq (bool): Print sequencing in SAM.
        cimpl (bool): Use the C++ implementation.
        force_restarted_search (bool): Force restarted search.
        prophyle_conf_string (str): ProPhyle configuration string.
    """

    _compile_prophyle_bin(parallel=True)
    index_fa = os.path.join(index_dir, 'index.fa')
    index_tree = os.path.join(index_dir, 'tree.nw')

    if k is None:
        k = pro.detect_k_from_index(index_dir)
        pro.message("Automatic detection of k-mer length: k={}".format(k))

    _test_tree(index_tree)

    if fq_pe_fn:
        pro.test_files(fq_fn, fq_pe_fn, allow_pipes=False)
    elif fq_fn != '-':
        pro.test_files(fq_fn, allow_pipes=False)

    pro.test_files(IND)

    pro.test_files(
        index_fa + '.bwt',
        #index_fa + '.pac',
        index_fa + '.sa',
        index_fa + '.ann',
        #index_fa + '.amb',
    )

    (bwt_s, sa_s) = pro.file_sizes(index_fa + '.bwt', index_fa + '.sa')
    if not abs(bwt_s - 2 * sa_s) < 1000:
        pro.error('Inconsistent index (SA vs. BWT)')
    #assert abs(bwt_s - 2 * pac_s) < 1000, 'Inconsistent index (PAC vs. BWT)'

    klcp_fn = "{}.{}.klcp".format(index_fa, k)
    if force_restarted_search:
        pro.message("Restarted search forced")
        use_rolling_window = False
    else:
        use_rolling_window = os.path.isfile(klcp_fn)
        if use_rolling_window:
            pro.message("k-LCP file found, going to use rolling window")
            pro.test_files(klcp_fn)
            (klcp_s, ) = pro.file_sizes(klcp_fn)
            if not abs(bwt_s - 4 * klcp_s) < 1000:
                pro.error('Inconsistent index (KLCP vs. BWT)')
        else:
            pro.message("k-LCP file not found, going to use restarted search")

    if cimpl:
        ASSIGN = C_ASSIGN
    else:
        ASSIGN = PY_ASSIGN

    if mimic_kraken:
        measure = "h1"
        tie_lca = True
        kmer_lca = True
        out_format = "kraken"

    cmd_assign = [ASSIGN]

    if not cimpl and prophyle_conf_string:
        cmd_assign += ['-c', prophyle_conf_string]

    cmd_assign += ['-m', measure, '-f', out_format]

    if annotate:
        cmd_assign += ['-A']

    if tie_lca:
        cmd_assign += ['-L']

    if kmer_lca:
        cmd_assign += ['-X']

    cmd_assign += [index_tree, k, '-']

    if fq_pe_fn:
        cmd_read = [READ, fq_fn, fq_pe_fn, '|']
        in_read = '-'
    else:
        cmd_read = []
        # fq_fn can be '-' as well
        in_read = fq_fn

    cmd_query = [
        IND, 'query', '-k', k, '-u' if use_rolling_window else '', '-b' if print_seq else '', index_fa, in_read, '|'
    ]

    command = cmd_read + cmd_query + cmd_assign
    pro.run_safe(command)


####################
# PROPHYLE ANALYZE #
####################


def prophyle_analyze(index_dir, out_prefix, input_fns, stats, in_format):

    cmd_analyze = [ANALYZE, '-s', stats, index_dir, out_prefix] + input_fns

    if in_format is not None:
        cmd_analyze += ['-f', in_format]

    pro.test_files(*filter(lambda x: x != "-", input_fns), test_nonzero=True)

    pro.run_safe(cmd_analyze)


######################
# PROPHYLE FOOTPRINT #
######################


def prophyle_footprint(index_dir):
    bwt_size = pro.file_sizes(os.path.join(index_dir, "index.fa.bwt"))[0]
    index_size = 2 * bwt_size
    print(pro.sizeof_fmt(index_size))


#####################
# PROPHYLE COMPRESS #
#####################


def prophyle_compress(index_dir, archive):
    _compile_prophyle_bin(parallel=True)
    tmp_dir = tempfile.mkdtemp()
    arcdir = index_dir.rstrip("/").split("/")[-1]
    tmp_arc_dir = os.path.join(tmp_dir, arcdir)

    # todo: should create a correct directory

    pro.message("Creating a temporary directory for files to compress")
    pro.makedirs(tmp_arc_dir)

    for x in FILES_TO_ARCHIVE:
        if x == "index.fa.bwt":
            continue
        pro.cp_to_dir(os.path.join(index_dir, x), tmp_arc_dir)

    bwt_fn_1 = os.path.join(index_dir, "index.fa.bwt")
    bwt_fn_2 = os.path.join(tmp_arc_dir, "index.fa.bwt")
    cmd = [IND, "debwtupdate", bwt_fn_1, bwt_fn_2]
    pro.run_safe(cmd)

    pro.message("Creating '{}'".format(archive))
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(tmp_arc_dir, arcname=arcdir)
    pro.message("File '{}' has been created".format(archive))


#######################
# PROPHYLE DECOMPRESS #
#######################


def prophyle_decompress(archive, output_dir, klcp):
    pro.test_files(archive)

    if not os.path.isdir(output_dir):
        pro.error("Directory '{}' does not exist.".format(output_dir))

    with tarfile.open(archive) as tar:
        names = tar.getnames()
        index_name = names[0]
        for x in FILES_TO_ARCHIVE:
            if not os.path.join(index_name, x) in names:
                pro.error("File '{}' is missing in the archive".format(x))

    index_dir = os.path.join(output_dir, index_name)

    index_exists = True
    for i in range(1, 7):
        fn = os.path.join(index_dir, ".complete.{}".format(i))
        if not os.path.isfile(fn):
            index_exists = False
            break
    if index_exists:
        pro.message("Index already exists")
        return

    _compile_prophyle_bin(parallel=True)

    pro.message("Decompressing core index files")
    cmd = ["tar", "xvf", archive, "-C", output_dir]
    pro.run_safe(cmd)
    fn = os.path.join(index_dir, ".complete.4")
    pro.rm(fn)

    pro.message("Reconstructing the index")
    pro.touch(os.path.join(index_dir, "index.fa"))
    pro.touch(os.path.join(index_dir, "index.fa.pac"))
    if klcp:
        config = pro.load_index_config(index_dir)
        cmd = [PROPHYLE, "index", "-k", config['k'], os.path.join(index_dir, "tree.nw"), index_dir]
    else:
        cmd = [PROPHYLE, "index", "-K", os.path.join(index_dir, "tree.nw"), index_dir]

    pro.run_safe(cmd)
    pro.message("Index reconstruction finished")


####################
# PROPHYLE COMPILE #
####################


def prophyle_compile(clean, parallel, force):
    _compile_prophyle_bin(clean=clean, parallel=parallel, force=force, silent=False)


########
# MAIN #
########


def parser():
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            if len(sys.argv) == 2:
                self.print_help()
                sys.exit(2)
            else:
                pro.error(message, 2)

    desc = """\
        Program: prophyle (phylogeny-based metagenomic classification)
        Version: {V}
        Authors: Karel Brinda, Kamil Salikhov, Simone Pignotti, Gregory Kucherov
        Contact: kbrinda@hsph.harvard.edu

        Usage:   prophyle <command> [options]
        """.format(V=version.VERSION)
    parser = MyParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(desc))

    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='ProPhyle {} (rev {}, commit {})'.format(version.VERSION, version.REVCOUNT, version.SHORTHASH),
    )

    _add_configuration_parameter(parser, visible=False)

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

    _add_configuration_parameter(parser_download)

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
        help='skip k-LCP construction (then restarted search only)',
    )

    parser_index.add_argument(
        '-T',
        dest='keep_tmp_files',
        action='store_true',
        help='keep temporary files from k-mer propagation',
    )

    parser_index.add_argument(
        '-A',
        help='autocomplete tree (names of internal nodes and FASTA paths)',
        dest='autocomplete',
        action='store_true',
    )

    parser_index.add_argument(
        '-R',
        help='switch propagation off (only re-assemble leaves)',
        dest='nonprop',
        action='store_true',
    )

    _add_configuration_parameter(parser_index)

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
        help='first file with reads in FASTA/FASTQ (- for standard input)',
    )

    parser_classify.add_argument(
        'reads_pe',
        metavar='<reads2.fq>',
        type=str,
        help='second file with reads in FASTA/FASTQ',
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
        '-m',
        dest='measure',
        choices=['h1', 'c1', 'h2', 'c2'],
        help='measure: h1=hit count, c1=coverage, h2=norm.hit count, c2=norm.coverage [{}]'.format(DEFAULT_MEASURE),
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
        '-P',
        dest='print_seq',
        action='store_true',
        help='incorporate sequences and qualities into SAM records',
    )

    parser_classify.add_argument(
        '-A',
        dest='annotate',
        action='store_true',
        help='annotate assignments (using tax. information from NHX)',
    )

    parser_classify.add_argument(
        '-L',
        dest='tie_lca',
        action='store_true',
        help='replace read assignments by their LCA',
    )

    parser_classify.add_argument(
        '-X',
        dest='kmer_lca',
        action='store_true',
        help='replace k-mer matches by their LCA',
    )

    parser_classify.add_argument(
        '-M',
        dest='mimic',
        action='store_true',
        help='mimic Kraken (equivalent to  "-m h1 -f kraken -L -X")',
    )

    parser_classify.add_argument(
        '-C',
        dest='cimpl',
        action='store_true',
        help='use C++ impl. of the assignment algorithm (experimental)',
        #help=argparse.SUPPRESS,
    )

    parser_classify.add_argument(
        '-K',
        dest='force_restarted_search',
        action='store_true',
        help='force restarted search mode',
    )

    _add_configuration_parameter(parser_classify)

    ##########

    parser_analyze = subparsers.add_parser(
        'analyze',
        help='analyze results (experimental)',
        formatter_class=fc,
    )

    parser_analyze.add_argument(
        'index_dir', metavar='{index_dir, tree.nw}', type=str, help='index directory or phylogenetic tree'
    )

    parser_analyze.add_argument(
        'out_prefix',
        metavar='<out.pref>',
        type=str,
        help="output prefix",
    )

    parser_analyze.add_argument(
        'input_fns',
        metavar='<classified.bam>',
        type=str,
        nargs='+',
        default=None,
        help="classified reads (use '-' for stdin)",
    )

    parser_analyze.add_argument(
        '-s', metavar=ANALYZE_STATS, type=str, dest='stats', choices=ANALYZE_STATS, default=ANALYZE_STATS[0],
        help="""statistics to use for the computation of histograms:
                w (default) => weighted assignments;
                u => unique assignments, non-weighted;
                wl => weighted assignments, propagated to leaves;
                ul => unique assignments, propagated to leaves."""
    )

    parser_analyze.add_argument(
        '-f', metavar=ANALYZE_IN_FMTS, type=str, dest='in_format', choices=ANALYZE_IN_FMTS, default=None,
        help="""Input format of assignments [auto]"""
    )

    _add_configuration_parameter(parser_analyze)

    ##########

    parser_footprint = subparsers.add_parser(
        'footprint',
        help='estimate memory footprint',
        formatter_class=fc,
    )

    parser_footprint.add_argument(
        'index_dir',
        metavar='<index.dir>',
        type=str,
        help='index directory',
    )

    _add_configuration_parameter(parser_footprint)

    ##########

    parser_compress = subparsers.add_parser(
        'compress',
        help='compress a ProPhyle index',
        formatter_class=fc,
    )

    parser_compress.add_argument(
        'index_dir',
        metavar='<index.dir>',
        type=str,
        help='index directory',
    )

    parser_compress.add_argument(
        'archive',
        metavar='<archive.tar.gz>',
        type=str,
        default=None,
        nargs="?",
        help='output archive [<index.dir>.tar.gz]',
    )

    _add_configuration_parameter(parser_compress)

    ##########

    parser_decompress = subparsers.add_parser(
        'decompress',
        help='decompress a compressed ProPhyle index',
        formatter_class=fc,
    )

    parser_decompress.add_argument(
        'archive',
        metavar='<archive.tar.gz>',
        type=str,
        help='output archive',
    )

    parser_decompress.add_argument(
        'output_dir',
        metavar='<output.dir>',
        type=str,
        nargs="?",
        default="./",
        help='output directory [./]',
    )

    parser_decompress.add_argument(
        '-K',
        dest='klcp',
        action='store_false',
        help='skip k-LCP construction',
    )

    _add_configuration_parameter(parser_decompress)

    ##########

    parser_compile = subparsers.add_parser(
        'compile',
        help='compile auxiliary ProPhyle programs',
        formatter_class=fc,
    )

    parser_compile.add_argument(
        '-C',
        dest='clean',
        action='store_true',
        help='clean files instead of compiling',
    )

    parser_compile.add_argument(
        '-F',
        dest='force',
        action='store_true',
        help='force recompilation',
    )

    parser_compile.add_argument(
        '-P',
        dest='parallel',
        action='store_true',
        help='run compilation in parallel',
    )

    _add_configuration_parameter(parser_compile)

    ##########

    return parser


def main():
    try:
        par = parser()
        args = par.parse_args()
        subcommand = args.subcommand

        global CONFIG
        prophyle_conf_string = pro.load_prophyle_conf(CONFIG, args.config)

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
                autocomplete=args.autocomplete,
                nonprop=args.nonprop,
            )
            pro.message('Index construction finished')
            pro.close_log()

        elif subcommand == "classify":
            # if args.log_fn is None:
            #    args.log_fn = os.path.join(args.index_dir, "log.txt")

            pro.open_log(args.log_fn)
            pro.message('Classification started')
            prophyle_classify(
                index_dir=args.index_dir,
                fq_fn=args.reads,
                fq_pe_fn=args.reads_pe,
                k=args.k,
                out_format=args.oform,
                mimic_kraken=args.mimic,
                measure=args.measure,
                tie_lca=args.tie_lca,
                kmer_lca=args.kmer_lca,
                annotate=args.annotate,
                print_seq=args.print_seq,
                cimpl=args.cimpl,
                force_restarted_search=args.force_restarted_search,
                prophyle_conf_string=prophyle_conf_string,  # already preprocessed
            )
            pro.message('Classification finished')
            pro.close_log()

        elif subcommand == "analyze":

            prophyle_analyze(
                index_dir=args.index_dir,
                out_prefix=args.out_prefix,
                input_fns=args.input_fns,
                stats=args.stats,
                in_format=args.in_format,
            )

        elif subcommand == "footprint":

            prophyle_footprint(index_dir=args.index_dir, )

        elif subcommand == "compress":

            if args.archive is None:
                archive = args.index_dir.rstrip("/") + ".tar.gz"
            else:
                archive = args.archive

            prophyle_compress(
                index_dir=args.index_dir,
                archive=archive,
            )

        elif subcommand == "decompress":

            prophyle_decompress(
                archive=args.archive,
                output_dir=args.output_dir,
                klcp=args.klcp,
            )

        elif subcommand == "compile":

            prophyle_compile(
                clean=args.clean,
                parallel=args.parallel,
                force=args.force,
            )

        else:
            msg_lns = par.format_help().split("\n")[2:]
            msg_lns = [x for x in msg_lns if x.find("optional arguments") == -1 and x.find("--") == -1]
            msg = "\n".join(msg_lns)
            msg = msg.replace("\n\n", '\n').replace("subcommands:\n", "Command:\n").replace("Usage", "\nUsage")
            msg = msg.replace("\n    compress", "\n\n    compress")
            print(file=sys.stderr)
            print(msg, file=sys.stderr)
            sys.exit(2)

    except BrokenPipeError:
        # pipe error (e.g., when head is used)
        sys.stderr.close()
        sys.stdout.close()
        exit(0)

    except KeyboardInterrupt:
        pro.error("Keyboard interrupt")

    finally:
        sys.stdout.flush()
        sys.stderr.flush()


if __name__ == "__main__":
    main()
