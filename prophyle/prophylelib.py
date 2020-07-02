#! /usr/bin/env python3
"""Auxiliary ProPhyle functions for running shell commands, file manipulation and tree manipulation.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import datetime
import ete3
import json
import os
import pathlib
import psutil
import shutil
import subprocess
import sys
import tempfile
import time
import gzip

###########
# LOGGING #
###########

log_file = None


def open_gzip(fn):
    """Open a file, possibly compressed with gzip.

    Args:
        fn (str): File name.
    """
    magic = b'\x1f\x8b\x08'
    l = len(magic)
    with open(fn, 'rb') as f:
        file_start = f.read(l)
        f.seek(0)
    # check if the file is compressed
    if file_start.startswith(magic):
        return gzip.open(fn, 'rt')
    # not compressed
    return open(fn, 'rt')


def open_log(fn):
    """Open a log file.

    Args:
        fn (str): File name.
    """

    global log_file
    if fn is not None:
        d = os.path.dirname(fn)
        if d != "":
            makedirs(d)
        log_file = open(fn, "a+")


def close_log():
    """Close a log file.
    """

    global log_file
    if log_file is not None:
        try:
            log_file.flush()
        finally:
            log_file.close()


def error(*msg, error_code=1):
    print('ProPhyle Error:', *msg, file=sys.stderr)
    sys.stdout.flush()
    sys.stderr.flush()
    close_log()
    sys.exit(error_code)


def message(*msg, subprogram='', upper=False, only_log=False):
    """Print a ProPhyle message to stderr.

    Args:
        *msg: Message.
        subprogram (str): Subprogram.
        upper (bool): Transform text to upper cases.
        only_log (bool): Don't print the message to screen (i.e., it would be only logged).
    """

    global log_file

    dt = datetime.datetime.now()
    fdt = dt.strftime("%Y-%m-%d %H:%M:%S")

    if upper:
        msg = map(str, msg)
        msg = map(str.upper, msg)

    log_line = '[prophyle{}] {} {}'.format(subprogram, fdt, " ".join(msg))

    if not only_log:
        print(log_line, file=sys.stderr)
    if log_file is not None:
        log_file.write(log_line)
        log_file.write("\n")
        log_file.flush()


###################
# TREE OPERATIONS #
###################


def load_nhx_tree(nhx_fn, validate=True):
    """Load a ProPhyle NHX tree.

    Args:
        nhx_fn (str): Newick/NHX tree.
        validate (bool): Validate the tree.
    """

    tree = ete3.Tree(nhx_fn, format=1)

    if validate:
        validate_prophyle_nhx_tree(tree)

    return tree


def save_nhx_tree(tree, nhx_fn):
    """Save a ProPhyle NHX tree.

    Args:
        tree (ete3.Tree): Ete3 tree.
        nhx_fn (str): Name of the file.
    """

    assert isinstance(tree, ete3.Tree)

    # make saving newick reproducible
    features = set()
    for n in tree.traverse():
        features |= n.features

    # otherwise some names stored twice – also as a special attribute
    features.remove("name")

    tree.write(
        features=sorted(features),
        format=1,
        format_root_node=True,
        outfile=nhx_fn,
    )


def validate_prophyle_nhx_tree(tree, verbose=True, throw_exceptions=True, output_fo=sys.stderr):
    """Validate an ETE3 tree with respect to ProPhyle requirements.

    Args:
        tree (ete3.Tree): Ete3 tree.
        verbose (bool): Verbose mode.
        throw_exceptions (bool): Throw an exception if the tree is not valid.
        output_fo (file): Output file object.
    """

    assert isinstance(tree, ete3.Tree), tree

    error = False

    existing_names_set = set()

    names_with_separator = []

    without_name = []
    empty_name = []
    forbidden_name = []

    duplicates = []

    for i, node in enumerate(tree.traverse("postorder")):
        noname = True
        try:
            nname = node.name
            noname = False
        except AttributeError:
            without_name.append((i, None))

        if not noname:
            if nname == '':
                empty_name.append((i, nname))
                error = True
            if nname in set(['A', '0']):
                forbidden_name.append((i, nname))
                error = True
            if nname in existing_names_set:
                duplicates.append((i, nname))
                error = True
            if "@" in nname:
                names_with_separator.append((i, nname))
                error = True
            existing_names_set.add(nname)

    def _format_node_list(node_list):
        return ", ".join(map(str, node_list))

    def _error_report(node_list, message):
        if len(node_list) > 0:
            print(
                "   * {} node(s) {}: {}".format(
                    len(node_list),
                    message,
                    _format_node_list(node_list),
                ), file=output_fo
            )

    if verbose:
        if error:
            print("Errors:".format(), file=output_fo)

        _error_report(without_name, "without name")
        _error_report(forbidden_name, "forbidden_name")
        _error_report(empty_name, "with empty name")
        _error_report(duplicates, "with a duplicate name")
        _error_report(names_with_separator, "with a name containing '@'")

    if throw_exceptions:
        if error:
            error("Invalid tree. The format of the tree is not correct. See the messages above.")

    if error:
        return False
    else:
        return True


def minimal_subtree(tree):
    """Take a tree and compute its minimal subtree. All singletons are removed and the highest branching node is
    used as the new root.

    Args:
        tree (ete3.Tree): Phylogenetic tree.
    """
    tree_copy = tree.copy()

    for n in tree_copy.traverse():
        if len(n.children) == 1:
            n.delete()

    new_root = tree_copy
    while len(new_root.children) == 1:
        new_root = new_root.children[0]

    new_tree = new_root.detach()
    return new_tree


def has_attribute(tree, attribute):
    for n in tree.traverse():
        if hasattr(n, attribute):
            return True
    return False


#####################
# FILE MANIPULATION #
#####################


def file_sizes(*fns):
    """Get file sizes in Bytes.

    Args:
        fns (str): File names.

    Returns:
        tuple(int): File sizes.
    """
    return tuple([os.stat(fn).st_size for fn in fns])


def sizeof_fmt(bs, suffix='B'):
    """Convert bytes to a human readable string.

    Based on https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size

    Args:
        bs (int): File size in bytes.

    Args:
        str: Human readable string.
    """

    num = bs
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def touch(*fns):
    """Touch files.

    Args:
        *fns: Files.
    """
    for fn in fns:
        if os.path.exists(fn):
            os.utime(fn, None)
        else:
            with open(fn, 'a'):
                pass


def rm(*fns):
    """Remove files (might not exists).

    Args:
        *fns: Files.
    """
    for fn in fns:
        try:
            os.remove(fn)
        except FileNotFoundError:
            pass


def cp_to_file(fn0, fn):
    """Copy file to file.

    Args:
        fn0 (str): Source file.
        fn (str): Target file.
    """

    # keep rewriting attributes
    shutil.copyfile(fn0, fn)


def cp_to_dir(fn0, d):
    """Copy file to dir.

    Args:
        fn0 (str): Source file.
        d (str): Target dir.
    """

    # keep rewriting attributes
    shutil.copy(fn0, d)


def makedirs(*ds):
    """Make dirs recursively (no error if existing, make parent directories as needed).

    Args:
        *ds: Dirs to create.
    """
    for d in ds:
        if not os.path.isdir(d):
            cmd = ['mkdir', '-p', d]
            run_safe(cmd, silent=True)


def existing_and_newer(fn0, fn):
    """Test if file fn exists and is newer than fn0. Exit if fn0 does not exist.

    Args:
        fn0 (str): Old file.
        fn (str): New file (to be generated from fn0).
    """

    if not os.path.isfile(fn0):
        error("Dependency '{}' does not exist".format(fn0))

    if not os.path.isfile(fn):
        return False

    if os.path.getmtime(fn0) <= os.path.getmtime(fn):
        return True
    else:
        return False


def existing_and_newer_list(fn0_l, fn):
    """Test if file fn exists and is newer than all files in fn0_l. Exit if some fn0 file does not exist.

    Args:
        fn0_l (list of str): Old files list.
        fn (str): New file (to be generated from fn0_l).
    """

    rs = [existing_and_newer(fn0, fn) for fn0 in fn0_l]
    some_false = False in rs
    return not some_false


def test_files(*fns, test_nonzero=False, allow_pipes=False):
    """Test if given files exist, and possibly if they are non-empty. If not, stop the program.

    Args:
        *fns: Files.
        test_nonzero (bool): Test if files have size greater than zero.
        allow_pipes (bool): Allow Unix pipes as input and don't test size.
    """

    for fn in fns:
        is_file = os.path.isfile(fn)
        is_pipe = pathlib.Path(fn).is_fifo()
        if allow_pipes:
            if not is_file or is_pipe:
                error('File "{}" does not exist.'.format(fn))
        else:
            if is_pipe:
                if not is_file:
                    error('File "{}" is a process substitution or a device.'.format(fn))
            else:
                if not is_file:
                    error('File "{}" does not exist.'.format(fn))

        if test_nonzero and not allow_pipes:
            if not file_sizes(fn)[0]:
                error('File "{}" has size 0.'.format(fn))


########
# MISC #
########


def run_safe(command, output_fn=None, output_fo=None, err_msg=None, thr_exc=True, silent=False):
    """Run a shell command safely.

    Args:
        command (list of str): Command to execute.
        output_fn (str): Name of a file for storing the output.
        output_fo (fileobject): Output file object. If both params are None, the standard output is used.
        err_msg (str): Error message if the command fails.
        thr_exc (bool): Through exception if the command fails. error_msg or thr_exc must be set.
        silent (bool): Silent mode (print messages only if the command fails).
    """

    assert output_fn is None or output_fo is None
    assert err_msg is not None or thr_exc
    assert len(command) > 0

    command_safe = []

    for part in command:
        part = str(part)
        if " " in part:
            part = '"{}"'.format(part)
        command_safe.append(part)

    command_str = " ".join(command_safe)

    if not silent:
        message("Shell command:", command_str)

    if output_fn is None:
        if output_fo is None:
            out_fo = sys.stdout
        else:
            out_fo = output_fo
    else:
        out_fo = open(output_fn, "w+")

    if out_fo == sys.stdout:
        p = subprocess.Popen("/bin/bash -e -o pipefail -c '{}'".format(command_str), shell=True)
    else:
        p = subprocess.Popen("/bin/bash -e -o pipefail -c '{}'".format(command_str), shell=True, stdout=out_fo)

    ps_p = psutil.Process(p.pid)

    max_rss = 0
    error_code = None
    while error_code is None:
        try:
            max_rss = max(max_rss, ps_p.memory_info().rss)
        except (psutil.NoSuchProcess, psutil.ZombieProcess, psutil.AccessDenied, OSError, IOError):
            pass
        #except psutil.NoSuchProcess as e:
        #    print("[prophylelib] Warning: psutil - NoSuchProcess (pid: {}, name: {}, msg: {})".format(e.pid, e.name, e.msg), file=sys.stderr)

        # wait 0.2 s
        time.sleep(0.2)
        error_code = p.poll()

    out_fo.flush()

    mem_mb = round(max_rss / (1024 * 1024.0), 1)

    if output_fn is not None:
        out_fo.close()

    if error_code == 0 or error_code == 141:
        if not silent:
            message("Finished ({} MB used): {}".format(mem_mb, command_str))
    else:
        message("Unfinished, an error occurred (error code {}, {} MB used): {}".format(error_code, mem_mb, command_str))

        if err_msg is not None:
            print('Error: {}'.format(err_msg), file=sys.stderr)

        if thr_exc:
            error("A command failed, see messages above.")

        sys.exit(1)


def save_index_config(index_dir, data):
    """Save a configuration dictionary (index.json in the index directory).

    Args:
        index_dir (str): Index directory.
    """

    with open(os.path.join(index_dir, 'index.json'), "w+") as data_file:
        json.dump(data, data_file, indent=4)


def load_index_config(index_dir):
    """Load a configuration dictionary (index.json in the index directory).

    Args:
        index_dir (str): Index directory.
    """
    with open(os.path.join(index_dir, 'index.json')) as data_file:
        data = json.load(data_file)
    return data


def detect_k_from_index(index_dir):
    """Detect k-mer size from a ProPhyle index.

    Args:
        index_dir (str): Index directory.
    """

    config = load_index_config(index_dir)
    return config['k']


def lower_nonsigleton(node):
    while len(node.children) == 1:
        node = node.children[0]
    return node


def load_prophyle_conf(globconf, json_strs):
    """Loads configuration from a JSON string or a file.

    The strings are first merged. If the string corresponds to a file,
    it is loaded. Otherwise it is interpreted as a json string
    (possibly without wrapping braces) and stored into a temporary
    file, whose name is then returned so that it can be passed to
    other programs.

    Params:
        glob_conf (dict): Global configuration dictionary.
        json_strs (str): File name with the dictionary.

    Returns:
        json_fn (str): Returns a file name of a JSON with the same information.
    """

    assert type(globconf) == dict

    if len(json_strs) == 0:
        return None

    string = (" ".join(json_strs)).strip()

    if os.path.isfile(string):
        # load JSON from a file & return the same filename
        json_filename = string
        with open(json_filename) as f:
            data = json.load(f)
        prophyle_conf_tmp_filename = json_filename

    else:
        # load JSON from the string & create a tmp file
        json_string = string
        if json_string[0] != "{" and json_string[-1] != "}":
            json_string = "".join(["{", json_string, "}"])
        data = json.loads(json_string)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            prophyle_conf_tmp_filename = f.name
            json.dump(data, f, indent=3)

    assert type(data) == dict, "The provided configuration is not a dictionary"
    globconf.update(data)

    try:
        if globconf['PRINT_CONF']:
            print("Configuration:", globconf, file=sys.stderr)
    except KeyError:
        pass

    return prophyle_conf_tmp_filename


## None if root
#def upper_nonsigleton(node):
#   while len(node.children)==1:
#       node = node.children[0]
#   return node

if __name__ == "__main__":
    sys.exit(1)
