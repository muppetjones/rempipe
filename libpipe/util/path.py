'''
Created on May 16, 2014

@author: sbush
@contact: stephen.james.bush at gmail
'''

import os.path
import re
import sys

import logging
log = logging.getLogger(__name__)
# startup logger


def isValidPath(path_str):
    '''Test for invalid characters in path name'''
    invalid_char_str = r'[\?\%\*\:\|\"\'\>\<]'
    if re.search(invalid_char_str, path_str):
        return False
    return True


def protect(path_str, abspath=True):
    '''Protect a given path string

    Checks the following for a given path:
        * Ensures no invalid characters (?%*:|"'><).
        * Removes trailing slash (if directory)
        * Expands any environment variables (and user)
        * Normalizes path components
        * Converts path component to absolute path
            -- good for ensuring file correctness
            -- bad for cross-system functionality

    Arguments:
        path_str

    Returns:
        A cleaned, normalized absolute path

    Raise:
        AssertionError: If contains invalid characters
    '''

    # check for valid path
    assert isValidPath(path_str)

    # remove trailing slash (if there)
    val_path = path_str.rstrip(os.path.sep)

    # expand variables in path (if any)
    val_path = os.path.expanduser(val_path)
    val_path = os.path.expandvars(val_path)

    # ensure absolute path
    # -- don't do it if no directory is included
    # otherwise it will add path to curr directory
    if abspath and os.path.dirname(val_path):
        # val_path = os.path.normpath(val_path)  # not needed--handled by abs
        val_path = os.path.abspath(val_path)

    return val_path


def common_directory(path_list):
    '''Return string for common directory for a given list of paths

    NOTE: os.path.dirname will clip a trailing slash and return the rest of the
          path, meaning `dirname` works fine.
    '''

    common = os.path.commonprefix(path_list)
    return os.path.dirname(common)


def unique_suffix(file_list, strict=True, suffix_terminator='_'):
    '''Return a list of path strings minus the common prefix (and extension)

    Arguments:
        file_list: Files with a common path
        strict (bool, default=True): If true, will eliminate all common
            characters. If false, the suffix will contain common characters
            back to the suffix terminator.
        suffix_terminator (char, default='_'): The character that denotes
            the end of the prefix.
    Return:
        A list of unique suffixes (missing file extensions)

    Examples:
        file_list = ['~/data/sample/S1.txt',
                     '~/data/sample/S2.txt',
                     '~/data/sample/S3.txt',
                     '~/data/sample/S4.txt',
                    ]
        suffixes_strict = unique_suffix(file_list)
        suffixes_relax = unique_suffix(
            file_list, strict=False,
            suffix_terminator='/')

        print(", ".join(suffixes_strict))
        1, 2, 3, 4

        print(", ".join(suffixes_relax))
        S1, S2, S3, S4
    '''

    common_prefix = os.path.commonprefix(file_list)
    idx_slash = common_prefix.rfind(os.sep)
    idx_uscr = common_prefix.rfind(suffix_terminator)
    if idx_slash <= idx_uscr:
        common_prefix = common_prefix[:idx_uscr + 1]
    len_prefix = len(common_prefix)

    return [os.path.splitext(f[len_prefix:])[0] for f in file_list]


def avoid_overwrite(filename, max_attempts=20, allow_empty=True):
    '''Return an overwrite-safe filename

    NOTE: If max number of attempts is hit, will return filename of
          oldest file.

    Arguments:
        filename (string)
        max_attempts (int): Number of filenames to try before failing
        allow_empty (bool): If True, will overwrite empty files
    Return:
        A file with os.path.exist == False
    '''

    protected_name = protect(filename)
    if not os.path.exists(protected_name):
        return protected_name

    check_names = ["{0}.{1}".format(protected_name, i)
                   for i in range(max_attempts)]
    for check in check_names:

        if not os.path.exists(check):
            return check

        elif allow_empty and os.path.isfile(check):
            try:
                file_size = os.path.getsize(check)
                if file_size == 0:
                    return check
            except OSError:
                pass
        else:
            pass

    # still here, so we didn't find anything
    # just use the oldest file
    age_sorted = sorted(check_names, key=os.path.getctime)
    return age_sorted[0]


def walk_safe_generator(dir_str, level=None):
    '''walk over the given path...return a dict with 'dir' and 'file' attr

    Note: This function used to take a function for filtering
        files as an argument. This feature was removed and the
        responsibility pushed to the caller. The reasoning is simple:
        Users most likely want to match using the full path, especially
        if looking for directories. Futhermore, it didn't make sense to
        call the same filter function multiple times (once for directories
        and once for files per iteration) if the user was only interested
        in either files or directories.

    Arguments:
        dir_str (path string): Path to search
        level (int): The maximum depth to descend; 0 is current directory only.
            Default: None (no maximum).

    Returns:
        A dict containing a list of directories (dir) and files (file).

    Raises:
        AssertionError: Invalid path (bad name or does not exist)
    '''

    # ------------------------ Ensure we were given a valid directory to search
    try:
        dir_str = protect(dir_str)
    except AssertionError:
        raise

    assert os.path.isdir(dir_str), "Not a valid dir: %s" % dir_str

    # ---------------------------------------------------- prep for depth check
    num_sep = dir_str.count(os.path.sep)

    for root, d, f in os.walk(dir_str, topdown=True):

        # ----------------------------------------------- check depth recursion
        if level is not None:
            num_sep_this = root.count(os.path.sep)
            if num_sep + level <= num_sep_this:
                del d[:]
                continue

        yield root, d, f  # makes this a generator function!

    return


def walk_safe(dir_str, level=None):
    '''walk over the given path...return a dict with 'dir' and 'file' attr

    Note: Do not use on large directories. You will regret it. Instead,
        use 'walk_safe_generator' directory.

    Arguments:
        dir_str (path string): Path to search
        level (int): The maximum depth to descend. Default: None

    Returns:
        A dict containing a list of directories (dir) and files (file).

    Raises:
        AssertionError: Invalid path (bad name or does not exist)

    See Also:
        walk_safe_generator
        walk_file
        walk_dir
    '''

    contents = {'dir': [], 'file': []}
    for root, d, f in walk_safe_generator(dir_str, level=level):

        d_path = [os.path.join(root, x) for x in d]
        f_path = [os.path.join(root, x) for x in f]

        contents['dir'].extend(d_path)
        contents['file'].extend(f_path)

    return contents


def walk_file(directory, level=None, extension=[], pattern=[]):
    '''Return a list of files in a given directory

    Note: 'patterns' uses regular expressions for pattern matching.
        Patterns are checked against the full path.

    Arguments:
        directory (string): The path to search
        level (int): The depth to recurse. Default: None (no limit).
        extension (string or list): Desired file extensions.
        pattern (string or list): Regex pattern string to match files against.

    Returns:
        A list of matching files found in the given directory

    Raises:
        None

    Examples:
        # Find all text files in a given directory
        files = walk_file('~/my_dir', extensions='txt')

        # Find all bowtie2 index files for a particular genome
        files = walk_file('~/genomes', patterns='GRCh37_p13.+bt2.*')
    '''

    files = walk_safe(directory, level=level)['file']

    # ---------------------------- create lambda functions for validating files
    if not extension:
        extension_match = lambda f: True
    else:
        extensions = (extension,) if isinstance(extension, str) \
            else tuple(extension)
        extension_match = lambda f: f.endswith(tuple(extensions))

    if not pattern:
        pattern_match = lambda f: True
    else:
        # combine all patterns into a single regular expression
        patterns = (pattern,) if isinstance(pattern, str) \
            else tuple(pattern)
        pattern_str = r'(?:' + r'|'.join(tuple(patterns)) + r')'
        pattern_re = re.compile(pattern_str, re.I)
        pattern_match = lambda f: True if pattern_re.search(f) else False

    # ----------------------------------------------------- find matching files
    matched_files = [f for f in files
                     if extension_match(f) and pattern_match(f)]

    return matched_files


def walk_dir(directory, level=None, pattern=[]):
    '''Return a list of directories in a given directory

    Note: 'patterns' uses regular expressions for pattern matching.
        Patterns are checked against the full path.

    Arguments:
        directory (string): The path to search
        level (int): The depth to recurse. Default: None (no limit).
        pattern (string or list): Regex pattern to match dirs against.

    Returns:
        A list of matching dirs found in the given directory

    Raises:
        None

    See Also:
        walk_safe
        walk_file
    '''

    dirs = walk_safe(directory, level=level)['dir']

    # ----------------------------- create lambda functions for validating dirs
    if not pattern:
        pattern_match = lambda f: True
    else:
        # combine all patterns into a single regular expression
        patterns = (pattern,) if isinstance(pattern, str) \
            else tuple(pattern)
        pattern_str = r'(?:' + r'|'.join(tuple(patterns)) + r')'
        pattern_re = re.compile(pattern_str, re.I)
        pattern_match = lambda f: True if pattern_re.search(f) else False

    # ------------------------------------------------------ find matching dirs
    matched_dirs = [d for d in dirs
                    if pattern_match(d)]

    return matched_dirs


def makedirs(dir_str, mode_oct=0o755):
    '''
    Safely and gracefully make a directory
    Raises OSError if bad path string or if error creating directory
    '''

    # ensure kosher path
    try:
        dir_str = protect(dir_str)
    except AssertionError as e:
        raise OSError(str(e))

    # attempt to make dir
    # raise unless directory already exists
    try:
        os.makedirs(dir_str, mode=mode_oct)
    except OSError:
        if not os.path.isdir(dir_str):
            raise

    return True
