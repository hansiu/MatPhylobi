import subprocess
import logging
import sys
import os
from Bio import Entrez
import datetime

try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import ast

logger = logging.getLogger('MatPhylobi')
logging.basicConfig(stream=sys.stdout, level=logging.INFO)


# -------------- Helpers --------------


def check(file, action, newname=''):
    """Checks for existence, creation etc of files. Created just for ease of use"""
    if isinstance(file, str) and isinstance(action, str):
        if not os.path.exists(file):
            if action == 'e':
                logger.critical('Non-existing path specified: ' + file)
                return False
            elif action == 'm':
                os.makedirs(file)
            elif action == 'e2':
                return False
        elif action == 'e2':
            return True
        elif action == 'r':
            os.remove(file)
        elif action == 'l':
            return [f for f in os.listdir(file)]
        elif action == 'n':
            os.rename(file, newname)
        return True


def check_blast():
    try:
        blast_command = subprocess.check_output('blastn -version', shell=True)
    except subprocess.CalledProcessError:
        logger.error("MatPhylobi requires blastn 2.5.0+ or later. "
                     "Try downloading blast+ from 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST'")
        sys.exit(1)

    blast_ver = blast_command.decode(sys.stdout.encoding)
    blast_ver = (blast_ver.split('\n')[0]).split(': ')[1].replace('\r', '').replace('+', '')
    if tuple(map(int, blast_ver.split('.'))) < (2, 5, 0):
        logger.error("MatPhylobi requires blastn 2.5.0+ or later. blastn " + str(blast_ver) + " detected.")
        sys.exit(1)


def validate_analyze_args(args, parser):
    """
    validate and prepare arguments
    """

    def check_most_important(gene_names, sequences, org_included, org_included_file, parser):
        if not (bool(gene_names) and bool(sequences) and (
                    bool(org_included) or bool(org_included_file))):
            parser.error('--gene_names, --sequences and --org_included(_file) must be provided together')
        if not len(sequences) == len(gene_names):
            parser.error('--sequences and --gene_names should be of same length')

    def check_gene_duplicates(gene_names, parser):
        from collections import Counter

        dup_gene_names = [item for item, count in Counter(gene_names).items() if count > 1]
        if len(dup_gene_names) != 0:
            parser.error('Some gene_names are duplicated in the input. Please correct that.')

    def merge_its(gene_names, sequences):
        allowed_its_names = ['ITS', 'ITS1', 'ITS2']
        its_found = [i for i in range(len(gene_names)) if gene_names[i].upper() in allowed_its_names]

        if its_found:
            itss = []

            for its_ind in its_found[::-1]:  # going from back so that we don't change index order
                itss.extend(sequences[its_ind])
                del gene_names[its_ind]
                del sequences[its_ind]

            gene_names.append('ITS')
            sequences.append(itss)

        return gene_names, sequences

    def prepare_org_lists(org_included_file, org_included, org_excluded_file, org_excluded):
        if org_included_file:
            with open(org_included_file) as orgfile:
                for line in orgfile:
                    org_included.append(line.strip('\n'))
        if org_excluded_file:
            with open(org_excluded_file) as exclfile:
                for line in exclfile:
                    org_excluded.append(line.strip('\n'))
        return list(set(args.org_included)), list(set(args.org_excluded))

    check_most_important(args.gene_names, args.sequences, args.org_included, args.org_included_file, parser)
    check_gene_duplicates(args.gene_names, parser)

    args.sequences = [seq.split(",") for seq in args.sequences]

    args.gene_names, args.sequences = merge_its(args.gene_names, args.sequences)

    args.org_included, args.org_excluded = prepare_org_lists(args.org_included_file, args.org_included,
                                                             args.org_excluded_file, args.org_excluded)

    args.blacklist = os.path.abspath(args.blacklist) if args.blacklist else ''
    args.synonyms = os.path.abspath(args.synonyms) if args.synonyms else ''

    return args


def prepare(args, parser):
    Entrez.email = 'matphylobi@gmail.com'
    Entrez.tool = 'MatPhylobi'
    if args.api_key:
        Entrez.api_key = args.api_key

    args.folder = os.path.abspath(args.folder)
    check(args.folder, 'm')
    logger.info('Dir: ' + args.folder)

    args.today = str(datetime.datetime.now())[:10].replace('-', '/')
    logger.info('Date: ' + args.today)

    if args.action == 'analyze':
        args = validate_analyze_args(args, parser)
    elif args.action == 'update':
        if not check(args.folder + '/config.ini', 'e'):
            clean_up(args.debug, args.folder, args.action, 1)

    return args


def clean_up(debug, folder, action, signal):
    if not debug:

        from shutil import rmtree
        if action.startswith('u'):
            rmtree(folder + '/update/wf', ignore_errors=True)
        else:
            rmtree(folder + '/normal/wf', ignore_errors=True)
    if signal:
        logger.warning('EXITED WITH ERRORS')
    else:
        logger.info('FINISHED')
    sys.exit(signal)


def makeconfig(gene_names, sequences, organisms, len_threshold='50:5000', its='None', query_cover='None',
               identity='None', distance=2, subsp='False', excluded='None', remote='True', blacklist='', synonyms='',
               mdate='1000/01/01', finished='False', update='False', date='', folder='.'):
    config = configparser.ConfigParser()
    config.add_section('StartValues')
    config.set('StartValues', 'Included_Organisms', str(organisms))
    config.set('StartValues', 'Excluded_Organisms', str(excluded))
    config.set('StartValues', 'Gene_Names', str(gene_names))
    config.set('StartValues', 'Sequences', str(sequences))
    config.set('StartValues', 'Length_Threshold', len_threshold)
    config.set('StartValues', 'Its_Joiner', str(its))
    config.set('StartValues', 'Query_Cover', str(query_cover))
    config.set('StartValues', 'Identity', str(identity))
    config.set('StartValues', 'String_Distance', str(distance))
    config.set('StartValues', 'Subspecies', str(subsp))
    config.set('StartValues', 'Remote_Blast', str(remote))
    config.set('StartValues', 'Blacklist', blacklist)
    config.set('StartValues', 'Synonyms', synonyms)
    config.add_section('Finish')
    config.set('Finish', 'Finished', str(finished))
    config.set('Finish', 'Date', date)
    config.set('Finish', 'Updated', str(update))
    config.set('Finish', 'Update_Date', mdate)

    with open(folder + '/config.ini', 'w') as configfile:
        config.write(configfile)


def readconfig(file='/config.ini', folder=''):
    config = configparser.ConfigParser()
    config.read(folder + file)
    genes = ast.literal_eval(config.get('StartValues', 'Gene_Names'))
    seqs = ast.literal_eval(config.get('StartValues', 'Sequences'))
    org = ast.literal_eval(config.get('StartValues', 'Included_Organisms'))
    excl = ast.literal_eval(config.get('StartValues', 'Excluded_Organisms'))
    lthres = config.get('StartValues', 'Length_Threshold')
    its = ast.literal_eval(config.get('StartValues', 'Its_Joiner'))
    qc = ast.literal_eval(config.get('StartValues', 'Query_Cover'))
    idy = ast.literal_eval(config.get('StartValues', 'Identity'))
    dist = config.getint('StartValues', 'String_Distance')
    subsp = ast.literal_eval(config.get('StartValues', 'Subspecies'))
    rem = config.getboolean('StartValues', 'Remote_Blast')
    blk = config.get('StartValues', 'Blacklist')
    syn = config.get('StartValues', 'Synonyms')
    if config.has_section('Finish'):
        date = config.get('Finish', 'Date')
        fin = config.getboolean('Finish', 'Finished')
        upd = config.getboolean('Finish', 'Updated')
        mdate = config.get('Finish', 'Update_Date')
        return genes, seqs, org, lthres, its, qc, idy, dist, subsp, excl, rem, blk, syn, mdate, fin, upd, date
    else:
        return genes, seqs, org, lthres, its, qc, idy, dist, subsp, excl, rem, blk, syn


def dameraulevenshtein(seq1, seq2):
    """
    Author: Michael Homer (adapted to work with Python3)
    Date: Sunday, April 26th, 2009
    License: MIT

    Calculate the Damerau-Levenshtein distance between sequences.

    This distance is the number of additions, deletions, substitutions,
    and transpositions needed to transform the first sequence into the
    second. Although generally used with strings, any sequences of
    comparable objects will work.

    Transpositions are exchanges of *consecutive* characters; all other
    operations are self-explanatory.

    This implementation is O(N*M) time and O(M) space, for N and M the
    lengths of the two sequences.

    >>> dameraulevenshtein('ba', 'abc')
    2
    >>> dameraulevenshtein('fee', 'deed')
    2

    It works with arbitrary sequences too:
    >>> dameraulevenshtein('abcd', ['b', 'a', 'c', 'd', 'e'])
    2
    """
    # codesnippet:D0DE4716-B6E6-4161-9219-2903BF8F547F
    # Conceptually, this is based on a len(seq1) + 1 * len(seq2) + 1 matrix.
    # However, only the current and two previous rows are needed at once,
    # so we only store those.
    oneago = None
    thisrow = list(range(1, len(seq2) + 1)) + [0]
    for x in range(len(seq1)):
        # Python lists wrap around for negative indices, so put the
        # leftmost column at the *end* of the list. This matches with
        # the zero-indexed strings and saves extra calculation.
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(seq2) + [x + 1]
        for y in range(len(seq2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (seq1[x] != seq2[y])
            thisrow[y] = min(delcost, addcost, subcost)
            # This block deals with transpositions
            if x > 0 and y > 0 and seq1[x] == seq2[y - 1] and seq1[x - 1] == seq2[y] and seq1[x] != seq2[y]:
                thisrow[y] = min(thisrow[y], twoago[y - 2] + 1)
    return thisrow[len(seq2) - 1]
