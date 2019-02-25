"""
Argument parsing, executing the right action and main function of MatPhylobi
"""

import os
import sys
import argparse
import logging
import pkg_resources
from .classes import Run
from .utils import prepare, check, readconfig, makeconfig, check_blast, clean_up

logger = logging.getLogger('MatPhylobi')
logging.basicConfig(stream=sys.stdout, level=logging.INFO)


# -------------- Execute runs from args --------------


def exec_update(args):
    """
    execute update run
    """
    r = Run('u', args.folder, args.debug)
    r.start()


def exec_from_inputfile(args):
    """
    execute from provided config file
    """
    args.path = os.path.abspath(args.path)
    if not check(args.path, 'e'):
        clean_up(args.debug, args.folder, args.action, 1)

    logger.info("You are using the inputfile. All parameters other than folder, API key and debug will be ignored")
    try:
        startargs = readconfig(args.path)
        makeconfig(*startargs[:13], date=args.today, folder=args.folder)

        r = Run('n', args.folder, args.debug)
        r.start()

    except TypeError:
        logger.critical("Wrong data format. Check the documentation")
        clean_up(args.debug, args.folder, args.action, 1)


def exec_from_args(args):
    """
    execute from arguments provided
    """
    outfolder = args.folder + '/normal/'
    check(outfolder, 'm')

    makeconfig(str(args.gene_names), str(args.sequences), str(args.org_included),
               len_threshold=args.len_threshold,
               its=str(args.its), query_cover=str(args.query_cover), identity=str(args.identity),
               distance=str(args.string_distance), subsp=str(args.subsp), excluded=str(args.org_excluded),
               remote=str(args.remote_blast), folder=args.folder, date=args.today, blacklist=args.blacklist,
               synonyms=args.synonyms)

    r = Run('n', args.folder, args.debug)
    r.start()


# -------------- Argument parsing --------------

def add_analyze_arguments(parser_analyze):
    parser_analyze.add_argument('-g', '--gene_names', nargs='+',
                                help='List of marker names in the same order as corresponding accesion IDs '
                                     'specified by -s argument. Example: ITS rps16.', type=str, required=True)
    parser_analyze.add_argument('-s', '--sequences', nargs='+',
                                help='List of accesion IDs for each marker in the same order as corresponding '
                                     'gene names specified by -g argument. Multiple sequences for given marker '
                                     'should be separated by comma. Example: FJ415117.1 KJ832104.1,GU395133.1.',
                                required=True)
    parser_analyze.add_argument('-o', '--org_included', help='List of taxa included in analysis.',
                                nargs='*', default=[], type=str)
    parser_analyze.add_argument('--org_included_file', help='Path to file with list of taxa included in analysis.',
                                type=str)
    parser_analyze.add_argument('-i', '--its',
                                help='Join conspecific ITS1 & ITS2 into contig. '
                                     'There are three depth levels of algorithm. '
                                     '-i: match pairs using both specimen voucher and the authorship; '
                                     '-ii: match pairs using specimen voucher only; '
                                     '-iii: check only species name and loosen all other cirteria.',
                                action='count')
    parser_analyze.add_argument('-l', '--len_threshold',
                                help='Sequence length threshold in min:max format. Default = 50:5000.',
                                type=str,
                                default='50:5000')
    parser_analyze.add_argument('-d', '--string_distance',
                                help='Maximum number of misspellings allowed for two specimen vouchers '
                                     'still being considered the same. Default = 2.',
                                type=int,
                                default=2)
    parser_analyze.add_argument('-q', '--query_cover', nargs='*',
                                help='Minimal query cover for sequences found in BLAST search. Integer in range 0-100. '
                                     'Default = 0 for all markers. Threshold for each marker can be set individually. '
                                     'Example for two markers: 10 25.',
                                type=int)
    parser_analyze.add_argument('-y', '--identity', nargs='*',
                                help='Minimal identity for sequences found in BLAST search. Integer in range 0-100. '
                                     'Default = 0 for all markers. Threshold for each marker can be set individually. '
                                     'Example for two markers: 10 25.',
                                type=int)
    parser_analyze.add_argument('--org_excluded', help='List of taxa excluded from analysis.', nargs='*', default=[],
                                type=str)
    parser_analyze.add_argument('--org_excluded_file',
                                help='Path to file with list of taxa excluded from analysis.',
                                type=str)
    parser_analyze.add_argument('-f', '--folder',
                                help='Folder for output files. Default = present directory.',
                                type=str, default='./')
    parser_analyze.add_argument('-b', '--remote_blast', help='Use if you do not have downloaded BLAST databases.',
                                action='store_true')
    parser_analyze.add_argument('-p', '--subsp', help='Prevent pooling sequences of subspecies.', action='store_true')
    parser_analyze.add_argument('--blacklist',
                                help='Path to file with list of blacklisted sequences. See docs for details.',
                                type=str, default='')
    parser_analyze.add_argument('--synonyms',
                                help='Path to file with list of synonyms. See docs for details.',
                                type=str, default='')
    parser_analyze.add_argument('--debug',
                                help='Do not remove temporary files.',
                                action='store_true')
    parser_analyze.add_argument('--api_key', help='Personal API key for Biopython Entrez allowing for '
                                                  'maximum of 10 requests per second.')

    return parser_analyze


def add_inputfile_arguments(parser_inputfile):
    parser_inputfile.add_argument('-t', '--path',
                                  help='Path to input file with all parameters of analysis '
                                       'instead of providing them on command line.', type=str, required=True)
    parser_inputfile.add_argument('-f', '--folder',
                                  help='Folder for output files. Default = present directory.',
                                  type=str, default='./')
    parser_inputfile.add_argument('--debug',
                                  help='Do not remove temporary files.',
                                  action='store_true')
    parser_inputfile.add_argument('--api_key', help='Personal API key for Biopython Entrez allowing for '
                                                    'maximum of 10 requests per second.')
    return parser_inputfile


def add_update_arguments(parser_update):
    parser_update.add_argument('-f', '--folder',
                               help='Folder containing results of analysis meant to be updated '
                                    'Default = present directory.',
                               type=str, default='./')
    parser_update.add_argument('--debug',
                               help='Do not remove temporary files.',
                               action='store_true')
    parser_update.add_argument('--api_key', help='Personal API key for Biopython Entrez allowing for '
                                                 'maximum of 10 requests per second.')
    return parser_update


def create_parser():
    parser = argparse.ArgumentParser(prog='MatPhylobi', description='MatPhylobi - a tool for automatic construction of '
                                                                    'molecular data matrices for phylogenetic inference '
                                                                    'based on GenBank records.')

    # Add version argument
    ver = pkg_resources.require("MatPhylobi")[0].version
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=ver))  # TODO add also citation to the version

    # Prepare subparsers for different actions
    actions = parser.add_subparsers(title="program action", help="either run an analysis or update an older run",
                                    dest='action')
    actions.required = True

    parser_analyze = actions.add_parser('analyze')
    parser_analyze = add_analyze_arguments(parser_analyze)
    parser_analyze.set_defaults(func=exec_from_args)

    parser_inputfile = actions.add_parser('inputfile')
    parser_inputfile = add_inputfile_arguments(parser_inputfile)
    parser_inputfile.set_defaults(func=exec_from_inputfile)

    parser_update = actions.add_parser('update')
    parser_update = add_update_arguments(parser_update)
    parser_update.set_defaults(func=exec_update)

    return parser


def main():
    check_blast()

    if len(sys.argv) == 1:
        sys.argv.append('-h')

    # Parse the arguments
    parser = create_parser()
    args, _ = parser.parse_known_args()
    args = prepare(args, parser)

    # Act accordingly to the user-chosen action
    if args.action == 'inputfile':
        exec_from_inputfile(args)
    elif args.action == 'analyze':
        exec_from_args(args)
    elif args.action == 'update':
        exec_update(args)

    clean_up(args.debug, args.folder, args.action, 0)
