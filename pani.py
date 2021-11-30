#!/usr/bin/env python3
"""Script for calculating parallel ANI"""

__author__ = "Lavanya Rishishwar"
__copyright__ = "Copyright 2021, Lavanya Rishishwar"
__credits__ = ["Lavanya Rishishwar"]
__license__ = "GPL"
__version__ = "1.1"
__maintainer__ = "Lavanya Rishishwar"
__email__ = "lavanyarishishwar@gmail.com"
__status__ = "Development"
__title__ = "pani.py"

from argparse import ArgumentParser, HelpFormatter
import glob
import tqdm
import sys
import os
import string
import random
import logging
import subprocess
import multiprocessing

VERSION = __version__
PROGRAM_NAME = __title__


class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


py_version = sys.version_info
if py_version[0] < 3 or py_version[1] < 4:
    sys.exit(f"Error: {__title__} requires Python version 3.4+ to work. Please install a newer version of Python.")


# Function definition
def run_cmd_shell(command_str):
    try:
        output = subprocess.check_output(command_str, encoding="utf-8", shell=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Encountered an error executing the command: ")
        logging.error(f"Error details:")
        logging.error(f"Exit code={e.returncode}")
        logging.error(f"Error message={e.output}")
        sys.exit(1)
    return output


def get_ani(comp_file):
    ref, query = comp_file
    prefix = ''.join(random.choice(string.ascii_letters) for i in range(10))
    subprocess.check_output("ls")
    cmd = f"dnadiff -p temp-{prefix} {ref} {query} 2>&1 > temp-{prefix}.log"
    run_cmd_shell(cmd)
    cmd = f"head -19 temp-{prefix}.report | tail -1 | awk '{{print $2}}'"
    ani = float(run_cmd_shell(cmd).strip())
    cmd = f"rm temp-{prefix}.*"
    run_cmd_shell(cmd)
    logging.debug(f"{ref} vs {query} = {ani}")
    return ani


def get_files(folder1: str, ext: str, distance: bool, folder2: str = None):

    run_list = []
    ani_mat = {}
    file_names = {"rows": [], "cols": []}
    file_bnames = {"rows": [], "cols": []}
    
    if folder2 is None:
        # Matrix run mode
        file_list = glob.glob(os.path.join(folder1, f"*.{ext}"))
        if len(file_list) == 0:
            logging.critical(f"{Colors.FAIL}Error: Can't seem to find any files with extension {ext} in" +
                             f" the directory {folder1}! Exiting. {Colors.ENDC}")
            sys.exit(1)

        logging.info(f"Found {len(file_list)} files.")
        # Init
        file_names["rows"] = file_list
        file_names["cols"] = file_list
        for ref_i in range(len(file_list)):
            ani_mat[file_list[ref_i]] = {}
            if distance:
                ani_mat[file_list[ref_i]][file_list[ref_i]] = 0
            else:
                ani_mat[file_list[ref_i]][file_list[ref_i]] = 100
            this_name = os.path.splitext(os.path.basename(file_list[ref_i]))[0]
            
            file_bnames["rows"].append(this_name)
            file_bnames["cols"].append(this_name)
            
            for query_i in range(ref_i + 1, len(file_list)):
                run_list.append([file_list[ref_i], file_list[query_i]])

    else:
        # Ref run mode
        ref_list = glob.glob(os.path.join(folder1, f"*.{ext}"))
        query_list = glob.glob(os.path.join(folder2, f"*.{ext}"))
        if len(ref_list) == 0:
            logging.critical(f"{Colors.FAIL}Error: Can't seem to find any reference files with extension {ext} in" +
                             f" the directory {folder1}! Exiting. {Colors.ENDC}")
            sys.exit(1)

        if len(query_list) == 0:
            logging.critical(f"{Colors.FAIL}Error: Can't seem to find any query files with extension {ext} in" +
                             f" the directory {folder1}! Exiting. {Colors.ENDC}")
            sys.exit(1)

        logging.info(f"Found {len(ref_list)} ref x {len(query_list)} query files.")

        file_names["rows"] = query_list
        file_names["cols"] = ref_list
        for query_i in range(len(query_list)):
            ani_mat[query_list[query_i]] = {}
            file_bnames["rows"].append(os.path.splitext(os.path.basename(query_list[query_i]))[0])
        for ref_i in range(len(ref_list)):
            file_bnames["cols"].append(os.path.splitext(os.path.basename(ref_list[ref_i]))[0])

        for query_i in range(len(query_list)):
            for ref_i in range(len(ref_list)):
                run_list.append([query_list[query_i], ref_list[ref_i]])

    return(run_list, file_bnames, ani_mat, file_names)


# Main
if __name__ == "__main__":
    parser = ArgumentParser(prog=PROGRAM_NAME, add_help=False, description=f'''
                                {PROGRAM_NAME} - Parallel script for computing pairwise ANIm values.
                                ''',
                            formatter_class=lambda prog: HelpFormatter(prog, width=120, max_help_position=120))

    parser.add_argument('--help', '-h', '--h', action='help', help='Show this help message and exit.')

    subparsers = parser.add_subparsers(title=f"{PROGRAM_NAME} commands")
    matrix_parser = subparsers.add_parser('matrix', help='All-against-all comparison',
                                         formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                    max_help_position=120))
    matrix_parser.set_defaults(sub_command='matrix')
    matrix_parser.add_argument("-f", "--folder", help="Folder in which input files are located", default="./", required=False,
                        type=str, dest="folder")
    matrix_parser.add_argument("-o", "--out", help="Output file name", default="ani.tsv", required=False, type=str, dest="out")
    matrix_parser.add_argument("-t", "--threads", help="How many threads to use?", default=10, required=False, type=int,
                        dest="threads")
    matrix_parser.add_argument("-e", "--ext", help="File extension to look for", default="fasta", required=False, type=str,
                        dest="ext")
    matrix_parser.add_argument("-d", "--distance", help="Calculate distance instead of similarity", action='store_true')
    matrix_parser.add_argument("-l", "--logfile", help="Log file to save logs in", required=False, type=str,
                        dest="logfile", default="run.log")


    ref_parser = subparsers.add_parser('ref', help='Comparison against reference',
                                         formatter_class=lambda prog: HelpFormatter(prog, width=120,
                                                                                    max_help_position=120))
    ref_parser.set_defaults(sub_command='ref')
    ref_parser.add_argument("-r", "--ref", help="Folder which contains reference files", default="./", required=False,
                        type=str, dest="reference")
    ref_parser.add_argument("-q", "--query", help="Folder which contains query files", default="./", required=False,
                        type=str, dest="query")
    ref_parser.add_argument("-o", "--out", help="Output file name", default="ani.tsv", required=False, type=str, dest="out")
    ref_parser.add_argument("-t", "--threads", help="How many threads to use?", default=10, required=False, type=int,
                        dest="threads")
    ref_parser.add_argument("-e", "--ext", help="File extension to look for", default="fasta", required=False, type=str,
                        dest="ext")
    ref_parser.add_argument("-d", "--distance", help="Calculate distance instead of similarity", action='store_true')
    ref_parser.add_argument("-l", "--logfile", help="Log file to save logs in", required=False, type=str,
                        dest="logfile", default="run.log")


    args, unknown_arguments = parser.parse_known_args()

    if len(unknown_arguments) > 0:
        print(Colors.HEADER)
        print("User specificed unknown arguments: " + " ".join([str(x) for x in unknown_arguments]))
        print("Please see the correct usage below:")
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    if "sub_command" not in args:
        print(Colors.HEADER)
        parser.print_help()
        print(Colors.ENDC)
        sys.exit(0)

    # TODO: Check if directory exist
    # TODO: Check if file list is non-zero
    # TODO: Check if user supplied extension without starting .

    logging.basicConfig(filename=args.logfile, filemode='w', level=logging.INFO, format='[%(asctime)s] %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter(fmt=f"[%(asctime)s] %(message)s", datefmt="%Y-%m-%d %I:%M:%S %p")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

    if args.sub_command == 'matrix':
        run_list, file_bnames, ani_mat, file_names = get_files(folder1=args.folder, ext=args.ext, distance=args.distance)
    elif args.sub_command == 'ref':
        run_list, file_bnames, ani_mat, file_names = get_files(folder1=args.reference, ext=args.ext, distance=args.distance, folder2=args.query)

    logging.info(f"Need to compute {len(run_list)} comparisons.")
    with multiprocessing.Pool(args.threads) as p:
        ani = list(tqdm.tqdm(p.imap(get_ani, run_list), total=len(run_list)))

    logging.info(f"Refactoring results.")
    for i in range(len(ani)):
        if args.distance:
            ani[i] = 100 - ani[i]
        ani_mat[run_list[i][0]][run_list[i][1]] = ani[i]
        if args.sub_command == 'matrix':
            ani_mat[run_list[i][1]][run_list[i][0]] = ani[i]

    logging.info(f"Printing results to file {args.out}")
    with open(args.out, "w") as f:
        f.write("\t" + "\t".join(file_bnames["cols"]) + "\n")
        for row_i in range(len(file_bnames["rows"])):
            f.write(file_bnames["rows"][row_i])
            for col_i in range(len(file_bnames["cols"])):
                f.write(f"\t{ani_mat[file_names['rows'][row_i]][file_names['cols'][col_i]]}")
            f.write("\n")
