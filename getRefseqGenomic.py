#!/usr/bin/env python2
"""
NAME: getKrakenFasta.py
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.0.1   2016/12/06    Initial version.

LICENCE
=======
2016, copyright Sebastian Schmeier (s.schmeier@gmail.com), sschmeier.com

template version: 1.6 (2016/11/09)
"""
from timeit import default_timer as timer
from multiprocessing import Pool
import sys
import os
import os.path
import argparse
import csv
import gzip
import bz2
import zipfile
import urllib
import hashlib
import subprocess
import time


__version__ = '0.0.1'
__date__ = '2016/12/06'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'


def parse_cmdline():
    """ Parse command-line args. """
    description = 'Download fasta-genomic sequences from ncbi using rsync.'
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument('-b',
        '--branch',
        dest='str_branch',
        metavar='BRANCH',
        type=str,
        default="bacteria,viral,fungi,protozoa,archaea",
        help='Branches of organisms to download separated by comma, e.g. bacteria,fungi, etc. [default=" bacteria,viral,fungi,protozoa,archaea"]')
    
    parser.add_argument('-l',
        '--level',
        dest='str_level',
        metavar='LEVEL',
        type=str,
        default="Complete Genome",
        help='Assembly - level of genomic sequences to include, separated by comma. For example: Chromosome,Contig,Scaffold. [default="Complete Genome"]')

    parser.add_argument('-a',
        '--assembly',
        dest='assemblystats',
        default=False,
        action='store_true',
        help='Print assembly stats for branches and exits.')

    group1 = parser.add_argument_group('Threading',
                                       'Multithreading arguments:')

    group1.add_argument(
        '-p',
        '--processes',
        metavar='INT',
        type=int,
        dest='process_number',
        default=1,
        help=
        'Number of sub-processes (concurrent downloads) to use.'+\
        ' It is only logical to not give more processes'+\
        ' than cpus/cores are available. [default: 1]')

    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args, parser


def load_file(filename):
    """ LOADING FILES """
    if filename in ['-', 'stdin']:
        filehandle = sys.stdin
    elif filename.split('.')[-1] == 'gz':
        filehandle = gzip.open(filename)
    elif filename.split('.')[-1] == 'bz2':
        filehandle = bz2.BZFile(filename)
    elif filename.split('.')[-1] == 'zip':
        filehandle = zipfile.Zipfile(filename)
    else:
        filehandle = open(filename)
    return filehandle


def my_func(args):
    """
    THIS IS THE ACCTUAL WORKFUNCTION THAT HAS TO BE EXECUTED MULTPLE TIMES.
    This function could be distributed to the cores requested.
    # do stuff
    Here we download a file, get a status and adjust fasta-header
    return (args, res)
    """
    fname = args[0]
    dnlurl = args[1]
    dest_dir = args[2]
    rsync_cmd = "rsync --times --copy-links --partial -aq %s %s"
    retcode = subprocess.call(rsync_cmd % (dnlurl, dest_dir), shell=True)
    return (args, retcode)


def parse_assemblyfile(branch, genomictypes=["Complete Genome"], dest_dir='genomes'):
    fname = 'genomes/refseq/%s/assembly_summary.txt' % branch
    url = 'rsync://ftp.ncbi.nlm.nih.gov/%s' %(fname)

    rsync_cmd = "rsync -R --times --copy-links --partial -aq %s %s"
    retcode = subprocess.call(rsync_cmd % (url, dest_dir), shell=True)
   
    jobs = []
    # read file, extract ftp paths and download each file
    oR = csv.reader(load_file(fname), delimiter = '\t')
    d = {}
    for a in oR:
        
        try:
            assembly_level =  a[11]
        except:
            continue

        version_status =  a[10]
        if version_status != 'latest':
            continue
    
        if assembly_level == 'assembly_level':
            continue

        d[assembly_level] = d.get(assembly_level, 0) + 1

        if assembly_level in genomictypes:
            ftp_path = a[19]
            name     = os.path.basename(ftp_path) + '_genomic.fna.gz'
            dnlurl   = os.path.join(ftp_path, name)
            dnlurl = dnlurl.replace('ftp://', 'rsync://')
            jobs.append((name, dnlurl, 'genomes/refseq/%s' % branch, branch))
    return jobs, retcode, d


def main():
    """ The main function. """
    args, parser = parse_cmdline()

    process_number = args.process_number
    if process_number < 1:
        parser.error('-p has to be > 0: EXIT.')

    branches = [s.strip() for s in args.str_branch.split(',')]
    types = [s.strip() for s in args.str_level.split(',')]

    job_list = []
    for branch in branches:
        job_list_branch, retcode, dStats = parse_assemblyfile(branch, types, 'genomes')
        job_list += job_list_branch
        if args.assemblystats:
            status = dStats.keys()
            status.sort()
            sys.stdout.write('Branch: %s\n'%branch)
            for s in status:
                sys.stdout.write('%s\t%i\n' %(s,dStats[s]))

    # exit if only stats should be displayed
    if args.assemblystats:
        return
        
    #-------------------------------------------------------------------------
    # MULTITHREADING
    #-------------------------------------------------------------------------
    start_time = timer()  # very crude timing
    # create pool of workers ---------------------
    pool = Pool(processes = process_number)
    result_list = pool.map_async(my_func, job_list, chunksize=1)  # chunksize=1 for correct progressbar
    pool.close()  # No more work

    jobs_total = len(job_list)
    progress_bar_length = 50
    #==============================
    while not result_list.ready():
        num_not_done = result_list._number_left
        num_done = jobs_total - num_not_done
        num_bar_done = num_done * progress_bar_length / jobs_total
        bar_str = ('=' * num_bar_done).ljust(progress_bar_length)
        percent = int(num_done * 100 / jobs_total)
        sys.stderr.write("JOBS (%s): [%s] (%s) %s%%\r" % (str(num_not_done).rjust(len(str(jobs_total))),
                                                          bar_str,
                                                          str(num_done).rjust(len(str(jobs_total))),
                                                          str(percent).rjust(3)))

        sys.stderr.flush()
        time.sleep(2)  # wait a bit: here we test every 2 sec
    # Finish the progress bar
    bar_str = '=' * progress_bar_length
    sys.stderr.write("JOBS (%s): [%s] (%i) 100%%\n" % ('0'.rjust(len(str(jobs_total))),
                                                       bar_str,
                                                       jobs_total))
    #result_list = result_list.get()
    end_time = timer()
    sys.stderr.write('PROCESS-TIME: %.1f sec\nDONE.\n\n' % (end_time - start_time))
    #-------------------------------------------------------------------------


    return


if __name__ == '__main__':
    sys.exit(main())

