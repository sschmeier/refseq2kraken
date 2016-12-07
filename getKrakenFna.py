#!/usr/bin/env python2
"""
NAME: getKrakenFna.py
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
from Bio import SeqIO # non-standard lib BioPython
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
import time


__version__ = '0.0.1'
__date__ = '2016/12/06'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'


def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = "Process fasta-genomic sequences from NCBI-refseq for inclusion in a KrakenDB. A branch's assembly_summary.txt as well as all of its *_genomic.fna.gz files are assumed to be in a sub-directory called like the branch, e.g. ./genomes/refseq/bacteria"
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument(
        dest='str_kraken',
        metavar='KrakenDB-DIR',
        type=str,
        help='Directory name for the processed fasta-files for kraken. Will be created if not found.')

    parser.add_argument('-b',
        '--branch',
        dest='str_branch',
        metavar='BRANCH',
        type=str,
        default="bacteria,viral,fungi,protozoa,archaea",
        help='Branches of organisms to download separated by comma, e.g. bacteria, fungi, etc. [default=" bacteria,viral,fungi,protozoa,archaea"]')
    
    parser.add_argument('-l',
        '--level',
        dest='str_level',
        metavar='LEVEL',
        type=str,
        default="Complete Genome",
        help='Assembly - level of genomic sequences to include, separated by comma. For example: Chromosome, Contig, Scaffold. [default="Complete Genome"]')
    parser.add_argument('-d',
        '--dir',
        dest='str_dir',
        metavar='DIRECTORY',
        type=str,
        default="./genomes/refseq/",
        help='Base directory for refseq fasta-files. Here, we assume sub-directories for branches, e.g. bacteria etc. [default="./genomes/refseq/"]')

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
        'Number of concurrent sub-processes to use.'+\
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


def new_file(outfile_name):
    # create outfile object
    if not outfile_name:
        outfileobj = sys.stdout
    elif outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(outfile_name, 'wb')
    else:
        outfileobj = open(outfile_name, 'w')
    return outfileobj


def my_func(args):
    """
    THIS IS THE ACCTUAL WORKFUNCTION THAT HAS TO BE EXECUTED MULTPLE TIMES.
    This function could be distributed to the cores requested.
    # do stuff
    args = (taxid, infile-path, outfile-path)
    
    return (args, res)
    """
    taxid = args[0]
    infilepath = args[1]
    outfilename = args[2]
    if not os.path.isfile(infilepath):
        sys.stderr.write('%s not found. SKIP\n'%(infilepath))
        return (args, 0)
    else:
        fasta_records = SeqIO.parse(load_file(infilepath), "fasta")
        outfile = new_file(outfilename)
    
    # here we iteracte over all records and change the header appropriately
    # >seq1|kraken:taxid|12345 original stuff
    for record in fasta_records:
        record.id = '%s|kraken:taxid|%s'%(record.id, taxid)
        SeqIO.write(record, outfile, "fasta")
        
    outfile.close()
    return (args, 1)


def parse_assemblyfile(branch, genomictypes=["Complete Genome"], dirpath='./genomes/refseq/', krakendir='./kraken'):
    basedir = os.path.join(dirpath, branch)
    fname = 'assembly_summary.txt'
    krakendir = os.path.join(krakendir, branch)

    if not os.path.isfile(os.path.join(basedir,fname)):
        sys.stderr.write("ERROR: '%s' not for branch '%s' at %s\nUse 'python getRefseqGenomic.py -b %s -p 4' to download fasta-sequences and assembly-summary.txt.\nEXIT\n\n"
                             % (fname, branch, os.path.join(basedir,fname), branch))
        sys.exit(1)
    else:
        jobs = []
        # read file, extract ftp paths and download each file
        oR = csv.reader(load_file(os.path.join(basedir,fname)), delimiter = '\t')
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
                name     = os.path.basename(a[19]) + '_genomic.fna.gz'
                filepath = os.path.join(basedir, name)
                taxid    = a[5]

                #fnameTax = name.replace('.fna.gz', '.tax.fna')
                fnameTax = name.replace('.fna.gz', '.tax.fna.gz')  # store gziped files
                jobs.append((taxid, filepath, os.path.join(krakendir, fnameTax)))
        
    return jobs, d


def main():
    """ The main function. """
    args, parser = parse_cmdline()

    branches = [s.strip() for s in args.str_branch.split(',')]
    types = [s.strip() for s in args.str_level.split(',')]
    dirpath = args.str_dir
    process_number = args.process_number
    if process_number < 1:
        parser.error('-p has to be > 0: EXIT.')

    job_list = []
    for branch in branches:
        job_list_br, dStats = parse_assemblyfile(branch, types, dirpath, args.str_kraken)
        job_list += job_list_br
        if args.assemblystats:
            status = dStats.keys()
            status.sort()
            sys.stdout.write('Branch: %s\n'%branch)
            for s in status:
                sys.stdout.write('%s\t%i\n' %(s,dStats[s]))
        else:
            # prepare directory structure for results
            if not os.path.exists(os.path.join(args.str_kraken, branch)):
                sys.stdout.write('Make directory for kraken-files: %s\n'%(os.path.join(args.str_kraken, branch)))
                os.makedirs(os.path.join(args.str_kraken, branch))

    # exit if only stats should be displayed
    if args.assemblystats:
        return

    #-------------------------------------------------------------------------
    # MULTITHREADING
    #-------------------------------------------------------------------------
    # For timing
    start_time = timer()  # very crude timing
    # create pool of workers ---------------------
    pool = Pool(processes=process_number)

    # "chunksize"" usually only makes a noticable performance
    # difference for very large iterables
    # Here I set it to 1 to get the progress bar working nicly
    # Otherwise it will not give the correct number of processes left
    # to process but rather the chunksize number.
    chunksize = 1
    result_list = pool.map_async(my_func, job_list, chunksize=chunksize)
    pool.close()  # No more work

    jobs_total = len(job_list)
    # Progress bar
    #==============================
    # This can be changed to make progressbar bigger or smaller
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
        time.sleep(1)  # wait a bit: here we test every sec
    # Finish the progress bar
    bar_str = '=' * progress_bar_length
    sys.stderr.write("JOBS (%s): [%s] (%i) 100%%\n" % ('0'.rjust(len(str(jobs_total))),
                                                       bar_str,
                                                       jobs_total))
    #result_list = result_list.get()
    # --------------------------------------------

    end_time = timer()
    sys.stderr.write('PROCESS-TIME: %.1f sec\nDONE.\n\n' % (end_time - start_time))

    
    return


if __name__ == '__main__':
    sys.exit(main())

