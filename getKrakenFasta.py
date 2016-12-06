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


__version__ = '0.0.1'
__date__ = '2016/12/06'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'


class cd:
    """
    Context manager for changing the current working directory
    and changing back to previous path when context is exited.
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = 'Download fasta-genomic sequences from ncbi and for building KrakenDB.'
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
        help='Branches of organisms to download separated by comma, e.g. bacteria, fungi, etc. [default=" bacteria,viral,fungi,protozoa,archaea"]')
    
    parser.add_argument('-t',
        '--type',
        dest='str_type',
        metavar='TYPE',
        type=str,
        default="Complete Genome,Chromosome,Contig,Scaffold",
        help='Type of genomic sequences to include, separated by comma. For example: Chromosome, Contig, Scaffold. [default="Complete Genome,Chromosome,Contig,Scaffold"]')

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
    Here we download a file, get a status and adjust fasta-header
    return (args, res)
    """
    taxid = args[0]
    fname = args[1]
    dnlurl = args[2]
    org = args[3]
    md5sum_url = args[4]

    fnameTax = fname.replace('.fna.gz', '.tax.fna')
    with cd(org):
        
        code = 0
        if not os.path.isfile(fnameTax):
            code = 1
            if not os.path.isfile(fname):
                urllib.urlretrieve(dnlurl, fname)
            
            outfile = new_file(fnameTax)
            listfasta = load_file(fname).readlines()
            desc = listfasta[0].split(' ')
            listfasta[0] = '%s|kraken:taxid|%s %s'%(desc[0], taxid, ' '.join(desc[1:]))
            for s in listfasta:
                outfile.write(s)
            outfile.close()

            # remove original file
            os.remove(fname)
        
    return (args, code)


def parse_assemblyfile(branch, genomictypes=["Complete Genome"]):
    fname = 'assembly_summary.txt'
    url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/%s/%s' %(branch, fname)

    # create dir if not available
    if not os.path.exists(branch):
        os.makedirs(branch)

    jobs = []
    # change into that dir
    with cd(branch):
        if not os.path.isfile(fname):
            sys.stdout.write("Downloading assembly-summary-file: %s\n"%url)
            urllib.urlretrieve(url, fname)
            
        # read file, extract ftp paths and download each file
        oR = csv.reader(load_file(fname), delimiter = '\t')
        for a in oR:
            try:
                version_status =  a[11]
            except:
                continue
            
            if version_status in genomictypes:
                ftp_path = a[19]
                md5sum   = os.path.join(ftp_path, 'md5checksums.txt')  # should I ever implement md5sum check
                name     = os.path.basename(ftp_path) + '_genomic.fna.gz'
                dnlurl   = os.path.join(ftp_path, name)
                taxid    = a[5]
                jobs.append((taxid, name, dnlurl, branch, md5sum))
        
    return jobs


def main():
    """ The main function. """
    args, parser = parse_cmdline()

    branches = [s.strip() for s in args.str_branch.split(',')]
    types = [s.strip() for s in args.str_type.split(',')]
    
    for branch in branches:
        job_list = parse_assemblyfile(branch, types)
        # Could be parallised now, however, ftp connections might time out
        # if to many concurrent downloads are run
        # currently do not require speed.
        
        num_jobs=len(job_list)
        i = 0
        for job in job_list:
            i+=1
            res = my_func(job)
            sys.stdout.write('%i/%i: %s,%i\n'%(i, num_jobs, res[0][1], res[1]))
            sys.stdout.flush()
        
    return


if __name__ == '__main__':
    sys.exit(main())

