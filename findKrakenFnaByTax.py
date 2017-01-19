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
    description = "Convenience script. Find processed fasta-genomic sequences basenames of files based on taxonomy ids."
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument(
        dest='str_assembly',
        metavar='assembly_summary.txt',
        type=str,
        help='Assembly summary file for the branch of the tax-ids.')

    parser.add_argument(
        dest='str_file',
        metavar='FILE',
        type=str,
        help='One column file containing the taxonomy id ["-" for reading from stdin]')

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


def parse_assemblyfile(taxa, fname):

    if not os.path.isfile(fname):
        sys.stderr.write("ERROR: '%s' not found.\nEXIT\n\n"
                             % (fname))
        sys.exit(1)
    else:
        jobs = {}
        # read file, extract ftp paths and download each file
        oR = csv.reader(load_file(fname), delimiter = '\t')
        # for each line in assembly_summary.txt
        for a in oR:
            if a[0][0] == '#':
                continue
            
            taxid = a[5]
            
            if taxid not in taxa:
                continue

            version_status =  a[10]
            basename     = os.path.basename(a[19]) + '_genomic.tax.fna'
            assembly_level =  a[11]
            if taxid not in jobs:
                jobs[taxid] = set()
                
            jobs[taxid].add((basename,
                            assembly_level,
                            version_status))
        
    return jobs

    
def main():
    """ The main function. """
    args, parser = parse_cmdline()
    file_assembly = args.str_assembly

    # load tax ids to add to files
    file_in = load_file(args.str_file)
    reader = csv.reader(file_in, delimiter = '\t')
    dTax = {}
    for a in reader:
        dTax[a[0]] = None

    job_list = []
    taxa = dTax.keys()
    jobs_br = parse_assemblyfile(taxa, file_assembly)
    job_list += jobs_br
       
    for taxid in taxa:
        if taxid not in jobs_br:
            sys.stdout.write('%s\tn/a\tn/a\tn/a\n' %( taxid ))
        else:
            for t in jobs_br[taxid]:
                asem_stat = t[1]
                v_stat = t[2]
                filepath = t[0]
                   
                sys.stdout.write('%s\t%s\t%s\t%s\n' %( taxid, asem_stat, v_stat, filepath ))
                
    return


if __name__ == '__main__':
    sys.exit(main())

