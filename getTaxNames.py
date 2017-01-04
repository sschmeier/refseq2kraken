#!/usr/bin/env python
"""
NAME: getTaxnames.py
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

v0.1   2017/01/04    Initial version.

LICENCE
=======
2016, copyright Sebastian Schmeier (s.schmeier@gmail.com), http://sschmeier.com

template version: 1.5 (2016/09/09)
"""
from timeit import default_timer as timer
from signal import signal, SIGPIPE, SIG_DFL
import sys
import os
import os.path
import argparse
import csv
import collections
import gzip
import bz2
import zipfile
import time

# increase the csv field size
csv.field_size_limit(sys.maxsize)

# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
signal(SIGPIPE, SIG_DFL)

__version__ = 'v0.1'
__date__ = '2017/01/04'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'


def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = 'Attach scientific names to Kraken result tax-ids.'
    version = 'version %s, date %s' % (__version__, __date__)
    epilog = 'Copyright %s (%s)' % (__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='%s' % (version))

    parser.add_argument(
        'str_names',
        metavar='FILE',
        help=
        'Kraken taxonomy names.dmp file.')
    parser.add_argument(
        'str_file',
        metavar='FILE',
        help=
        'Kraken result file. [if set to "-" or "stdin" reads from standard in]')
   
   
   
    parser.add_argument('-o',
                        '--out',
                        metavar='STRING',
                        dest='outfile_name',
                        default=None,
                        help='Out-file. [default: "stdout"]')

    parser.add_argument(
        '--eval',
        action='store_true',
        dest='eval',
        default=False,
        help='Calc sen/spec. EXPERIMENTAL [default: False]')

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


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()

    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wb')
    else:
        outfileobj = open(args.outfile_name, 'w')

    # build the tax dict
    oF = load_file(args.str_names)
    reader = csv.reader(oF, delimiter = '\t')
    dict_tax = {}
    for a in reader:
        name_type = '_'.join(a[6].split(' '))
        if name_type not in dict_tax:
            dict_tax[name_type] = {}

        tax = int(a[0])
        dict_tax[name_type][tax] = dict_tax[name_type].get(tax, []) + [a[2]]
    oF.close()

    # load results
    oF = load_file(args.str_file)
    reader = csv.reader(oF, delimiter ='\t')

    list_keys = ['scientific_name',
                     'authority',
                     'synonym',
                     'type_material',
                     'genbank_common_name',
                     'equivalent_name',
                     'genbank_synonym',
                     'blast_name']
    
    iNum = 0
    iNotC = 0
    for a in reader:
        iNum += 1
        if a[0] == 'C':
            name = ''
            name_test = ''

            tax = int(a[2])
        
            for key in list_keys:
                if tax in dict_tax[key]:
                    name = '%s:%s'%(key,'|'.join(dict_tax[key][tax]))
                    break
            
            if args.eval:
                tax_test = int(a[1].split('|')[-1].strip())

                for key in list_keys:
                    if tax_test in dict_tax[key]:
                        name_test = '%s:%s'%(key,'|'.join(dict_tax[key][tax_test]))
                        break
                    
            outfileobj.write('%s\t%s\t%s\n' % ('\t'.join(a), name, name_test))
                
        else:
            outfileobj.write('%s\t''\t''\n' % ('\t'.join(a)))        
            iNotC += 1
           
    outfileobj.close()
    return


if __name__ == '__main__':
    sys.exit(main())

