# Download refseq-genomic data and prepare it for [Kraken](https://ccb.jhu.edu/software/kraken/)

Similar approach as in:
[http://www.opiniomics.org/building-a-kraken-database-with-new-ftp-structure-and-no-gi-numbers/](http://www.opiniomics.org/building-a-kraken-database-with-new-ftp-structure-and-no-gi-numbers/)

I borrowed freely, however, this is a `Python2` rewrite.

Its improved in two ways:

1. I split the tasks of  data download and kraken data preparation.
2. Both separate tasks allow multithreading to be used.

## Requirements

1. Python2
2. BioPython
3. Linux shell with rsync and Git
4. ~10GB of space for compressed ncbi refseq fasta-files
5. ~40GB of space for processed uncompressed kraken-readable fasta-files
6. ~130GB if a complete Kraken database is build without restricting its size (e.g. with --max-db-size 20)

## Download refseq genomic fasta-data via rsync (getRefseqGenomic.py)

This script will retrieve genomic data from refseq via rsync. It saves on downloads as only
files that updated or are new will be downloaded in sub-sequent runs.

**Warning! Using this script will make one rsync call to the ftp-server from ncbi per file you 
want to download. In case of bacteria and all assembly levels, this will result in ~70000 ftp-server 
accesses. There might be a limit on what ncbi allows in terms of connections to their ftp-server. 
If you overdo it, ncbi might take action against you. Use this script at your own risk.**


```bash
usage: getRefseqGenomic.py [-h] [--version] [-b BRANCH] [-l LEVEL] [-a]
                           [-p INT]

Download fasta-genomic sequences from ncbi using rsync.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -b BRANCH, --branch BRANCH
                        Branches of organisms to download separated by comma,
                        e.g. bacteria, fungi, etc. [default="
                        bacteria,viral,fungi,protozoa,archaea"]
  -l LEVEL, --level LEVEL
                        Assembly - level of genomic sequences to include,
                        separated by comma. For example: Chromosome, Contig,
                        Scaffold. [default="Complete Genome"]
  -a, --assembly        Print assembly stats for branches and exits.

Threading:
  Multithreading arguments:

  -p INT, --processes INT
                        Number of sub-processes (concurrent downloads) to use.
                        It is only logical to not give more processes than
                        cpus/cores are available. [default: 1]

Copyright Sebastian Schmeier (s.schmeier@gmail.com)
```

The simplest call would be:

```bash
# all default branches + all default "assembly_level"s of assemblies.
python getRefseqGenomic.py -p 8

# for only one branch only
python getRefseqGenomic.py -b archaea -p 8
```

This would download "Complete Genome" assembly genomic sequences
for "bacteria, viral, fungi, protozoa, archaea". However, this might not be what
you necessarily want . Instead one can use the following to download "Contig" and "Scaffold" assemblies for archaea:

```bash
python getRefseqGenomic.py -b archaea -l "Contig,Scaffold" -p 8
```

You can check the stats on the assembly levels with `-a`. This will only
download the assembly_summary.txt for the branch(es) and counts the different "assembly_level"s
available.

```bash
python getRefseqGenomic.py -a -b protozoa
Branch: protozoa
Chromosome      25
Complete Genome 2
Contig  2
Scaffold        45
```

Should you at a later stage re-run the command, `rsync` makes sure to only
download changed files (**Attention:** in terms of filesize, not content).

## Convert fasta-headers to work with Kraken (getKrakenFna.py)

This script will take fasta-files and create new "uncompressed" (Kraken
needs uncompressed files) fasta-files with each header changed to a form:
`>seq1|kraken:taxid|12345 blah`. This allows to use the new ncbi files (without
GI identifiers) with `Kraken`. Requieres the third-party `BioPython` lib.


```bash
usage: getKrakenFna.py [-h] [--version] [-b BRANCH] [-l LEVEL] [-d DIRECTORY]
                       [-a] [-p INT]
                       KrakenDB-DIR

Process fasta-genomic sequences from NCBI-refseq for inclusion in a KrakenDB.
A branch's assembly_summary.txt as well as all of its *_genomic.fna.gz files
are assumed to be in a sub-directory called like the branch, e.g.
./genomes/refseq/bacteria

positional arguments:
  KrakenDB-DIR          Directory name for the processed fasta-files for
                        kraken. Will be created if not found.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -b BRANCH, --branch BRANCH
                        Branches of organisms to download separated by comma,
                        e.g. bacteria, fungi, etc. [default="
                        bacteria,viral,fungi,protozoa,archaea"]
  -l LEVEL, --level LEVEL
                        Assembly - level of genomic sequences to include,
                        separated by comma. For example: Chromosome, Contig,
                        Scaffold. [default="Complete Genome"]
  -d DIRECTORY, --dir DIRECTORY
                        Base directory for refseq fasta-files. Here, we assume
                        sub-directories for branches, e.g. bacteria etc.
                        [default="./genomes/refseq/"]
  -a, --assembly        Print assembly stats for branches and exits.

Threading:
  Multithreading arguments:

  -p INT, --processes INT
                        Number of concurrent sub-processes to use. It is only
                        logical to not give more processes than cpus/cores are
                        available. [default: 1]

Copyright Sebastian Schmeier (s.schmeier@gmail.com)
```

Usage:

```bash
python getKrakenFna.py -b archaea -p 8 kraken_201612
# will create a dir kraken_201612/archaea with converted files.
```

## Putting it all together

```bash
# clone repo
git clone https://github.com/sschmeier/refseq2kraken.git kraken
cd kraken

# Download refseq => here only "Complete Genome"
# assemblies, e.g. the default
python getRefseqGenomic.py -p 8

# convert to kraken format => again only "Complete Genome"
# assemblies here, e.g. the default
python getKrakenFna.py -p 8 kraken_201612

# install kraken if you must
conda create -n kraken kraken-all
source activate kraken

# build a new minikraken database 
# download taxonomy
kraken-build --download-taxonomy --db kraken-db-bva_201612

# for each branch, add all fna in the directory to the database
for dir in bacteria viral archaea ; do
    find kraken_201612/$dir/ -name '*.fna' -print0 | xargs -0 -I{} -n1 -P8 kraken-build --add-to-library {} --db kraken-db-bva_201612;
done

# build the actual database restrict its size here to 20GB
kraken-build --build --db kraken-db-bva_201612 --max-db-size 20

# remove intermediate files
kraken-build --clean --db kraken-db-bva_201612

# classify seqs
kraken --db kraken-db-bva_201612 test_seqs/bva_test.fa > test_seqs/bva-results.txt

# get report of all taxa found in sample
kraken-report --show-zeros --db kraken-db-bva_201612 test_seqs/bva-results.txt | sort -n -k 5 | gzip > test_seqs/bva-results-report.txt.gz

# ONLY for testing the 
# attach taxonomy names of classification resutls + test seqs original tax names (use without --eval for non-test case)
# python getTaxNames.py kraken-db-bva_201612/taxonomy/names.dmp test_seqs/bva-results.txt --eval > test_seqs/bva-results-names.txt
```
