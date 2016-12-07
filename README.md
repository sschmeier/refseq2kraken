# Download refseq-genomic data and prepare it for Kraken

Similar approach as in:
[http://www.opiniomics.org/building-a-kraken-database-with-new-ftp-structure-and-no-gi-numbers/](http://www.opiniomics.org/building-a-kraken-database-with-new-ftp-structure-and-no-gi-numbers/)

I borrowed freely, however, this is a `Python2` rewrite.

I improved in two ways:

1. I split the tasks of  data download and kraken data preparation.
2. Both separate tasks allow multithreading to be used.


##1. Download refseq genomic fasta-data via rsync (getRefseqGenomic.py)

This script will retrieve genomic data from refseq via rsync. It saves on downloads as only
files that updated or are new will be downloaded in sub-sequent runs.

```bash
usage: getRefseqGenomic.py [-h] [--version] [-b BRANCH] [-t TYPE] [-a]
                           [-p INT]

Download fasta-genomic sequences from ncbi using rsync.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -b BRANCH, --branch BRANCH
                        Branches of organisms to download separated by comma,
                        e.g. bacteria, fungi, etc. [default="
                        bacteria,viral,fungi,protozoa,archaea"]
  -t TYPE, --type TYPE  Type of genomic sequences to include, separated by
                        comma. For example: Chromosome, Contig, Scaffold.
                        [default="Complete Genome,Chromosome,Contig,Scaffold"]
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
# all default branches + all default "version_status" of assemblies.
python getRefseqGenomic.py -p 8

# for only one branch only
python getRefseqGenomic.py -b archaea -p 8
```

This would download "Complete Genome, Chromosome, Contig, Scaffold" genomic sequences
for "bacteria, viral, fungi, protozoa, archaea". However, this might not be what
you necessarily want, e.g. for bacteria the download becomes huge. Instead one
can use the following to only download "Compel Genome" for bacteria:

```bash
python getRefseqGenomic.py -b bacteria -t "Complete Genome" -p 8
```

You can check the stats on the assembly versions with `-a`. This will only
download the assembly_summary.txt for the branch(es) and counts the different "version_status"
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
download changed files (in terms of filesize).


##2. Convert fasta-headers to work with Kraken (getKrakenFna.py)
This script will take fasta-files and create a new "uncompressed" (as Kraken
needs uncompressed files) fasta-files with each header changed to a form:
`>seq1|kraken:taxid|12345 blah`. This allows to use the new ncbi files (without
GI identifiers) with `Kraken`. Requieres the third-party `BioPython` lib.


```bash
usage: getKrakenFna.py [-h] [--version] [-b BRANCH] [-t TYPE] [-d DIRECTORY]
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
  -t TYPE, --type TYPE  Type of genomic sequences to include, separated by
                        comma. For example: Chromosome, Contig, Scaffold.
                        [default="Complete Genome,Chromosome,Contig,Scaffold"]
  -d DIRECTORY, --dir DIRECTORY
                        Base directory for refseq fasta-files. Here, we assume
                        sub-directories for branches, e.g. bacteria etc.
                        [default="./genomes/refseq/"]
  -a, --assembly        Print assembly stats for branches and exit.

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
# Download refseq
python getRefseqGenomic.py -b fungi -p 8
python getRefseqGenomic.py -b archaea -p 8
python getRefseqGenomic.py -b viral -p 8
python getRefseqGenomic.py -b protozoa -p 8
python getRefseqGenomic.py -b bacteria -t "Complete Genome" -p 8

# convert to kraken format
python getKrakenFna.py -b fungi -p 8 kraken_201612
python getKrakenFna.py -b viral -p 8 kraken_201612
python getKrakenFna.py -b archaea -p 8 kraken_201612
python getKrakenFna.py -b protozoa -p 8 kraken_201612
python getKrakenFna.py -b bacteria -p 8 -t "Complete Genome" kraken_201612

# install kraken if you must
conda create -n kraken kraken-all
source activate kraken

# build a new database 
# download taxonomy
kraken-build --download-taxonomy --db kraken-db-bvapf_201612

# for each branch, add all fna in the directory to the database
for dir in bacteria viral archaea protozoa fungi; do
        for fna in `ls kraken_201612/$dir/*.fna`; do
                kraken-build --add-to-library $fna --db kraken-db-bvapf_201612;
        done;
done

# build the actual database
kraken-build --build --db kraken-db-bvapf_201612
```
