# getKrakenFasta.py
## Download refseq fasta genomic data and convert header to work with Kraken

Similar approach as in:
[http://www.opiniomics.org/building-a-kraken-database-with-new-ftp-structure-and-no-gi-numbers/](http://www.opiniomics.org/building-a-kraken-database-with-new-ftp-structure-and-no-gi-numbers/)

I borrowed freely, however, this is a `Python2` rewrite. It does not require any
fancy third-party libraries and should work _out-of-the-box_.

```bash
$ python getKrakenFasta.py -h
usage: getKrakenFasta.py [-h] [--version] [-b BRANCH] [-t TYPE]

Download fasta-genomic sequences from ncbi and for building KrakenDB.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -b BRANCH, --branch BRANCH
                        Branches of organisms to download separated by comma,
                        e.g. bacteria, fungi, etc. [default="
                        bacteria,viral,fungi,protozoa,archaea"]
  -t TYPE, --type TYPE  Type of genomic sequences to include. For example:
                        Chromosome, Contig, Scaffold. [default="Complete
                        Genome"]

Copyright Sebastian Schmeier (s.schmeier@gmail.com)
```


```bash
$ python getKrakenFasta.py -b bacteria
$ python getKrakenFasta.py -b archaea
$ python getKrakenFasta.py -b viral
$ python getKrakenFasta.py -b protozoa
$ python getKrakenFasta.py -b fungi

# build a new database 
# download taxonomy
$ kraken-build --download-taxonomy --db kraken_bvfpa_201612

# for each branch, add all fna in the directory to the database
$ for dir in fungi protozoa archaea viral bacteria; do
        for fna in `ls $dir/*.fna`; do
                kraken-build --add-to-library $fna --db kraken_bvfpa_201612;
        done;
done

# build the actual database
$ kraken-build --build --db kraken_bvfpa_201612
```
