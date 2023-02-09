# HMSSS2
HMSSS version 2

# HMS-S-S
HMS-S-S: a tool for the identification of sulfur metabolism-related genes and analysis of operon structures in genome and metagenome assemblies

## Installing HMSSS on Linux
You can install HMSSS by downloading it directly from GitHub in compiled or non-compiled form.

1. Download the latest release from github

2. In a terminal, 'cd' to the downloaded package

3. Extract the files

4. Test you can run by './dist/HMSSS -h' or 'python HMSSS/HMSSS.py -h'

5. Installation of required external programs HMSSS depends on:

    5.1 HMMER3 package depends on the hmmer3 package, which can be downloaded from hmmer.org
    
    5.2 Prodigal for translation of nucleotide fasta github.com/hyattpd/Prodigal
  
6. That's it! You can now run HMSSS on a directory of protein sequence fasta files with gff files or nucleotide fasta files

## Running HMSSS

To run HMSSS on your own data type:

./HMSSS -f Directory

Replace "Directory" with the directory containing your input fasta files, with one file per species. The names of the files should match the identifiers of the genomes, since the name of the files will later be used for identification. When using protein fasta files with gff files, the names of the related files should be the same except for the file extension. HMSSS will look for input fasta files with any of the following filename extensions:

* .fna
* .fna.gz
* .faa
* .faa.gz
* .fasta

In case of any .fna file extension HMSSS will try to transcripe to protein fasta via prodigal. 
HMSSS also comes with several options which are scribed in the help accessed by `-h` and in the following:

### Extending the HMM library and gene cluster patterns
* `-l` sets the HMM library. By default this is set to the library in the source folder which includes the sulfur related HMMs. However this library can either be extended by or changed to any other HMM library compatible with the HMMER3 package. 
* `-t` sets the threshold file. Specific threshold scores for each HMM in the library are located here in a tab separated file. Each name is assigned to threshold score. In case of extended libraries the threshold scores should be set here.
* `-p` sets the syntenic gene patterns to be detected. All genes are listed in a tab separated file. Each line corresponds to one syntenic cluster. Collinearity is defined by the order of appearence. The first word in the line defines the name of the whole genecluster and is used as keyword.
*  `-c` sets the number of CPUs to be used by HMMER
*  `-nt` sets the number of nucleotides between two genes to be considered as syntenic. The distance is calculated between the closest ends of two genes.
*  `-mc` sets the minimal fraction of the gene cluster to be present to assign a keyword. By default 0.5. 
*  `-t`
*  `-t`
*  `-t`


























