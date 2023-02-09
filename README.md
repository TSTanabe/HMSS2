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
