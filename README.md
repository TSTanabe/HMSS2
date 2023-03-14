# HMSSS2
HMSSS version 2

# HMS-S-S
HMS-S-S: a tool for the identification of sulfur metabolism-related genes and analysis of operon structures in genome and metagenome assemblies. It searches fasta files for sulfur metabolism associated proteins using hidden markov models and defined threshold scores. Furthermore, the genes of the detected proteins are analyzed for their position in the genome. The detected gene clusters are then named with a keyword if it is a known pattern of a gene cluster. HMSSS can also be extended with other compatible HMMs.

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
* `-p` sets the syntenic gene patterns to be detected. All genes are listed in a tab separated file. Each line corresponds to one syntenic cluster. Collinearity is defined by the order of appearence. The first word in the line defines the name of the whole genecluster and is used as keyword. A gencluster can be given several different keywords, but not the same one more than once.
*  `-c` sets the number of CPUs to be used by HMMER
*  `-nt` sets the number of nucleotides between two genes to be considered as syntenic. The distance is calculated between the closest ends of two genes.
*  `-mc` sets the minimal fraction of the gene cluster to be present to assign a keyword. if the match between the defined gene cluster pattern and the examined gene cluster is greater than this threshold, the corresponding keyword is given to the gene cluster. 
### Work step regulation
*  `-redo_csb` start at the collinear synthenic block prediction. This will only include data already stored in the database
*  `-redo_csb_naming`start at the synthenic block naming. This function restarts the comparison of gene clusters in the database with the entered gene cluster naming patterns. Old naming patterns are not overwritten.
*  `-redo_search` FASTA files wihich have already been searched and have an entry in the database will not be ignored but searched again.
*  `-redo_taxonomy` make the taxonomic assignment again.


### Result files and output
Results are stored in a local database which can be accessed to retrieve different results of interest. The local database can also be extended by later searched. If not defined before the search HMSSS will create a new local database for each run.
*  `-r` sets the directory for all results to be stored.
*  `-db` sets the database to be created/extended or from which results should be retrieved
*  `-gtdb` sets the path to a metadata file from the GTDB. This is required if it is desired to use the taxonomic information from GTDB.

The output from a database requires the `-db` opotion to define the database from which the desired output is taken. As the output is normally a set of sequences from a certain protein, possibly with a defined genomic vicinity or from specified taxonomic group there are several options to limit the number of retrieved sequences.
Limiting output to certain genomes:

*  `-dll` sets the level of taxonomy. If output should be limit to a group of organisms sharing the same taxonomic group this option sets the level of taxonomy between superkingdom and species.
*  `-dlt` sets the name of the taxonomic group. Together with `-dll` this defines the taxon for which results should be fetched
*  `-dlp` limits retrieved results to organisms which encode for the specified protein
*  `-dlk` limits retrieved results to organisms which encode for a genecluster with the specified keyword

Sequences for proteins can be retrieved and written to fasta files with the following commands:

*  `-fl` & `-ft` retrieve all sequences from all organisms of this taxon.`-fl` specifies taxonomic hirarchy level,`-ft` specifies the name of the taxon
*  `-fd` fetch sequences for proteins with the given domain. If several domains are desired these should be separated by whitespace characters. These will be handeled as connected by an logical OR, which means any protein matching one of the given domains will be retrieved.
*  `-fk` fetch sequences for proteins from geneclusters with the given keyword. If several keywords are desired these should be separated by whitespace characters. The connection between the keywords can be set to AND or OR by the `-kc` option. With AND proteins from geneclusters matching all given keywords will be fetched. Otherwise proteins from geneclusters matching any of the given keyword are retrieved. 
*  `-fd` & `-fk` in combination retrieves all sequences matching both, the given keyword(s) and the given domain(s).

The output always includes several files with reports of the written sequences and some subsets for proteins. Each file starts with a short summary of the given command, followed by the name of the protein whose sequences were written to the file. Proteins which have more than one domain previously detected by HMSSS are separately written to files and names by the all detected domains. Furthermore subsets prepared whch include sequences from genomes only encoding for a single ortholog. These are marked by the extension `singleton`. Second subset includes all remaining sequences and is marked by the extions `doublicate`. The metadata of the output sequences are presented in a report in a text file in a short summary. Columns of this file report the proteinID, domains, domains scores, domain coordinates in the proteins sequence, as well as the contig, gene stat and end, strand and locustag. If a genecluster is present the keyword , completeness of the gepattern and collinearity is also added to this list.   The columns of this file described the following data:
        
        Protein attributes:
        proteinID get_domains domain_scores domain_coordinates
        gene_contig gene_start gene_end gene_strand gene_locustag
        Cluster attributes:
        keyword completeness csb
        Taxonomy lineage 
        
        Reports in a textfile a abbreviated summary of the protein information with cluster and lineage



Information about the presence of given proteins and/or keywords in a taxon or species can be retrieved and written to tsv files. This also includes iTol dataset compatible files but sequences will not be retrieved:

*  `-dfd` retrieve presence/absence in the genome for given protein
*  `-dfk` retrieve presence/absence in the genome for given keyword
*  `-dff` retrieve presence/absence in the genome for given proteins with two or more domains in a genome
*  `-dfpk` retrieve presence/absence in the genome for given protein inside a gene cluster with the given keyword(s)
*  `-dmc` minimum similarity between pattern and gene cluster to be considered present
*  `-dcb` retrieve presence/absence in the genome for given keyword

The output contains the number for the presence of the desired proteine/keywords at each taxonomy level in absolute and relative values, each normalized to the number of genomes in the given taxonomy level. An iTol binary dataset is also output, with the specified names consisting of the genome identifiers and the taxonomic line

### Processing result files
Protein sequences are written to files with identifiers retrieved from the local database. These can be directly used. However if desired HMSSS also comes with some functions to alter these files.


*  `-merge_fasta` Merges one or more files with .faa extension into a single file without doublicates.
*  `-filter_fasta` retrieve presence/absence in the genome for given keyword
*  `-filter_limits` retrieve presence/absence in the genome for given protein
*  `-concat_alignment` retrieve presence/absence in the genome for given keyword
*  `-add_taxonomy` retrieve presence/absence in the genome for given keyword
*  `-add_genomic_context` retrieve presence/absence in the genome for given keyword
*  `-create_gene_cluster_dataset` retrieve presence/absence in the genome for given keyword
*  `-create_type_range_dataset` retrieve presence/absence in the genome for given keyword

















