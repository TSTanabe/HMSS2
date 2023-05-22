#!/usr/bin/python


import sqlite3
from . import myUtil
#import re


"""
Einbindung des assembly perl scripts in python, soll dabei helfen die taxdump funktion von bioperl zu nutzen die bei biopython nicht verfügbar ist. Somit soll verhindert werden dass exzessiv viele Anfragen an den NCBI Server gesendet werden und lokal gearbeitet werden kann.

Python übernimmt nur die funktion eine Textdatei tsv auszuwerten und in die datenbank zu schreiben.
Perl übernimmt die Funktion alle assembly statistic files auszuwerten und gesondert in eine tsv datei zu schreiben zusammen mit der phylogenie aus taxdump. Diese Datei bleibt auch erhalten damit der nutzer auch eine weitere Auswertung vornehmen kann mittels R oder ähnlichen programmen.
dem Perl script muss daher ein directory mit assembly statistic files übergeben werden. diese funktion sollte von dem übergeordneten prozess vollzogen werden, da auch die fasta files in einem directory liegen müssen.
"""

def taxdump_perl_script(perlscript,directory,filepath,taxdump):
    """
    3.10.22
    
    Args:
        perlscript  executable compiled perl script
        directory   Directory of the assembly statistic files
        filepath    Path to Outputfile
        taxdump     location of the taxdump
        
    Output:
        TSV-file    Columns are: filename superkingdom clade phylum class order family genus
                                 species strain taxid biosample bioproject genbank refseq
                                 
    Warning:
        The perl script requires the NCBI taxdump files in the source folder 
        cwd/source/taxdump [..]
    """
    myUtil.command(f'{perlscript} -d {directory} -f {filepath} -t {taxdump}')
    return

def parse_gtdb_metadata(directory,queued_genomeIDs,filepath):
    """
    12.10.22
    
    Args:
        directory   Directory of the metadata file
        filepath    Path to Outputfile
        
    Output:
        TSV-file    Columns are: filename superkingdom clade phylum class order family genus
                                 species strain taxid biosample bioproject genbank refseq completeness contamination
                                 
    """
    genomeIDs = dict.fromkeys(queued_genomeIDs, 1)

    #file öffnen
    writer = open(filepath,"a")
    with open (directory,"r") as reader:
        
        header = reader.readline()
        headers = header.split("\t")
        #header_dict = {k: v for v, k in enumerate(headers)} #create header => index dictionary
        
    #accession	ambiguous_bases	checkm_completeness	checkm_contamination	checkm_marker_count	checkm_marker_lineage	checkm_marker_set_count	checkm_strain_heterogeneity	coding_bases	coding_density	contig_count	gc_count	gc_percentage	genome_size	gtdb_genome_representative	gtdb_representative	gtdb_taxonomy	gtdb_type_designation	gtdb_type_designation_sources	gtdb_type_species_of_genus	l50_contigs	l50_scaffolds	longest_contig	longest_scaffold	lsu_23s_contig_len	lsu_23s_count	lsu_23s_length	lsu_23s_query_id	lsu_5s_contig_len	lsu_5s_count	lsu_5s_length	lsu_5s_query_id	lsu_silva_23s_blast_align_len	lsu_silva_23s_blast_bitscore	lsu_silva_23s_blast_evalue	lsu_silva_23s_blast_perc_identity	lsu_silva_23s_blast_subject_id	lsu_silva_23s_taxonomy	mean_contig_length	mean_scaffold_length	mimag_high_quality	mimag_low_quality	mimag_medium_quality	n50_contigs	n50_scaffolds	ncbi_assembly_level	ncbi_assembly_name	ncbi_assembly_type	ncbi_bioproject	ncbi_biosample	ncbi_contig_count	ncbi_contig_n50	ncbi_country	ncbi_date	ncbi_genbank_assembly_accession	ncbi_genome_category	ncbi_genome_representation	ncbi_isolate	ncbi_isolation_source	ncbi_lat_lon	ncbi_molecule_count	ncbi_ncrna_count	ncbi_organism_name	ncbi_protein_count	ncbi_refseq_category	ncbi_rrna_count	ncbi_scaffold_count	ncbi_scaffold_l50	ncbi_scaffold_n50	ncbi_scaffold_n75	ncbi_scaffold_n90	ncbi_seq_rel_date	ncbi_spanned_gaps	ncbi_species_taxid	ncbi_ssu_count	ncbi_strain_identifiers	ncbi_submitter	ncbi_taxid	ncbi_taxonomy	ncbi_taxonomy_unfiltered	ncbi_total_gap_length	ncbi_total_length	ncbi_translation_table	ncbi_trna_count	ncbi_type_material_designation	ncbi_ungapped_length	ncbi_unspanned_gaps	ncbi_wgs_master	protein_count	scaffold_count	ssu_contig_len	ssu_count	ssu_gg_blast_align_len	ssu_gg_blast_bitscore	ssu_gg_blast_evalue	ssu_gg_blast_perc_identity	ssu_gg_blast_subject_id	ssu_gg_taxonomy	ssu_length	ssu_query_id	ssu_silva_blast_align_len	ssu_silva_blast_bitscore	ssu_silva_blast_evalue	ssu_silva_blast_perc_identity	ssu_silva_blast_subject_id	ssu_silva_taxonomy	total_gap_length	trna_aa_count	trna_count	trna_selenocysteine_count
        accession_idx = headers.index('accession')
        taxonomy_idx = headers.index('gtdb_taxonomy')
        completeness_idx = headers.index('checkm_completeness')
        contamination_idx = headers.index('checkm_contamination')
        ncbi_taxid_idx = headers.index('ncbi_taxid')
        ncbi_bioproject_idx = headers.index('ncbi_bioproject')
        ncbi_biosample_idx = headers.index('ncbi_biosample')
        ncbi_genbank_idx = headers.index('ncbi_genbank_assembly_accession')
        ncbi_typestrain_idx = headers.index('gtdb_type_designation')
        ncbi_strain_idx = headers.index('ncbi_strain_identifiers')
        
        for line in reader.readlines():
           
            listing = ()
            line.strip("\n")
            line = line.split("\t")    
            lineage = []
            typestrain = "0"
            
            genomeID = myUtil.getGenomeID(line[accession_idx])
            if genomeID in genomeIDs:
                #parse gtdb taxnonomy field
                gtdb_taxonomy = line[taxonomy_idx]
                gtdb_taxonomy.strip("\n")
                gtdb_taxonomy = gtdb_taxonomy.split(";")
                #domain phylum class order family genus species
                for field in gtdb_taxonomy:
                    lineage.append(field[3:])
                
                if line[ncbi_typestrain_idx] == "type strain of species":
                    typestrain = "1"
                    
                #filename superkingdom clade phylum class order family genus species strain taxid biosample bioproject genbank refseq completeness contamination typestrain
                listing = (line[accession_idx],lineage[0],"",lineage[1],lineage[2],lineage[3],lineage[4],lineage[5],lineage[6],\
                           line[ncbi_strain_idx],line[ncbi_taxid_idx],line[ncbi_biosample_idx],line[ncbi_bioproject_idx],line[ncbi_genbank_idx],"",\
                           line[completeness_idx],line[contamination_idx],typestrain)
                string = "\t".join(listing)
                writer.write(string+"\n")
    writer.close()
    return

def parse_custom_metadata(directory,queued_genomeIDs,filepath):
    """
    23.05.23
    
    Args:
        directory   Directory of the metadata file
        filepath    Path to Outputfile
        
    Output:
        TSV-file    Columns are: filename superkingdom clade phylum class order family genus
                                 species strain taxid biosample bioproject genbank refseq completeness contamination
                                 
    """
    genomeIDs = dict.fromkeys(queued_genomeIDs, 1)

    #file öffnen
    writer = open(filepath,"a")
    with open (directory,"r",newline="") as reader:
        
        header = reader.readline()
        header = header.replace("\n","")
        headers = header.split("\t")
        #header_dict = {k: v for v, k in enumerate(headers)} #create header => index dictionary
    #genomeID superkingdom clade phylum class order family genus species strain taxid biosample bioproject genbank refseq completeness contamination typestrain
        accession_idx = headers.index('genomeID')
        taxonomy_superkingdom = headers.index('superkingdom')
        taxonomy_clade = headers.index('clade')
        taxonomy_phylum = headers.index('phylum')
        taxonomy_class = headers.index('class')
        taxonomy_order = headers.index('order')
        taxonomy_family = headers.index('family')
        taxonomy_genus = headers.index('genus')
        taxonomy_species = headers.index('species')
        completeness_idx = headers.index('completeness')
        contamination_idx = headers.index('contamination')
        ncbi_taxid_idx = headers.index('taxid')
        ncbi_bioproject_idx = headers.index('bioproject')
        ncbi_biosample_idx = headers.index('biosample')
        ncbi_genbank_idx = headers.index('genbank')
        ncbi_typestrain_idx = headers.index('typestrain')
        ncbi_strain_idx = headers.index('strain')
        
        for line in reader.readlines():
           
            listing = ()
            line = line.replace("\n","")
            line = line.split("\t")    
            lineage = []
            typestrain = "0"
            
            genomeID = myUtil.getGenomeID(line[accession_idx])
            if genomeID in genomeIDs:
                ##parse gtdb taxnonomy field
                #gtdb_taxonomy = line[taxonomy_idx]
                #gtdb_taxonomy.strip("\n")
                #gtdb_taxonomy = gtdb_taxonomy.split(";")
                #domain phylum class order family genus species
                #for field in gtdb_taxonomy:
                #    lineage.append(field[3:])
                
                if line[ncbi_typestrain_idx]:
                    typestrain = "1"
  
                #TSV-file Columns are: filename => 0 superkingdom => 1 clade => 2 phylum => 3
                #                      class => 4 order => 5 family => 6 genus => 7 
                #                      species => 8 strain => 9 taxid => 10 biosample => 11
                #                      bioproject => 12 genbank => 13 refseq => 14 completeness => 15 contamination => 16 typestrain => 17

                    
                #filename superkingdom clade phylum class order family genus species strain taxid biosample bioproject genbank refseq completeness contamination typestrain
                listing = (line[accession_idx],line[taxonomy_superkingdom],line[taxonomy_clade],line[taxonomy_phylum],line[taxonomy_class],line[taxonomy_order],line[taxonomy_family],line[taxonomy_genus],line[taxonomy_species],\
                           line[ncbi_strain_idx],line[ncbi_taxid_idx],line[ncbi_biosample_idx],line[ncbi_bioproject_idx],line[ncbi_genbank_idx],"",\
                           line[completeness_idx],line[contamination_idx],typestrain)
                string = "\t".join(listing)
                writer.write(string+"\n")
    writer.close()
    return

def insert_database_assembly_statistics(database,filepath,genomeIDs):
    """
    3.10.22
    
    Args:
        database   Directory of the assembly statistic files
        filepath   TSV-file   Columns are: filename superkingdom clade phylum class
                                           order family genus species strain taxid 
                                           biosample bioproject genbank refseq
        genomeIDs   Set of genomeIds to be altered in the database
    Warning:
        The perl script requires the NCBI taxdump files in the source folder 
        cwd/source/taxdump [..]
    """
    
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("""PRAGMA foreign_keys = ON;""")
        try:
            os = open(filepath, "r")
        except:
            print("WARNING: No taxonomy found to assign. Taxonomy list in the database directory is not present")
            return
        with open(filepath, "r") as reader:
            reader.readline()
            for line in reader.readlines():
                line.strip("\n") #remove all newline characters, hopefully also all wrong formats from the string
                ar = line.split("\t") #array with tsv fields
                
                #TSV-file Columns are: filename => 0 superkingdom => 1 clade => 2 phylum => 3
                #                      class => 4 order => 5 family => 6 genus => 7 
                #                      species => 8 strain => 9 taxid => 10 biosample => 11
                #                      bioproject => 12 genbank => 13 refseq => 14 completeness => 15 contamination => 16 typestrain => 17
                #print(ar[1],ar[2],ar[3],ar[4],ar[5],ar[6],ar[7],ar[8],ar[9],ar[10],ar[11],ar[12],ar[13],ar[15],ar[16],ar[17])
                
                genomeID = myUtil.getGenomeID(ar[0])
                cur.execute("""UPDATE Genomes SET Superkingdom = ?, Clade = ?, Phylum = ?, Class = ?, Ordnung = ?, Family = ?, Genus = ?, Species = ?, Strain =?, NCBITaxon = ?, NCBIBiosample = ?, NCBIBioproject = ?, NCBIAssembly = ?, completeness = ?, contamination = ?, TypeStrain = ? WHERE genomeID = ?; """,(ar[1],ar[2],ar[3],ar[4],ar[5],ar[6],ar[7],ar[8],ar[9],ar[10],ar[11],ar[12],ar[13],ar[15],ar[16],ar[17],genomeID))
            
        con.commit()
        
    return 

#parse_gtdb_metadata("ptest/ar53_metadata_r207.tsv",("GB_GCA_002502215.1","GB_GCA_002502415.1"),"ptest/Output.txt")

