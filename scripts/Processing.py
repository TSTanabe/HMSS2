#!/usr/bin/python
import os               # For file and directory operations
import subprocess       # For running external commands (e.g., MAFFT, trimAl)
import shutil           # For finding executables in the system
import sys              # For getting system-specific parameters and functions
import csv              # For reading and writing CSV (TSV) files
import sqlite3          # For interacting with SQLite databases

from Bio import SeqIO    # For reading and writing sequence files (FASTA format)

from . import myUtil     # Custom module for utility functions (assumed to be in the same package)
from . import ParseReports  # Custom module for parsing reports (assumed to be in the same package)
from . import Csb_finder    # Custom module for finding conserved sequence blocks (assumed to be in the same package)



#########################################################################
####################### MAIN OUTPUT ROUTINE #############################
#########################################################################


##############################################################
##########    Fasta File Postprocessing ON COMMAND    ########
##############################################################

def merge_fasta(options):
    """
    Concatenates all .faa files in a directory, ensuring all identifiers are unique using Biopython.

    Args:
        directory: Directory containing .faa files.
        output_file: Name of the output file for the concatenated sequences.
    """
    directory = options.merge_fasta
    output_file = options.merge_fasta + "/concat.faa"
    # Dictionary to store unique sequences by identifier
    unique_sequences = {}

    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".faa"):
            file_path = os.path.join(directory, filename)
            
            # Parse the .faa file using SeqIO
            for record in SeqIO.parse(file_path, "fasta"):
                # Add the record to the dictionary if the identifier is not already present
                if record.id not in unique_sequences:
                    unique_sequences[record.id] = record

    # Write all unique sequences to the output file using SeqIO
    with open(output_file, 'w') as outfile:
        SeqIO.write(unique_sequences.values(), outfile, "fasta")

    print(f"Concatenated file created: {output_file}")


def filter_length_fasta(filepath,minimal,maximal):
    #22.11.22
    output = filepath+"_filter_"+str(minimal)+str(maximal)
    filter_writer = open(output,"w")

    for record in SeqIO.parse(filepath, "fasta"):
        sequence_length = len(str(record.seq))
        if sequence_length > minimal and sequence_length < maximal:
            SeqIO.write(record,filter_writer,"fasta")
    
    
    filter_writer.close()    

    return output























####################################################
##########    Concat alignments routines    ########
####################################################

#concats .fasta_aln files in a directory. if none are found
#it aligns with mafft and removes gaps with trimal at 95 % 
            
            
def concat_alignments(options):
    filepaths, output, trim_output = get_alignments(options)
    
    concat_dict = {}
    identifier_dict = {}
    alignment_lengths = {}
    removed_doublicates = set()
    
    # Specify the desired order of filepaths manually, if needed
    # Example: filepaths_sorted = ["domain1.fasta", "domain2.fasta", "domain3.fasta"]
    # Ensure that filepaths_sorted matches the order in which you want to concatenate the domains
    filepaths_sorted = sorted(filepaths, key=lambda x: x)  # Replace with custom order if needed

    print("Concatenation order of files:", filepaths_sorted)

    # Initialize concat_dict and identifier_dict (headers for the final print)
    for filepath in filepaths_sorted:
        genomeIDs_in_file = set()
        for record in SeqIO.parse(filepath, "fasta"):
            genomeID = record.id.split("-")[0]
            
            #Skip and remove the doublicates
            if genomeID in removed_doublicates:
                continue
            if genomeID in genomeIDs_in_file:
                #if double remove from the process and notify
                concat_dict.pop(genomeID)
                identifier_dict.pop(genomeID)
                removed_doublicates.add(genomeID)
                print(f"WARNING: Encountered doublicates of {genomeID} in {filepath}. Removed from concatenation.")
                continue
            else:
                genomeIDs_in_file.add(genomeID)
            
            concat_dict[genomeID] = ""    
                
            try:
                dom_type = record.description.split(' ')[1]
                dom_type = dom_type.split('_')[-1]
            except:
                dom_type = ''

            if genomeID in identifier_dict.keys():
                description = identifier_dict[genomeID]
                identifier_dict[genomeID] = description + '_' + dom_type
            else:
                identifier_dict[genomeID] = record.id + '_' + dom_type
                
            alignment_lengths[filepath] = len(str(record.seq))

    # Concatenate sequences in the specified order
    for filepath in filepaths_sorted:
        present_set = set()
        for record in SeqIO.parse(filepath, "fasta"):
            genomeID = record.id.split("-")[0]
            if genomeID in concat_dict:
                concat_dict[genomeID] += str(record.seq)
                present_set.add(genomeID)
            
        absent_genomes = set(concat_dict.keys()).difference(present_set)
        length = alignment_lengths[filepath]
        
        #all genomes that were not included get gaps instead or get removed
        for absent in absent_genomes:
            if options.gaps:
                print(f"WARNING: Missing sequence for {absent} in {filepath}. Concatenating gaps instead.")
                concat_dict[absent] += length * "-"
            else:
                print(f"WARNING: Missing sequence for {absent} in {filepath}. Removing concatenated sequence.")
                concat_dict.pop(absent)
    
    # Write the concatenated sequences to the output file
    with open(output, "w") as concat_file:
        for genomeID, sequence in concat_dict.items():
            header = identifier_dict[genomeID]
            concat_file.write(f">{header}\n")
            concat_file.write(f"{sequence}\n")

    print(f"\nWrote concatenated alignments to: {output}")

    # Trim the concatenated alignment using trimAl
    remove_gaps_with_trimal(output, trim_output)

    return
    
def get_alignments(options):
    """
    02.11.22
        Args:
            filepaths with the files to be processed
            directory for the concatenated seqs file
        Alignment files should be concated if header is the same            
    """

    alignment_files = myUtil.getAllFiles(options.concat_alignment,".fasta_aln")
    output = options.concat_alignment+"/concat.fasta_aln"
    trim_output = options.concat_alignment+"/trimmed_concat.fasta_aln"
    
    #If no alignments were found check for faa files that can be aligned
    if not alignment_files:
        print(f"WARNING: There were no alignments with .fasta_aln ending found in {options.concat_alignment}")
        fasta_files = myUtil.getAllFiles(options.concat_alignment,".faa")
        if not fasta_files:
            print(f"ERROR: There were no .faa fasta files in {options.concat_alignment}")
            return
        #align each .faa file
        for fasta_file in fasta_files:
            input_dir = os.path.dirname(fasta_file)
            base_name = os.path.splitext(os.path.basename(fasta_file))[0] # Extract the filename without the .faa extension
            output_fasta = os.path.join(input_dir, f"{base_name}.fasta_aln") # Create the output filename with .aln extension in the same directory as the input file
            print(fasta_file,output_fasta)
            #Align with default mafft
            align_fasta_with_mafft(fasta_file, output_fasta)


    alignment_files = myUtil.getAllFiles(options.concat_alignment,".fasta_aln")            
    return alignment_files,output,trim_output




def align_fasta_with_mafft(input_fasta, output_fasta):
    """
    Aligns a given .faa FASTA file using MAFFT.
    
    Args:
        input_fasta: Path to the input FASTA (.faa) file.
        output_fasta: Path to the output aligned FASTA file.
    """
    # Find MAFFT executable
    mafft = myUtil.find_executable("mafft")  # Ensure this function finds the MAFFT executable

    # Run MAFFT alignment
    try:
        with open(output_fasta, "w") as output_file:
            # Pass the file object to stdout and stderr
            subprocess.run([mafft, "--thread", "2", "--auto", input_fasta], stdout=output_file, stderr=subprocess.PIPE, check=True)
        print(f"Alignment complete: {input_fasta} -> {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during MAFFT alignment: {e.stderr.decode('utf-8')}")
        raise
    except FileNotFoundError:
        print("MAFFT executable not found. Please make sure MAFFT is installed and in your system's PATH.")



def remove_gaps_with_trimal(input_fasta, output_alignment, gap_threshold=0.95):
    """
    Remove columns with gaps using trimAl based on the specified gap threshold.

    Args:
        input_fasta: Path to the input FASTA file (aligned).
        output_alignment: Path to the output trimmed alignment file.
        gap_threshold: Proportion of gaps allowed in a column (default: 0.95).
    """
    trimal = myUtil.find_executable("trimal")
    try:
        # Run the trimAl command with the gap threshold
        subprocess.run([
            trimal, 
            "-in", input_fasta, 
            "-out", output_alignment, 
            "-gt", str(gap_threshold),
            "-keepheader"
        ], check=True)
        #print(f"Trimming complete: {input_fasta} -> {output_alignment}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during trimming with trimAl: {e.stderr.decode('utf-8')}")
        raise












###################################################################
##########    Add taxonomy to fasta/alignments routines    ########
###################################################################

def taxonomy_comprehension(options):
    if os.path.isdir(options.add_taxonomy): 
        #if directory was provided
        #name all sequences and concat all files
        fasta_files = myUtil.getAllFiles(options.add_taxonomy,".faa")
        for fasta in fasta_files:
            add_taxonomy(options.database_directory,fasta)

        fasta_files = myUtil.getAllFiles(options.add_taxonomy,".fasta_aln")
        for fasta in fasta_files:
            add_taxonomy(options.database_directory,fasta)        
            
    elif os.path.isfile(options.add_taxonomy):
        #if single file was provided
        #name all sequences
        add_taxonomy(options.database_directory,options.add_taxonomy,options.dataset_divide_sign)
            

def add_taxonomy(database, filepath, trennzeichen=';'):
    """
    12.03.23
    Taxonomy lineage and domain type should be added to a FASTA file and replace the
    previous headers. 'trennzeichen' is the character used to separate different
    pieces of information in the header for better readability.
    
    IMPORTANT: NAMING OF DOUBLICATES HAS TO BE THE SAME AS FOR iTol PROTEIN DATASETS
    """
    print("Collecting genome identifiers")

    # Collect all genomeIDs from the FASTA file
    genomeID_list = []
    record_dict = {}
    sequences = []

    # First pass: collect all genomeIDs and sequences
    for record in SeqIO.parse(filepath, "fasta"):
        genomeID, proteinID = record.id.split("-", maxsplit=1)
        genomeID = myUtil.getGenomeID(genomeID)
        genomeID_list.append(genomeID)
        
        # Check if genome has multiple proteins in the FASTA
        if genomeID not in record_dict:
            record_dict[genomeID] = 1
        else:
            record_dict[genomeID] += 1
        
        sequences.append((record, genomeID))  # Store record and genomeID for later use

    # Fetch all taxonomic data in a single query
    taxon_dict = {}
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        placeholder = ",".join("?" for _ in genomeID_list)
        query = f"""SELECT genomeID, Superkingdom, Phylum, Class, Ordnung, Family, Genus, Species 
                    FROM Genomes 
                    WHERE genomeID IN ({placeholder})"""
        cur.execute(query, genomeID_list)
        rows = cur.fetchall()
        for row in rows:
            taxon_dict[row[0]] = row  # Store taxonomy data keyed by genomeID

    # Second pass: write sequences with updated headers
    with open(filepath + ".taxon_names", "w") as writer:
        for record, genomeID in sequences:
            if genomeID in taxon_dict:
                lineage = myUtil.taxonomy_lineage(taxon_dict[genomeID], trennzeichen)
                dataset_range_line = lineage
                if record_dict[genomeID] > 1:
                    dataset_range_line += f"_{record_dict[genomeID]}"
            else:
                dataset_range_line = record.id  # Use original ID if taxonomy is missing

            try:
                dom_type = record.description.split(' ')[1]
                dom_type = f" {dom_type} "
            except IndexError:
                dom_type = ''
            
            # Write the updated sequence header and sequence
            writer.write(f">{dataset_range_line}{dom_type}\n")
            writer.write(str(record.seq) + "\n")

    print(f"Taxonomy-added FASTA file written to {filepath}.taxon_names")
    return
    
    
##########################################################################
##########    Add genomic context to fasta/alignments routines    ########
##########################################################################
    
    
def add_genomic_context(database, filepath):
    """
    22.02.23
    This routine takes a FASTA file, iterates through the sequences, and
    writes down the specific genomic context in which this was found together with
    some genomic features.
    """
    cp_dict = {}
    print("Adding genomic context")

    # Collect protein IDs from the FASTA file
    protein_ids = [record.id.split('-', 1)[-1] for record in SeqIO.parse(filepath, "fasta")]

    # Fetch cluster IDs for all protein IDs in a single batch query
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = f"SELECT proteinID, clusterID FROM Proteins WHERE proteinID IN ({','.join(['?']*len(protein_ids))})"
        cur.execute(query, protein_ids)
        results = cur.fetchall()
        cp_dict = {protein_id: cluster_id for protein_id, cluster_id in results}

        cur.execute("""PRAGMA foreign_keys = ON;""")

    # Fetch all required information for the context in one go
    context_data = {}
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cluster_ids = list(cp_dict.values())
        if cluster_ids:
            query = f"""SELECT DISTINCT Proteins.proteinID, Proteins.genomeID, Proteins.clusterID, contig,
                        start, end, strand, sequence, domain, domStart, domEnd, score, dom_count 
                        FROM Proteins 
                        LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID 
                        WHERE Proteins.clusterID IN ({','.join(['?']*len(cluster_ids))})"""
            cur.execute(query, cluster_ids)
            rows = cur.fetchall()

            # Organize the fetched data
            for row in rows:
                proteinID, genomeID, clusterID, contig, start, end, strand, sequence, domain, domStart, domEnd, score, dom_count = row
                if proteinID not in context_data:
                    context_data[proteinID] = {
                        "genomeID": genomeID,
                        "clusterID": clusterID,
                        "contig": contig,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "sequence": sequence,
                        "domains": []
                    }
                if domain:
                    context_data[proteinID]["domains"].append((domain, domStart, domEnd, score))

    # Write the genomic context to the output file
    with open(filepath + "_gene_vicinity", "w") as writer:
        writer.write("Index\tProteinID\tDomain(s)\tHit_score\thit_align\tcontig\tstart\tend\tstrand\tSuperkingdom\tClade\tPhylum\tClass\tOrdnung\tFamily\tGenus\tSpecies\n")

        # Iterate over protein IDs and write the context data
        for proteinID, clusterID in cp_dict.items():
            protein_data = context_data.get(proteinID, {})
            genomeID = protein_data.get("genomeID", "")
            domains = ";".join([f"{d[0]}({d[1]}-{d[2]}):{d[3]}" for d in protein_data.get("domains", [])])
            taxon_data = fetch_taxonomy(database, genomeID)  # Fetch taxonomy data for genomeID

            writer.write(f"0\t{proteinID}\t{domains}\t{protein_data.get('sequence', '')}\t{protein_data.get('contig', '')}\t{protein_data.get('start', '')}\t{protein_data.get('end', '')}\t{protein_data.get('strand', '')}\t")
            writer.write("\t".join(taxon_data) + "\n")

    # Sorting the gene vicinity data
    sort_gene_vicinity(filepath + "_gene_vicinity", filepath + "_gene_vicinity_sorted")

    print("Genomic context added and saved.")


def fetch_taxonomy(database, genomeID):
    """
    Fetch taxonomy data for a given genomeID.
    """
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT Superkingdom, Clade, Phylum, Class, Ordnung, Family, Genus, Species FROM Genomes WHERE genomeID = ?", [genomeID])
        row = cur.fetchone()
        return row if row else [""] * 8  # Return empty fields if no taxonomy data is found


def sort_gene_vicinity(input_file, output_file):
    """
    Sort gene vicinity data by taxonomic ranks and save to a new file.
    """
    sort_headers = ['Superkingdom', 'Phylum', 'Class', 'Ordnung', 'Family', 'Genus', 'Species']

    # Read the input file and store the rows as dictionaries
    with open(input_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = [row for row in reader]

    # Sort the rows hierarchically based on the sort columns
    sorted_rows = sorted(rows, key=lambda x: tuple(x[header] for header in sort_headers))

    # Write the sorted rows to the output file
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(sorted_rows)

def deprecated_add_genomic_context(database,filepath):
    """
    22.02.23
        Args:
            database for the sequence file 
            filepath with the file to be processed
        This routine takes a fasta file, iterates through the sequences and
        writes down the specific genomic context in which this was found together with
        some genomic features
        
        nicht schön aber funktioniert vielleicht so
    """
    cp_dict = {}
    print("Adding genomic context")
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = "SELECT proteinID, clusterID FROM Proteins WHERE proteinID IN ({})"

        # Group protein IDs into batches of 1000 for efficient querying
        batch_size = 1000
        protein_ids = [f"'{record.id.split('-', 1)[-1]}'" for record in SeqIO.parse(filepath, "fasta")]

        for i in range(0, len(protein_ids), batch_size):
            batch = ','.join(protein_ids[i:i+batch_size])
            cur.execute(query.format(batch))
            results = cur.fetchall()
            for protein_id, cluster_id in results:
                cp_dict[protein_id] = cluster_id

        
        
        cur.execute("""PRAGMA foreign_keys = ON;""")
        with open(filepath+"_gene_vicinity", "a") as writer:
            writer.write("Index\tProteinID\tDomain(s)\tHit_score\thit_align\tcontig\tstart\tend\tstrand\tSuperkingdom\tClade\tPhylum\tClass\tOrdnung\tFamily\tGenus\tSpecies")
                    
        print(f"\tAssigning genetic environment:")
        for current_proteinID,current_clusterID in cp_dict.items():
            print(f"\tFetched: {current_proteinID} {current_clusterID}",end="\r")
            protein_dict = dict()
            fusion_protIDs = dict()
            query = "SELECT DISTINCT Proteins.proteinID,Proteins.genomeID,Proteins.clusterID,contig,start,end,strand,sequence,domain,domStart,domEnd,score,dom_count from Proteins LEFT JOIN Domains ON Proteins.proteinID = Domains.proteinID WHERE Proteins.clusterID = ?"
            cur.execute(query,[current_clusterID,])
            #Fill the protein dict
            
            for row in cur:
                # 0 => proteinID, 1 => genomeID, 2 => clusterID, 3 => contig,
                # 4 => start, 5 => end, 6 => strand, 7 => sequence,
                # 8 => domain, 9 => domStart, 10 => domEnd, 11 => score,
                if row[0] in protein_dict.keys():
                    protein = protein_dict[row[0]]
                    protein.add_domain(row[8],row[9],row[10],row[11])
                else:
                    protein = ParseReports.Protein(row[0],row[8],row[9],row[10],row[11])
                    protein.set_genomeID(row[1])
                    protein.set_clusterID(row[2])
                    protein.set_gene_contig(row[3])
                    protein.set_gene_start(row[4])
                    protein.set_gene_end(row[5])
                    protein.set_gene_strand(row[6])
                    protein_dict[row[0]] = protein

                if row[12] > 1:
                    fusion_protIDs[row[0]] = 1
                    
            #Also collect domains from fusion proteins
            
            if fusion_protIDs:
                count = 1        
                query = """SELECT DISTINCT proteinID,domain,domStart,domEnd,score FROM Domains WHERE proteinID = ? """
                add_array = []
                proteins_array = []
                args = []
                                    
                for proteinID in fusion_protIDs.keys():
                    arguments = args.copy()
                    arguments.insert(0,proteinID)
                    cur.execute(query,arguments)
                    for row in cur:
                        print(f"\tFetched fused domains: {count}",end="\r")
                        count = count + 1
                        protein = protein_dict[row[0]]
                        protein.add_domain(row[1],row[2],row[3],row[4])
                        
                print("")
            #an diesem punkt haben wir ein vollständiges unsortiertes protein_dict aber keine txonomy information
            # taxonomie info holen
            taxon_dict = dict()
            query = "SELECT genomeID,Superkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species FROM Genomes WHERE genomeID = ? "
            genomeID = current_proteinID.split('-')[0]
            cur.execute(query,[genomeID,])
            for row in cur:
                taxon_dict[row[0]] = row[1:]
            
            #create indices for the
            proteinID_list = sorted(protein_dict, key=lambda x: protein_dict[x].gene_start)
            try:
                index = proteinID_list.index(current_proteinID) * -1
            except:
                index = 0    
            with open(filepath+"_gene_vicinity", "a") as writer:
                for proteinID in proteinID_list:
                    #Write each line includes all information about one protein
                    protein = protein_dict[proteinID]
                    proteinlist = protein.get_protein_list()  #list representation of a protein object
            
                    out = str(index) + '\t' +'\t'.join(proteinlist) + '\t'.join(taxon_dict[genomeID]) + '\n'
                    index = index + 1
                    writer.write(out)
            
    
    # Specify input and output file names
    input_file = filepath+"_gene_vicinity"
    output_file = filepath+"_gene_vicinity_sorted"

    # Specify the headers of the columns to sort onSuperkingdom,Clade,Phylum,Class,Ordnung,Family,Genus,Species 
    sort_headers = ['Superkingdom', 'Phylum', 'Class', 'Ordnung', 'Family', 'Genus', 'Species']

    # Read in the input file and store the rows as a list of dictionaries
    with open(input_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = [row for row in reader]

    # Sort the rows hierarchically based on the sort columns
    sorted_rows = sorted(rows, key=lambda x: (x[sort_headers[0]], x[sort_headers[1]],x[sort_headers[2]],x[sort_headers[3]],x[sort_headers[4]],x[sort_headers[5]]))

    # Write the sorted rows to the output file
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        for row in sorted_rows:
            writer.writerow(row)



def help_routine_select_random_marker_set(database,filepath):

    print("Collecting genome identifier")
    level_dict = {}
    genomeID_dict = {}
    
        
    print(f"Fetching taxonomic data from {database}")
    #For performance reasons this is written here without a specific subroutine but execute shall be the same select as in Datasets
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        #cur.execute("""SELECT DISTINCT Phylum FROM Genomes WHERE NOT Phylum LIKE "-" AND NOT Phylum LIKE "_" OR NOT Phylum LIKE "1" OR NOT Phylum LIKE "0";""")
        phyla = ["Abyssubacteria","Acidobacteriota","Actinobacteriota","Aerophobota","Aquificota","Armatimonadota","Atribacterota","Aureabacteria","Bacteroidota","Bdellovibrionota","Bipolaricaulota","Caldisericota","Calditrichota","Calescibacterota","Campylobacterota","Chlamydiota","Chloroflexota","Chrysiogenota","Cloacimonadota","Coprothermobacterota","Cyanobacteria","Deferribacterota","Deinococcota","Delongbacteria","Dependentiae","Desantisbacteria","Desulfobacterota","Dictyoglomota","Dormibacterota","Edwardsbacteria","Eisenbacteria","Elusimicrobiota","Eremiobacterota","Fermentibacterota","Fibrobacterota","Firestonebacteria","Firmicutes","Fusobacteriota","Gemmatimonadota","Goldbacteria","Hydrogenedentota","Krumholzibacteriota","Latescibacterota","Lindowbacteria","Margulisbacteria","Marinisomatota","Mcinerneyibacteriota","Methylomirabilota","Moduliflexota","Muirbacteria","Myxococcota","Nitrospinota","Nitrospirota","Omnitrophota","Patescibacteria","Planctomycetota","Poribacteria","Proteobacteria","Ratteibacteria","Riflebacteria","Schekmanbacteria","Spirochaetota","Sumerlaeota","Synergistota","Tectomicrobia","Thermodesulfobiota","Thermosulfidibacterota","Thermotogota","Verrucomicrobiota","WOR-3","Wallbacteria","Zixibacteria"]
        
        for row in phyla:
            level_dict[row] = 1
                
        
        for level in level_dict.keys():
            cur.execute("""SELECT genomeID from Genomes WHERE Phylum = ? LIMIT 5""",[level])        
            for row in cur:
                genomeID_dict[row[0]] = 1
                        
    writer = open(filepath+".limited_phyla","w")
    for record in SeqIO.parse(filepath, "fasta"):
        genomeID = record.id.split("_Bacteria")[0]
        #genomeID = myUtil.getGenomeID(genomeID)
        if genomeID in genomeID_dict.keys():
            writer.write(">"+genomeID+"\n")
            writer.write(str(record.seq) + "\n")
        
    
    
    
    writer.close()


