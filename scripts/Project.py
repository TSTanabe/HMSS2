#!/usr/bin/python

import os
from datetime import datetime

def prepare_result_space(options,project="project"):
    """
        6.10.22
        Args:
            directory   directory for the result files to be stored
        Return:
            path        directory path for the result files
            continue    boolean if old project should be continues
        Options taken are
        Flowchart for this routine
            Checks if user input directory is already present
             if not 
                try to create this one
             if yes 
                check if database is already present
                 if no
                    create new one
                 if yes
                    continue writing on this existing database
            when database is created follow with the other report files        
            
    """
    
    logfile = 0
    now = datetime.now()
    timestamp = str(datetime.timestamp(now))
    
    if isProjectFolder(options): #check if existing project folder was given
        
        #print(options.result_files_directory)
        #print(options.fasta_initial_hit_directory)
        
        #Define directories from the existing project
        options.database_directory = options.result_files_directory+"/database.db"
        options.fasta_initial_hit_directory = options.result_files_directory+"/Hit_list"
        options.fasta_output_directory = options.result_files_directory+"/Sequences"
        options.Csb_directory = options.result_files_directory+"/Collinear_syntenic_blocks"
        
        
        #Define files from the existing project
        options.divergent_output_file = options.result_files_directory+"/div_output_file.faa"
        options.csb_output_file = options.Csb_directory+"/Csb_output.txt"
        options.gene_clusters_file = options.Csb_directory+"/All_gene_clusters.txt"
        
        if os.path.isfile(options.Csb_directory+"/csb_instances.json"):
            options.data_computed_Instances_json = options.Csb_directory+"/csb_instances.json"
        else:
            options.data_computed_Instances_json = None

        
        
        

    else: #if new project is created
        if not os.path.isdir(options.result_files_directory):
            try:
                os.mkdir(options.result_files_directory)
            except:
                raise Exception(f"\nERROR: No writing rights.")
        #Creates an new project and overwrites the result file directory with the project folder. All
        #other subfolders and files are stored in this project folder
        options.result_files_directory = create_project(options.result_files_directory,project)        

        options.database_directory = options.result_files_directory+"/database.db"
        

        options.fasta_initial_hit_directory = options.result_files_directory+"/Hit_list"
        os.mkdir(options.fasta_initial_hit_directory)
        options.fasta_output_directory = options.result_files_directory+"/Sequences"
        os.mkdir(options.fasta_output_directory)
        options.Csb_directory = options.result_files_directory+"/Collinear_syntenic_blocks"
        os.mkdir(options.Csb_directory)
        
        
        options.divergent_output_file = options.result_files_directory+"/div_output_file.faa"
        options.csb_output_file = options.Csb_directory+"/Csb_output.txt"
        options.gene_clusters_file = options.Csb_directory+"/All_gene_clusters.txt"
        options.data_computed_Instances_json = options.Csb_directory+"/csb_instances.json"

        


def create_project(directory, projectname="project"):
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d_%H-%M-%S")  # z. B. "2025-04-16_14-53-21"
    directory = os.path.join(directory, f"{timestamp}_{projectname}")
    
    try:
        os.mkdir(directory)
    except Exception:
        raise Exception("\nERROR: No writing rights.")
    
    return directory
    

def isProjectFolder(options):
    
    #Database may be provided, then consider this as a result directory
    if options.database_directory and os.path.isfile(options.database_directory):
        options.result_files_directory = os.path.dirname(options.database_directory)
    
    #Confirm that in the directory is everything
    try:
        if not os.path.isfile(options.result_files_directory+"/database.db"):
            print("No database file found. Creating new project folder.")
            return 0
        if not os.path.isdir(options.result_files_directory+"/Sequences"):
            print("No sequence directory found. Creating new project folder.")
            return 0
        if not os.path.isdir(options.result_files_directory+"/Hit_list"):
            print("No hit directory found. Creating new project folder.")
            return 0
        if not os.path.isdir(options.result_files_directory+"/Collinear_syntenic_blocks"):
            print("No collinear syntenic block directory found. Creating new project folder.")
            return 0
    except:
        raise Exception()
    return 1
    
    
def any_process_args_provided(args, default_values):
    for arg, default in default_values.items():
        if getattr(args, arg) != default:
            #print(f"{arg} was not {default} but {getattr(args, arg)}")
            return True
    return False
    

   
     






