#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use File::Copy;
use File::Find;
use Getopt::Std;
use Cwd qw(getcwd);
use Error qw(:try);

use Bio::DB::Taxonomy;



my %opts;
getopts('d:f:t:', \%opts);


##################################################################
################ MAIN BODY #######################################
##################################################################
#
#
#

#Globale Variablen
my $directory = $opts{d}; #Directory of the assembly statistic files
my $taxdump = $opts{t};
my $output_file = $opts{f}; #Path to the output file
my %assembly_to_taxid;
my %assembly_stats;
my %taxonid;

&getGenomesAssembly();
&getGenomesAssembly2();


##################################################################
################ MAIN BODY #######################################
##################################################################
#
#
#

sub getGenomesAssembly{
##
#Collect assembly statistics and taxonid
#Input directory with genome assembly statistics
##
my $assemblyFiles = &getallFiles($directory,"assembly_(report|stats).txt");

	for my $file (@$assemblyFiles){
	#print "Assembly stats from $file\n";
	my ($genomeID,$Biosample,$Bioproject,$Refseq,$Genbank,$Taxid) = (0,0,0,0,0,0,0);


	open my $fh, "<",$file or die "$!";
		while(<$fh>){
			#print $_,"\n";
			if($_ =~ /Taxid:\s+(\d+)\s/){$Taxid=$1;}
			elsif($_ =~ /BioSample:\s+(SAMN\d+)\s/){$Biosample=$1;}
			elsif($_ =~ /BioProject:\s+(\w{5}\d+)\s/){$Bioproject=$1;}
			#elsif($_ =~ /Organism name:\s+(.*)\n/){$Organism=$1;}
			elsif($_ =~ /GenBank assembly accession:\s+(GCA\S+)\s/){$Genbank=$1;}
			elsif($_ =~ /RefSeq assembly accession:\s+(GCF\S+)\s/){$Refseq=$1;}
			elsif($Taxid and $Biosample and $Bioproject and $Refseq){last;}
		}
	close $fh;

	$taxonid{$Taxid} = "1";
	$assembly_to_taxid{$file} = $Taxid;
	$assembly_stats{$file} = "$Taxid\t$Biosample\t$Bioproject\t$Genbank\t$Refseq";
	}



##
#Collect the lineage with the taxonID
##
my ($num,@taxonid);

#print "Starting taxonomy search\n";
my $cwd = getcwd;
my $taxonDB = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => "$taxdump", -nodesfile => "$taxdump/nodes.dmp" , -namesfile => "$taxdump/names.dmp");
#$num = $taxonDB->get_num_taxa();
#print "There are $num taxa stores in the TaxonDB $taxonDB\n";
while (my ($key, $value) = each(%taxonid)) {

	my %lineage;
	my $taxonid = $key;
	my $taxon;
	#print "Retrieving lineage for taxon id $taxonid\n";
	#Retrieve lineage of Phylogeny with rank from flatfileDB
		do{
		
		$taxon = $taxonDB->get_taxon(-taxonid => $taxonid);
			if(defined $taxon){
			$lineage{$taxon->rank} = $taxon->scientific_name;
			$taxonid = $taxon->parent_id;
			}else{
			$taxonid{$key} = "\t\t\t\t\t\t\t\t\n"; # If taxid id completely undefined
			next;
			}
		}until($taxon->rank eq "superkingdom");

	#Write correct ranks into Database
	#	strain species genus family order class phylum kingdom superkingdom
	
	
	my $species = $lineage{'species'};
	unless($lineage{"strain"}){
	$lineage{"strain"} = "";
	}elsif($species){
	$lineage{"strain"} =~ s/$species //g;
	}
	#print "Concating lineage information\n\n";
	my @ranks = ("superkingdom","clade","phylum","class","order","family","genus","species","strain");
	my $string = "";
	for my $i (0..$#ranks){
		my $rank = $ranks[$i];
		if(defined $lineage{$rank}){
			#print "$rank defined \n";
			$string .= $lineage{$rank};
			$string .= "\t";
		}else{
			#print "$rank undefiend \n";
			$string .= "\t";
		}
	}
	$taxonid{$key} = $string
}
#print "Finished taxon assignment\n";

##
#Write results to file
##

open(my $fh2, '>>', $output_file);
while (my ($file, $Taxid) = each(%assembly_to_taxid)) {

print $fh2 $file."\t".$taxonid{$Taxid}.$assembly_stats{$file}."\t\t\t\t\n";                 
#filename superkingdom clade phylum class order family genus species strain taxid biosample bioproject genbank refseq completeness contamination typestrain

}
close $fh2;

}





sub getGenomesAssembly2{
##
#Collect assembly statistics and taxonid
#Input directory with genome assembly statistics
##
my $assemblyFiles = &getallFiles($directory,"data_summary.tsv");

	for my $file (@$assemblyFiles){
	#print "Assembly stats from $file\n";
	my ($genomeID,$Biosample,$Bioproject,$Refseq,$Genbank,$Taxid) = (0,0,0,0,0,0,0);


	
	open my $fh, '<', $file or die "Could not open file '$file' $!";

	my @headers = split(/\t/, <$fh>);

	while (my $line = <$fh>) {
	  chomp $line;
	  my %data;
	  my @values = split(/\t/, $line);
	  for (my $i = 0; $i < scalar(@headers); $i++) {
	    $data{$headers[$i]} = $values[$i];
	  }
	$Taxid = $data{"Taxonomy id"};
	$genomeID = $data{"Assembly Accession"};
	$taxonid{$Taxid} = "1";
	$assembly_to_taxid{$genomeID} = $Taxid;
	#$assembly_stats{$file} = "$Taxid\t".$data{"BioSample"}."\t".$data{"BioProject"}."\t$genomeID\t$genomeID";
	$assembly_stats{$genomeID} = "$Taxid\t \t \t$genomeID\t$genomeID";
	}

	close $fh;
	}

##
#Collect the lineage with the taxonID
##
my ($num,@taxonid);

#print "Starting taxonomy search\n";
my $cwd = getcwd;
my $taxonDB = Bio::DB::Taxonomy->new(-source => 'flatfile', -directory => "$taxdump", -nodesfile => "$taxdump/nodes.dmp" , -namesfile => "$taxdump/names.dmp");
#$num = $taxonDB->get_num_taxa();
#print "There are $num taxa stores in the TaxonDB $taxonDB\n";
while (my ($key, $value) = each(%taxonid)) {

	my %lineage;
	my $taxonid = $key;
	my $taxon;
	#print "Retrieving lineage for taxon id $taxonid\n";
	#Retrieve lineage of Phylogeny with rank from flatfileDB
		do{
		
		$taxon = $taxonDB->get_taxon(-taxonid => $taxonid);
			if(defined $taxon){
			$lineage{$taxon->rank} = $taxon->scientific_name;
			$taxonid = $taxon->parent_id;
			}else{
			$taxonid{$key} = "\t\t\t\t\t\t\t\t\n"; # If taxid id completely undefined
			next;
			}
		}until($taxon->rank eq "superkingdom");

	#Write correct ranks into Database
	#	strain species genus family order class phylum kingdom superkingdom
	
	
	my $species = $lineage{'species'};
	unless($lineage{"strain"}){
	$lineage{"strain"} = "";
	}elsif($species){
	$lineage{"strain"} =~ s/$species //g;
	}
	#print "Concating lineage information\n\n";
	my @ranks = ("superkingdom","clade","phylum","class","order","family","genus","species","strain");
	my $string = "";
	for my $i (0..$#ranks){
		my $rank = $ranks[$i];
		if(defined $lineage{$rank}){
			#print "$rank defined \n";
			$string .= $lineage{$rank};
			$string .= "\t";
		}else{
			#print "$rank undefiend \n";
			$string .= "\t";
		}
	}
	#print($string,"\n");
	$taxonid{$key} = $string
}
#print "Finished taxon assignment\n";

##
#Write results to file
##

open(my $fh2, '>>', $output_file);
while (my ($file, $Taxid) = each(%assembly_to_taxid)) {

print $fh2 $file."\t".$taxonid{$Taxid}.$assembly_stats{$file}."\t\t\t\t\n";                 
#filename superkingdom clade phylum class order family genus species strain taxid biosample bioproject genbank refseq completeness contamination typestrain

}
close $fh2;

}










sub getallFiles{
#Get Directory and Pattern for Files to search for recursively
my $directory = shift;
my @patterns = @_;
my @files;
my $pattern = join("|",@patterns);
	find ( sub {
		
		return unless -f;       #Must be a file
		return unless /$pattern/;  #Must match regex
		push @files, $File::Find::name;
	}, $directory );
	
return \@files;
}
