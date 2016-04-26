#!/usr/bin/perl

#Matthew Orton

#The purpose of this program is to parse accession numbers in a txt file specified by the user and access genbank for these numbers
#and their associated species and sequence information. It will also output in two separate files the relevant binomial name and dna sequence
#or binomial name and amino acid sequence of a gene specified by the user

use strict;
use warnings;
#Using the Bioperl SeqIO module for outputting to fasta format
use Bio::SeqIO;
#Also using Bioperl DB::GenBank module for accessing of the genbank files we need
use Bio::DB::GenBank;
#To ouput the file path where the perl file is located (should be the directory where the input file was run)
use Cwd 'abs_path';
#Creating a new genbank object $db_o
my $db_o = Bio::DB::GenBank->new;

#Defining both of our command line arguments: filename and genename
my ($filename, $genename) = @ARGV;

#Checking to make sure both a filename and a gene name have been entered as command line arguments
#If one of these has not been entered then an error will be presented to the user
if (not defined $filename) {
  die "Error! Please enter a file name\n";
}
if (not defined $genename) {
  die "Error! Please enter a gene name\n";
}

#Opening the file defined by filename and pushing each line of that file to an array called @filelines
#Will present an error if the file cannot be opened
open(F, "<", $filename) or die("Error! Cannot open the file: $!\n");
my @filelines = ();
while(<F>) { chomp; push(@filelines, $_); }

#Defining the names of the files for both dna and amino acid sequences
#Inserting the variable name $genename instead of hardcoding a genename into the output
#Should also output to the current directory
my $outfileDNA = "dna_matthew_$genename.fa";
my $outfileAA = "aa_matthew_$genename.fa";

#Defining the output file in a new SeqIO object for the DNA sequences and specifying fasta format
#Will output to current working directory
my $seqoutDNA = Bio::SeqIO->new(
                             -file => ">$outfileDNA",
                             -format => 'Fasta',
                             );
#Defining the outfile file in a new SeqIO object for the amino acid sequences and specifying fasta format
#Will output to current working directory
my $seqoutAA= Bio::SeqIO->new(
                             -file => ">$outfileAA",
                             -format => 'Fasta',
                             );
#Starting with a for loop to interate through each accession id within the input file
#for loop will iterate as long as less than scalar of our filelines array
my $i;
for($i=0; $i<scalar(@filelines); $i++){
    #Defining our first seq object which will use the Genbank package to allow us to access
    #the associated genbank data for each one of the accession numbers in our file
    my $seq_o = $db_o->get_Seq_by_acc($filelines[$i]);
    #Creating a species object to grab the species data for the associated gene
    my $species_o = $seq_o->species;
    #Grabbing the specific node_name of this species object to give us our binomial name
    my $species_name = $species_o->node_name;
    #Replacing spaces with underscores to the binomial names using the substitute operator
    $species_name =~s/ /_/g;
    #Grabbing SeqFeatures of our seq_o object and storing in a variable called $object
    for my $object ($seq_o->get_SeqFeatures) {
        #Making sure the primary tag equals CDS
        if ($object->primary_tag eq "CDS") {
            #Then checking if a gene tag is contained within the CDS of the given accession number
            if ($object->has_tag('gene')) {
                #If there is a valid gene tag, then we define a new variable known as $gene that will store the values of this tag
                for my $gene ($object->get_tag_values('gene')){
                    #If the gene tag values match with the specified genename in the command line argument, then further information will be gathered
                    if ($gene eq $genename){
                      #To grab the DNA sequence I chose to use the spliced_seq function in bioperl
                      #This function will also grab the dna sequence of the specified gene
                      #It will also grab the sequence in the correct orientation if the sequence is located on the reverse strand
                      #I tested to make sure of this by using the gene F19K19.1 of accession AC011808 which is -1 strand gene
                        for my $dnaSeq ($object->spliced_seq->seq){
                            my $seq_o2 = Bio::Seq->new (-id => $species_name,
                                                        -seq => $dnaSeq);
                            $seqoutDNA->write_seq($seq_o2);
                        }
                        #For amino acid sequences, first we check for the translation tag,
                        if ($object->has_tag('translation')) {
                            #If our specified gene has a translation tag for the current accession number,
                            #then the associated amino acid sequence will be added to $aaSeq
                            for my $aaSeq ($object->get_tag_values('translation')){
                                #Once again, another seq object is created,
                                #This time with our species name and amino acid sequence
                                #However, this amino acid sequence will be written to a different file for amino acids only
                                my $seq_o3 = Bio::Seq->new (-id => $species_name,
                                                            -seq => $aaSeq);
                                $seqoutAA->write_seq($seq_o3);
                            }
                        }
                    }
                }
            }
        }
    }
}

#Outputting a message if the seqOutDNA varibale is defined
#If it is, then will output a message indicating that the what the file is and the path of where it is located
if ($seqoutDNA) {
  my $abs_path = abs_path($filename);
  $abs_path =~ s/$filename//;
  print "The file $outfileDNA has now been deposited to the path: $abs_path\n";
}

#Outputting a message if the seqOutAA variable is defined
#If it is, then will output a message indicating that the what the file is and the path of where it is located
if ($seqoutDNA) {
  my $abs_path = abs_path($filename);
  $abs_path =~ s/$filename//;
  print "The file $outfileAA has now been deposited to the path: $abs_path\n";
}

