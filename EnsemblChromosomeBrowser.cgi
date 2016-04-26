#!/usr/bin/perl

#Authored by Matthew Orton

#Note that you will need the use of a lib file with all of the relevant Ensembl modules to use this script.

#This program will take a chromosomal region specified by the user for the species Equus Caballus (horse) and
#assuming it passes a validation check will search for the relevant genes relating to ths region on Ensembl. If genes
#are found then it will show a table with important details on these genes including id, start pos, end pos, strand, length, description,
#external id, status and number of transcripts.
#It will also display a graphical view of the region specified by the user taken from the Ensembl website.
#In the graphical view protein coding genes will appear in red while non protein coding will appear in black.

use warnings;
use strict;
#using CGI module for the HTML code
use CGI qw/:standard/;
#using the EnsEMBL::Registry module to gain access to the registry
use Bio::EnsEMBL::Registry;
#For the graphical view we need these two packages
use Bio::Graphics;
use Bio::SeqFeature::Generic;

#starting a new registry variable for Ensembl
my $registry = 'Bio::EnsEMBL::Registry'; 
$registry->load_registry_from_db( 
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

#Defining start position, end position and chromosome chosen by user outside of
#any loop so they can be used in whichever subroutine or loop they need to be
my $startPosition;
my $endPosition;
my $chromosome;

#Grabbing adaptor for entire horse genome
my $horse_adaptor = $registry->get_adaptor( 'Equus Caballus', 'Core', 'Slice' );

#Then for each horse chromosome grabbing the slice specific to that chromosome
#first making an array with chromosomes of the horse genome
my @chr = qw/1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 X MT/;

my $i;
my @horseSlice;
#pushing horse adaptors for each specified chromosome to the horseslice array with a for loop
for($i=0; $i<scalar(@chr); $i++){
    push @horseSlice, $horse_adaptor->fetch_by_region( 'chromosome', "$chr[$i]");
}
#then determining end positions from these slices with another for loop
my @endPos;
for($i=0; $i<scalar(@horseSlice); $i++){
    push @endPos, $horseSlice[$i]->end();
}

#Defining a hash for all of the chromosomes and their end positions, this will allow me to use the end positions 
#as a validation check depending on which chromosome was chosen to make sure the end position is not greater than this 
my %horseChrEnd = ("1"  => "$endPos[0]",  "2"  => "$endPos[1]",  "3"  => "$endPos[2]",  "4" =>  "$endPos[3]",  "5"  => "$endPos[4]",
                   "6"  => "$endPos[5]",  "7"  => "$endPos[6]",  "8"  => "$endPos[7]",  "9" =>  "$endPos[8]",  "10" => "$endPos[9]",
                   "11" => "$endPos[10]", "12" => "$endPos[11]", "13" => "$endPos[12]", "14" => "$endPos[13]", "15" => "$endPos[14]",
                   "16" => "$endPos[15]", "17" => "$endPos[16]", "18" => "$endPos[17]", "19" => "$endPos[18]", "20" => "$endPos[19]",
                   "21" => "$endPos[20]", "22" => "$endPos[21]", "23" => "$endPos[22]", "24" => "$endPos[23]", "25" => "$endPos[24]",
                   "26" => "$endPos[25]", "27" => "$endPos[26]", "28" => "$endPos[27]", "29" => "$endPos[28]", "30" => "$endPos[29]",
                   "31" => "$endPos[30]", "X"  => "$endPos[31]", "MT" => "$endPos[32]");

#Then doing the validation checks for the chromosomal region entered
#The program must meet all of these checks to proceed on to the gene table view and graphical representation view
if (param()) {
    my @searchErrors;
    $startPosition = param("startposition");
    $endPosition = param("endposition");
    $chromosome = param("horsechr");
    if ($startPosition < 1) {
        push @searchErrors, "Error! Start position must be greater or equal to 1.";
    }
    #start cannot be greater or equal to end
    if ($startPosition >= $endPosition) {
        push @searchErrors, "Error! Start position must be less than end position.";
    }
    #end cannot be less than or equal to start
    if ($endPosition <= $startPosition) {
        push @searchErrors, "Error! End position must be greater than start position.";
    }
    #length must be smaller than 10e7
    if (($endPosition - $startPosition) > 10000000) {
        push @searchErrors, "Error! Chromosomal region must be less than or equal to 10e7 (10000000 bp) in length.";
    }
    #length must be greater than 1e3
    if (($endPosition - $startPosition) < 1000) {
        push @searchErrors, "Error! Chromosomal region must be greater than or equal to 1e3 (1000 bp) in length.";
    }
    #Limiting to integer values
    if ($startPosition !~ /^[0-9]{1,10}$/ && $endPosition !~ /^[0-9]{4,10}$/) {
        push @searchErrors, "Error! Chromosomal regions specified must be integer values.";
    }
    #Lastly if the endposition specified is greater than the length of that chromosome than an error will be shown
    if ($endPosition > $horseChrEnd{$chromosome}) {
        push @searchErrors, "Error! End position must be less than or equal to $horseChrEnd{$chromosome} for chromosome $chromosome.";
    }
    #If there are errors, then the search form will reappear with the errors contained in the searchErrors array
    #Errors will appear in red above the form elements to emphasize what went wrong and what the user has to fix.
    if (@searchErrors) {
        print "Content-type: text/html\n\n";
        print horsetop_html1();
        foreach (@searchErrors) {
            print "<font color=red>$_<br></font>";
        }
        print horsesearchform2();
        print horsebottom_html();
    #However if all validation checks are passed than a new page will appear with a table containing all of the genes contained
    #within the chromosome region specified by the user as well as columns for each attribute of the gene
    #below this table a graphical view of the chromosomal region will also appear
    } else {
        #Fetching entire region that was specified by the user, define this as regionChosen
        my $regionChosen = $horse_adaptor->fetch_by_region( 'chromosome', "$chromosome", $startPosition, $endPosition );
        
        #First we can do the graphical representation to keep it outside of the while loops below
        #Within the graphical view, protein coding genes will appear in red font while non protein coding will appear in black font
        #genes will also display based on their respective strands in the chromosome region
        #The graphical view will show exactly where each gene is located in relation to the region specified
        
        #First defining length for the graphic
        my $lengthGraphic = $endPosition - $startPosition;
        #Then we have to define a panel to contain the graphic, to do this we have to specify the start and end positions set by the user
        my $chrHorsePanel = Bio::Graphics::Panel->new(-length => $lengthGraphic, -width  => 800, -pad_left=>100, -pad_right=>100,
        -start=>$startPosition,-end=>$endPosition);
        #We also need a scale for the track of the graphic so this is what trackScale is defined as
        #again defined by start and end positions
        my $trackScale = Bio::SeqFeature::Generic->new(-start => $startPosition,-end => $endPosition);
        #Defining the appearance and properties of the track that will scale based on the region specified
        $chrHorsePanel->add_track($trackScale,-glyph => 'arrow',-tick => 3,-fgcolor => 'black',-double  => 2);
        #Then we can define our chrTrack for each gene within a foreach loop so each gene gets its own graphic
        foreach my $genesFound (@{$regionChosen->get_all_Genes()}) {
            #Each item is grabbed by its respective subroutine
            my $idString = idstring($genesFound);
            #start
            my $startPos = startstring($genesFound);
            #end
            my $endPos = endstring($genesFound);
            #strand
            my $strand = $genesFound->strand();
            #biotype for protein coding vs non
            my $biotype = $genesFound->biotype();
        #As mentioned, protein coding = red while non coding = black for font to emphasiz the difference, we can do this with an if loop
        #if biotype eq protein coding than red, if not then black
        my $chrTrack;
        if ($biotype eq "protein_coding") {
            $chrTrack = $chrHorsePanel->add_track(-glyph => 'transcript2',-label => 1, -fontcolor => 'red', 
            -bgcolor => 'blue', -description => "$biotype", -stranded => '1');
        } else {
            $chrTrack = $chrHorsePanel->add_track(-glyph => 'transcript2',-label => 1, -fontcolor => 'black', 
            -bgcolor => 'blue', -description => "$biotype", -stranded => '1'); 
        }
        #then we create an object to represent each gene feature for each gene, this includes, id, start, end and strand attributes
        my $horseGeneFeat = Bio::SeqFeature::Generic->new(-display_name => "$idString ($startPos-$endPos)", -start => $startPos, -end => $endPos, -strand => $strand);
        #horse genes found in the specified region are then added to the panel
        $chrTrack->add_feature($horseGeneFeat);
        }
        #Table printing
        #Grabbing all relevant genes to that position again for the table
        my $genesFound2 = $regionChosen->get_all_Genes();
        #dereference to check if genes exist in next step
        my @genesFound2 = @$genesFound2;
        #If no genes are present then will print top heading, table headings, message and link back to search form
        if (!@genesFound2) {
             print "Content-type: text/html\n\n";
             print horsetop_html3();
             print horseviewtable();
             print nogenes();
             print horselink();
             print horsebottom_html(); 
        } else {
            #Then printing HTML with viewtable form and graphic, report header is shown first, then table, then graphic, then link
            print "Content-type: text/html\n\n";
            print horsetop_html2();
            print horseviewtable();
            #Now for the table
            while ( my $horseGenes = shift @{$genesFound2} ) {
                #Grabbing all relevant gene data with subroutines created and $genesFound2
                #id
                my $idString = idstring($horseGenes);
                #start
                my $startPos = startstring($horseGenes);
                #end
                my $endPos = endstring($horseGenes);
                #strand
                my $strand = strandstring($horseGenes);
                #Also displaying strand as + or -
                if ($strand == 1) {
                    $strand = "+";
                } else {
                    $strand = "-";
                }
                #Length of gene
                my $length = ($endPos - $startPos);
                #Common identifier or external name for the gene
                my $external = $horseGenes->external_name();
                #biotype of the gene
                my $biotype = $horseGenes->biotype();
                #description of the gene
                my $desc = $horseGenes->description();
                #status of the gene
                my $status = $horseGenes->status();
              
                #Also grabbing transcript related data based on region chosen
                my $horseTranscripts = $horseGenes->get_all_Transcripts();
                while ( my $horseTranscript = shift @{$horseTranscripts} ) {
                    #Defining an array with all gene ids, a gene id will appear multiple times if it has more than
                    #one transcript
                    my @idString;
                    #Then using a hash I can count how many times a gene is duplicated based on how many transcripts it has
                    push @idString, $idString;
                    my %idCounts;
                    $idCounts{$_}++ for @idString;
                    #Then will actually print the table with all the relevant gene information
                    #Gene ids will also have the functionality of being a link to the Ensembl website where that specific id can be found and all of its available information
                    print "<tr><td><a href='http://uswest.ensembl.org/Gene/Summary?db=core;g=$idString' target='_blank'>$idString</td><td>$startPos</td><td>$endPos</td><td>$strand</td><td>$length</td><td>$desc</td><td>$external</td><td>$biotype</td><td>$status</td><td>$idCounts{$idString}</td></tr>";
                }
            }
            #Creating an actual png file of the chrPanel and all associated genes in this panel
            open FH, ">horseGenes.png" or die $!;
            print FH $chrHorsePanel->png;
            close FH;
            #Finally printing the graphical representation in another form below the table
            print horseviewgraphic();
            #then below this we print our link back to the searchform if the user wants to return to this
            print horselink();
            print horsebottom_html();
        }
    }
} else {
#Default form and HTML that presents itself before any user input is done
print "Content-type: text/html\n\n";
print horsetop_html1();
print horsesearchform();
print horsebottom_html();
}

#Subroutines for grabbing some of the associated gene data from the chromosome region specified
#Each item is its own subroutine so that each item can be individually inserted to the table
#Grabbing the gene id:
sub idstring {
    my $featureId = shift;
    my $id = $featureId->stable_id();
    return ( $id );
}

#Start position of gene
sub startstring {
    my $featureStart = shift;
    my $start = $featureStart->seq_region_start();
    return ( $start );
}

#End position of gene
sub endstring {
    my $featureEnd = shift;
    my $end = $featureEnd->seq_region_end();
    return ( $end );
}

#Strand of the gene
sub strandstring {
    my $featureStrand = shift;
    my $strand = $featureStrand->strand();
    return ( $strand );
}

#Subroutine section for all associated printed html
#top_html1 will have the main headers for the search form
#this also includes a title for the browser
sub horsetop_html1 {
   return<<TOP;
   <head>
        <title>Welcome to the Equus Caballus (Horse) Chromsome Browser</title>
            <h2>Welcome to the Equus Caballus (Horse) Chromsome Browser</h2>
                <h3>Search a Chromosome Region of the Equus Caballus Genome:</h3>
   </head>
TOP
}

#html top for the table and graphic page
#this also includes a title for the browser
sub horsetop_html2 {
   return<<TOP;
   <head>
        <title>Report for Chromosome Region $startPosition to $endPosition of Equus Caballus Chromosome $chromosome</title>
            <h4>Report for Chromosome Region $startPosition to $endPosition of Equus Caballus Chromosome $chromosome</h4>
            <h4>Showing All Genes</h4>
   </head>
TOP
}

#html top for view page if no genes since it wouldnt say: showing all genes
#this also includes a title for the browser
sub horsetop_html3 {
   return<<TOP;
   <head>
        <title>Report for Chromosome Region $startPosition to $endPosition of Equus Caballus Chromosome $chromosome</title>
            <h4>Report for Chromosome Region $startPosition to $endPosition of Equus Caballus Chromosome $chromosome</h4>
   </head>
TOP
}

#This sub will present the default form with all of the required fields for the chromosomal position entry
#this form will also show a default chromosome and start/end positions so that a user can simply search those to check functionality
#if errors are present, form will not reset but will instead retain users entries -see searchform2
#Default is MT chromosome from positions 1-6000
sub horsesearchform {
    return<<FORM;
    <form action="$0" method ="post"><br>
    <table border="2" width="95%">
        <tr><td align="center">Start Chromosome Position:</td><td align="center"><input type="text" name="startposition" value="1" style="width: 64px"></td>
        <td align="center">End Chromosome Position:</td><td align="center"><input type="text" name="endposition" value="6000" style="width: 65px"></td>
        <td align="center">Select a Chromosome: <select align="center" name="horsechr", size="1">
        <option value="1">1</option><option  value="2">2</option><option  value="3">3</option><option  value="4">4</option>
        <option value="5">5</option><option  value="6">6</option><option  value="7">7</option><option  value="8">8</option>
        <option value="9">9</option><option  value="10">10</option><option value="11">11</option><option value="12">12</option>
        <option value="13">13</option><option value="14">14</option><option value="15">15</option><option value="16">16</option>
        <option value="17">17</option><option value="18">18</option><option value="19">19</option><option value="20">20</option>
        <option value="21">21</option><option value="22">22</option><option value="23">23</option><option value="24">24</option>
        <option value="25">25</option><option value="26">26</option><option value="27">27</option><option value="28">28</option>
        <option value="29">29</option><option value="30">30</option><option value="31">31</option><option value="X">X</option>
        <option value="MT" selected>MT</option>
        </select></td>
        <td align="center" colspan="2.5"><input type="submit"></td></tr>
    </table>
    <h4>General Information</h4>
        <ul>
        <li>To use this tool, simply input a region of a chromosome you are interested in using the text fields and drop down menu.</li>
        <li>A new page is then generated showing a table with all of the genes associated with this region and their details.</li>
        <li>The tool will also generate a graphical view of the of the chromosomal region showing all genes in that region.</li>
        <li>The horse genome build used for this tool is from the Ensembl API and is called EquCab2.</li>
        <li>This build is a Whole Genome Shotgun (WGS) assembly at 6.79x and was sequenced from a female thoroughbred named "Twilight".</li>

        </ul>
    <h4>Tips on searching chromosome regions:</h4>
        <ul>
        <li>Chromosome regions must be within the start and end positions of a specified chromosome.</li>
        <li>Chromosome regions must be between 1e3 (1000 bp) and 10e7 (10000000 bp) in length and positions must be integers.</li>
        <li>Also keep in mind that the mitochondrial chromosome is very small: 16,660 bp if searching it.</li>
        </ul>
    <h4>Example regions with at least 3 genes:</h4>                                  
        <ul>
        <li>Chromosome MT between positions 1 and 6000.</li>
        <li>Chromosome 29 between positions 270 to 975000.</li>
        <li>Chromosome 16 between positions 125000 and 500000.</li>
        </ul>
    <h4>Example region with no genes:</h4>
        <ul>
        <li>Chromosome 1 between positions 450 and 4500.</li>
        </ul>
    </form>
FORM
}
#Searchform2 will be presented if there is an error with the user input to the form
#the user input will be retained in this form if there is an error since we use param to grab the inputted values
sub horsesearchform2 {
    my $s = param('startposition');
    my $e = param('endposition');
    my $h = param('horsechr');
    return<<FORM;
    <form action="$0" method ="post"><br>
    <table border="2" width="95%">
        <tr><td align="center">Start Chromosome Position:</td><td align="center"><input type="text" name="startposition" value="$s" style="width: 64px"></td>
        <td align="center">End Chromosome Position:</td><td align="center"><input type="text" name="endposition" value="$e" style="width: 65px"></td>
        <td align="center">Select a Chromosome: <select align="center" name="horsechr", size="1">
        <option value="$h" selected>$h</option>
        <option value="1">1</option><option  value="2">2</option><option  value="3">3</option><option  value="4">4</option>
        <option value="5">5</option><option  value="6">6</option><option  value="7">7</option><option  value="8">8</option>
        <option value="9">9</option><option  value="10">10</option><option value="11">11</option><option value="12">12</option>
        <option value="13">13</option><option value="14">14</option><option value="15">15</option><option value="16">16</option>
        <option value="17">17</option><option value="18">18</option><option value="19">19</option><option value="20">20</option>
        <option value="21">21</option><option value="22">22</option><option value="23">23</option><option value="24">24</option>
        <option value="25">25</option><option value="26">26</option><option value="27">27</option><option value="28">28</option>
        <option value="29">29</option><option value="30">30</option><option value="31">31</option><option value="X">X</option>
        <option value="MT">MT</option>
        </select></td>
        <td align="center" colspan="2.5"><input type="submit"></td></tr>
    </table>
    <h4>General Information</h4>
        <ul>
        <li>To use this tool, simply input a region of a chromosome you are interested in using the text fields and drop down menu.</li>
        <li>A new page is then generated showing a table with all of the genes associated with this region and their details.</li>
        <li>The tool will also generate a graphical view of the of the chromosomal region showing all genes in that region.</li>
        <li>The horse genome build used for this tool is from the Ensembl API and is called EquCab2.</li>
        <li>This build is a Whole Genome Shotgun (WGS) assembly at 6.79x and was sequenced from a female thoroughbred named "Twilight".</li>
        </ul>
    <h4>Tips on searching chromosome regions:</h4>
        <ul>
        <li>Chromosome regions must be within the start and end positions of a specified chromosome.</li>
        <li>Chromosome regions must be between 1e3 (1000 bp) and 10e7 (10000000 bp) in length and positions must be integers.</li>
        <li>Also keep in mind that the mitochondrial chromosome is very small: 16,660 bp if searching it.</li>
        </ul>
    <h4>Example regions with at least 3 genes:</h4>
        <ul>
        <li>Chromosome MT between positions 1 and 6000.</li>
        <li>Chromosome 29 between positions 270 to 975000.</li>
        <li>Chromosome 16 between positions 125000 and 500000.</li>
        </ul>
    <h4>Example region with no genes:</h4>
        <ul>
        <li>Chromosome 1 between positions 450 and 4500.</li>
        </ul>
    </form>
FORM
}
    
#subroutine for the printing of the table
sub horseviewtable{
    return<<FORM;
    <table border="1" CELLSPACING="0" width="80%">
        <tr><th>Gene id</th>
            <th>Start</th>
            <th>End</th>
            <th>Strand</th>
            <th>Length</th>
            <th>Description</th>
            <th>External Name</th>
            <th>Gene Type</th>
            <th>Status</th>
            <th>Number of Transcripts</th></tr>
FORM
}

#separate html for the graphic in order to print below the table and also print a margin inbetween
sub horseviewgraphic{
    return<<FORM;
    </table>
    <br />
    <table border="1" width="90%" style="margin-top:-2;">
        <tr align=center><td><img src='horseGenes.png'/></td></tr>
    </table>
FORM
}

#table shown if no genes present
sub nogenes{
    return<<FORM;
    <table border="1" width="80%">
        <tr align=center><td>No genes could be found.</td></tr>
    </table>
FORM
}

#will return you to the searchform page
sub horselink{
    return<<FORM;
        <table border="0" width="90%">
            <tr align=center><td><b><a href="http://zenit.senecac.on.ca/~bif724_161a17/matthew_a3.cgi">New Search</a></b></td></tr>
        </table>
FORM
}

#bottom html, same for both the search form and table/graphic form
sub horsebottom_html {
    return<<BOTTOM;
BOTTOM
}
1;


