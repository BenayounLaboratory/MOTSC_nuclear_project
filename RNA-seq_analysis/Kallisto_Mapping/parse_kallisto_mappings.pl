#! /usr/bin/perl

use warnings;
use strict;

#2016-09-22
# parse kallisto maps into one file

unless (@ARGV >= 2) {
	print "parse_kallisto_mappings.pl <outfile name> <kallistofile1> ... <kallistofileN>";

}

my $outname = shift @ARGV;

# to get gene length on the first file
my $start = 1;

my %outresults = ();

# prepare header
my $header = "GeneSymbol\tENS_GID\tENS_TID\teff_length\t";
$header .= join("\t",@ARGV);
$header .= "\n";

# loop on files
foreach my $file (@ARGV) {
	
	open(FILE,$file) or die "Could not open $file: $!\n";
	
	# skip first line
	my $headerline = <FILE>;
	
	if ($start == 1) {
		
		while (my $line = <FILE>) {
			#target_id	length	eff_length	est_counts	tpm
			#STPG1|ENSG00000001460|ENST00000003583	2544	2394.73	0	0
			
			my @linedata = get_line_data($line);
			#print "@linedata\n:";

			push(@{$outresults{$linedata[0]}}, ($linedata[2],$linedata[3]) );
		}
		
		$start = 0;
	
		
	} else {
	
		while (my $line = <FILE>) {
			my @linedata = get_line_data($line);
			push(@{$outresults{$linedata[0]}}, $linedata[3] );
		}
		
	}
	
	close FILE;
}



#### output results
open(OUT,'>',$outname) or die "Could not open $outname: $!\n";

print OUT $header;

foreach my $gene (sort keys %outresults) {
	
	#print $gene."\n";

	my @geneinfo =  split(/\|/, $gene);
	
	#print @geneinfo."\n";
	
	my $outline = $geneinfo[0]."\t".$geneinfo[1]."\t".$geneinfo[2]."\t".join("\t",@{$outresults{$gene}})."\n";
	print OUT $outline;
}


close OUT;

exit;

###########################################################
# SUBROUTINES
###########################################################

###########################################################
# a subroutine that separates fields from a data line and
# returns them in an array

sub get_line_data {

    my $line = $_[0];
    
    chomp $line;  

    my @linedata = split(/\t/, $line);
        
    return @linedata;
}
