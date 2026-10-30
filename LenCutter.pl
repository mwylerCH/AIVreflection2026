#!/usr/bin/perl

# Author Michele Wyler, IVI Mittelhäusern, Switzerland


use warnings;
use strict;
use English;
use File::Temp qw/ tempdir /;
use Getopt::Long 'HelpMessage';
use File::Basename;
use Cwd;
use POSIX;

# script to clean GISAID Influenza fastas
# length filtering +/- 20%

my $FASTA = $ARGV[0];

### Reference  -------------------------------------------------------------

my $TOLLERANCE = 0.2;

my %REF = (
'PB2' => 2340,
'PB1' => 2340,
'PA' => 2230,
'HA' => 1770,
'NP' => 1560,
'NA' => 1460,
'MP' => 1030,
'NS' => 890,
);


### Read Fasta ------------------------------------------------------------

# read in fasta file, check header and put into temp folder
open(IN, $FASTA) or die "can't open $FASTA";
my %seqs = ();
my $header = '';

while (my $line = <IN>){
    chomp $line;
	$line =~ s/^\s*$//;
    if ($line =~ m/^>(.*)$/){
        $header = $line;
       # Check if segments are annotated (as '|1|)
        if ($header !~ m/\|[[:digit:]]\|/){
                print "How are each segment defined? Not with '|number|'.\n";
                exit;
        }
    } else {
        $seqs{"$header"} .= $line;
    }
}
close (IN);



### Header Check  ------------------------------------------------------------
my @regex = keys %REF;

# test first if header name of input contains '|PB|'
while ( my ($header,$sequence) = each %seqs ) {
	if( $header !~ m/\|PB2\|/ ){
		if( $header !~ m/\|PB1\|/ ){
			if( $header !~ m/\|PA\|/ ){
				if( $header !~ m/\|HA\|/ ){
					if( $header !~ m/\|NP\|/ ){
						if( $header !~ m/\|NA\|/ ){
							if( $header !~ m/\|MP\|/ ){
								if( $header !~ m/NS/ ){
									print "Use '|' as delimiter (for PB2, PB3,etc) in $header\n";
									exit;								
								}
							}
						}
					}
				}
			}
		}
	}

}

### Filtering  ------------------------------------------------------------
	
while ( my ($header,$sequence) = each %seqs ) {
        # get length of sequence
        my $FASTAlen = length($seqs{"$header"});
	# Test Len
	while (my ($NAME, $LEN) = each %REF ){
		# test which reference should be used
		if ($header =~ m/\|$NAME\|/) {
			# test the length
	               if ($FASTAlen >= 0.8*$LEN && $FASTAlen <= 1.2*$LEN) {
				# print out new fasta
	                	print "$header\n$sequence\n";
	               }

		}
	}

}
