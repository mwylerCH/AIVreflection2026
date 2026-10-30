#!/usr/bin/perl

use warnings;
use strict;
use English;
use File::Temp qw/ tempdir /;
use Getopt::Long 'HelpMessage';
use File::Basename;
use Cwd;
use POSIX;

# screen fastas for strange content

my $FASTA = $ARGV[0];
#my $FASTA="/home/mwyler/tempInfluenza/Chicken/bigH5/pezzetti/ultimi/dieci/xal";


## Read fasta --------------------------------
# unzip if needed
if ($FASTA =~ m/.gz$/) {
    open(IN, "gunzip -c $FASTA |") or die "can't open $FASTA";
}
else {
    open(IN, $FASTA) or die "can't open $FASTA";
}

## read input sequence inside a hash
my %seqs = ();
my $header = '';

while (my $line = <IN>){
    # Strip any trailing CR and LF line endings safely
    $line =~ s/[\r\n]+$//;
    # Skip empty lines
    $line =~ s/^\s*$//;
    
    if ($line =~ m/^>(.*)$/){
        # test if not duplicated (remove from hash if already present)
        if (exists($seqs{"$line"})) {
                delete($seqs{"$line"});
        }
        $header = $line;
    } else {
        # Strip internal/trailing whitespaces from the sequence line
        $line =~ s/\s+//g;
        $seqs{"$header"} .= $line;
    }
}
close (IN);

## control sequence ------------------------
# and print out

# check only IUPAC code (Commas removed so they are not treated as literal characters)
while ( my ($header,$sequence) = each %seqs ) {
        if ($sequence =~ m/[^ACGTURYSWKMBDHVN]/i ){
                delete($seqs{"$header"});
        } else {
                print "$header\n$sequence\n";
        }
}
