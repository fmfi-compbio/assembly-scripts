#! /usr/bin/perl -w

use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin qw($Bin);  #$Bin is the directory with the script
use lib $Bin;        #add bin to the library path
use shared;
use POSIX;

my $USAGE = "$0 [options] k < fasta > kmers

Print canonical versions of all k-mers in a set fo sequences.
Argument k is the length of the k-mer, preferable odd number.
K-mers containing N are skipped, only A,C,G,T,N are allowed in sequences.

-c  consider sequences in fasta as circular (e.g. mtDNA or tandem repeat motifs)

Output format

default: just list of k-mers
-e  extended: canonical k-mer, original kmer, chromosome, start pos
-b  write kmers in bed format with middle position, valid only for linear
";

my %Options;
getopts('cbe', \%Options);

die $USAGE  unless @ARGV==1;
my ($K) = @ARGV;
die $USAGE unless $K=~/^[0-9]+$/ && $K>0;
die $USAGE if exists $Options{'e'} && exists $Options{'b'};

while(1) {
    my ($name, $seq) = read_fasta(\*STDIN);
    last unless defined $name;

    $name = normalize_fasta_name($name);
    $name =~ s/>//;
    $seq = uc $$seq;
    my $n = length($seq);
    die unless $seq=~/^[ACGTN]+$/;

    if(exists $Options{'c'}) {
	if($n < $K) {
	    my $copies = $K/$n + 1;
	    $seq = $seq x $copies;
	    die "bad copies: $seq $copies $K" if length($seq) < $K;
	}
	$seq .= substr($seq, 0, $K-1);
    } else {
	$n-=$K-1;   # how far we print k-mers
    }

    for(my $i=0; $i<$n; $i++) {
	my $kmer=substr($seq, $i, $K);
	die unless length($kmer)==$K;
	next if $kmer=~/N/;  # skip Ns
	# find alphabetically smaller from kmer and its reverse complement
	my $kmer2 = reverse($kmer);
	$kmer2=~tr/ACGT/TGCA/;
	if($kmer le $kmer2) {
	    $kmer2=$kmer;
	}
	# $kmer is original, $kmer2 is canonical
	if(exists $Options{'e'}) {
	    print join("\t", $kmer2, $kmer, $name, $i), "\n";
	} elsif(exists $Options{'b'}) {
	    my $middle = POSIX::floor($i+$K/2);
	    print join("\t", $name, $middle, $middle+1, $kmer2), "\n";
	} else {
	    print $kmer2, "\n";
	}
    }
    
    
}
