#! /usr/bin/perl -w

use strict;

my $USAGE = "$0 fasta contig start end [strand [name]] > output.fasta

If start and/or end are negative, change them to length+x+1,
i.e. -1 means end of sequence.
 "; 

die $USAGE unless @ARGV>=4 && @ARGV<=6;

my ($fasta, $contig, $start, $end, $strand, $name) = @ARGV;
die unless -r $fasta;

if(!defined $strand) { $strand = '+'; }
die "Bad strand '$strand'" unless $strand=~/^[+-]$/;

if($end<0 || $start<0) {
    my $run = "faSize -detailed \"$fasta\"";
    my $res = `$run`;
    my @lines = split "\n", $res;
    my $len;
    foreach my $line (@lines) {
	my @parts = split " ", $line;
	die "bad format $line" unless @parts==2;
	if($parts[0] eq $contig) {
	    die if defined $len;
	    $len = $parts[1];
	}
    }
    die "finding sequence size failed" unless defined $len && $len =~ /^\d+$/;
    if($end<0) {
	$end = $len+$end+1;
    }
    if($start<0) {
	$start = $len+$start+1;
    }
}


if(!defined $name) { 
    $name = "$contig:$start-$end"; 
    if($strand eq '-') { $name .= "($strand)"; }
}

my $bed = join("\t", $contig, $start, $end, $name, 0, $strand);

my $run = "echo \"$bed\" | bedtools getfasta -fi \"$fasta\" -bed stdin -fo stdout -s -name";
my $ret = system($run);
exit $ret;
