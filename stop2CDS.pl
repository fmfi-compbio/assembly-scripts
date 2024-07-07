#! /usr/bin/perl -w

use strict;

my $USAGE = "$0 in.gtf > out.gtf

Take exons in input, convert them to CDS,
omit original CDS and copy frame from CDS/stop_codon to exons
with the same 5' end (wrt strand of the gene).

Assumes no UTRs are included in exons.
";



die $USAGE unless @ARGV==1;
my $in;
open $in, "<", $ARGV[0] or die;
my $hash = {};
my @lines;
while(my $line = <$in>) {
    chomp $line;
    my @parts = split "\t", $line;
    if($parts[2] eq "CDS" || $parts[2] eq "stop_codon") {
	my $key = get_key(\@parts);
	die if exists $hash->{$key} && $hash->{$key} ne $parts[7];
	$hash->{$key} = $parts[7];
    }
    if($parts[2] ne "CDS") {
	push @lines, \@parts;
    }
}

foreach my $parts (@lines) {
    if($parts->[2] eq "exon") {
	$parts->[2] = "CDS";
	die unless $parts->[7] eq ".";
	my $key = get_key($parts);
	die unless exists $hash->{$key};
	$parts->[7] = $hash->{$key};
    }
    
    print join("\t", @$parts), "\n";
}
    

sub get_key {
    my ($parts) = @_;

    if($parts->[6] eq "+") {
	return sprintf("%s %d +", $parts->[0], $parts->[3]);
    } else {
	return sprintf("%s %d -", $parts->[0], $parts->[4]);
    }
}
