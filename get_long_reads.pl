#! /usr/bin/perl -w

use strict;
use Getopt::Std;
use Data::Dumper;

my $USAGE = "$0 [-v] [-f] [-r] contig coord sign max_gap min_aln min_free len_free len_aln < psl > bed

-f Input is paf format, not psl

-v Input is paf.view, not psl

-r reverse strand in bed file, useful in telomeres

Take psl with reads (or contigs) aligned to a reference assembly.
Get bed file for extracting reads aligning close to place suggested 
by contig and coord (0-based).

If sign is -, we care about reads aligning to the right of coord,
which however extend, perhaps unaligned, to the left. 
These can be used to build a new consensus of what is to the left of coord.

max_gap: how far from coord on the right can the aln start
(aln must overlap (coord, coord+max_gap)

min_aln: how far from coord should aln extend - typically min_aln > max_gap

min_free: minimum length of read beyond coord (if aln extends to the
left of coord, we do not extract exact coordinates but assume 
approximately even distribution of indels)

len_free: maximum length extracted from reads to the left of coord
(or less, if read shorter). Typically len_free >= min_free

len_aln: maximum length extracted from reads to the right of coord
(or less if read shorter). Typically min_aln > len_aln to anchor the read.

If sign is +, \"left\" and \"right\" are interchanged.


$0 -f tig00000001 811281 - 500 2000 1000 1000 1000 < nanopore-mapped.paf > reads.bed 2> reads.paf

    ";

my %Options;
getopts('fvr', \%Options);

die $USAGE  unless @ARGV==8;
my ($contig, $coord, $sign, $max_gap, $min_aln, 
    $min_free, $len_free, $len_aln) = @ARGV;

die unless $sign=~/^[+-]$/;

my $coord_int = make_int($coord, 0, $max_gap, $sign);
my $far_coord = move_left($coord, -$min_aln, $sign);
my $far_int = [$far_coord, $far_coord+1];
while(my $line = <STDIN>) {

    my $rec;

    if(exists $Options{'f'}) {
	$rec = parse_paf_line($line);
    } elsif(exists $Options{'v'}) {
	$rec = parse_pafview_line($line);
    } else {
	$rec = parse_psl_line($line);
    }
	
    next if ! defined $rec;
    
    # check contig
    next unless $rec->{'tName'} eq $contig;
    #interval covered by aln
    my $aln_int = [$rec->{'tStart'},$rec->{'tEnd'}];

    # check that overlaps close to $coord
    next unless overlaps($coord_int, $aln_int);
    # check that extends far enough
    next unless overlaps($far_int, $aln_int);

    my $strand = $rec->{'strand'};
    die unless $strand=~/^[+-]$/;
    my $r_len = $rec->{'qSize'};

    # first $coord relative to aln as a fraction of aln len
    my $frac = ($coord-$aln_int->[0])/($aln_int->[1]-$aln_int->[0]);

    # interval of whole read
    my $r_int = [0, $r_len];    
    # get aln int in read
    # for - strand, relative to opposite strand (exchange 0 and $r_len)
    my $ra_int = [$rec->{'qStart'},$rec->{'qEnd'}];
    if($strand eq "-") {
	$ra_int = [$r_len-$ra_int->[1], $r_len-$ra_int->[0]];
    }
    
    # estimate position of $coord inside read    
    my $coord_r;
    $coord_r = $ra_int->[0]+$frac*($ra_int->[1]-$ra_int->[0]);

    # min_free int should be overlapped
    my $free_coord = move_left($coord_r, $min_free, $sign);
    my $free_int = [$free_coord, $free_coord+1];
    next unless overlaps($free_int, $r_int);

    # get interval to be taken from read
    my $ideal_int = make_int($coord_r,  $len_free, $len_aln, $sign);
    my $result = intersection($ideal_int, $r_int);

    # find how much in first and second part
    my $len1 = $coord_r-$result->[0];
    my $len2 = $result->[1]-$coord_r;
    # turn interval for - reads
    if($strand eq "-") {
	$result = [$r_len-$result->[1], $r_len-$result->[0]];
    }

    if(exists $Options{'r'}) {
	die unless $strand eq "+" || $strand eq "-";
	$strand = ($strand eq "+") ? "-" : "+";
    }
    
    # print bed 
    printf "%s\t%d\t%d\t%s_%d_%d\t0\t%s\n",
	$rec->{'qName'}, $result->[0], $result->[1],
	$rec->{'qName'}, $len1, $len2, $strand;
    # print psl
    print STDERR $line;
}

exit 0;


sub parse_format_line {
    my ($line, $keys) = @_;

    my @parts = split " ", $line;
    return undef if @parts==0;

    die "Bad line $line" unless scalar(@parts)>=scalar(@$keys);
    
    my %rec;
    @rec{@$keys} = @parts;
    return \%rec;
}

sub parse_psl_line {
    my ($line) = @_;
    
    my @psl_keys = qw/matches misMatches repMatches nCount qNumInsert
    qBaseInsert tNumInsert tBaseInsert strand qName qSize qStart qEnd
    tName tSize tStart tEnd blockCount blockSizes qStarts tStarts/; 

    return parse_format_line($line, \@psl_keys);
}

sub parse_paf_line {
    my ($line) = @_;

    my @paf_keys = qw/qName qSize qStart qEnd strand tName tSize
    tStart tEnd matches alnLen quality/;

    return parse_format_line($line, \@paf_keys);
}

sub parse_pafview_line {
    my ($line) = @_;

    my @pafview_keys = qw/matches strand qName qSize qStart qEnd tName tSize
    tStart tEnd pid/;

    return parse_format_line($line, \@pafview_keys);
}


sub move_left {
    my ($x, $dist, $sign) = @_;
    return (($sign eq "-") ? $x-$dist : $x+$dist);
}

sub make_int {
    my ($x, $len_left, $len_right, $sign) = @_;
    my $s = ($sign eq "-") ? $x-$len_left : $x-$len_right;
    my $e = ($sign eq "-") ? $x+$len_right : $x+$len_left;
    
    return [$s, $e];
}

sub min {
    my ($x, $y) = @_;
    return (($x<$y) ? $x : $y);
}

sub max {
    my ($x, $y) = @_;
    return (($x>$y) ? $x : $y);
}

sub intersection {
    my ($int1, $int2) = @_;
    my $s = max($int1->[0], $int2->[0]);
    my $e = min($int1->[1], $int2->[1]);
    die unless $s<$e;
    return [$s, $e];
}

sub overlaps {
    my ($int1, $int2) = @_;
    return $int1->[0] <= $int2->[1]
	&& $int1->[1] >= $int2->[0];
}

