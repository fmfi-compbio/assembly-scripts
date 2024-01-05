#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use POSIX;

my $MAX_LEN = 2000000;        # max. length of one region
my $MIN_SEP = 1000;           # min separation between gene and cut
my $END_SEP = 5000;           # max. amount of intergenic at sequence ends
my $GENE_CUT_PENALTY = 3;     # penalty for cutting inside gene (integer, >=1)

my $USAGE = "Usage: $0 <input> <output>

Deletes overlapping genes, alternative splicing.
Output only CDS, start_codon and stop_codon.

Also print list of intervals of length at most $MAX_LEN
that if possible does not cut any gene or closer than
$MIN_SEP to any gene.

";

die $USAGE unless scalar @ARGV==2;
my $input = $ARGV[0];
my $output = $ARGV[1];

open IN, "<$input" or die "Cannot open $input";

my %intervals;  #$intervals{seq name}{transcript id}{from/to/weight}
                # weight is the length of CDS
my @sequences;

#read input
while (my $line = <IN>) {
    $line =~ s/\s+$//;
    my @parts = split '\t', $line;
    next if ! @parts;
    if (scalar @parts<9) {
	warn "Omitting incomplete line";
	next;
    }
    die "No transcript_id in a line '$line'" 
	unless $parts[8] =~ /\btranscript_id\s+\"([^\"]*)\";/;
    
    my $id = $1;
    my $seq = $parts[0];
    my $from = $parts[3];
    my $to = $parts[4];
    if ($from>$to) {
	warn "from>to: $line ";
	next;
    }
    my $len = $to - $from + 1;

    next unless $parts[2] eq 'CDS' || $parts[2] eq 'stop_codon';

    if(!exists $intervals{$seq}) {
	push @sequences, $seq;
    }
    
    if(!exists $intervals{$seq}{$id}) {
	$intervals{$seq}{$id} = {from=>$from, to=>$to, weight=>$len};
    }
    else {
	$intervals{$seq}{$id}{from} = min($intervals{$seq}{$id}{from}, $from);
	$intervals{$seq}{$id}{to} = max($intervals{$seq}{$id}{to}, $to);
	$intervals{$seq}{$id}{weight} += $len;
    }
}   
close IN;

#compute
my %chosen;  #$chosen{seq}{transcript id} = 1 if chosen
foreach my $seq (@sequences) {
    my $res = dyn_prog_deoverlap($intervals{$seq});
    $chosen{$seq} = $res;

    my @cuts = dyn_prog_cut($intervals{$seq});
    foreach my $elem (@cuts) {
	print $seq, " ", $elem->{from}, " ", $elem->{to}, "\n";
    }
}

#write output
open IN, "<$input" or die "Cannot open $input";
open OUT, ">$output" or die "Cannot open $output";

# read file again and filter
while (my $line = <IN>) {
    my @parts = split '\t', $line;

    if (scalar @parts<9) {
	next;
    }
    die "$line " unless $parts[8] =~ /\btranscript_id\s+\"([^\"]*)\";/;
    
    my $id = $1;
    my $seq = $parts[0];

    next unless $parts[2] eq 'CDS' 
	|| $parts[2] eq 'start_codon' || $parts[2] eq 'stop_codon';

    if(exists $chosen{$seq}{$id}) {
	print OUT $line;
    }
}   
close IN;
close OUT;
exit 0;


#------------------------------------
sub dyn_prog_cut
{
    my ($intervals) = @_;

    #add $MIN_SEP to both sides
    my @expand;
    foreach my $int (values %$intervals) {
	my $from = max(1, $int->{from}-$MIN_SEP);
	my $to = $int->{to}+$MIN_SEP;
	if($from<=$to) {
	    push @expand, {from=>$from, to=>$to};
	}
    }
    @expand = sort {$a->{from} <=> $b->{from}} @expand;

    #compute union of intervals
    my @union;
    my $n = scalar @expand;
    my $i=0;
    while($i<$n) {
	my $from = $expand[$i]{from};
	my $to = $expand[$i]{to};

	$i++;
	while($i<$n && $expand[$i]{from}<=$to) {
	    $to = max($to, $expand[$i]{to});
	    $i++;
	}
	push @union, {from=>$from, to=>$to};
    }

    #handle empty union - one interval
    my $m = scalar @union;
    if($m==0) {
	return ({from=>1, to=>$MAX_LEN});
    }

    # extent of the interval
    my $start = max(1, $union[0]{from} + $MIN_SEP - $END_SEP);
    my $end = $union[$m-1]{to} - $MIN_SEP + $END_SEP;

    # trivial case with one interval
    if($end-$start+1<$MAX_LEN) {
	return ({from=>$start, to=>$end});
    }
    
    my $max_cost = POSIX::floor(($end-$start)/$MAX_LEN+2)*$GENE_CUT_PENALTY;

    #furthest point we can get with $i intervals
    my @extent = (($start-1)) x ($max_cost+1);  
    #is it inside?
    my @inside = (0) x ($max_cost+1);

    $i = 0;
    while($extent[$i]<$end) {
	die unless $i+$GENE_CUT_PENALTY<=$max_cost;

	my $next = $extent[$i]+$MAX_LEN;
	my $last_free = find_last_free(\@union, $next);
	if($last_free == $next) {   # next is not inside interval
	    store(\$extent[$i+1], \$inside[$i+1], $next, 0);
	}
	else {
	    my $cost = $i+$GENE_CUT_PENALTY;
	    store(\$extent[$cost], \$inside[$cost], $next, 1);
	    store(\$extent[$i+1], \$inside[$i+1], $last_free, 0);
	}
	$i++;
    }

    $extent[$i] = $end;

    my @endpoints;
    while($i>=0) {
	push @endpoints, $extent[$i];
	if($inside[$i]) {
	    $i -= $GENE_CUT_PENALTY;
	}
	else {
	    $i--;
	}
    }

    my @parts;
    my $last = pop @endpoints;
    while(defined(my $next = pop @endpoints)) {
	push @parts, {from=>$last+1, to=>$next};
	$last = $next;
    }

    return @parts;
}

#------------------------------------
sub find_last_free
{
    #find the first place before or at $pos that is outside 
    #all intervals, or at right boundary. Intervals are sorted and disjoint.
    
    my ($intervals, $pos) = @_;

    # bin. search to find first interval that ends > $pos
    my $from = 0;
    my $to = scalar(@$intervals) -1;
    
    # we are after last interval
    if($intervals->[$to]{to} < $pos) { return $pos; }

    while($from<$to) {
	my $m = POSIX::floor(($from+$to)/2);
	if($intervals->[$m]{to} > $pos) {
	    $to = $m;
	}
	else {
	    $from = $m+1;
	}
    }

    #check
    die unless $from==$to && $from>=0 && $from<scalar(@$intervals)
	&& $intervals->[$from]{to}>$pos 
	&& ($from==0 || $intervals->[$from-1]{to}<=$pos);
	

    if($intervals->[$from]{from}<=$pos) {   #in an interval
	return $intervals->[$from]{from}-1;
    }
    else {                                  #out of an interval
	return $pos; 
    }
}


#------------------------------------
sub store
{
    my ($val, $back, $newval, $newback) = @_;

    #compare $newval with $$val 
    #and if bigger, replace $$val with $newval and $$back with $newback

    if($$val<$newval) {
	$$val = $newval;
	$$back = $newback;
    }
}

#------------------------------------
sub dyn_prog_deoverlap
{
    my ($intervals) = @_;


    my $n = scalar keys %$intervals;
    if($n==0) { return {}; }

    #sort by interval end to index array
    my @idx = keys %$intervals;
    @idx = sort {$intervals->{$a}{to} <=> $intervals->{$b}{to}} @idx;
    
    #@a: dynamic programming array - best solution 
    #ending at $intervals->{$idx[$i]}{to}
    #@prev: previous value to look at 
    #@max: index of max. position so far
    my @a; 
    my @prev;
    my @max;
    
    #start with a single interval
    push @a, $intervals->{$idx[0]}{weight};
    push @prev, -1;
    push @max, 0;
    
    foreach my $i (1..$n-1) {
	
	my $curr_from = $intervals->{$idx[$i]}{from};

	#find the last ending before $curr_from by bin.search
	#(if any) 
	my $pos = bin_search(\@idx, $intervals, $curr_from);
	
	if($pos<0) {
	    push @a, $intervals->{$idx[$i]}{weight};
	    push @prev, -1;
	}
	else {
	
	    my $max_pos = $max[$pos];
	    
	    my $sol = $a[$max_pos] + $intervals->{$idx[$i]}{weight};
	    
	    push @a, $sol;
	    push @prev, $max_pos;
	}
	 
	if($a[$i] > $a[$max[$i-1]]) {
	    push @max, $i;
	}
	else {
	    push @max, $max[$i-1];
	}
    }
    
    #find max. value in the array
    my $max=0;
    foreach my $i (1..$n-1) {
	if($a[$i]>$a[$max]) {
	    $max = $i;
	}
    }
    
    #gather hash of chosen id's
    my %chosen;
    my $pos = $max[$n-1];
    while($pos>=0) {
	die unless $max[$pos]==$pos;
	$chosen{$idx[$pos]} = 1;
	$pos = $prev[$pos];
    }
     
    return \%chosen;
}


sub bin_search
{
    #find the last interval ending before $pos by bin.search
    my ($idx, $intervals, $pos) = @_;

    my $n = scalar(@$idx);
    if($intervals->{$idx->[0]}{to}>=$pos) { return -1; }

    my $from=0; my $to=$n-1;
    while($from<$to) {
	my $m = POSIX::floor(($from+$to+1)/2);
	if($intervals->{$idx->[$m]}{to}<$pos) {
	    $from = $m;
	}
	    else {
		$to = $m-1;
	    }
    }
    my $res = $from;   # result of bin. search
    
    die unless $from == $to && $res<$n && $res>=0;
    die unless $intervals->{$idx->[$res]}{to}<$pos;
    die unless $res==$n-1 || $intervals->{$idx->[$res+1]}{to}>=$pos;
    
    return $res;
}

#---------------------
sub max 
{
    my ($a, $b) = @_;
    return $a if($a>$b);
    return $b;
}
sub min 
{
    my ($a, $b) = @_;
    return $a if($a<$b);
    return $b;
}
