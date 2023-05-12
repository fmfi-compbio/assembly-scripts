#! /usr/bin/perl -w

use strict;
use Getopt::Std;
use File::Temp qw/tempdir/;

# import shared.pm from script directory
use FindBin qw($Bin);  #$Bin is the directory with the script
use lib $Bin;        #add bin to the library path
use shared;


my $USAGE = "
$0 [options] fasta_list region_list > output.fa

Region list is a file with 5 columns whitespace spearated
0: old_contig_name, 1:start, 2:end, 3:strand, 4:new_name
If a new contig consists of several parts, they have the same name in the last column and should be adjacent
The file may contain also empty lines and comment lines starting with #

Beware: currently no checks if old contigs exists, are unique and sufficiently long (TODO)

Options
-b filename   Output filename with join positions
";

my %Options;
getopts('b:', \%Options);

# check arguments
die $USAGE unless @ARGV==2;
my ($fasta_list, $region_list) = @ARGV;
my $bed_file;
if (exists $Options{'b'}) {
    $bed_file = $Options{'b'}
}
foreach my $f ($fasta_list, $region_list) {
    die "cannot find '$f'" unless -r $f;
}

my $dir = tempdir( CLEANUP => 1 );
print STDERR "Temp dir $dir\n";

# reformat region_list to a proper bed file
# check that names are unique, add suffix _CONT to each contig continuation
my %names;
my $prev_name;
my $in;
open $in, "<", $region_list or die;
my $out;
open $out, ">", "$dir/regions.bed" or die;
while(my $line = <$in>) {
    my @parts = split " ", $line;
    next if @parts==0 || $line =~ /^\#/;   # skip empty lines and comments

    die "Bad number of columns in line '$line'" unless @parts==5;
    die "Bad start in line '$line'" unless $parts[1]=~/^\d+$/;
    die "Bad end in line '$line'" unless $parts[2]=~/^\d+$/;
    die "Bad strand in line '$line'" unless $parts[3]=~/^[+-]$/;
    my $name = $parts[4];
    die "new contig names cannot end in _CONT" if $name =~ /_CONT$/;
    if(! defined $prev_name || $prev_name ne $name) {
        die "Duplicated name $name" if exists $names{$name};
	$names{$name} = 1;
	$prev_name = $name;	
    } else {
	$name = $name . "_CONT";
    }
    print $out join("\t", @parts[0..2], $name, 0, $parts[3]), "\n";    
}
close($out) or die;
close($in) or die;

# read the list of fasta files, check that they exist
my @fasta;
open $in, "<", $fasta_list or die;
while(my $line = <$in>) {
    $line =~ s/\s+$//;
    $line =~ s/^\s+//;
    die "missing fasta $line" unless -r $line;
    push @fasta, $line;
}
close $in or die;

# create concatenated fasta
my_run("cat " . join(" ", @fasta) . " > $dir/tmp1.fa");

# get regions
my_run("bedtools getfasta -fi $dir/tmp1.fa -bed $dir/regions.bed -s -name | perl -ne 's/\\([+-]\\)\$//; print' >$dir/tmp2.fa");
# remove fasta header from _CONT regions
my_run("grep -v '>.*_CONT\$' $dir/tmp2.fa > $dir/tmp3.fa");
my_run("fastareformat $dir/tmp3.fa");  # to stdout
unlink("$dir/tmp1.fa.fai") or warn "fai not found";