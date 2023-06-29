#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Std;
use Bio::Tools::CodonTable;
use FindBin qw($Bin);  #$Bin is the directory with the script
use lib $Bin;        #add bin to the library path
use shared;

my $USAGE =
"$0 [options] <input fasta> <input gtf> <min. support> <min. length> <transcript output fasta> <protein output fasta>

Read CDS coordinates from gtf and print transcript fasta 
to standard output. 

Min. support is in percent, use -10 to get all.
Min. length is in nucleotides.

Options:
-s  do not prepend sequence name in front of transcript_id
    (by default sequence name added if id does not start with it)

-S  do not add support data to sequence name

-g number genetic code (default 1)

-c do not add * for the final stop codon
";

my %Options;
getopts('sSg:c', \%Options);
if(!exists $Options{'g'}) {
    $Options{'g'}=1;
}

die $USAGE unless @ARGV == 6;
my ($input_fasta, $input_gtf, $min_support, $min_length,
    $output_rna, $output_prot) = @ARGV;

my %genes;         #$genes{$gene_name}[$exon_num]{gene/from/to/strand}
my %contig2genes;  #$contig2genes{$contig}{$gene_name} = 1 if gene on contig

#read gff and store coordinates
open IN, "<$input_gtf" or die;
while(my $line = <IN>) {
    $line =~ s/\s+$//;
    next if ($line eq '' || $line =~ /^#/);
    my ($contig, $source, $feature, $from, $to, $score, $strand, 
	$frame, $rest) = split '\t', $line;
    warn "Bad line 'line'" unless defined $rest;

    next unless $feature eq 'CDS' || $feature eq 'stop_codon';

    #split the rest, store in a hash
    my %parts;
    foreach my $part (split ';', $rest) {
	die unless $part =~ /^\s*(\S+)\s+\"(.*)\"\s*$/;
	my ($first, $second) = ($1, $2);    
	$parts{$first} = $second;
    }

    die unless $from <= $to;

    my $support = -1;
    if(exists $parts{'transcript_support'} && !exists $Options{'S'}) {
	$support = ($parts{'transcript_support'} =~ tr/X/X/);
    }
	
    my $gene_name = $parts{'transcript_id'};
    if(index($gene_name, $contig)!=0 && ! exists $Options{'s'}) {
	$gene_name = $contig . '-' . $gene_name;
    }

    push @{$genes{$gene_name}}, {'gene'=>$gene_name, 
				 'from'=>$from, 'to'=>$to, 
				 'frame'=>$frame,
				 'strand'=>$strand, 'support'=>$support};
    $contig2genes{$contig}{$gene_name} = 1;

}
close IN;


#read fasta and print proteins
open IN, "<$input_fasta" or die;
open OUT_RNA, ">$output_rna" or die;
open OUT_PROT, ">$output_prot" or die;
while(1) {
   my ($name, $seq) = read_fasta(\*IN);

   last if !defined($name);
   $name = normalize_fasta_name($name);

   foreach my $gene (keys %{$contig2genes{$name}}) {
       #print STDERR $gene, " ";
       my ($transcript, $support, $frame) 
	   = get_transcript($genes{$gene}, $seq);

       my $name = '>' . $gene;
       if(defined $support && $support >=0) {
	   $name .= sprintf " support %d-%d%%", $support*10, $support*10+10;
       }

       if(defined($transcript) 
	  && length($transcript)>=$min_length
	  && $support*10>=$min_support) {
	   write_fasta(\*OUT_RNA, $name, \$transcript);

	   my $to_translate = substr($transcript, $frame);

	   my $protein = translate($to_translate, $Options{'g'});
	   if(exists $Options{'c'}) {  # skip last * for stop codon if present and -c option chosen
	       $protein =~ s/\*$//;
	   }
	   
	   write_fasta(\*OUT_PROT, $name, \$protein);
       }
   }
}
close IN;
close OUT_RNA;
close OUT_PROT;


exit 0;

############################
sub get_transcript {
    my ($gene, $seq) = @_;

    my $n = scalar @$gene;
    my $len = length($$seq);
    #print STDERR $n, "\n";

    #check strand and support
    my $strand = $gene->[0]{'strand'};
    my $support = -1;
    foreach my $exon (0..$n-1) {
	die unless $strand eq $gene->[$exon]{'strand'};
	if(exists $gene->[$exon]{'support'} &&  $gene->[$exon]{'support'} >= 0) {
	    if($support < 0) {
		$support = $gene->[$exon]{'support'};
	    }
	    die "support not the same in all exons"  unless $support eq $gene->[$exon]{'support'};
	}
    }

    if($strand eq '+') {
	#sort in increasing order
	@$gene = sort {$a->{'from'} <=> $b->{'from'}} @$gene;
    }
    else {
	@$gene = sort {$b->{'from'} <=> $a->{'from'}} @$gene;
    }

    my $result  = '';
    foreach my $exon_idx (0..$n-1) {
	my $exon = $gene->[$exon_idx];
	
	if($exon->{'from'} < 1 || $exon->{'to'} > $len) {
	    print STDERR "Skipping transcript extending beyond sequence "
                . $exon->{'gene'} . " (seq len $len)\n";
	    return;
	}
	my $exon_len = $exon->{'to'} - $exon->{'from'} + 1;

	my $exon_seq = substr($$seq, $exon->{'from'}-1, $exon_len);
	if($strand eq '-') {
	    $exon_seq = reverse_seq($exon_seq);
	}
	$result .= $exon_seq;
    }

    my $frame = $gene->[0]{'frame'};

    return ($result, $support, $frame);
}


####################################################
sub reverse_seq
{
   my ($seq) = @_;

   $seq =~
       tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;

   return scalar reverse $seq;
}


####################################################
sub translate
{
    my ($seq, $code) = @_;
    my $CodonTable = Bio::Tools::CodonTable->new( -id => $code);
    my $result = $CodonTable->translate($seq);

    return $result;
}

