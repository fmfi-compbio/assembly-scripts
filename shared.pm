package shared;

use strict;
use POSIX;

BEGIN {
    use Exporter   ();
    our (@ISA, @EXPORT, @EXPORT_OK);
    @ISA = qw(Exporter);
    # symbols to export by default
    @EXPORT = qw(my_run read_fasta normalize_fasta_name write_fasta
                 );
}


############################
sub my_run
{
    my ($run, $die) = @_;
    if(!defined($die)) { $die = 1; }

    my $short = substr($run, 0, 20);

    print STDERR $run, "\n";
    my $res = system("bash", "-c", $run);
    if($res<0) {
        die "Error in program '$short...' '$!'";
    }
    if($? && $die) {
        my $exit  = $? >> 8;
        my $signal  = $? & 127;
        my $core = $? & 128;

        die "Error in program '$short...' "
            . "(exit: $exit, signal: $signal, dumped: $core)\n\n ";
    }
}

############################
sub read_fasta {
    # Read one fasta sequence from the fasta file
    # Return undef, if no sequence found; 
    #        name and reference to sequence otherwise.
    
    my ($input) = @_;
    
    # read name of fasta sequence
    my $name;
    while(1) {
        my $line = <$input>;
        
        # end of file
        return unless defined $line;

        # skip empty lines
        next if $line =~ /^\s*$/;

        # parse the name
        $line =~ s/\s+$//;
        if($line =~ /^>/) {
            $name = $line; 
            last;
        }
        else { die "Incorrect fasta file '$line'."; }
    }

    # read fasta sequence
    my $sequence = "";
    while(1) {
        my $file_pos = tell($input);
        my $line = <$input>;
        
        # end of file
        last unless defined $line; 

        # if there is a beginning of a new sequence
        if($line =~ /^>/) {
            seek($input, $file_pos, 0);
            last;
        }

        # remove all whitespaces from line and append it to the sequence
        $line =~ s/\s//g;
        $sequence .= $line;
    }
    
    return ($name, \$sequence);
}


############################
sub normalize_fasta_name
{
    my ($name) = @_;

    $name =~ s/\s+$//;          # remove trailing white space
    $name =~ s/^>\s*//;         # remove > at the beginning of line
    $name =~ s/^(\S+)\s.*$/$1/; # remove everything after first space

    #if it has gi, take only number
    if ( $name =~ /gi\|(.*)\|/ ) {
        $name = $1;
    }

    return $name;
}


############################
sub write_fasta {
    my ($file, $name, $seq) = @_;

    print $file $name, "\n";
    my $n = length($$seq);
    my $linelen = 60;
    my $i=0;
    while($i<$n) {
        print $file (substr($$seq, $i, $linelen), "\n");
        $i+=$linelen;
    }
}



1;
