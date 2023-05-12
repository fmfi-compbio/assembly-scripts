package shared;

use strict;
use POSIX;

BEGIN {
    use Exporter   ();
    our (@ISA, @EXPORT, @EXPORT_OK);
    @ISA = qw(Exporter);
    # symbols to export by default
    @EXPORT = qw(my_run 
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


1;
