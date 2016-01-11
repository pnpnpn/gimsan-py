#!/usr/bin/perl -w

use strict;
use Benchmark qw(:all);



my @paramlst = (
    '../../combo3_1485-001/combo3_1485-001-001.implant.seq -l 13 -t 50 -L 100 -s 123 -em 50 -entsamp -best_ilr -markov 0 -r 1 -print_runs -json',
    '../../combo3_1485-001/combo3_1485-001-001.implant.seq -l 13 -t 30 -L 30 -em 0 -s 123 -gibbsamp -p 0.05 -best_ent -markov 3 -ds -r 1 -print_runs -zoops 0.2 -json',
);

if (scalar(@ARGV) == 0) {
    printf STDERR ("Execute with parameter [0 .. %d]\n", scalar(@paramlst) -1 );
    for(my $i = 0; $i < scalar(@paramlst); $i++) {
        printf STDERR ("[$i] " . $paramlst[$i] . "\n");
    }
    exit(1);
}


foreach my $i (@ARGV) {
	if ($i !~ /^\d+$/) {
		next;
	}
	my $istr = sprintf("%02d", $i);
	(my $outfile = $0) =~ s/\.pl/_$istr\.out/;
	$outfile =~ s/run\_//;
	(my $errfile = $outfile) =~ s/\.out/\.err/;
	my $exec = '../gibbsmarkov.out';
	my $param = $paramlst[$i];

	print STDERR "Running job $i ...\n";
	my $start = new Benchmark();
	my $cmd = ("$exec $param 1>$outfile 2>$errfile");
    print STDERR "$cmd\n";
    system($cmd);
	my $end = new Benchmark();
	my $diff = timediff($end, $start);

	open(my $fh, ">>$errfile");
	printf $fh ("Benchmark: %s\n", timestr($diff, 'all'));
	close($fh);
}
