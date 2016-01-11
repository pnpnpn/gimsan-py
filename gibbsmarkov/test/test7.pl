#!/usr/bin/perl -w

use strict;
use Benchmark qw(:all);

(my $outfile = $0) =~ s/\.pl/\.out/;
(my $errfile = $outfile) =~ s/\.out/\.err/;
my $exec = '../gibbsmarkov.out';
my $param = ' -l 13 -t 50 -L 50 -em 0 -s 123 -gibbsamp -best_clr -markov 3 -ds -r 1 -p 0.1 -print_runs';

my $start = new Benchmark();
system("$exec ../../combo3_1485-001/combo3_1485-001-001.implant.seq $param 1>$outfile 2>$errfile");
my $end = new Benchmark();
my $diff = timediff($end, $start);

open(my $outfh, ">>$outfile");
printf $outfh ("Benchmark: %s\n", timestr($diff, 'all'));
close($outfh);
