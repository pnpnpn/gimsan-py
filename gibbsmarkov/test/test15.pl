#!/usr/bin/perl -w

#This script tests ZOOPS

use strict;
use Benchmark qw(:all);

(my $outfile = $0) =~ s/\.pl/\.out/;
(my $errfile = $outfile) =~ s/\.out/\.err/;
my $exec = '../gibbsmarkov.out';
my $param = ' -l 13 -t 100 -L 5 -em 0 -s 123 -gibbsamp -p 0.05 -best_ent -markov 3 -ds -r 1 -print_runs -zoops 0.2';

my $start = new Benchmark();
system("$exec ../../combo3_1485-001/combo3_1485-001-001_with_N_826.implant.seq $param 1>$outfile 2>$errfile");
my $end = new Benchmark();
my $diff = timediff($end, $start);

open(my $outfh, ">>$outfile");
printf $outfh ("Benchmark: %s\n", timestr($diff, 'all'));
close($outfh);
