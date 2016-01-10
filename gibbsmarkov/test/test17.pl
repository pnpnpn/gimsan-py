#!/usr/bin/perl -w

#This script tests whether ZOOPS mode will skip sequence that is shorter than the given span

use strict;
use Benchmark qw(:all);

(my $outfile = $0) =~ s/\.pl/\.out/;
(my $errfile = $outfile) =~ s/\.out/\.err/;
my $exec = '../gibbsmarkov.out';
my $param = ' -l 13 -t 50 -L 5 -em 0 -s 123 -gibbsamp -p 0.05 -best_ent -markov 5 -ds -zoops 0.2';

my $start = new Benchmark();
system("$exec input/short_seq1.fa $param 1>$outfile 2>$errfile");
my $end = new Benchmark();
my $diff = timediff($end, $start);

open(my $outfh, ">>$outfile");
printf $outfh ("Benchmark: %s\n", timestr($diff, 'all'));
close($outfh);
