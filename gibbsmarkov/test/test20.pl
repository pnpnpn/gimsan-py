#!/usr/bin/perl -w

#This script tests psp

use strict;
use Benchmark qw(:all);

(my $outfile = $0) =~ s/\.pl/\.out/;
(my $errfile = $outfile) =~ s/\.out/\.err/;
my $exec = '../gibbsmarkov.out';
my $param = ' -l 15 -t 5 -F 50 -em 0 -s 123 -gibbsamp -p 0.05 -ps 0.00 -best_ent -markov 0 -zoops 0.2 -bfile input/NOT_IN_PROMOTER_1000_1000.fa';
#my $param = ' -l 8 -cput 300 -L 200 -em 0 -s 123 -gibbsamp -p 0.05 -best_ent -markov 5 -ds -zoops 0.2 -psp input/ABF1_YPD.psp_c';

my $start = new Benchmark();
system("$exec input/combo3_1485-001-001.implant.seq $param 1>$outfile 2>$errfile");
my $end = new Benchmark();
my $diff = timediff($end, $start);

open(my $outfh, ">>$outfile");
printf $outfh ("Benchmark: %s\n", timestr($diff, 'all'));
close($outfh);
