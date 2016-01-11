#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import sys

import time_benchmark
import os
import json
import numpy as np
import random
from argparse_plus import ArgumentParserPlus

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from gibbslib.simple_logging import *
from gibbslib.markov_seq_sampler import MarkovSeqSampler 
from gibbslib.window_seq_sampler import WindowSeqSampler

class NullsetGenerator():
    def __init__(self, seq_sampler, template_fasta):
        import logging
        self.logger = logging.getLogger(self.__class__.__name__)

        self.seq_sampler = seq_sampler
        self.template_fasta = template_fasta

    def generate_single_sample(self, iter_idx):
        record_lst = []
        self.seq_sampler.reset_state()
        for i,template_record in enumerate(self.template_fasta):
            new_seq_res = self.seq_sampler.generate_seq(template_record)
            new_description = json.dumps(new_seq_res['description'], sort_keys=True)
            record = SeqRecord(
                    new_seq_res['seq'], 
                    id="null-%05d-%05d" % (iter_idx, i),
                    description=new_description,
                    )
            record_lst.append(record)
        return record_lst

    def generate_samples(self, iters, output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        
        for i in range(1, iters+1):
            null_sample = self.generate_single_sample(i)
            filename = os.path.join(output_dir, "null-%05d.fa" % i)
            SeqIO.write(null_sample, filename, "fasta")


def open_fasta_file(fasta_file):
    """
    Return sequence record list
    """
    with open(fasta_file, "rU") as fh:
        return list(SeqIO.parse(fh, "fasta"))
    raise Exception('Invalid FASTA file')

def run_generation(
        rand_obj,
        nullset_mode, 
        template_file, 
        genome_file, 
        iters, 
        outdir, 
        window_sampler_config,
        verbose):
    markov_order = 3
    template_fasta = open_fasta_file(template_file) 

    if genome_file and nullset_mode == 'window':
        genome_fasta = open_fasta_file(genome_file)
        seq_sampler = WindowSeqSampler(genome_fasta, window_sampler_config, rand_obj, verbose=verbose)
        if not seq_sampler.is_genome_data_sufficient(template_fasta):
            raise Exception("Insufficient genomic data. Consider: \n[1] Use a genome file with more data, " \
                    "or \n[2] Do not specify a genomic file (GIMSAN would use Markov nullset sampling instead). " \
                    "\nWe recommend method [1] if possible.")
    else:
        seq_sampler = MarkovSeqSampler(markov_order, template_fasta, rand_obj, verbose=verbose)


    nullset_generator = NullsetGenerator(seq_sampler, template_fasta)
    nullset_generator.generate_samples(iters, outdir)


if __name__ == '__main__':
    benchmark = time_benchmark.Benchmark()

    #defaults
    description = """
Stand-alone version of sample_from_GC_content_range. Construct null datasets based on GC-content percentages
This is particularly created for the GIMSAN application.

Default window size: 100 bp
Default GC-content bin size: 5%%
If --genome is not specified, Markov model is used to generate the nullset.
""" 


    epilog = """
Examples:

%(prog)s --template=testdata/ABF1_YPD_mod.fasta --outdir=testnull --mode=markov
""" 
    #iters, genome, output-dir, binsize, wndsize

    argp = ArgumentParserPlus(description=description, epilog=epilog)
    argp.add_argument('--mode', required=True, choices=['markov', 'window'], help="Markov or window-sampling")
    argp.add_argument('--template', required=True, help="template FASTA file")
    argp.add_argument('--genome', help="genome FASTA file")
    argp.add_argument('--iters', help="", type=int, default=10)
    argp.add_argument('--outdir', help="", required=True)
    argp.add_argument('-v', '--verbose', action='store_true')
    args = argp.parse_args()

    import logging
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    simple_stream_logging(log_level)
    rand_obj = random.Random()
    run_generation(
            rand_obj, 
            args.mode, 
            args.template, 
            args.genome, 
            args.iters, 
            args.outdir, 
            {}, 
            args.verbose)
    benchmark.print_time(sys.stderr)

