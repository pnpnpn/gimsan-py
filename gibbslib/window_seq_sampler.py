#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import sys
import logging
import re

from .dna_markov_model import DnaMarkovModel, DnaUtil
from .seq_sampler import SeqSampler
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC


class WindowSeqSampler(SeqSampler):
    def __init__(self, fasta, cfg, rand_obj, verbose = False):
        """
        fasta - list of SeqRecord (large genome file)
        rand_obj - Random number generator
        """
        super(WindowSeqSampler, self).__init__()

        #1. open genome file as one single string (removing gaps)
        #2. bin by GC content based on parameter (wndsize and numbins)
        self.logger = logging.getLogger(self.__class__.__name__)
        self.verbose = verbose

        self.rand_obj = rand_obj
        self.sampler_name = 'window'

        self.wndsize = int(cfg['window_size'])
        self.numbins = int(cfg['num_gc_bins'])

        self.genome_seq = None
        self.bin_left_lst = None #range of each bin (left-endpoint only)
        self.bin_to_poslst = None #bin-index to a list of positions (for sampling)

        self.convert_to_single_genome_seq(fasta)
        self.bin_genome_seq()

        self.reset_state()

    def reset_state(self):
        self.used_wnd_set = set()

    def is_genome_data_sufficient(self, template_fasta):
        #ensure genomic sequence is at least XX times larger than template
        multiplier_threshold = 30
        template_total_len = sum([len(record.seq) for record in template_fasta])
        if template_total_len * multiplier_threshold < len(self.genome_seq):
            return True
        return False

    def sample_wnd_pos_from_bin(self, bin_idx):
        count = 0 
        while(count < 999999):
            pos = self.rand_obj.choice(self.bin_to_poslst[bin_idx])
            if pos not in self.used_wnd_set:
                self.used_wnd_set.add(pos)
                return pos
            count += 1
        raise Exception('Failed to find non-used window when generating nullset')
 
    def sample_window_from_bin(self, bin_idx):
        pos = self.sample_wnd_pos_from_bin(bin_idx)
        return self.genome_seq[pos:pos+self.wndsize]

    def generate_ACGT_seq(self, gapless_template_seq):
        pos_to_gc = self.compute_pos_to_gc(gapless_template_seq, complete_only=False)
        null_ACGT_seq = Seq('', self.alphabet)
        for pos, gc in pos_to_gc.items():
            bin_idx = self.determine_bin_from_gc_content(gc)
            null_wnd = self.sample_window_from_bin(bin_idx)
            null_ACGT_seq += null_wnd

        #truncate right-end window (since the most right-end template window is probably incomplete)
        null_ACGT_seq = null_ACGT_seq[:len(gapless_template_seq)]

        if len(null_ACGT_seq) != len(gapless_template_seq):
            raise Exception( 'Null sequence generated has inconsistent length null=%d and template=%d' % (
                        len(null_ACGT_seq), len(gapless_template_seq)))

        self.logger.debug('Window-sampling GC comparison template=%.1f and null=%.1f' % (
            GC(gapless_template_seq), GC(null_ACGT_seq)))


        return null_ACGT_seq

    def determine_bin_from_gc_content(self, gc_content):
        """
        Naive slow method (consider using binary tree search)
        """
        for bin_idx, left_gc in enumerate(self.bin_left_lst):
            if bin_idx == self.numbins - 1:
                right_gc = 101.0 
            else:
                right_gc = self.bin_left_lst[bin_idx + 1]

            if left_gc <= gc_content and gc_content < right_gc:
                return bin_idx

        raise Exception('Cannot find bin from GC-content: %.1f' % gc_content)

    def compute_bin_to_poslst(self, pos_to_gc):
        #compute bin to position list
        self.bin_to_poslst = []
        for i in range(self.numbins):
            self.bin_to_poslst.append([])
        for pos, gc in pos_to_gc.items():
            bin_idx = self.determine_bin_from_gc_content(gc)
            self.bin_to_poslst[bin_idx].append(pos)

        for bin_idx, poslst in enumerate(self.bin_to_poslst):
            if bin_idx == 0:
                bin_left = 0.0
            else:
                bin_left = self.bin_left_lst[bin_idx]

            if bin_idx == self.numbins - 1:
                bin_right = 100.0
            else:
                bin_right = self.bin_left_lst[bin_idx+1]
            self.logger.info('Number of windows for bin %2d with [%4.1f, %4.1f): %d' % (bin_idx, bin_left, bin_right, len(poslst)))

    def compute_pos_to_gc(self, seq, complete_only=True):
        pos_to_gc = {} #position to GC-content 

        for i in range(0, len(seq), self.wndsize):
            wnd = seq[i:i+self.wndsize]
            if len(wnd) == self.wndsize or not complete_only:
                #only allow for full wndsize windows (the last window is drop if it is not full)
                pos_to_gc[i] = GC(seq[i:i+self.wndsize])
        return pos_to_gc

    def bin_genome_seq(self):
        pos_to_gc = self.compute_pos_to_gc(self.genome_seq)
        self.logger.info("Number of null windows: %d" % len(pos_to_gc))
        self.compute_bin_range(pos_to_gc.values())
        self.compute_bin_to_poslst(pos_to_gc)

    def compute_bin_range(self, gc_lst):
        gc_lst.sort()

        approx_binsize = float(len(gc_lst)) / self.numbins

        self.bin_left_lst = [] #left-endpoint of bin range
        left_idx = None
        left_approx_idx = None 
        for i in range(self.numbins):
            if i == 0:
                left_approx_idx = 0
            else:
                left_approx_idx += approx_binsize
            left_idx = int(left_approx_idx)
            self.bin_left_lst.append(gc_lst[left_idx])

        self.bin_left_lst[0] = -1e-10 #force the first bin to have left-endpoint as 0.0
        self.logger.info('GC-content bins: %s' % self.bin_left_lst)

    def convert_to_single_genome_seq(self, fasta):
        self.genome_seq = '' . join([str(s.seq).upper() for s in fasta])

        self.logger.info('Length of genome sequence with ambiguous DNA: %d' % len(self.genome_seq))

        #remove non-ACGT
        self.genome_seq = re.sub(r'[^ACGT]', '', self.genome_seq)
        self.genome_seq = Seq(self.genome_seq, self.alphabet)
        self.logger.info('Length of genome sequence ACGT-only: %d' % len(self.genome_seq))

