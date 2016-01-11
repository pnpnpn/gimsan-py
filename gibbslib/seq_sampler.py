#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import sys
import logging
import re

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC


class SeqSampler(object):
    def __init__(self):
        self.alphabet = IUPAC.unambiguous_dna
        self.rand_obj = None

    def generate_seq(self, template_record):
        """
        Returns a hash with 'seq', etc
        """
        raise NotImplementedError()

    def reset_state(self):
        raise NotImplementedError()


    def fake_generate_seq(self, template_record):
        seqstr = '' . join([self.rand_obj.choice(['A','C','G','T']) for i in range(len(template_record))]) 
        seq = Seq(seqstr, self.alphabet)

        description = {
                'sampler' : 'fake (%s)' % self.sampler_name,
                'template' : template_record.description,
                }
        return { 'seq' : seq, 'description' : description}

    def temporarily_remove_gaps(self, template_seq):
        rec = re.compile(r"[^ACGTacgt]")
        gap_lst = [(m.start(), m.group()) for m in rec.finditer(str(template_seq))]
        gapless_template = Seq(rec.sub('', str(template_seq)), template_seq.alphabet)
        return {
                'gapless_template' : gapless_template,
                'gap_data' : {'list' : gap_lst},
                }

    def insert_gaps(self, ACGT_seq, gap_data):
        gap_lst = gap_data['list']
        new_seq = Seq('', self.alphabet)
        ACGT_pointer = 0
        for (gap_pos, gap_char) in gap_lst:
            if len(new_seq) < gap_pos:
                ACGT_right_pos = gap_pos - len(new_seq) + ACGT_pointer #excluding this pos
                new_seq += ACGT_seq[ACGT_pointer:ACGT_right_pos]
                ACGT_pointer = ACGT_right_pos
            new_seq += gap_char
            if len(new_seq)-1 != gap_pos:
                raise Exception('Inconsistency when adding gap')
        new_seq += ACGT_seq[ACGT_pointer:]
        return new_seq
 
    def generate_seq(self, template_record):
        null_seq = self.generate_gapped_seq(template_record.seq)
        description = {
                'sampler' : self.sampler_name,
                'template' : template_record.description,
                'null_GC' : round(GC(null_seq),1),
                'template_GC' : round(GC(template_record.seq),1),
                }
        return { 'seq' : null_seq, 'description' : description}

    def generate_gapped_seq(self, template_seq):
        res = self.temporarily_remove_gaps(template_seq)
        gapless_template = res['gapless_template']
        gap_data = res['gap_data']

        null_ACGT_seq = self.generate_ACGT_seq(gapless_template)
        null_seq = self.insert_gaps(null_ACGT_seq, gap_data)

        if len(null_seq) != len(template_seq):
            self.logger.info("null: " + null_seq)
            self.logger.info("template: " + template_seq)
            raise Exception('Inconsistent null=%d and template=%d sequence length' % (len(null_seq), len(template_seq)))
        return null_seq





