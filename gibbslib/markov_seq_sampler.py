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


class MarkovSeqSampler(SeqSampler):
    def __init__(self, markov_order, fasta, rand_obj, verbose=False):
        """
        fasta - list of SeqRecord
        rand_obj - Random number generator
        """
        super(MarkovSeqSampler, self).__init__()

        self.logger = logging.getLogger(self.__class__.__name__)
        self.sampler_name = 'markov'
        self.rand_obj = rand_obj

        self.dna_util = DnaUtil()
        self.markov = DnaMarkovModel(markov_order, rand_obj)
        self.markov.train(fasta)

        if verbose:
            self.markov.print_model(self.markov.cntmat, sys.stderr)


    def generate_ACGT_seq(self, gapless_template_seq):
        return Seq(self.markov.sample_seq(len(gapless_template_seq)), self.alphabet)

    def reset_state(self):
        #No state for markov chain generation
        pass



