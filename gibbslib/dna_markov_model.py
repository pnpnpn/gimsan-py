#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import numpy as np

class DnaUtil:
    def __init__(self):
        self.nt2int_dct = {'A':0, 'C':1, 'G':2, 'T':3}
        self.int2nt_arr = ['A', 'C', 'G', 'T', 'N']
        self.acgt_set = set(['A', 'C', 'G', 'T'])

    def is_acgt_only(self, seq):
        return all(nt in self.acgt_set for nt in seq)

    def nt2int(self, nt):
        #return '4' for unknown/gap
        return self.nt2int_dct.get(nt, 4)

    def seq_to_intarr(self, seq):
        return [self.nt2int(nt) for nt in seq]

    def intarr_to_seq(self, intarr):
        return '' . join([self.int2nt(i) for i in intarr])

    def int2nt(self, i):
        return self.int2nt_arr[i]

class DnaMarkovModel():
    def __init__(self, markov_order, rand_obj):
        import logging
        self.nalphas = 4
        self.logger = logging.getLogger(self.__class__.__name__)
        self.order = markov_order
        self.dna_util = DnaUtil()
        self.rand_obj = rand_obj

        self.reset_kmer_dct()

        #count matrix
        self.cntmat = np.ones([self.nalphas] * (self.order + 1), dtype = int) 
        self.freqmat = {}

        #self.freqmat[ord][tuple(coord)]

    def print_model(self, mat, fptr, val_fmt="%s"):
        fmt = "%s  " + val_fmt
        for tup,val in np.ndenumerate(mat):
            kmer = self.dna_util.intarr_to_seq(tup)
            print(fmt % (kmer, val), file=fptr)

    def reset_kmer_dct(self):
        self.kmer_dct = {}

    def train(self, seq_record_lst):
        for record in seq_record_lst:
            self.logger.debug("record-length: %d" % len(record))
            self.tally_kmers_from_seq(record)
        self.push_kmers_to_countmat()
        self.compute_freqmat()

    def sample_from_prob_arr(self, prob_arr):
        """
        Sample based on probability array. This returns the index
        """
        tot = 0.0
        rand = self.rand_obj.random()
        for i, p in enumerate(prob_arr):
            tot += p
            if tot >= rand:
                return i
        return len(prob_arr)-1

    def sample_seq(self, seqlen):
        intarr = [] #convert to seq later
        for i in range(seqlen):
            ord = min(self.order, len(intarr))
            coord = intarr[-ord:]
            prob_arr = self.freqmat[ord][tuple(coord)]
            intarr.append(self.sample_from_prob_arr(prob_arr))
            #self.logger.info(intarr)
        return self.dna_util.intarr_to_seq(intarr)

    def compute_freqmat(self):
        """
        Compute model of order 1, 2, etc
        """
        cmat = self.cntmat.copy()

        for ord in range(self.order, -1, -1):
            res = self.compute_freqmat_helper(cmat, ord)
            self.freqmat[ord] = res['fmat']
            cmat = res['cmat_minus']

        self.logger.debug(self.freqmat)

    def generate_all_permuations(self, ord):
        """
        Generate all permutations for alphabet
        e.g. for ord=3, we have (0,0,0), (0,0,1), ...
        """
        lim = self.nalphas ** ord
        for i in range(lim):
            coord = [None] * ord
            for k in range(ord):
                coord[k] = i % self.nalphas
                i //= self.nalphas
            yield tuple(coord)

    def compute_freqmat_helper(self, cmat, ord):
        assert(cmat.ndim == ord + 1)
        cmat_minus = np.sum(cmat, ord)
        fmat = cmat.astype(np.float) #makes a copy

        for coord in self.generate_all_permuations(ord):
            #self.logger.info("fmat.ndim=%d  cmat_minus.ndim=%d" % (fmat.ndim, cmat_minus.ndim))
            fmat[coord] /= cmat_minus[coord]
        return {
                'cmat_minus' : cmat_minus, 
                'fmat' : fmat,
                }

    def push_kmers_to_countmat(self):
        for kmer,cnt in self.kmer_dct.items():
            intarr = self.dna_util.seq_to_intarr(kmer)
            assert(len(intarr) == self.order + 1) 
            self.cntmat[tuple(intarr)] += cnt

        #reset
        self.reset_kmer_dct()
        self.logger.debug(self.cntmat)

    def tally_kmers_from_seq(self, seq_record):
        seq = seq_record.seq.upper()
        kmer_len = self.order + 1
        for i in range(len(seq) - kmer_len + 1):
            kmer = seq[i:i+kmer_len]
            assert(len(kmer) == kmer_len)
            if not self.dna_util.is_acgt_only(kmer):
                continue
            self.kmer_dct[kmer] = self.kmer_dct.get(kmer, 0) + 1
        self.logger.debug('kmer_dct size: %d' % len(self.kmer_dct))


