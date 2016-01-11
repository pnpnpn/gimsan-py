#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import os
import re
import sys
import time
import ConfigParser

from gibbslib.simple_logging import *
from gibbslib.gimsan_exception import *

class GimsanJob(object):
    def __init__(self, outdir, config):
        self.config = config

        self.gimsan_home= os.path.expanduser(self.config.get('common', 'gimsan_home')) 
        self.outdir = outdir
        self.num_proc = self.config.getint('common', 'num_proc')
        self.nullset_size = self.config.getint('nullset', 'size')
        self.nullset_outdir = os.path.join(self.outdir, 'nullset/')

        if self.config.has_option('common', 'genome'):
            self.genome_filename = self.config.get('common', 'genome')
        else:
            self.genome_filename = None

        self.check_path()
        
        self.width_lst = [int(w.strip()) for w in self.config.get('common', 'widths').split(",")]

    def check_path(self):
        if not os.path.isdir(self.gimsan_home):
            raise Exception('Directory gimsan_home does not exist: %s' % self.gimsan_home)

    def get_gibbsmarkov_meta(self, w):
        finder_id = 'gm_width%03d' % w
        return {
                'finder_id' : finder_id,
                'width' : w,
                'finder_outdir' : os.path.join(self.outdir, finder_id),
                }

    def get_finders(self):
        self.finder_lst = []
        for w in self.width_lst:
            gibbsmarkov_meta = self.get_gibbsmarkov_meta(w)
            finder_id = gibbsmarkov_meta['finder_id']
            finder_outdir = gibbsmarkov_meta['finder_outdir']
            finder = GibbsMarkovFinder(
                    finder_id, 
                    w, 
                    self.gimsan_home, 
                    finder_outdir, 
                    self.genome_filename, 
                    self.config)
            self.finder_lst.append(finder)



class GibbsMarkovFinder():
    def __init__(self, id, width, gimsan_home, finder_outdir, genome_filename, config):
        self.id = id
        self.width = width
        self.gimsan_home = gimsan_home 
        self.finder_outdir = finder_outdir
        self.genome_filename = genome_filename
        self.config = config
        self.markov_order  =  self.config.getint('gibbsmarkov', 'markov')

        #self.gm_exec = os.path.join(self.gimsan_home, config.get('gibbsmarkov', 'path'))
        self.gm_exec = os.path.join(gimsan_home, 'gibbsmarkov/gibbsmarkov.out')
        if not os.path.isfile(self.gm_exec):
            raise Exception('GibbsMarkov executable missing: %s' % self.gm_exec)

    def get_stop_crit(self):
        string = self.config.get('gibbsmarkov', 'stop_crit')
        cols = string.split(':')
        if len(cols) != 2:
            raise InvalidConfigParamError('Invalid stop_crit: %s' % string)
        number = int(cols[1])
        if cols[0] == 'cput':
            return '-cput %d' % number
        elif cols[0] == 'cycle':
            return '-t %d' % number
        else:
            raise InvalidConfigParamError('Invalid stop_crit: %s' % string)

    def get_rand_seed(self):
        if self.config.has_option('gibbsmarkov', 'randseed'):
            return '-s %d' % self.config.getint('gibbsmarkov', 'randseed')
        else:
            return ''

    def get_per_seq_model(self):
        string = self.config.get('gibbsmarkov', 'per_seq_model')
        rec = re.compile(r'zoops:([\d\.]+)')
        if string == 'oops':
            return '-oops'
        elif rec.search(string):
            m = rec.search(string)
            zoops_weight = float(m.group(1))
            if zoops_weight >= 1.0 or zoops_weight <= 0.0:
                raise InvalidConfigParamError('The zoops_weight parameter must be (0.0, 1.0): %f' % zoops_weight)
            return '-zoops %.2f' % zoops_weight
        else:
            raise InvalidConfigParamError('Invalid per_seq_model: %s' % string)

    def get_rapid_conv(self):
        return '-L %d' % self.config.getint('gibbsmarkov', 'rapid_conv')

    def get_double_strand(self):
        return self.config.getint('gibbsmarkov', 'double_strand') == 1

    def get_markov(self):
        return '-markov %d' % self.markov_order

    def get_genome(self):
        return '-bfile %s' % self.genome_filename if self.genome_filename else ''

    def get_command(self, fasta_file):
        param_dct = {
                'rapid_conv' : self.get_rapid_conv(),
                'ds' : '-ds' if self.get_double_strand() else '',
                'markov': self.get_markov(),
                'per_seq_model' : self.get_per_seq_model(),
                'stop_crit' : self.get_stop_crit(),
                'bfile' : self.get_genome(),
                'rand_seed' : self.get_rand_seed(),
                }

        param_base = "-gibbsamp -best_clr -json -p 0.05 -em 0 %s" % (' ' . join(param_dct.values()),)
        cmd = '%s %s -l %d %s' % (self.gm_exec, fasta_file, self.width, param_base)
        return cmd

    def get_finder_outdir(self):
        return self.finder_outdir



