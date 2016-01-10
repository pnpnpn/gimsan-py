#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import time_benchmark
import os
import sys
import json
import re
import time
import numpy as np
import ConfigParser
import subprocess

from argparse_plus import ArgumentParserPlus
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gibbslib.batch_experiment import *
from gibbslib.simple_logging import *
from gibbslib.gimsan_base import GimsanJob
from gibbslib.gimsan_exception import *
import gibbslib.nullset_generator as nullset_generator

if sys.version_info < (2,7,0):
    sys.stderr.write("You need Python 2.7 or later to run this script\n")
    exit(1)

def cmd_exec_func(job):
    import logging
    logging.info('%s : %s' % (job['job_id'], job['cmd']))
    time_wall_start = time.time()
    ret_code = subprocess.call(job['cmd'], shell=True)
    diff_time_wall = (time.time() - time_wall_start) / 60.0
    logging.info('%s DONE (wall time %.3fm)' % (job['job_id'], diff_time_wall))

class GimsanSubmitJob(GimsanJob):
    def __init__(
            self, 
            name, 
            template_file, 
            outdir, 
            config, 
            is_overwrite=False, 
            skip_null=False, 
            dryrun=False, 
            verbose=False):
        super(GimsanSubmitJob, self).__init__(outdir, config)

        self.name = name
        self.logger = logging.getLogger(self.__class__.__name__)
        self.template_file = template_file
        self.verbose = verbose
        self.dryrun = dryrun
        self.is_overwrite = is_overwrite
        self.skip_null = skip_null

        self.check_submit_path()
        self.setup_finders() #for each width

    def check_submit_path(self):
        if not os.path.isdir(self.outdir):
            os.mkdir(self.outdir)
        elif not self.is_overwrite:
            raise AlreadyExistOutputDirError('Directory already exist: %s' % self.outdir)

        if not os.path.isfile(self.template_file):
            raise MissingFileError('Missing template file: %s' % self.template_file)

    def setup_finders(self):
        self.get_finders()
        for finder in self.finder_lst:
            finder_outdir = finder.finder_outdir
            if not os.path.isdir(finder_outdir):
                os.mkdir(finder_outdir)
            elif not self.is_overwrite:
                raise AlreadyExistOutputDirError('Directory already exist: %s' % finder_outdir)

    def generate_nullset_helper(self):
        window_sampler_config = None
        if self.config.has_option('common', 'genome'):
            genome_file = os.path.expanduser(self.config.get('common', 'genome'))
            mode = 'window' 
            window_sampler_config = dict(self.config.items('window_sampling'))
        else:
            genome_file = None
            mode = 'markov'

        import random
        if self.config.has_option('nullset', 'randseed'):
            nullset_randseed = self.config.getint('nullset', 'randseed')
        else:
            nullset_randseed = random.randint(1, 2**31-1)
        self.logger.info('Nullset random seed: %d' % nullset_randseed)
        nullset_rand_obj = random.Random(nullset_randseed)

        nullset_generator.run_generation(
                nullset_rand_obj,
                mode, 
                self.template_file, 
                genome_file,
                self.nullset_size,
                self.nullset_outdir,
                window_sampler_config,
                self.verbose)


    def batch_execute(self):
        job_lst = []

        for finder in self.finder_lst:
            for i in range(0, self.nullset_size+1):
                #cmd = "sleep 5 && echo %d" % idx
                if i == 0:
                    fasta_file = self.template_file
                    stdout_file = os.path.join(finder.get_finder_outdir(), 'motif-%05d.stdout' % i)
                    stderr_file = os.path.join(finder.get_finder_outdir(), 'motif-%05d.stderr' % i)
                else:
                    fasta_file = os.path.join(self.nullset_outdir, 'null-%05d.fa' % i)
                    stdout_file = os.path.join(finder.get_finder_outdir(), 'null-%05d.stdout' % i)
                    stderr_file = os.path.join(finder.get_finder_outdir(), 'null-%05d.stderr' % i)
                cmd = "%s 1>%s 2>%s" % (finder.get_command(fasta_file), stdout_file, stderr_file)
                job = {
                        'cmd' : cmd,
                        'job_id' : '%s|%s|idx=%05d' % (self.name, finder.id, i),
                        }

                job_lst.append(job)

        if self.dryrun:
            for job in job_lst:
                self.logger.info(job['cmd'])
        else:
            import multiprocessing
            pool = multiprocessing.Pool(processes=self.num_proc)
            pool.map(cmd_exec_func, job_lst)

    def submit_job(self):
        if not self.skip_null:
            self.generate_nullset_helper()
        self.batch_execute()

if __name__ == '__main__':
    benchmark = time_benchmark.Benchmark()

    #defaults
    description = """
Submit GIMSAN job
""" 


    epilog = """
Examples:

%(prog)s --conf=conf_examples/test_window_sampling.cfg -v
%(prog)s --conf=conf_examples/test_window_sampling.cfg -v --overwrite 
%(prog)s --conf=conf_examples/test_window_sampling.cfg -v --overwrite --dryrun --skip-null
""" 

    argp = ArgumentParserPlus(description=description, epilog=epilog)
    argp.add_argument('--conf', required=True, help="", dest="conf_file")
    argp.add_argument('--overwrite', action="store_true", help="") 
    argp.add_argument('--skip-null', action="store_true", help="") 
    argp.add_argument('--dryrun', action="store_true", help="")
    argp.add_argument('-v', '--verbose', action='store_true')
    args = argp.parse_args()

    import logging
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    simple_stream_logging(log_level)

    args.conf_file = os.path.expanduser(args.conf_file)

    batch_exp = BatchExperiment(args.conf_file, args.overwrite)
    batch_exp.mkdir_submit_directories()
    batch_exp.copy_conf_file()

    for exp in batch_exp.experiments:
        gsj = GimsanSubmitJob(
                exp['name'],
                exp['fasta_file'], 
                exp['outdir'],
                batch_exp.config, 
                is_overwrite = args.overwrite, 
                skip_null = args.skip_null, 
                dryrun = args.dryrun,
                verbose = args.verbose)
        gsj.submit_job()
    batch_exp.flag_success()

    benchmark.print_time(sys.stderr)

