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

from tornado import template
from argparse_plus import ArgumentParserPlus
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Motif import Motif

from gibbslib.batch_experiment import *
from gibbslib.simple_logging import *
from gibbslib.gimsan_base import GimsanJob
from gibbslib.gimsan_exception import *



class FinderResult():
    def __init__(self, finder_meta, nullset_size, statsig_dir, outdir, gimsan_home):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.finder_id = finder_meta['finder_id']
        self.finder_id_alt = self.finder_id.replace('_', '-')
        self.finder_outdir = finder_meta['finder_outdir'] #per-finder
        self.width = finder_meta['width']
        self.nullset_size = nullset_size
        self.statsig_dir = statsig_dir
        self.outdir = outdir #per-template

        self.load_template_json()
        self.pval_r_path = os.path.join(gimsan_home, 'misc', "conf_pval_only.R")

        #results
        self.pvalue_comment = None
        self.kmer_lst = None
        self.kmer_filename = None
        self.logowidth = 5

        self.weblogo_basename = None
        self.weblogo_revcompl_basename = None

        self.coldep_outfile = None
        self.coldep_num_pairs = None

    def construct_weblogo(self, weblogo_filename, weblogo_revcompl_filename):
        self.weblogo_basename = os.path.basename(weblogo_filename)
        self.weblogo_revcompl_basename = os.path.basename(weblogo_revcompl_filename)

        motif = Motif(alphabet=IUPAC.unambiguous_dna)
        for kmer in self.kmer_lst:
            motif.add_instance(Seq(kmer, motif.alphabet))

        logowidth_normal = self.construct_weblogo_helper(weblogo_filename, motif)

        #reverse complement
        motif_revcompl = motif.reverse_complement()
        logowidth_revcompl = self.construct_weblogo_helper(weblogo_revcompl_filename, motif_revcompl)

        self.logowidth = max(self.logowidth, logowidth_normal, logowidth_revcompl)

    def construct_weblogo_helper(self, weblogo_filename, motif):
        logowidth = (20.0/45.0) * self.width + 2;
        motif.weblogo(weblogo_filename, logowidth=logowidth)
        #return logowidth * 50 #width to pixel conversion
        return logowidth

    def load_json(self, filename):
        json_dct = {}
        with open(filename, 'rb') as fh:
            json_dct = json.loads(fh.read())

        if json_dct is None or "kmers" not in json_dct:
            raise InvalidMotifJsonFileError("File does not seem to be a valid json file for motif results: %s" % json_filename)
        return json_dct

    def get_pvalue_comment_from_rout(self):
        rout_filename = self.get_rout_filename()
        with open(rout_filename, 'rb') as fh:
            for ln in fh:
                ln = unicode(ln, 'utf-8')
                if 'MLE of the p-value' in ln:
                    self.pvalue_comment = ln.strip()
                    break

        if not self.pvalue_comment:
            raise ParsingMotifResultError('Cannot find P-value comment: %s' % rout_filename)

    def load_template_json(self):
        json_filename = os.path.join(self.finder_outdir, 'motif-00000.stdout')
        self.template_json = self.load_json(json_filename)

    def extract_and_write_kmers(self):
        self.kmer_lst = [l[1] for l in self.template_json['kmers'].values()]
        self.kmer_filename = os.path.join(self.statsig_dir, '%s.kmers' % self.finder_id)
        self.write_kmers_file(self.kmer_lst, self.kmer_filename)

    def write_kmers_file(self, curr_kmer_lst, curr_kmer_filename):
        self.logger.info('Writing kmers: %s' % curr_kmer_filename)
        with open(curr_kmer_filename, 'wb') as fh:
            for kmer in curr_kmer_lst:
                print('>FASTA header', file=fh)
                print(kmer, file=fh)

    def get_rscript_text(self, template_score, null_scores_path):
        params = {
                'pval_r_path' : self.pval_r_path,
                'null_scores_path' : null_scores_path,
                'template_score' : template_score,
                }

        rscript_text = """
source("%(pval_r_path)s")
library(MASS)
sample<-scan("%(null_scores_path)s")
getConfPvalLat(%(template_score)s, sample, conf=0.1, mins=7, maxs=200)
""" % params

        return rscript_text

    def get_rout_filename(self):
        return os.path.join(self.statsig_dir, '%s.r_out' % self.finder_id)

    def extract_and_write_scores(self):
        score_lst = [None] * (self.nullset_size + 1)
        for i in range(self.nullset_size+1):
            if i == 0:
                json_filename = os.path.join(self.finder_outdir, 'motif-00000.stdout')
            else:
                json_filename = os.path.join(self.finder_outdir, 'null-%05d.stdout' % i )
            json_dct = self.load_json(json_filename)
            score_lst[i] = json_dct['score_ranking_runs']

        #write gm_width008.scores file
        nullscores_filename = os.path.join(self.statsig_dir, '%s.scores' % self.finder_id)
        self.logger.info('Writing null scores: %s' % nullscores_filename)
        with open(nullscores_filename, 'wb') as fh:
            print('\n' . join([str(s) for s in score_lst[1:]]), file=fh)

        #write R file
        rscript_text = self.get_rscript_text(score_lst[0], nullscores_filename)
        rscript_filename = os.path.join(self.statsig_dir, '%s.R' % self.finder_id)
        self.logger.info('Writing R script: %s' % rscript_filename)
        with open(rscript_filename, 'wb') as fh:
            print(rscript_text, file=fh)
        return rscript_filename


class GimsanResultManager(GimsanJob):
    def __init__(self, name, template_file, outdir, config, conf_file, is_overwrite=False, dryrun=False, verbose=False):
        super(GimsanResultManager, self).__init__(outdir, config)

        self.logger = logging.getLogger(self.__class__.__name__)

        self.name = name
        self.template_file = template_file
        self.outdir = outdir
        self.conf_file = conf_file
        self.verbose = verbose
        self.dryrun = dryrun
        self.is_overwrite = is_overwrite

        self.css_outdir = os.path.join(self.outdir, 'css')
        self.js_outdir = os.path.join(self.outdir, 'js')
        self.statsig_dir = os.path.join(self.outdir, 'statsig')
        self.r_path = os.path.expanduser(self.config.get('result', 'r_path'))
        self.check_result_path()
        self.get_finders()

        #column-dependency
        self.column_dependency_exec = os.path.join(self.gimsan_home, 'column_dependency_app/column_dependency.out')
        if not os.path.isfile(self.column_dependency_exec):
            raise Exception('Column-Dependency executable missing: %s' % self.column_dependency_exec)

    def check_result_path(self):
        if not os.path.isdir(self.outdir):
            raise MissingDirError('Missing output directory: %s' % self.outdir)

        if not os.path.isdir(self.statsig_dir):
            os.mkdir(self.statsig_dir)
        elif not self.is_overwrite:
            raise AlreadyExistOutputDirError('Directory already exist: %s' % self.statsig_dir)

        if not os.path.isdir(self.css_outdir):
            os.mkdir(self.css_outdir)
        elif not self.is_overwrite:
            raise AlreadyExistOutputDirError('Directory already exist: %s' % self.css_outdir)


        if not os.path.isdir(self.js_outdir):
            os.mkdir(self.js_outdir)
        elif not self.is_overwrite:
            raise AlreadyExistOutputDirError('Directory already exist: %s' % self.js_outdir)

    def get_all_finder_meta(self):
        lst = []
        for width in self.width_lst:
            lst.append(self.get_gibbsmarkov_meta(width))
        return lst


    def generate_finder_result_list(self):
        finder_meta_lst = self.get_all_finder_meta()

        finder_res_lst = []
        for finder_meta in finder_meta_lst:
            finder_res_lst.append(FinderResult(finder_meta, self.nullset_size, self.statsig_dir, self.outdir, self.gimsan_home))

        rscript_jobs = []
        for finder_res in finder_res_lst:
            rscript_filename = finder_res.extract_and_write_scores()
            cmd = "%s -f %s &>%s" % (self.r_path, rscript_filename, finder_res.get_rout_filename())
            job = {
                    'cmd' : cmd,
                    'job_id' : rscript_filename,
                    }
            rscript_jobs.append(job)

        #run R in parallel
        if self.dryrun:
            for job in rscript_jobs:
                self.logger.info(job['cmd'])
        else:
            import multiprocessing
            pool = multiprocessing.Pool(processes=self.num_proc)
            pool.map(subprocess_exec_func, rscript_jobs)

        #pvalue
        for finder_res in finder_res_lst:
            finder_res.get_pvalue_comment_from_rout()

        #weblogo
        img_dir = os.path.join(self.outdir, 'images')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        for finder_res in finder_res_lst:
            finder_res.extract_and_write_kmers()
            weblogo_filename = os.path.join(img_dir, '%s.png' % finder_res.finder_id)
            weblogo_revcompl_filename = os.path.join(img_dir, '%s.revcompl.png' % finder_res.finder_id)
            finder_res.construct_weblogo(weblogo_filename, weblogo_revcompl_filename)

        #column dependency
        self.compute_column_dependency(finder_res_lst)

        return finder_res_lst

    def compute_column_dependency(self, finder_res_lst):
        coldep_dir = os.path.join(self.outdir, 'coldep')
        if not os.path.isdir(coldep_dir):
            os.mkdir(coldep_dir)

        if self.config.has_option('column_dependency', 'randseed'):
            randseed_param = '-s %d' % self.config.getint('column_dependency', 'randseed')
        else:
            randseed_param = ''

        job_lst = []
        for finder_res in finder_res_lst:
            coldep_fileroot = '%s.coldep' % finder_res.finder_id
            stdout_fn = coldep_fileroot + ".txt"
            stderr_fn = coldep_fileroot + ".stderr"
            finder_res.coldep_outfile = os.path.join('coldep', stdout_fn)
            cmd = "%s -fsa %s %s 1>%s 2>%s" % (
                    self.column_dependency_exec,
                    finder_res.kmer_filename,
                    randseed_param,
                    os.path.join(coldep_dir, stdout_fn),
                    os.path.join(coldep_dir, stderr_fn))
            job = {
                    'cmd' : cmd,
                    'job_id' : coldep_fileroot,
                    }
            job_lst.append(job)

        #run R in parallel
        if self.dryrun:
            for job in job_lst:
                self.logger.info(job['cmd'])
        else:
            import multiprocessing
            pool = multiprocessing.Pool(processes=self.num_proc)
            pool.map(subprocess_exec_func, job_lst)

        for finder_res in finder_res_lst:
            full_path = os.path.join(self.outdir, finder_res.coldep_outfile)
            with open(full_path, 'rb') as fh:
                for ln in fh:
                    m = re.search(r'for statistically significant pairs \((\d+) pairs\)', ln)
                    if m:
                        finder_res.coldep_num_pairs = int(m.group(1))
                        break

            if finder_res.coldep_num_pairs is None:
                raise Exception('Unable to find statistically significant pairs from %s' % full_path)

            if finder_res.coldep_num_pairs == 0:
                finder_res.coldep_btn_style = 'btn-default'
            else:
                finder_res.coldep_btn_style = 'btn-success'


    def generate_html(self):
        finder_res_lst = self.generate_finder_result_list()
        gm_finder0 = self.finder_lst[0]

        tp = None
        out_tp_file = os.path.join(self.gimsan_home, 'misc/output_template.html')
        with open(out_tp_file) as fh:
            tp = template.Template(fh.read())

        if tp is None:
            raise MissingFileError('Unable to generate HTML from template: %s' % out_tp_file)
        output_html = tp.generate(
                experiment_name = self.name,
                config_filename = os.path.join('../meta', os.path.basename(self.conf_file)),
                fsa_filename = os.path.basename(self.template_file),
                nullset_size = self.nullset_size,
                per_seq_model_comment = gm_finder0.get_per_seq_model(),
                stop_crit_comment = gm_finder0.get_stop_crit(),
                rapid_conv = gm_finder0.get_rapid_conv(),
                double_strand_comment = 'yes' if gm_finder0.get_double_strand() else 'no',
                markov_order = gm_finder0.markov_order,
                genomic_file_comment = self.genome_filename if self.genome_filename else 'input FASTA file',
                finder_res_lst = finder_res_lst,
                )

        output_html_file = os.path.join(self.outdir, 'output.html')
        self.logger.info('Writing HTML file to: %s' % output_html_file)
        with open(output_html_file, 'wb') as fh:
            print(output_html, file=fh)

        self.copy_html_assets()

    def copy_html_assets(self):
        lst = [
                (os.path.join(self.gimsan_home, 'misc', 'css', "bootstrap.min.css"), self.css_outdir),
                (os.path.join(self.gimsan_home, 'misc', 'js', "bootstrap.min.js"), self.js_outdir),
                (os.path.join(self.gimsan_home, 'misc', 'js', "jquery-1.10.2.min.js"), self.js_outdir),
                ]
        for l in lst:
            os.system('cp -v %s %s' % (l[0], l[1]))

def subprocess_exec_func(job):
    import logging
    logging.info('(%s): %s' % (job['job_id'], job['cmd']))
    ret_code = subprocess.call(job['cmd'], shell=True)


if __name__ == '__main__':
    benchmark = time_benchmark.Benchmark()

    #defaults
    description = """
Generate GIMSAN result
"""

    epilog = """
Examples:

%(prog)s --dir=testout -v
"""

    argp = ArgumentParserPlus(description=description, epilog=epilog)
    argp.add_argument('--dir', required=True, help="main output directory used with gimsan_submit.py")
    argp.add_argument('--overwrite', action="store_true", help="")
    argp.add_argument('--dryrun', action="store_true", help="")
    argp.add_argument('-v', '--verbose', action='store_true')
    args = argp.parse_args()

    import logging
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    simple_stream_logging(log_level)
    args.dir = os.path.expanduser(args.dir)

    conf_file = BatchExperiment.get_conf_file(args.dir)
    batch_exp = BatchExperiment(conf_file, args.overwrite)

    for exp in batch_exp.experiments:
        gr_manager = GimsanResultManager(
                exp['name'],
                exp['fasta_file'],
                exp['outdir'],
                batch_exp.config,
                batch_exp.conf_file,
                is_overwrite = args.overwrite,
                dryrun = args.dryrun,
                verbose = args.verbose)
        gr_manager.generate_html()

    benchmark.print_time(sys.stderr)
