#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import re
import os
import sys
import json
import time
import ConfigParser

from gibbslib.simple_logging import *
from gibbslib.gimsan_exception import *

class BatchExperiment(object):
    def __init__(self, conf_file, is_overwrite=False):

        if not os.path.isfile(conf_file):
            raise Exception('Missing config file: %s' % conf_file)

        config = ConfigParser.SafeConfigParser()
        config.optionxform=str #preserves upper/lower case for option names
        config.read([conf_file])

        self.conf_file = conf_file
        self.config = config
        self.main_output_dir = os.path.expanduser(self.config.get('common', 'main_output_dir'))

        self.meta_dir = os.path.join(self.main_output_dir, 'meta/')
        self.success_file = os.path.join(self.meta_dir, 'SUCCESS.json')
        self.is_overwrite = is_overwrite


        experiment_dct = {}
        for name in self.config.options('fasta'):
            if re.search(r'\W', name):
                raise InvalidConfigParamError('Name of experiment in [fasta] cannot be non-word char: %s' % name)
            filename = os.path.expanduser(self.config.get('fasta', name))
            experiment_dct[name] = {
                    'name' : name,
                    'fasta_file' : filename,
                    'outdir' : os.path.join(self.main_output_dir, name),
                    }
        self.experiments = sorted(experiment_dct.values(), key=lambda d : d['name'])

    def mkdir_submit_directories(self):
        if not os.path.isdir(self.main_output_dir):
            os.mkdir(self.main_output_dir)
        elif not self.is_overwrite:
            raise AlreadyExistOutputDirError('Directory already exist: %s' % self.main_output_dir)
 
    @staticmethod
    def get_conf_file(main_output_dir):
        if not os.path.isdir(main_output_dir):
            raise Exception('Incorrect result directory to obtain results: %s' % main_output_dir)
            
        success_file = os.path.join(main_output_dir, 'meta', 'SUCCESS.json')

        json_dct = None
        with open(success_file, 'rb') as fh:
            json_dct = json.loads(fh.read())
        if json_dct is None:
            raise Exception('Cannot load json file %s' % success_file)
        
        return os.path.join(main_output_dir, 'meta', json_dct['conf_file'])


    def copy_conf_file(self):
        if not os.path.isdir(self.meta_dir):
            os.mkdir(self.meta_dir)
        elif not self.is_overwrite:
            raise AlreadyExistOutputDirError('Directory already exist: %s' % self.meta_dir)
 
        os.system('cp -v %s %s' % (self.conf_file, self.meta_dir))

    def flag_success(self):
        if os.path.isfile(self.success_file) and not self.is_overwrite:
            raise Exception('Success flag already exist: %s' % self.success_file)

        json_dct = {
            'status' : 'success',
            'conf_file' : os.path.basename(self.conf_file),
            }
        with open(self.success_file, 'wb') as fptr:
            print(json.dumps(json_dct, sort_keys=True), file=fptr)







            

