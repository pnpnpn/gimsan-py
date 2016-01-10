#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

import time_benchmark
import os
import sys
import json

from argparse import ArgumentParser


if __name__ == '__main__':
    benchmark = time_benchmark.Benchmark()

    description = """
""" 


    epilog = """
Examples:
""" 
    #iters, genome, output-dir, binsize, wndsize

    argp = ArgumentParser(description=description, epilog=epilog)
    argp.add_argument('--json', required=True)
    argp.add_argument('-v', '--verbose', action='store_true')
    args = argp.parse_args()

    import logging
    if args.verbose:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO


    with open(args.json, 'rb') as fh:
        dct = json.loads(fh.read())
        print(json.dumps(dct, sort_keys=True, indent=2))

    benchmark.print_time(sys.stderr)

