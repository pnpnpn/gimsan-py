#!/usr/bin/env python

import sys
import logging

def simple_stream_logging(log_level):
    handler = logging.StreamHandler(sys.stderr)
    simple_logging(log_level, handler)

def simple_logging(log_level, handler):
    handler.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s', "%Y%m%d %H:%M:%S")
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)
    logging.getLogger().setLevel(log_level)

