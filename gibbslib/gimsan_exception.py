#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

class AlreadyExistOutputDirError(Exception):
    pass

class InvalidMotifJsonFileError(Exception):
    pass

class ParsingMotifResultError(Exception):
    pass

class InvalidConfigParamError (Exception):
    pass

class MissingDirError(Exception):
    pass


class MissingFileError(Exception):
    pass

