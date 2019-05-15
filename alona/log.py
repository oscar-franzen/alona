"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Details:
 https://alona.panglaodb.se/

 Contact:
 Oscar Franzen, <p.oscar.franzen@gmail.com>
"""

import sys
import logging

from .colors import (RED, GREEN, END)

def init_logging(loglevel, logfile):
    if loglevel == 'regular':
        _ll = logging.INFO
    elif loglevel == 'debug':
        _ll = logging.DEBUG

    logging.basicConfig(filename=logfile,
                        level=_ll,
                        format='%(asctime)s %(message)s')

    console = logging.StreamHandler()
    # Logging messages which are less severe than level will be ignored.
    # Level         Numeric value
    # CRITICAL      50
    # ERROR         40
    # WARNING       30
    # INFO          20
    # DEBUG         10
    # NOTSET        0
    console.setLevel(_ll)
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def log_info(msg):
    logging.info('%s%s%s', GREEN, msg, END)

def log_error(abobj=None, msg='error'):
    logging.error('%s%s%s', RED, msg, END)

    if abobj is not None:
        abobj.cleanup()

    sys.exit(1)
