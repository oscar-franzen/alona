import sys
import logging

from .colors import (red,green,end)

def init_logging(loglevel):
    if loglevel == 'regular':
        _ll = logging.INFO
    elif loglevel == 'debug':
        _ll = logging.DEBUG
        
    logging.basicConfig(filename='alona.log',
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

    return

def log_info(msg):
    logging.info('%s%s%s' % (green,msg,end))
    
def log_error(ab=None,msg='error'):
    logging.error('%s%s%s' % (red,msg,end))
    
    if ab != None:
        ab.cleanup()
        
    sys.exit(1)
