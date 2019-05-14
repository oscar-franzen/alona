import logging
from .colors import (RED, END)
from ._version import __version__

def show_logo(nologo):
    if not nologo:
        print('''
  %s##%s   #       ####  #    #   ##   
 %s#  #%s  #      #    # ##   #  #  #  
%s#    #%s #      #    # # #  # #    # 
%s######%s #      #    # #  # # ###### 
%s#    #%s #      #    # #   ## #    # 
%s#    #%s ######  ####  #    # #    #  v.%s
''' % (RED, END, RED, END, RED, END, RED, END, RED, END, RED, END, __version__))
    else:
        logging.debug('hiding logo')
