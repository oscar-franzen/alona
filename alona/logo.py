import logging

import alona
from .constants import (RED, END)

def show_logo(nologo):
    if not nologo:
        print('''
  %s##%s   #       ####  #    #   ##   
 %s#  #%s  #      #    # ##   #  #  #  
%s#    #%s #      #    # # #  # #    # 
%s######%s #      #    # #  # # ###### 
%s#    #%s #      #    # #   ## #    # 
%s#    #%s ######  ####  #    # #    #  v.%s
''' % (RED, END, RED, END, RED, END, RED, END, RED, END, RED, END, alona.__version__))
    else:
        logging.debug('hiding logo')
