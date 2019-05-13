import logging
from .colors import (red,
                     end)
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
''' % (red,end,red,end,red,end,red,end,red,end,red,end,__version__))
    else:
        logging.debug('hiding logo')
