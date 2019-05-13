import logging
from .colors import (red,
                     end)

def show_logo(nologo):
    if not nologo:
        print('''
  %s##%s   #       ####  #    #   ##   
 %s#  #%s  #      #    # ##   #  #  #  
%s#    #%s #      #    # # #  # #    # 
%s######%s #      #    # #  # # ###### 
%s#    #%s #      #    # #   ## #    # 
%s#    #%s ######  ####  #    # #    # 
''' % (red,end,red,end,red,end,red,end,red,end,red,end))
    else:
        logging.debug('hiding logo')
