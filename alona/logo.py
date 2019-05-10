import logging

def show_logo(nologo):
    if not nologo:
        print('''
  ##   #       ####  #    #   ##   
 #  #  #      #    # ##   #  #  #  
#    # #      #    # # #  # #    # 
###### #      #    # #  # # ###### 
#    # #      #    # #   ## #    # 
#    # ######  ####  #    # #    # 
''')
    else:
        logging.debug('hiding logo')
