"""
 alona

 Description:
 An analysis pipeline for scRNA-seq data.

 How to use:
 https://github.com/oscar-franzen/alona/

 Contact:
 Oscar Franzen <p.oscar.franzen@gmail.com>
"""

import sys

from .celltypes import AlonaCellTypePred

class AlonaFindmarkers(AlonaCellTypePred):
    """
    Find markers class.
    """

    def __init__(self):
        super().__init__()

    def findMarkers(self):
        print(self.data_norm)
        print(self.leiden_cl)
        print('hello from findMarkers()')
