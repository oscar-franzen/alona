```
  ##   #       ####  #    #   ##   
 #  #  #      #    # ##   #  #  #    A pipeline for cell type prediction from
#    # #      #    # # #  # #    #        single cell RNA sequencing data.
###### #      #    # #  # # ###### 
#    # #      #    # #   ## #    #           "Life is better at the beach"
#    # ######  ####  #    # #    #          
```

# Description
`alona` is a Python-based software pipeline for analysis of single cell RNA sequencing data. `alona` also exists as a parallel cloud-based service [0].

[0] http://alona.panglaodb.se/

# Installation
### Requirements
* Linux (alona should work on MacOS too, but it is untested)
* Python >= 3.6

### From GitHub
```bash
# Clone the repository
git clone https://github.com/oscar-franzen/alona/

# Enter the directory
cd alona

# Install the package
pip3 install .
```

# Usage
```bash
python3 -m alona \
        --mrnafull \
        --dark_bg \
        --hvg 2000 \
        --leiden_res 0.1 \
        --output test \
        --loglevel debug \
        --header yes \
        --minexpgenes 0.001 \
        --nomito input.mat
```

## Contact
* Oscar Franzen <p.oscar.franzen@gmail.com>

## Reference
A manuscript is in preparation.

# License
GPLv3
