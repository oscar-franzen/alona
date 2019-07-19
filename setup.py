# python3 setup.py sdist
import setuptools

from distutils.core import setup
from distutils.extension import Extension

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="alona-oscarf",
    version="0.1",
    author="Oscar Franzen",
    author_email="p.oscar.franzen@gmail.com",
    description="Analysis pipeline for single cell RNA sequencing data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/oscar-franzen/alona",
    packages=setuptools.find_packages(),
    install_requires=['click>=7.0', 'matplotlib>=3.0.3', 'numpy>=1.16.3',
                      'pandas>=0.24.2', 'scipy>=1.2.1', 'scikit-learn>=0.21.0',
                      'leidenalg>=0.7.0', 'umap-learn>=0.3.9', 'statsmodels>=0.9.0',
                      'python-igraph>=0.7.1', 'seaborn>=0.9.0'],
    include_package_data=True,
    python_requires='>=3.6',
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPLv3",
        "Operating System :: OS Independent",
    ],

    ext_modules=[Extension('annlib', ['./alona/ANN/brute.cpp',
                                 './alona/ANN/bd_tree.cpp',
                                 './alona/ANN/bd_search.cpp',
                                 './alona/ANN/bd_pr_search.cpp',
                                 './alona/ANN/bd_fix_rad_search.cpp',
                                 './alona/ANN/ANN.cpp',
                                 './alona/ANN/kd_tree.cpp',
                                 './alona/ANN/kd_split.cpp',
                                 './alona/ANN/kd_search.cpp',
                                 './alona/ANN/kd_pr_search.cpp',
                                 './alona/ANN/kd_fix_rad_search.cpp',
                                 './alona/ANN/kd_dump.cpp',
                                 './alona/ANN/kd_util.cpp',
                                 './alona/ANN/NN.cc'],
                          include_dirs=['./alona/ANN/', './alona/ANN/ANN/'])]
)
