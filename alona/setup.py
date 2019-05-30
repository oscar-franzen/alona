import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="alona-oscarf",
    version="0.0.1",
    author="Oscar Franzen",
    author_email="p.oscar.franzen@gmail.com",
    description="Prediction of cell types from scRNA-seq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/oscar-franzen/alona",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPLv3",
        "Operating System :: OS Independent",
    ],
)
