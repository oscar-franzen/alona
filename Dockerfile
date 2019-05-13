# alona Dockerfile

# Pull base images.
FROM ubuntu:18.04

MAINTAINER Oscar Franzen "p.oscar.franzen@gmail.com"

# Install some requirements
RUN apt-get update && \
    apt-get install -y tzdata && \
    apt-get install -y r-base python3 python3-pip && \
    apt-get install -y texlive-full && \
    apt-get install -y libxml2-dev libudunits2-dev libcairo2 libcairo2-dev && \
    apt-get install -y libcurl4-openssl-dev libssl-dev && \
    apt-get install -y libssh2-1-dev openssl curl

# R Packages
RUN R -e "install.packages(c('devtools','ggplot2','data.table'), \
                           dependencies=TRUE, \
                           repos='http://cran.rstudio.com/')"

# Python3 packages
RUN pip3 install click

# Define working directory.
WORKDIR /alona

RUN mkdir /alona/alona/
RUN mkdir /alona/alona/genome

# copy the package into the image
COPY ./alona/*.py /alona/alona/
COPY ./alona/genome/human_to_mouse_1_to_1_orthologs.tsv /alona/alona/genome/

ENV ALONA_INSIDE_DOCKER Y

ENTRYPOINT ["/alona/alona.py"]
