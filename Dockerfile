# alona Dockerfile

# Pull base images.
FROM ubuntu:18.04

MAINTAINER Oscar Franzen "p.oscar.franzen@gmail.com"

# Install some requirements
RUN apt-get update && \
    apt-get install -y tzdata && \
    apt-get install -y r-base python3 && \
    apt-get install -y texlive-full && \
    apt-get install -y libxml2-dev libudunits2-dev libcairo2 libcairo2-dev && \
    apt-get install -y libcurl4-openssl-dev libssl-dev && \
    apt-get install -y libssh2-1-dev openssl curl

# R Packages
RUN R -e "install.packages(c('devtools','ggplot2'), \
                           dependencies=TRUE, \
                           repos='http://cran.rstudio.com/')"

# Python3 packages
RUN pip3 install click

# Define working directory.
WORKDIR /alona

COPY alona.py /alona

ENV ALONA_INSIDE_DOCKER Y

ENTRYPOINT ["/alona/alona.py"]
