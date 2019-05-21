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
RUN pip3 install scanpy
RUN pip3 install pandas
RUN pip3 install numpy

# Define working directory.
WORKDIR /alona

RUN mkdir /alona/alona/
RUN mkdir /alona/alona/genome
RUN mkdir /alona/alona/irlbpy
RUN mkdir /alona/alona/ANN

# copy the package into the image
COPY ./alona/*.py /alona/alona/
COPY ./alona/irlbpy/*.py /alona/alona/irlbpy/
COPY ./alona/genome/* /alona/alona/genome/
COPY ./alona/ANN/ /alona/alona/ANN

RUN cd /alona/alona/ANN/ && /alona/alona/ANN/compile

ENV ALONA_INSIDE_DOCKER Y

ENTRYPOINT ["/alona/alona.py"]
