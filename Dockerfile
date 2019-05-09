# alona Dockerfile

# Pull base images.
FROM ubuntu:18.04

MAINTAINER Oscar Franzen "p.oscar.franzen@gmail.com"

# Install some requirements
RUN apt-get update && \
    apt-get install -y tzdata && \
    apt-get install -y r-base python3 && \
    apt-get install -y texlive-full
    
RUN R -e "install.packages(c('ggplot2'), \
                           dependencies=TRUE, \
                           repos='http://cran.rstudio.com/')"

# Define working directory.
WORKDIR /alona

COPY alona.py /alona

ENV ALONA_INSIDE_DOCKER Y

ENTRYPOINT ["/alona/alona.py"]
