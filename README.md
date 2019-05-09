# alona
alona is a software pipeline for analysis of single cell RNA sequencing data. It exists as a cloud-based service and as a standalone software, which runs inside a [Docker](https://en.wikipedia.org/wiki/Docker_(software)) container.

Cloud alona: http://alona.panglaodb.se/

## Installing the standalone
First install [Docker](https://en.wikipedia.org/wiki/Docker_(software)) by following instructions from https://docs.docker.com/install/

### Alternative 1: Download a pre-build Docker image

### Alternative 2: Create the Docker image from a `Dockerfile`
1. Create a new directory and place the `Dockerfile` in this directory.

2. To build the image issue this command. It will take a while to download all dependencies:

```bash
docker build --tag=alona .
```

3. Inspect that the image is there:
```bash
docker images
```

## Contact
* Oscar Franzen <p.oscar.franzen@gmail.com>

## Reference
A manuscript is in preparation.
