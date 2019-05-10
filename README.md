# alona
alona is a software pipeline for analysis of single cell RNA sequencing data. It exists as a cloud-based service and as a standalone software, which runs inside a [Docker](https://en.wikipedia.org/wiki/Docker_(software)) container.

Cloud based version of alona: http://alona.panglaodb.se/

## Installing the standalone
First install [Docker](https://en.wikipedia.org/wiki/Docker_(software)) by following instructions from https://docs.docker.com/install/

Check that docker is working:

```
docker
```

You should see:
```
Usage:	docker [OPTIONS] COMMAND

A self-sufficient runtime for containers

Options:
      --config string      Location of client config files (default "/home/rand/.docker")
  -D, --debug              Enable debug mode
  -H, --host list          Daemon socket(s) to connect to
  -l, --log-level string   Set the logging level ("debug"|"info"|"warn"|"error"|"fatal") (default "info")
      --tls                Use TLS; implied by --tlsverify
      --tlscacert string   Trust certs signed only by this CA (default "/home/rand/.docker/ca.pem")
      --tlscert string     Path to TLS certificate file (default "/home/rand/.docker/cert.pem")
      --tlskey string      Path to TLS key file (default "/home/rand/.docker/key.pem")
      --tlsverify          Use TLS and verify the remote
  -v, --version            Print version information and quit

Management Commands:
  builder     Manage builds
  config      Manage Docker configs
  container   Manage containers
  engine      Manage the docker engine
  image       Manage images
  network     Manage networks
  node        Manage Swarm nodes
  plugin      Manage plugins
  secret      Manage Docker secrets
  service     Manage services
  stack       Manage Docker stacks
  swarm       Manage Swarm
  system      Manage Docker
  trust       Manage trust on Docker images
  volume      Manage volumes
```

### First clone the repository

```bash
git clone https://github.com/oscar-franzen/alona
```

### Alternative 1: Download a pre-build Docker image

### Alternative 2: Create the Docker image from a `Dockerfile`
This alternative may be a little bit slower but should work equally well.

1. Create a new directory and place the `Dockerfile` in this directory.

```bash
cd alona
```

2. To build the image issue this command. It will take a while to download all dependencies:

```bash
docker build --tag=alona .
```

3. Inspect that the image is there:
```bash
docker images
```

## Usage
Create a new container using the installed image and enter it using Docker's interactive mode:
```bash
docker run -it --entrypoint /bin/bash alona

python3 -m alona
```

## Contact
* Oscar Franzen <p.oscar.franzen@gmail.com>

## Reference
A manuscript is in preparation.
