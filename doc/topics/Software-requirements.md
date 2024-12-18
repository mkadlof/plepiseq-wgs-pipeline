# Software requirements

This pipeline is designed to run on a general purpose GNU/Linux operating system. 
The pipeline is written in [Nextflow](https://www.nextflow.io/docs/latest/index.html), which is a language for writing bioinformatics pipelines. It is designed to be portable and scalable, and can be run on a variety of platforms, including local machines, clusters, and cloud computing environments.

Our pipeline is containerized using [Docker](https://www.docker.com/), which is a platform for developing, shipping, and running applications in containers. Containers allow a developer to package up an application with all the parts it needs, such as libraries and other dependencies, and ship it all out as one package. This makes it easy to deploy the application on any machine that supports Docker.

However, our containers are run differently than usual docker workflow. Containers are run by nextflow, instead of manual execution of docker. Nextflow take care of  mounting volumes and deciding which container should be run and when.

## Compatibility 
Pipeline was tested with following software versions:
 * **Operating system**:
   * Ubuntu 20.04.06 LTS
   * Debian 12.5
 * **Docker**:
   * 24.0.7
   * 26.0.2
 * **Nextflow**:
   * 24.04.4.5917

It is known that pipeline will not work with docker 20.10.5.
