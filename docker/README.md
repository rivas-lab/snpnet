# Docker image of `snpnet`

We make the Docker image available at [Docker Hub](https://hub.docker.com/repository/docker/yosuketanigawa/snpnet).

## How to use the Docker image?

### Docker

```{bash}
docker run yosuketanigawa/snpnet:latest R

> library(snpnet)
```

### Singularity

```{bash}
singularity pull docker://yosuketanigawa/snpnet:latest
singularity -s exec snpnet_latest.sif R

> library(snpnet)
```

## How to build the Docker image and push it to the Docker Hub?

```{bash}
docker build . -t yosuketanigawa/snpnet:latest
docker push yosuketanigawa/snpnet:latest
```

## Version log

- `v1.4.1`: (`yosuketanigawa/snpnet:v1.4.1`): We use the following versions of the software
  - `R` from `jupyter/r-notebook:2021-11-20`
  - `plink2`: Alpha 2.3 final (24 Jan 2020) Linux 64-bit Intel
  - `snpnet` version 1.4.1
