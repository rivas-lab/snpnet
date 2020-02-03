FROM jupyter/datascience-notebook:latest
# FROM artemklevtsov/r-alpine:latest
# FROM rstudio/r-base:3.6-centos7
# RUN yum -y install curl wget openssl

RUN R -e "install.packages(c('devtools', 'tidyverse', 'data.table', 'glmnet', 'ROCR', 'crochet', 'dotCall64'), repos = 'http://cran.us.r-project.org', dependencies=TRUE)"

# Copy additional files
WORKDIR /opt
COPY misc/zstd-1.4.4.tar.gz .
RUN tar -xzvf zstd-1.4.4.tar.gz
WORKDIR /opt/zstd-1.4.4
RUN make -j install PREFIX=/opt/conda

# Copy additional files
WORKDIR /opt
COPY misc/plink-ng-4a14c82.tar.gz .
# COPY misc/glmnetPlus-7446f97.tar.gz .
COPY misc/snpnet-v.0.2.0.tar.gz .
RUN tar -xzvf plink-ng-4a14c82.tar.gz

RUN R -e "install.packages(c('BGData'), repos = 'http://cran.us.r-project.org', dependencies=TRUE)"
RUN R -e "devtools::install_github('junyangq/glmnetPlus')"

RUN R -e "Sys.setenv(TAR = '/bin/tar'); install.packages(c('plink-ng-4a14c82/2.0/pgenlibr', 'snpnet-v.0.2.0.tar.gz'), repos = NULL, type='source');"

WORKDIR /usr/local/bin
COPY misc/plink2 .

WORKDIR /opt
COPY misc/snpnet_wrapper.sh .
COPY misc/snpnet_wrapper.R .
COPY misc/snpnet_misc.sh .
COPY misc/snpnet_misc.R .

# # RUN apt install curl
# # RUN R -e "install.packages(c('devtools', 'glmnet'), repos = 'http://cran.us.r-project.org', dependencies=TRUE)"
# # install other packages
# # RUN R -e "devtools::install_github('junyangq/glmnetPlus'); install.packages('snpnet-v.0.2.0.tar.gz', repos = NULL, type='source')"
# # RUN R -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('junyangq/glmnetPlus'); install.packages('glmnet_2.0-20.tar.gz', repos = NULL, type='source'); install.packages('snpnet-v.0.2.0.tar.gz', repos = NULL, type='source')"
# # RUN R -e "Sys.setenv(TAR = '/bin/tar'); install.packages('glmnet_2.0-20.tar.gz', repos = NULL, type='source'); install.packages('snpnet-v.0.2.0.tar.gz', repos = NULL, type='source')"

# # RUN apt-get install curl wget openssl
# RUN yum -y install curl wget openssl
# # RUN R -e "install.packages(c('tidyverse'), repos = 'http://cran.us.r-project.org', dependencies=TRUE)"

