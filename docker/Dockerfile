# FROM jupyter/r-notebook:latest
FROM jupyter/r-notebook:2021-11-20

# Install zstandard (zstd) vis conda
RUN conda install -c conda-forge zstandard

# Install R packages from CRAN
RUN R -e "install.packages(c('devtools', 'tidyverse', 'data.table', 'glmnet', 'ROCR', 'crochet', 'dotCall64'), repos = 'http://cran.us.r-project.org', dependencies=TRUE)"

# Install R packages from GitHub
RUN R -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('chrchang/plink-ng', subdir='2.0/pgenlibr');" \
&&  R -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('chrchang/plink-ng', subdir='2.0/cindex');" \
&&  R -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github(c('junyangq/glmnetPlus'))"

# Add plink2 binary
#
# https://www.cog-genomics.org/plink/2.0/
# Alpha 2.3 final (24 Jan 2020) Linux 64-bit Intel
ENV PATH="/opt/bin:${PATH}"
USER root
ADD https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip /opt/bin/plink2_linux_avx2.zip
RUN unzip /opt/bin/plink2_linux_avx2.zip -d /opt/bin \
&&  rm /opt/bin/plink2_linux_avx2.zip
USER $NB_UID

# Install the latest snpnet package
RUN R -e "Sys.setenv(TAR = '/bin/tar'); devtools::install_github('rivas-lab/snpnet@v1.4.1')"
