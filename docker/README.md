# Docker image

Local

```{bash}
docker build .
docker tag  ef83b7a0d37c  yosuketanigawa/201912_test_repo:20191221
docker push yosuketanigawa/201912_test_repo:20191221
```

Sherlock

```{bash}
cd /scratch/groups/mrivas/users/ytanigaw/simg
singularity pull docker://yosuketanigawa/201912_test_repo:20191221
```

```{bash}
docker run yosuketanigawa/201912_test_repo:20191221 /tmp/snpnet_wrapper.sh /Users/yosuketanigawa/repos/rivas-lab/snpnet/inst/extdata/sample /Users/yosuketanigawa/repos/rivas-lab/snpnet/inst/extdata/sample.phe QPHE gaussian test_QPHE
```

--> it seems like tidyverse is not installed..

--> Switch to Rstudio:R
