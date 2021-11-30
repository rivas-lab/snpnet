####################################################################
# A copy of psychometric::CI.Rsq() function written by Thomas D. Fletcher
#
# https://rdrr.io/cran/psychometric/man/CI.Rsq.html
# https://rdrr.io/cran/psychometric/src/R/CI.Rsq.R
#
# As this is the only function we use from psychometric package (GPL2),
# we place a copy here for now.
#
# 2021.11.29
####################################################################
"CI.Rsq" <-
function(rsq, n, k, level=.95)
 {
noma <- 1-level
sersq <- sqrt((4*rsq*(1-rsq)^2*(n-k-1)^2)/((n^2-1)*(n+3)))
zs <- - qnorm(noma/2)
mez <- zs*sersq
lcl <- rsq - mez
ucl <- rsq + mez
mat <- data.frame(Rsq = rsq, SErsq = sersq, LCL = lcl, UCL = ucl)
return(mat)
}
