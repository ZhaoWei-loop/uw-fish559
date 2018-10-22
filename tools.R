# library(devtools)
# https://github.com/kaskr/TMB_contrib_R

# devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
library("TMBhelper")

# devtools::install_github("kaskr/TMB_contrib_R/TMBdebug")
library("TMBdebug")

# devtools::install_github("kaskr/TMB_contrib_R/TMBphase")
libraru("TMBphase")

### Find 'string' in 'file' and read vector from next line
readVec <- function(string, file){
  txt <- readLines(file)
  skip <- match(string, txt)
  vec <- scan(file, quiet=TRUE, skip=skip, nlines=1)
  return(vec)
}

### Find 'string' in 'file' and read matrix with 'nrow' rows from next line
readMat <- function(string, file, nrow){
  txt <- readLines(file)
  skip <- match(string, txt)
  mat <- as.matrix(read.table(file, skip=skip, nrows=nrow))
  dimnames(mat) <- NULL
  return(mat)
}

## Function to read a basic AD Model Builder fit.
## Use for instance by:
## simple.fit <- readFit('c:/admb/examples/simple')
## Then the object 'simple.fit' is a list containing sub-objects
# 'names', 'est', 'std', 'cor', and 'cov' for all model
# parameters and sdreport quantities.
readFit <- function(file){
	ret <- list()
	parfile <- as.numeric(scan(paste(file,'.par', sep=''), what='', n=16, quiet=TRUE)[c(6,11,16)])
	ret$nopar <- as.integer(parfile[1])
	ret$nlogl <- parfile[2]
	ret$maxgrad <- parfile[3]
	file <- paste(file,'.cor', sep='')
	lin <- readLines(file)
	ret$npar <- length(lin)-2
	ret$logDetHess <- as.numeric(strsplit(lin[1], '=')[[1]][2])
	sublin <- lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!=''])
	ret$names <- unlist(lapply(sublin, function(x) x[2]))
	ret$est <- as.numeric(unlist(lapply(sublin, function(x) x[3])))
	ret$std <- as.numeric(unlist(lapply(sublin, function(x) x[4])))
	ret$cor <- matrix(NA, ret$npar, ret$npar)
	corvec <- unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)]))
	ret$cor[upper.tri(ret$cor, diag=TRUE)] <- as.numeric(corvec)
	ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)]
	ret$cov <- ret$cor*(ret$std%o%ret$std)
	return(ret)
}
