
# don't include *out components
#out <- out[-c(1,4,8,9)] 

gelmanRubin <- function( out,nChains,thin,start ){
  library(coda); library(rjags)
nChains <- 5
thin <- 1
start <- 1
ans <- vector("list", nChains)
for (ch in 1:nChains) {
  ans.ch <- vector("list", length(out))
  vnames.ch <- NULL
  for (i in seq(along = out)) {
    varname <- names(out)[[i]]
    d <- dim(out[[i]])
    if (length(d) < 3) {
      stop("Invalid dimensions for sampled output")
    }
    vardim <- d[1:(length(d) - 2)]
    nvar <- prod(vardim)
    niter <- d[length(d) - 1]
    nchain <- d[length(d)]
    values <- as.vector(out[[i]])
    var.i <- matrix(NA, nrow = niter, ncol = nvar)
    for (j in 1:nvar) {
      var.i[, j] <- values[j + (0:(niter - 1)) * nvar +
        (ch - 1) * niter * nvar]
    }
    vnames.ch <- c(vnames.ch, rjags:::coda.names(varname, vardim))
    ans.ch[[i]] <- var.i
  }
  ans.ch <- do.call("cbind", ans.ch)
  colnames(ans.ch) <- vnames.ch
  ans[[ch]] <- mcmc(ans.ch, start = start, thin = thin) } 

o <- mcmc.list(ans)
gD <- gelman.diag(o) 
return(gD)
}