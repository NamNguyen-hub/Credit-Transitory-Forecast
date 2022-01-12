library(matconv)

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

filepath1 = "one_sided_hp_filter_kalman.m"
out <- mat2r(inMat = filepath1)
out$rCode

write.table(out$rCode, "HP1sided.R", sep=',', quote=FALSE, row.names=FALSE, col.names=FALSE)



nargin <- function() {
  if(sys.nframe()<2) stop("must be called from inside a function")
  length(as.list(sys.call(-1)))-1
}