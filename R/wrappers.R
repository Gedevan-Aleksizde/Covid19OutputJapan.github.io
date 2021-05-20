# Matlab fundamental functions

zeros <- function(...){
  array(0, c(...))
}

ones <- function(...){
  array(1, c(...))
}
flip <- function(x){
  rev(x)
}

nnz <- function(x){
  sum(x != 0)
}

saveas <- function(path, dir = "", device = png, dev.args = NULL){
  do.call(dev.copy, append(dev.args, list(device = device, filename = file.path(dir, path))))
  dev.off()
}
