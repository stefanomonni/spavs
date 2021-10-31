sample2=function(x, size, replace = FALSE, prob = NULL){
  if (missing(size)){ size <- length(x)}
  x[sample.int(length(x), size, replace, prob)]
}
 

