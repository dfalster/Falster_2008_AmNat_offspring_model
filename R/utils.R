
to.pdf <- function(expr, filename, ..., verbose = TRUE) {
    if (!file.exists(dirname(filename))) 
        dir.create(dirname(filename), recursive = TRUE)
    if (verbose) 
        cat(sprintf("Creating %s\n", filename))
    pdf(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
} 
