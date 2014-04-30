
seq_log <- function(lower, upper, multiply) {
    i <- 1
    temp <- NULL
    while (i * lower < upper) {
        temp <- c(temp, lower * i)
        i <- i * multiply
    }
    return(temp)
}

fig.setup <- function(ylim, ylab, y_ax, y_ax2, xlim, xlab, x_ax, x_ax2, cex.A, cex.L, 
    log = "") {
    plot(1:2, 1:2, type = "n", log = log, axes = F, ann = F, xlim = xlim, ylim = ylim, 
        xaxs = "i", yaxs = "i", las = 1)
    
    axis(2, at = y_ax, labels = y_ax, las = 1, tck = 0.03, cex.axis = cex.A, adj = 0.5)
    axis(4, at = y_ax, labels = F, tck = 0.03)
    axis(2, at = y_ax2, labels = F, las = 1, tck = 0.015, cex.axis = cex.A, adj = 0.5)
    axis(4, at = y_ax2, labels = F, tck = 0.015)
    axis(1, at = x_ax, labels = x_ax, las = 1, tck = 0.03, cex.axis = cex.A)
    axis(3, at = x_ax, labels = F, las = 1, tck = 0.03)
    axis(1, at = x_ax2, labels = F, las = 1, tck = 0.015, cex.axis = cex.A, adj = 0.5)
    axis(3, at = x_ax2, labels = F, tck = 0.015)
    box()
    mtext(xlab, side = 1, line = 3, outer = F, at = NA, cex = cex.L)
    mtext(ylab, side = 2, line = 3, outer = F, at = NA, cex = cex.L)
}


to.pdf <- function(expr, filename, ..., verbose = TRUE) {
    if (!file.exists(dirname(filename))) 
        dir.create(dirname(filename), recursive = TRUE)
    if (verbose) 
        cat(sprintf("Creating %s\n", filename))
    pdf(filename, ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
} 
