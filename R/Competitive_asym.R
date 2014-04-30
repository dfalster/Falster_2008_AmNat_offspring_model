postscript("plots/z.eps", horizontal = FALSE, onefile = FALSE, height = 6, width = 6, 
    pointsize = 10, paper = "special")

c_LRI <- 0.927/2
s <- 0.944

X <- seq.log(0.001, 20, 1.01)
Y <- (s * c_LRI)^(1/(1 - q)) * exp(2 * (r - q) * (1 - r)/r/X/(q - 1))
Y.assym <- (s * c_LRI)^(1/(1 - q))
plot(X, Y, type = "l", log = "", ylim = c(0, 0.5), lty = "solid", lwd = 1, axes = F, 
    ann = F, xaxs = "i", yaxs = "i")
curve(0 * x + Y.assym, from = 0, to = max(X), lty = "dotted", add = T)

box()
axis(1, at = seq(0, 20, 2), labels = seq(0, 20, 2), las = 1, tck = 0.012)
axis(1, at = seq(0, 20, 1), labels = F, tck = 0.005)
mtext(expression(paste("Scaling constant, ", sigma[2])), side = 2, line = 3, outer = F, 
    at = NA, adj = 0.5, cex = 1.1)
mtext("Degree of competitive asymmetry, z", side = 1, line = 3, outer = F, at = NA, 
    adj = 0.5, cex = 1.1)

# add extra lines
G <- 0.8
points(X, (s * c_LRI)^(1/(1 - q)) * exp(2 * (r - q) * (1 - G * (r - q) - q)/G/r/X/(q - 
    1)), type = "l", lty = "dashed", lwd = 1)
G <- 1.2
points(X, (s * c_LRI)^(1/(1 - q)) * exp(2 * (r - q) * (1 - G * (r - q) - q)/G/r/X/(q - 
    1)), type = "l", lty = "dashed", lwd = 1)

dev.off() 
