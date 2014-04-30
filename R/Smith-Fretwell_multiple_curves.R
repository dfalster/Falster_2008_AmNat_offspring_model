
# -------------------------------------------------------------------------------
postscript("plots/SF-mulitple-curve.eps", horizontal = FALSE, onefile = FALSE, height = 6, 
    width = 6, pointsize = 10, paper = "special")
c_LRI <- 0.927/2
D <- 1
b_LRI <- 0.9519
s <- 0.944
xmin <- 10^-4
xmax <- min(Wa, 100)
par(pty = "s")
plot(1:2, 1:2, type = "n", log = "", ylim = c(0, 1), xlim = c(1e-04, 0.8), axes = F, 
    ann = F, xaxs = "i", yaxs = "i", family = "helvetica")
box()
mtext(expression(paste("Survival to maturity,  ", S[alpha], "(", W[paste(0, ",m")], 
    ",", W[paste(0, ",r*")], ")")), side = 2, line = 3, outer = F, at = NA, adj = 0.5, 
    cex = 1.1)
mtext(expression(paste("Offspring mass of mutant, ", W[paste(0, ",m")])), side = 1, 
    line = 3, outer = F, at = NA, adj = 0.5, cex = 1.1)


for (i in c(0.25, 1.2, 2)) {
    Wa <- i
    W0 <- ESS_analytic(Wa, LRI(Wa))
    X <- seq.log(0.001 * W0, 1, 1.01)
    Y <- s_mut(Wa, X, W0, Wa, LRI(Wa))
    points(X, Y, type = "l", lty = "solid", lwd = 2)
    axis(1, at = c(W0), labels = c(expression(paste(W[paste(0, ",r*")]))), las = 1, 
        tck = 0.012)
    
    # add cost curve
    S_RES <- s_mut(Wa, W0, W0, Wa, LRI(Wa))
    
    points(W0, S_RES, type = "p", pch = 20, cex = 1.5)
    curve(x * (S_RES/W0), from = 1e-06, to = 2 * W0, add = T, lty = "solid", lwd = 1.5)
    points(c(0, W0), c(S_RES, S_RES), type = "l", lty = "dotted", lwd = 1)
    points(c(W0, W0), c(0, S_RES), type = "l", lty = "dotted", lwd = 1)
}

dev.off() 
