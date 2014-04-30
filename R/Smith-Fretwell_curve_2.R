
# -------------------------------------------------------------------------------
# X axis real numbers
postscript("plots/SF curve 2.eps", horizontal = FALSE, onefile = FALSE, height = 6, 
    width = 6, pointsize = 10, paper = "special")
c_LRI <- 0.927/2
D <- 1
b_LRI <- 0.9519
s <- 0.944
Wa <- 1
xmin <- 10^-4
xmax <- min(Wa, 100)
W0_ESS <- ESS_analytic(Wa, LRI(Wa))
S_RES <- s_mut(Wa, W0_ESS, W0_ESS, Wa, LRI(Wa))

X <- seq.log(10^-3, 2, 1.001)
Y <- s_mut(Wa, W0_ESS * X, W0_ESS, Wa, LRI(Wa))
par(pty = "s")
plot(X * W0_ESS, Y, type = "l", log = "", ylim = c(0, 1), lty = "solid", lwd = 2, 
    axes = F, ann = F, xaxs = "i", yaxs = "i", family = "helvetica")
box()
axis(1, at = c(W0_ESS), labels = c(expression(paste(W[paste(0, ",r*")]))), las = 1, 
    tck = 0.012)
mtext(expression(paste("Survival to maturity,  ", S[alpha], "(", W[paste(0, ",m")], 
    ",", W[paste(0, ",r*")], ")")), side = 2, line = 3, outer = F, at = NA, adj = 0.5, 
    cex = 1.1)
mtext(expression(paste("Offspring mass of mutant, ", W[paste(0, ",m")])), side = 1, 
    line = 3, outer = F, at = NA, adj = 0.5, cex = 1.1)

# add cost curve
points(W0_ESS, S_RES, type = "p", pch = 20, cex = 1.5)
curve(x * (S_RES/W0_ESS), from = 0 * W0_ESS, to = 2 * W0_ESS, add = T, lty = "solid", 
    lwd = 1.5)
points(c(0, W0_ESS), c(S_RES, S_RES), type = "l", lty = "dotted", lwd = 1)
points(c(W0_ESS, W0_ESS), c(0, S_RES), type = "l", lty = "dotted", lwd = 1)

dev.off() 
