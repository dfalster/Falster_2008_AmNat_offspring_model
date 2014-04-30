# plot mut survival curve - smith fretwell curve at optimum

postscript("plots/SF curve 1.eps", horizontal = FALSE, onefile = FALSE, height = 6, 
    width = 6, pointsize = 10, paper = "special")
c_LRI <- 0.927/2
D <- 1
b_LRI <- 0.9519
s <- 0.944

Wa <- 1
xmin <- 10^-4
xmax <- min(Wa, 100)
W0_ESS <- ESS_analytic(Wa, LRI(Wa))

X <- seq.log(10^-3, 2, 1.001)
Y <- s_mut(Wa, W0_ESS * X, W0_ESS, Wa, LRI(Wa))
par(pty = "s")
plot(X, Y, type = "l", log = "", ylim = c(0, 1), lty = "solid", lwd = 2, axes = F, 
    ann = F, xaxs = "i", yaxs = "i", family = "helvetica")
axis(2, at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1), las = 1, tck = 0.012)
axis(2, at = seq(0, 1, 0.5), labels = F, tck = 0.005)
axis(4, at = seq(0, 1, 0.1), labels = F, tck = 0.012)
axis(4, at = seq(0, 1, 0.5), labels = F, tck = 0.005)
axis(1, at = seq(0, 2, 0.5), labels = seq(0, 2, 0.5), las = 1, tck = 0.012)
axis(1, at = seq(0, 2, 0.1), labels = F, tck = 0.005)
axis(3, at = seq(0, 2, 0.5), labels = F, las = 1, tck = 0.012)
axis(3, at = seq(0, 2, 0.1), labels = F, tck = 0.005)
mtext(expression(paste("Survival to maturity,  ", S[alpha], "(", W[paste(0, ",m")], 
    ",", W[paste(0, ",r*")], ")")), side = 2, line = 3, outer = F, at = NA, adj = 0.5, 
    cex = 1.1)
mtext(expression(paste("Offspring mass of mutanat relative to resdient, ", W[paste(0, 
    ",m")]/W[paste(0, ",r*")])), side = 1, line = 3, outer = F, at = NA, adj = 0.5, 
    cex = 1.1)

# add cost curve
S_RES <- s_mut(Wa, W0_ESS, W0_ESS, Wa, LRI(Wa))
points(1, S_RES, type = "p", pch = 20, cex = 1.5)
curve(x * (S_RES), from = 0, to = 2, add = T, lty = "solid", lwd = 1.5)
points(c(0, 1), c(S_RES, S_RES), type = "l", lty = "dotted", lwd = 1)
points(c(1, 1), c(0, S_RES), type = "l", lty = "dotted", lwd = 1)

# add new plot
points(X, Y, type = "l", lty = "solid", lwd = 2, col = "blue")
dev.off() 
