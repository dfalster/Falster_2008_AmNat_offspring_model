# MS FIG 1
postscript("plots/fig3.eps", horizontal = FALSE, onefile = FALSE, height = 6, width = 6, 
    pointsize = 10, paper = "special")
# LOAD MAMMAL DATA
X <- obs.mam$Wa
Y <- obs.mam$W0
Y_ax <- seq.log(10^-9, 10^3, 10)
X_ax <- seq.log(10^-7, 10^5, 10)
Y_lab <- c(expression(10^-9), expression(10^-8), expression(10^-7), expression(10^-6), 
    expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2), expression(10^-1), 
    expression(10^0), expression(10^1), expression(10^2), expression(10^3))
X_lab <- c(Y_lab[3:13], expression(10^4), expression(10^5))

par(pty = "s")
plot(X, Y, type = "p", log = "xy", axes = F, ann = F, ylim = c(min(Y_ax), max(Y_ax)), 
    xlim = c(min(X_ax), max(X_ax)), col = "red", pch = 19, cex = 1, xaxs = "i", yaxs = "i", 
    family = "helvetica")
axis(2, at = Y_ax, labels = Y_lab, las = 1, tck = 0.012, , cex.axis = 0.8)
axis(4, at = Y_ax, labels = F, tck = 0.012)
axis(1, at = X_ax, labels = X_lab, las = 1, tck = 0.012, cex.axis = 0.8)
axis(3, at = X_ax, labels = F, las = 1, tck = 0.012)
mtext(expression(paste("Offspring size (kg) ", " [log scale]")), side = 2, line = 3, 
    outer = F, at = NA, adj = 0.5, cex = 1.1)
mtext(expression(paste("Adult size (kg) ", " [log scale]")), side = 1, line = 3, 
    outer = F, at = NA, adj = 0.5, cex = 1.1)

# LOAD PLANT DATA
points(obs.plt$Wa, obs.plt$W0, type = "p", pch = 19, cex = 1, col = "blue")

# ------------------------------------------------------------ ADD FITTED CURVES
# - ver 2 --> use different values of D
X <- seq.log(10^-7, 10^5, 2)
c_LRI <- 0.927/2
b_mam <- 0.9519
b_plt <- 0.632
b_step <- (b_mam/b_plt)^(1/10)
s_mam <- 0.944
s_plt <- 0.00677/17
s_step <- (s_mam/s_plt)^(1/10)

# mammals
D <- 1
i <- 1
while (i < 4) {
    b_LRI <- b_mam * b_step^(i - 3)
    s <- s_mam * s_step^(i - 3)
    points(X, ESS_analytic(X, LRI(X)), type = "l", lty = "solid", lwd = 1, col = "red")
    i <- i + 1
}

# best fit
b_LRI <- b_mam
s <- s_mam
points(X, ESS_analytic(X, LRI(X)), type = "l", lty = "solid", lwd = 1)

# plants
D <- 2
i <- 1
N <- length(X)
while (i < 16) {
    b_LRI <- b_plt * b_step^(i - 8)
    s <- s_plt * s_step^(i - 8)
    X2 <- X[1:(N/2)]
    points(X2, ESS_analytic(X2, LRI(X2)), type = "l", lty = "solid", lwd = 1, col = "blue")
    X2 <- X[(N/2):(N + 1)]
    points(X2, ESS_analytic(X2, LRI(X2)), type = "l", lty = "solid", lwd = 1, col = "blue")
    i <- i + 1
}

# best fit
b_LRI <- b_plt
s <- s_plt
X2 <- X[1:(N/2)]
points(X2, ESS_analytic(X2, LRI(X2)), type = "l", lty = "solid", lwd = 1)
X2 <- X[(N/2):(N + 1)]
points(X2, ESS_analytic(X2, LRI(X2)), type = "l", lty = "solid", lwd = 1)
# 1:1
points(X, X, type = "l", lty = "solid", lwd = 2)

dev.off() 
