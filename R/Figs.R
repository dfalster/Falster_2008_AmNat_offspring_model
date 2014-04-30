# returns vector from lo to hi with multiplication steps of incr
seq.log <- function(lo, hi, incr) {
    temp <- NULL
    while (lo < hi) {
        temp <- c(temp, lo)
        lo <- lo * incr
    }
    if (max(temp) < hi)
        temp <- c(temp, hi)
    temp
}

load.data <- function(){
	mammals <- read.table("data/mammals.txt", header = TRUE, sep = "\t", strip.white = TRUE,
        blank.lines.skip = TRUE)
	plants <-  read.table("data/plants.txt", header = TRUE, sep = "\t", strip.white = TRUE,
        blank.lines.skip = TRUE)
	list(mammals=mammals, plants=plants)
}

Fig2 <- function() {

    data <- load.data()

    par(mfrow = c(1, 2))

    # 1. Annual reproductive output
    # MAMMALS
    i <- !is.na(data$mammals$Wa * data$mammals$A)
    X <- data$mammals$Wa[i]
    Y <- data$mammals$A[i]
    Y_ax <- seq.log(10^-6, 10^4, 10)
    X_ax <- seq.log(10^-6, 10^4, 10)
    X_lab <- Y_lab <- c(expression(10^-6), expression(10^-5), expression(10^-4), expression(10^-3),
        expression(10^-2), expression(10^-1), expression(10^0), expression(10^1),
        expression(10^2), expression(10^3), expression(10^4))

    par(pty = "s")
    plot(X, Y, type = "p", log = "xy", axes = FALSE, ann = FALSE, ylim = c(min(Y_ax),
        max(Y_ax)), xlim = c(min(X_ax), max(X_ax)), col = "darkgray", pch = 19, cex = 1.2,
        xaxs = "i", yaxs = "i", family = "helvetica")
    points(X, Y, type = "p", pch = 1, cex = 1.3, lwd = 0.25)

    axis(2, at = Y_ax, labels = Y_lab, las = 1, tck = 0.012, , cex.axis = 0.8)
    axis(4, at = Y_ax, labels = FALSE, tck = 0.012)
    axis(1, at = X_ax, labels = X_lab, las = 1, tck = 0.012, cex.axis = 0.8)
    axis(3, at = X_ax, labels = FALSE, las = 1, tck = 0.012)
    mtext("Annaul reproductive investment (kg)", side = 2, line = 3, outer = FALSE,
        at = NA, adj = 0.5, cex = 1.1)
    mtext("Adult size (kg)", side = 1, line = 3, outer = FALSE, at = NA, adj = 0.5,
        cex = 1.1)

    # fit slope
    fit <- slope.test(log10(Y), log10(X), test.value = 0.75, method = "SMA")
    Int <- 10^(mean(log10(Y)) - fit$b * mean(log10(X)))
    curve(Int * x^fit$b, min(X), max(X), add = TRUE, lwd = 1.5, col = "darkgray")

    # LOAD PLANT DATA
    i <- !is.na(data$plants$Wa * data$plants$RI)

    X <- data$plants$Wa[i]
    Y <- data$plants$RI[i]
    points(X, Y, type = "p", pch = 1, cex = 1.3, lwd = 0.25)
    fit <- slope.test(log10(Y), log10(X), test.value = 0.75, method = "SMA")
    Int <- 10^(mean(log10(Y)) - fit$b * mean(log10(X)))
    curve(Int * x^fit$b, min(X), max(X), add = TRUE, lwd = 1.5)


    # 3. Lifetime RI (mammals)
    i <- !is.na(data$mammals$Wa * data$mammals$E * data$mammals$A)
    X <- data$mammals$Wa[i]
    Y <- data$mammals$E[i] * data$mammals$A[i]

    par(pty = "s")
    plot(X, Y, type = "p", log = "xy", axes = FALSE, ann = FALSE, ylim = c(min(Y_ax),
        max(Y_ax)), xlim = c(min(X_ax), max(X_ax)), col = "darkgray", pch = 19, cex = 1.2,
        xaxs = "i", yaxs = "i", family = "helvetica")
    points(X, Y, type = "p", pch = 1, cex = 1.3, lwd = 0.25)
    axis(2, at = Y_ax, labels = Y_lab, las = 1, tck = 0.012, , cex.axis = 0.8)
    axis(4, at = Y_ax, labels = FALSE, tck = 0.012)
    axis(1, at = X_ax, labels = X_lab, las = 1, tck = 0.012, cex.axis = 0.8)
    axis(3, at = X_ax, labels = FALSE, las = 1, tck = 0.012)
    mtext("Lifetime reproductive investment (kg)", side = 2, line = 3, outer = FALSE,
        at = NA, adj = 0.5, cex = 1.1)
    mtext("Adult size (kg)", side = 1, line = 3, outer = FALSE, at = NA, adj = 0.5,
        cex = 1.1)

    # fit slope
    fit <- slope.test(log10(Y), log10(X), test.value = 1, method = "SMA")
    Int <- 10^(mean(log10(Y)) - fit$b * mean(log10(X)))
    curve(Int * x^fit$b, min(X), max(X), add = TRUE, lwd = 1.5, col = "darkgray")
}

Fig3 <- function() {

    # MAMMALS
    data <- load.data()
    i <- !is.na(data$mammals$Wa * data$mammals$W0)
    X <- data$mammals$Wa[i]
    Y <- data$mammals$W0[i]
    Y_ax <- seq.log(10^-9, 10^3, 10)
    X_ax <- seq.log(10^-7, 10^5, 10)
    Y_lab <- c(expression(10^-9), expression(10^-8), expression(10^-7), expression(10^-6),
        expression(10^-5), expression(10^-4), expression(10^-3), expression(10^-2),
        expression(10^-1), expression(10^0), expression(10^1), expression(10^2),
        expression(10^3))
    X_lab <- c(Y_lab[3:13], expression(10^4), expression(10^5))

    par(pty = "s")
    plot(X, Y, type = "p", log = "xy", axes = FALSE, ann = FALSE, ylim = c(min(Y_ax),
        max(Y_ax)), xlim = c(min(X_ax), max(X_ax)), col = "darkgray", pch = 19, cex = 1.2,
        xaxs = "i", yaxs = "i", family = "helvetica")
    points(X, Y, type = "p", pch = 1, cex = 1.3, lwd = 0.25)

    axis(2, at = Y_ax, labels = Y_lab, las = 1, tck = 0.012, , cex.axis = 0.8)
    axis(4, at = Y_ax, labels = FALSE, tck = 0.012)
    axis(1, at = X_ax, labels = X_lab, las = 1, tck = 0.012, cex.axis = 0.8)
    axis(3, at = X_ax, labels = FALSE, las = 1, tck = 0.012)
    mtext(expression(paste("Offspring size (kg) ", " [log scale]")), side = 2, line = 3,
        outer = FALSE, at = NA, adj = 0.5, cex = 1.1)
    mtext(expression(paste("Adult size (kg) ", " [log scale]")), side = 1, line = 3,
        outer = FALSE, at = NA, adj = 0.5, cex = 1.1)

    # ADD FITTED CURVES
	p=get.pars("mammals")
    points(X, ESS_analytic(X, LRI(X, p), p), type = "l", lty = "solid", lwd = 2, col = "darkgray")

    # LOAD PLANT DATA
    i <- !is.na(data$plants$Wa * data$plants$W0)
    X <- data$plants$Wa[i]
    Y <- data$plants$W0[i]
    points(X, Y, type = "p", pch = 1, cex = 1.3, lwd = 0.25)
	p=get.pars("plants")
     points(X, ESS_analytic(X, LRI(X, p), p), type = "l", lty = "solid", lwd = 2)
    # 1:1
    points(c(10^-7, 10^5), c(10^-7, 10^5), type = "l", lty = "dashed", lwd = 1.5)
}

Fig4 <- function() {
    Wa <- 1
    xmin <- 10^-4
    xmax <- min(Wa, 100)
	p=get.pars("mammals")

    W0_ESS <- ESS_analytic(Wa, LRI(Wa, p), p)
    S_RES <- s_mut(Wa, W0_ESS, W0_ESS, Wa, LRI(Wa, p), p)

    X <- seq.log(10^-3, 2, 1.001)
    Y <- s_mut(Wa, W0_ESS * X, W0_ESS, Wa, LRI(Wa, p), p)
    par(pty = "s")
    plot(X * W0_ESS, Y, type = "l", log = "", ylim = c(0, 1), lty = "solid", lwd = 2,
        axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i", family = "helvetica")
    box()
    axis(1, at = c(W0_ESS), labels = c(expression(paste(W[paste(0, ",r*")]))), las = 1,
        tck = 0.012)
    mtext(expression(paste("Survival to maturity,  ", S[alpha], "(", W[paste(0, ",m")],
        ",", W[paste(0, ",r*")], ")")), side = 2, line = 3, outer = FALSE, at = NA,
        adj = 0.5, cex = 1.1)
    mtext(expression(paste("Offspring mass of mutant, ", W[paste(0, ",m")])), side = 1,
        line = 3, outer = FALSE, at = NA, adj = 0.5, cex = 1.1)

    # add cost curve
    points(W0_ESS, S_RES, type = "p", pch = 20, cex = 1.5)
    curve(x * (S_RES/W0_ESS), from = 0 * W0_ESS, to = 2 * W0_ESS, add = TRUE, lty = "solid",
        lwd = 1.5)
    points(c(0, W0_ESS), c(S_RES, S_RES), type = "l", lty = "dotted", lwd = 1)
    points(c(W0_ESS, W0_ESS), c(0, S_RES), type = "l", lty = "dotted", lwd = 1)
}

Fig5 <- function() {

	p=get.pars("mammals")

    X <- seq.log(0.001, 20, 1.01)
    Y <- (p$s * p$c_LRI)^(1/(1 - p$q)) * exp(2 * (p$r - p$q) * (1 - p$r)/p$r/X/(p$q - 1))
    Y.assym <- (p$s * p$c_LRI)^(1/(1 - p$q))
    plot(X, Y, type = "l", log = "", ylim = c(0, 0.5), lty = "solid", lwd = 1, axes = FALSE,
        ann = FALSE, xaxs = "i", yaxs = "i")
    curve(0 * x + Y.assym, from = 0, to = max(X), lty = "dotted", add = TRUE)
    box()
    axis(1, at = seq(0, 20, 2), labels = seq(0, 20, 2), las = 1, tck = 0.012)
    axis(1, at = seq(0, 20, 1), labels = FALSE, tck = 0.005)
    mtext(expression(paste("Scaling constant, ", sigma[2])), side = 2, line = 3,
        outer = FALSE, at = NA, adj = 0.5, cex = 1.1)
    mtext("Degree of competitive asymmetry, z", side = 1, line = 3, outer = FALSE,
        at = NA, adj = 0.5, cex = 1.1)

    # add extra lines
    for(G in c(0.8, 1.2))
	    points(X, (p$s * p$c_LRI)^(1/(1 - p$q)) * exp(2 * (p$r - p$q) * (1 - G * (p$r - p$q) - p$q)/G/p$r/X/(p$q -
        1)), type = "l", lty = "dashed", lwd = 1)
}

Fig6 <- function() {
    # investigate change in slope as function B and competitive assymetry

    # optimal fit
    data <- load.data()

    Wa <- c(mean(data$mammals$Wa, na.rm=TRUE), max(data$mammals$Wa, na.rm=TRUE))

    plot(NA, type = "n", xlim = c(0, 1), ylim = c(0, 1), lwd = 1, axes = FALSE,
        ann = FALSE, xaxs = "i", yaxs = "i")

    box()
    axis(1, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, tck = 0.012)
    axis(1, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.005)
    axis(3, at = seq(0, 1, 0.2), labels = FALSE, tck = 0.012)

    axis(3, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.005)
    axis(2, at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), las = 1, tck = 0.012)
    axis(2, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.005)
    axis(4, at = seq(0, 1, 0.2), labels = FALSE, tck = 0.012)

    axis(4, at = seq(0, 1, 0.1), labels = FALSE, tck = 0.005)
    mtext("Slope with constant RGR", side = 1, line = 3)
    mtext("Slope with declining RGR", side = 2, line = 3)

	p <- get.pars("mammals")
	p$RGR_switch  = 0

    for (alph_k in c(64, 16, 4, 1, 0.25)) {
        X <- Y <- NULL
        i <- 1

        # set min value of b_LRI tried, root solving fails if too low and cutrade-off interacts with alph_k
        min.b <- 0.2
        if(alph_k==1) min.b <- 0.4
        if(alph_k==0.25) min.b <- 0.6

        for (b_LRI in seq(1, min.b, -0.1)) {
        	p$alph_k <- alph_k
        	p$b_LRI <- b_LRI
        	y <- ESS_numeric(Wa,  LRI(Wa, p), p$k, p)
            Y[i] <-  sign(cor(log10(y),log10(Wa))) * sqrt(var(log10(y))/var(log10(Wa))) # SMA slope = ratio of variance on logged data
            X[i] <- (b_LRI - p$q)/(1 - p$q)
            i <- i + 1
        }
        points(X, Y, type = "l")
    }
    points(c(0, 1), c(0, 1), type = "l", lty = "dashed")

}
