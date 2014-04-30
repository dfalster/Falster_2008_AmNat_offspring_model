# Dan Falster Functions for Optimal Offspring size model last modified Dec 4 2006

# Parameters
get.pars <- function(taxa="mammals"){

    if(taxa=="mammals"){
        D = 1
        b_LRI = 0.9519
        s = 0.944
    } else {     #plants
        D = 2
        b_LRI = 0.632
        s = 0.00677/17
    }
    list(q = 0.2,  # density independent mortality
        r = 0.75,  # thinning mortality
        k = 10^-8,
        alph_c = 2,
        alph_k = 3,  # for assymetric competition model
        c_LRI = 0.927/2,
        D=D,
        b_LRI=b_LRI,
        s=s,
        RGR_switch = 1,  # 1 = assume constant RGR
        g = 0.8 # defines fraction of max size for maturity (only needed if RGR_switch !=1)
 )
}

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

# offspring production - n estimates offspring production through scaling with
# body size, n_obs uses observed values for E and RI
n <- function(w0, lri) {
    lri/w0
}

# Offspring size-number tradeoff
RI <- function(wa, p) {
    p$c_RI * wa^p$b_RI
}

# Rerproductive investment as function adult size
RL <- function(wa, p) {
    p$c_RL * wa^p$b_RL
}

# Rerproductive lifespan as function adult size
LRI <- function(wa, p) {
    p$c_LRI * wa^p$b_LRI
}

# size at thinning onset as function offsprinf size, adult size, and lifetime
# reproductive investment LRI
wt <- function(w0, wa, lri, p) {
    (1/(p$s * lri) * wa^p$r * w0^(1 - p$q) * p$D^(p$k/w0 - p$q))^(1/(p$r - p$q))
}

# surival functions
s_est_phase <- function(w, w0, p) {
    p$s * (w/w0)^(-p$k/w0)
}

s_densindep_phase <- function(w, w0, p) {
    p$s * p$D^(p$q - p$k/w0) * (w0/w)^p$q
}

s_ST_bound <- function(w, w0, wa, lri, p) {
    (w/wa)^(-p$r)/n(w0, lri)
}

# survival to adulthood - for resident only
s_adult <- function(w0, wa, wt, p) {
    p$s * p$D^(p$q - p$k/w0) * wt^(p$r - p$q) * w0^p$q * wa^(-p$r)
}

# survival to size w of mutant in env of residents
s_mut <- function(w, w0.mut, w0.res, wa, lri, p) {
    # standradise array dimensions for delaing with plotting etc
    LEN <- max(length(w), length(w0.mut), length(w0.res), length(wa), length(lri))

    # print(LEN)

    if (length(w) == LEN)
        W <- w
    else
        W <- rep(w[1], LEN)

    if (length(w0.mut) == LEN)
        W0.MUT <- w0.mut
    else
        W0.MUT <- rep(w0.mut[1], LEN)

    if (length(w0.res) == LEN)
        W0.RES <- w0.res
    else
        W0.RES <- rep(w0.res[1], LEN)

    if (length(wa) == LEN)
        WA <- wa
    else
        WA <- rep(wa[1], LEN)

    if (length(lri) == LEN)
        LRI2 <- lri
    else
        LRI2 <- rep(lri[1], LEN)


    # solves for point of thinning onset
    WT <- wt(W0.RES, WA, LRI2, p)

    gam <- gamma(W0.MUT, W0.RES, WA, WT, p)  # accounts for changes in RGR with size

    WT_mut <- gam * WT * W0.MUT/W0.RES

    surv <- rep(0, times = LEN)

    for (i in 1:LEN) {
        # check to ensure offspring not larger than adult!
        if (W[i] > WA[i])
            surv[i] <- 0
        else if (W[i] <= p$D * W0.MUT[i]) { # establishment phase
            surv[i] <- s_est_phase(W[i], W0.MUT[i], p)
        } else { # excess juveniles - crash at onset of thinning
            if (WT_mut[i] < p$D * W0.MUT[i]) {
                s_crash <- s_ST_bound(p$D * W0.RES[i], W0.RES[i], WA[i], LRI2[i], p) /
                    s_est_phase(p$D * W0.RES[i], W0.RES[i], p)  # magnitude of crash
                surv[i] <- p$s *p$ D^(-p$k/W0.MUT[i]) * (s_crash/alpha_ij(W0.MUT[i], W0.RES[i], p)) *
                  (W[i]/p$D/W0.MUT[i])^(-p$r * alpha_ij(W0.MUT[i], W0.RES[i], p))
            } else if (WT_mut[i] > WA[i]) { # too few juveniles - thinning starts after Wa --> density indept mortality only
                surv[i] <- s_densindep_phase(W[i], W0.MUT[i], p)
            } else { # normal case - thinning starts between D*w0 & Wa
                if (W[i] <= WT_mut[i])
                  surv[i] <- s_densindep_phase(W[i], W0.MUT[i], p)
                else
                  surv[i] <- p$s * p$D^(-p$k/W0.MUT[i]) * (gam[i] * WT[i]/p$D/W0.RES[i])^(-p$q) *
                    (W[i] * W0.RES[i]/WT[i]/W0.MUT[i]/gam[i])^(-p$r * alpha_ij(gam[i] *
                      W0.MUT[i], W0.RES[i], p))
            }
        }
    }
    surv
}

# adjustment function for RGR
gamma <- function(w0.mut, w0.res, wa, wt.res, p) {
    gam <- (((wa/wt.res)^0.25 * (1 - (w0.res/w0.mut)^0.25) + (wa/w0.mut)^0.25 - p$g^0.25)/((wa/w0.res)^0.25 -
        p$g^0.25))^4

    if (p$RGR_switch == 1) {
        gam <- 0 * gam + 1
    }
    gam
}

# for calculation of ESS R0 function - fitness
R0 <- function(w0.mut, w0.res, wa, lri, p) {
    n(w0.mut, lri) * s_mut(wa, w0.mut, w0.res, wa, lri, p)
}

# selection gradient
dR0 <- function(w0, wa, lri, p) {
    DERIV <- 10^-6  # percentage difference when taking numerical derivatives

    (R0(w0 * (1 + DERIV), w0, wa, lri, p) - R0(w0 * (1 - DERIV), w0, wa, lri, p))/(w0 *
        2 * DERIV)
}

# solve for ESS with analyictal solution and adult size scaling
ESS_analytic <- function(wa, lri, p) {
    exp(2 * (p$r - p$q) * (1 - p$r)/p$r/p$alph_k/(p$q - 1)) * (p$s * lri/wa^p$q)^(1/(1 - p$q))
}

# solve for ESS by locating zero in selection gradient
ESS_numeric <- function(wa, lri, min, p) {

    Sol <- 0 * wa

    # do for all values of wa if a vector
    for (i in 1:length(wa))
        Sol[i] <- uniroot(dR0, wa[i] * c(min, 1), tol = 10^-9, wa = wa[i], lri = lri[i], p=p)$root

    Sol
}

# competition effect
alpha_ij <- function(w1, w2, p) {
    p$alph_c/(1 + (w1/w2)^p$alph_k)
}
