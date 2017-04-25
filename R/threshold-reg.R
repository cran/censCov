thlm <- function(formula, data, cens = NULL, subset, method = "cc", B = 0,
                 x.upplim = NULL, t0 = NULL, control = thlm.control()) {
    Call <- match.call()
    mnames <- c("", "formula", "data", "subset", "cens")
    cnames <- names(Call)
    cnames <- cnames[match(mnames, cnames, 0)]
    mcall <- Call[cnames]
    mcall[[1]] <- as.name("model.frame")
    obj <- eval(mcall, parent.frame())
    y <- as.numeric(model.extract(obj, "response"))
    cens <- as.numeric(model.extract(obj, "cens"))
    if (is.null(cens) | formula[[3]] == 1 | all(cens == 1)) {
        out <- lm(formula, data = obj)
        method <- "lm"
    } else {
        mat <- as.matrix(model.matrix(attr(obj, "terms"), obj))
        mat <- mat[, apply(mat, 2, function(x) !all(x == 1))]
        u <- mat[,1]
        z <- as.matrix(mat[,-1])
        if (method %in% c("cc", "complete cases")) {
            ## out <- lm(formula, subset = cens == 1)
            out <- cc.reg(y, u, cens, z)
            names(out$a2) <- names(out$a2.sd) <- colnames(mat)[-1]
            method <- "clm"
        }
        if (method %in% c("reverse", "rev", "reverse survival")) {
            out <- rev.surv.reg(y, u, cens, z)
            method <- "rev"
        }
        if (method %in% c("deletion", "deletion threshold", "dt")) {
            op2 <- par(mar = c(3.5, 3.5, 2.5, 2.5))
            out <- threshold.reg.m1(y, u, cens, z, x.upplim, t0, control)
            names(out$a2) <- names(out$a2.sd) <- colnames(mat)[-1]
            method <- "dt"
            out$a1.sd <- NA
            par(op2)
        }
        if (method %in% c("complete", "complete threshold", "ct")) {
            op2 <- par(mar = c(3.5, 3.5, 2.5, 2.5))
            out <- threshold.reg.m2(y, u, cens, z, x.upplim, t0, control)
            names(out$a2) <- names(out$a2.sd) <- colnames(mat)[-1]
            method <- "ct"
            out$a1.sd <- NA
            par(op2)
        }
        if (method == "all") {
            out <- NULL
            out$cc <- cc.reg(y, u, cens, z)
            names(out$cc$a2) <- names(out$cc$a2.sd) <- colnames(mat)[-1]
            out$rev <- rev.surv.reg(y, u, cens, z)
            op2 <- par(mfrow = c(2, 1), mar = c(3.5, 3.5, 2.5, 2.5))
            out$dt <- threshold.reg.m1(y, u, cens, z, x.upplim, t0, control)
            names(out$dt$a2) <- names(out$dt$a2.sd) <- colnames(mat)[-1]
            out$dt$a1.sd <- NA
            out$ct <- threshold.reg.m2(y, u, cens, z, x.upplim, t0, control)
            names(out$ct$a2) <- names(out$ct$a2.sd) <- colnames(mat)[-1]
            out$ct$a1.sd <- NA
            par(op2)
        }
    }
    if (B > 0 & method %in% c("dt", "ct", "all")) {
        a1.bp <- a1.bp2 <- rep(NA, B)
        k <- 1
        t0 <- out$threshold
        n <- length(y)
        while(k <= B) {
            bp.id <- sample(1:n, n, replace = TRUE)
            y.bp <- y[bp.id]
            u.bp <- u[bp.id]
            z.bp <- as.matrix(z[bp.id,])
            delta.bp <- cens[bp.id]
            temp <- diff(sort(u.bp))
            temp2 <- sort(temp[temp>0])
            u.bp <- u.bp + runif(n, 0, min(temp2) / n)
            if((sum(u.bp[delta.bp == 1] <= t0) > 0) & (sum(u.bp > t0) > 0)) {
                if (method == "dt") 
                    a1.bp[k] <- threshold.reg.m1(y.bp, u.bp, delta.bp, z.bp, x.upplim, t0)$a1
                if (method == "ct")
                    a1.bp[k] <- threshold.reg.m2(y.bp, u.bp, delta.bp, z.bp, x.upplim, t0)$a1
                k <- k + 1
            }
            if (method == "all" & (sum(u.bp[delta.bp == 1] <= out$dt$threshold) > 0) & (sum(u.bp > out$dt$threshold) > 0) &
                (sum(u.bp[delta.bp == 1] <= out$ct$threshold) > 0) & (sum(u.bp > out$ct$threshold) > 0)) {
                a1.bp[k] <- threshold.reg.m1(y.bp, u.bp, delta.bp, z.bp, x.upplim, out$dt$threshold)$a1
                a1.bp2[k] <- threshold.reg.m2(y.bp, u.bp, delta.bp, z.bp, x.upplim, out$ct$threshold)$a1
                k <- k + 1
            }
        }
        if (method == "all") {
            out$dt$a1.sd <- sd(a1.bp)
            out$ct$a1.sd <- sd(a1.bp2)
        } else { out$a1.sd <- sd(a1.bp) }
    }
    out$method <- method
    out$Call <- Call
    out$covNames <- names(obj)[-ncol(obj)]
    out$B <- B
    class(out) <- "thlm"
    return(out)
}

cc.reg <- function(y, u, delta, z = NULL) {
    fit <- lm(y ~ u + z, subset = delta == 1)
    alpha1.est <- fit$coef["u"]
    alpha1.sd <- coef(summary(fit))[2,2]
    alpha2.est <- fit$coef[-(1:2)]
    alpha2.sd <- coef(summary(fit))[-(1:2),2]
    power <- (((alpha1.est + alpha1.sd * qnorm(0.975))<0) ||
              (0<(alpha1.est - alpha1.sd * qnorm(0.975))))*1
    ## output results
    list(a1 = alpha1.est, a1.sd = alpha1.sd,
         a2 = alpha2.est, a2.sd = alpha2.sd,
         power = power)
}

rev.surv.reg <- function(y, u, delta, z = NULL) {
    if (is.null(z)) fit <- coxph(Surv(u, delta) ~ y)
    else fit <- coxph(Surv(u, delta) ~ y + z)
    alpha1.est <- fit$coef["y"]
    alpha1.sd <- coef(summary(fit))[1,3]
    alpha1.pval <- coef(summary(fit))[1,5]
    power <- (alpha1.pval < 0.05)*1
    list(a1 = alpha1.est, a1.sd = alpha1.sd, power = power)
}

thlm.control <- function(t0.interval = NULL, t0.plot = TRUE) {
    list(t0.interval = t0.interval, t0.plot = t0.plot)
}

opt.threshold.m1 <- function(y, u, delta, x.upplim = 1.75, t0.interval = NULL, t0.plot = FALSE) {
    n <- length(y)
    if (is.null(t0.interval)) t0.interval <- range(u)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    obj.m1 <- function(t0){
        ind <- sum(km.time <= t0)
        surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
            (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
        ## approx(km.time, km.est, t0)$y
        km.time.2 <- c(0, km.time[-length(km.time)])
        km.est.2 <- c(1, km.est[-length(km.est)]) 
        vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
            (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
        est.cond.x <- sum(vec.a) / surv.t0 + t0
        bias.correct <- est.cond.x - sum(u * delta * (u <= t0)) / sum(delta * (u <= t0))
        n1 <- sum((u <= t0 & delta == 1))
        n2 <- sum(u > t0) 
        return(abs(bias.correct) / (sqrt(1 / n1 + 1 / n2))) ## assuming equal variance    
    }
    t0.opt <- optimize(obj.m1, t0.interval, tol = 1e-5, maximum = TRUE)
    if (t0.plot) {
        t0.vec <- seq(t0.interval[1], t0.interval[2], length = 50)[2:49]
        t0.thr <- unlist(sapply(t0.vec, obj.m1))
        t0.vec <- c(t0.vec, t0.opt$maximum)
        t0.thr <- c(t0.thr, t0.opt$objective)
        t0.thr <- t0.thr[order(t0.vec)]
        t0.vec <- t0.vec[order(t0.vec)]
        plot(t0.vec, t0.thr, "l", lty = 2, lwd = 2, main = "", xlab = "", ylab = "")
        mtext(expression(bold("Threshold estimation with deletion threshold regression")), 3, line = .5, cex = 1.2)
        title(xlab = "Time", ylab = "Objective function", line = 2)
        fit <- loess.smooth(t0.vec, t0.thr, degree = 2)
        lines(fit$x, fit$y, lty = 1, lwd = 2)
        abline(v = t0.opt$maximum, lty = "dotted", lwd = 1.5)
        legend("topright", c("raw value", "smoothed"), lty = 1:2, bty = "n")
    }
    list(t0.opt = t0.opt$maximum, obj.val = t0.opt$objective)
}

threshold.reg.m1 <- function(y, u, delta, z = NULL, x.upplim = 1.75, t0 = NULL, control = thlm.control()) {
    if (is.null(x.upplim)) x.upplim <- max(u)
    if (is.null(t0)) t0 <- opt.threshold.m1(y, u, delta, x.upplim, control$t0.interval, control$t0.plot)$t0.opt
    n <- length(y)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    ind <- sum(km.time <= t0)
    surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
        (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
    km.time.2 <- c(0, km.time[-length(km.time)])
    km.est.2 <- c(1, km.est[-length(km.est)]) 
    vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
        (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
    est.cond.x <- sum(vec.a) / surv.t0 + t0
    delta2 <- 1 * (u > t0)
    fit.lm <- lm(y ~ delta2 + z, subset = ((u <= t0 & delta == 1) | (u > t0)))
    beta1.est <- coef(fit.lm)[2]
    beta1.sd <- coef(summary(fit.lm))[2,2]
    power <- as.numeric(1 * (abs(beta1.est) / beta1.sd > qnorm(.975)))
    beta2.est <- coef(fit.lm)[-(1:2)]
    beta2.sd <- coef(summary(fit.lm))[-(1:2), 2]
    bias.correct <- est.cond.x - sum(u * delta * (u <= t0)) / sum(delta * (u <= t0))
    alpha1.est <- beta1.est / bias.correct
    ## percentage of different portions
    percent.del <- sum((1 - delta) * (u <= t0))/n
    percent.below <- sum(delta * (u <= t0))/n
    percent.above <- sum(u > t0) / n
    list(a1 = alpha1.est, b1 = beta1.est, b1.sd = beta1.sd, power.b1 = power,
         a2 = beta2.est, a2.sd = beta2.sd, cond.x = est.cond.x, pdel = percent.del,
         pbelow = percent.below, pabove = percent.above, bias.cor = bias.correct, threshold = t0)
}

opt.threshold.m2 <- function(y, u, delta, x.upplim = 1.75, t0.interval = NULL, t0.plot = FALSE) {
    n <- length(y)
    if (is.null(t0.interval)) t0.interval <- range(u)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    obj.m2 <- function(t0) {
        ind <- sum(km.time <= t0)
        surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
            (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
        km.time.2 <- c(0, km.time[-length(km.time)])
        km.est.2 <- c(1, km.est[-length(km.est)])
        vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
            (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
        est.cond.x <- sum(vec.a) / surv.t0 + t0 
        vec.b <- (km.time - km.time.2) * (km.est.2 + km.est) / 2
        est.mean.x <- sum(vec.b)
        est.prob.u <- approx(u, ecdf(u)(u), t0)$y
        bias.correct <- (est.cond.x - est.mean.x) / est.prob.u
        n1 <- sum(u <= t0)
        n2 <- sum(u > t0)
        return(abs(bias.correct) / (sqrt(1/n1 + 1/n2)))
    }
    t0.opt <- optimize(obj.m2, t0.interval, tol = 1e-5, maximum = TRUE)
    if (t0.plot) {
        t0.vec <- seq(t0.interval[1], t0.interval[2], length = 50)[2:49]
        t0.thr <- unlist(sapply(t0.vec, obj.m2))
        t0.vec <- c(t0.vec, t0.opt$maximum)
        t0.thr <- c(t0.thr, t0.opt$objective)
        t0.thr <- t0.thr[order(t0.vec)]
        t0.vec <- t0.vec[order(t0.vec)]        
        plot(t0.vec, t0.thr, "l", lty = 2, lwd = 2, main = "", xlab = "", ylab = "")
        mtext(expression(bold("Threshold estimation with complete threshold regression")), 3, line = .5, cex = 1.2)
        title(xlab = "Time", ylab = "Objective function", line = 2)
        fit <- loess.smooth(t0.vec, t0.thr, degree = 2)
        lines(fit$x, fit$y, lty = 1, lwd = 2)
        abline(v = t0.opt$maximum, lty = "dotted", lwd = 1.5)
        legend("topright", c("raw value", "smoothed"), lty = 1:2, bty = "n")
    }
    list(t0.opt = t0.opt$maximum, obj.val = t0.opt$objective)
}


threshold.reg.m2 <- function(y, u, delta, z = NULL, x.upplim = 1.75, t0 = NULL, control = thlm.control()) {
    if (is.null(x.upplim)) x.upplim <- max(u)
    if (is.null(t0)) t0 <- opt.threshold.m2(y, u, delta, x.upplim, control$t0.interval, control$t0.plot)$t0.opt
    n <- length(y)
    km.fit <- survfit(Surv(u, delta) ~ 1)
    if (delta[which.max(u)]) {
        km.time <- summary(km.fit)$time
        km.est <- summary(km.fit)$surv
    } else {
        km.time <- c(summary(km.fit)$time, max(u))
        km.est <- c(summary(km.fit)$surv, 0)
    }
    km.time[which.max(km.time)] <- max(km.time, x.upplim)
    ind <- sum(km.time <= t0)
    surv.t0 <- km.est[ind+1] + (km.est[ind] - km.est[ind + 1]) *
        (km.time[ind + 1] - t0) / (km.time[ind + 1] - km.time[ind])
    km.time.2 <- c(0, km.time[-length(km.time)])
    km.est.2 <- c(1, km.est[-length(km.est)]) 
    vec.a <- (km.time > t0) * (km.time - pmax(km.time.2, t0)) *
        (pmin(km.est.2, surv.t0) + pmin(km.est, surv.t0)) / 2
    est.cond.x <- sum(vec.a) / surv.t0 + t0
    vec.b <- (km.time - km.time.2) * (km.est.2 + km.est) / 2
    est.mean.x <- sum(vec.b)
    est.prob.u <- approx(u, ecdf(u)(u), t0)$y
    delta2 <- 1 * (u > t0)
    fit.lm <- lm(y ~ delta2 + z)
    beta1.est <- coef(fit.lm)[2]
    beta1.sd <- coef(summary(fit.lm))[2,2]
    power <- as.numeric(1 * (abs(beta1.est) / beta1.sd > qnorm(.975)))
    beta2.est <- coef(fit.lm)[-(1:2)]
    beta2.sd <- coef(summary(fit.lm))[-(1:2), 2]
    bias.correct <- (est.cond.x - est.mean.x) / est.prob.u
    alpha1.est <- beta1.est / bias.correct
    percent.below <- sum(u <= t0) / n
    percent.above <- sum(u > t0) / n
    list(a1 = alpha1.est, b1 = beta1.est, b1.sd = beta1.sd, power.b1 = power,
         a2 = beta2.est, a2.sd = beta2.sd, cond.x = est.cond.x, mean.x = est.mean.x,
         pbelow = percent.below, pabove = percent.above, bias.cor = bias.correct, threshold = t0)
}
