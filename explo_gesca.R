library(gesca)


GSCA <- function (myModel, data, group.name = NULL, group.equal = NULL, 
          nbt = 100, itmax = 100, ceps = 1e-05, moption = 0, missingvalue = -9999) 
{
  if (!is.character(myModel) || length(myModel) > 1 || nchar(myModel) == 
      0) 
    stop("formula should be a string")
  if (!is.data.frame(data) && !is.matrix(data) || !nrow(data) || 
      !ncol(data)) 
    stop("data should be a non-empty matrix or dataframe")
  if (moption != 0) {
    if (moption != 1 && moption != 2 && moption != 3) 
      stop("moption should be 1, 2 or 3")
    cat("Missing value is defined to [", missingvalue, "]\n")
  }
  else moption = 0
  if (!is.null(group.name)) {
    orig.gname <- group.name
    if (length(group.name) > 1) 
      stop("group.name should be a column name or an index")
    if (is.na(match(group.name, colnames(data)))) {
      if ((group.name - as.integer(group.name)) != 0 || 
          ncol(data) < group.name || group.name < 0) 
        stop("group.name should be a column name or an index")
    }
    else group.name <- match(group.name, colnames(data))
    group.var <- as.character(data[, group.name])
    data = data[order(group.var), ]
  }
  else {
    group.var <- rep(1, nrow(data))
    orig.gname <- NULL
  }
  ret <- gesca:::parse.formula(myModel, m.data = data, group.equal = group.equal)
  if (!is.null(ret$W002)) 
    msg <- capture.output(res <- gsc.mg(ret$Z0, group.var, 
                                             ret$W00, ret$W002, ret$C00, ret$B00, ret$loadtype, 
                                             ret$loadtype2, ceps = ceps, nbt = nbt, itmax = itmax, 
                                             moption = moption, missingvalue = missingvalue))
  else msg <- capture.output(res <- gsc.mg(ret$Z0, group.var, 
                                            ret$W00, ret$C00, ret$B00, ret$loadtype, ceps = ceps, 
                                            nbt = nbt, itmax = itmax, moption = moption, missingvalue = missingvalue))
  res$niter <- res$niter
  res$eps <- ceps
  res$wname <- colnames(ret$Z0)
  if (!is.null(ret$W002)) 
    res$lname <- c(names(ret$loadtype), names(ret$loadtype2))
  else res$lname <- names(ret$loadtype)
  res$grp <- names(table(group.var))
  res$gname <- orig.gname
  res$B00 <- ret$B00
  res$C00 <- ret$C00
  res$W00 <- ret$W00
  res$nbt <- nbt
  if (!is.null(ret$W002)) {
    res$loadtype <- c(ret$loadtype, ret$loadtype2)
    res$W002 <- ret$W002
  }
  else res$loadtype <- ret$loadtype
  class(res) <- c("gesca", class(res))
  invisible(res)
}

#------------------------------------------------------------------------------------

gsc.mg <- function (z0, group_var, W00, C00, B00, loadtype = matrix(1, 
                                                          1, ncol(W00)), nbt = 100, 
          itmax = 100, ceps = 1e-05, moption = 0, 
          missingvalue = NULL) 
{
  if (moption == 1) {
    nrow <- nrow(z0)
    row_index <- rep(1, nrow)
    for (i in 1:nrow) {
      if (sum(which(z0[i, ] == missingvalue)) > 0) {
        row_index[i] = 0
      }
    }
    rindex <- which(row_index != 0)
    z0 <- z0[rindex, ]
    group_var <- group_var[rindex]
  }
  nobs_tot <- nrow(z0)
  nvar <- ncol(z0)
  ng <- length(unique(group_var))
  nobs_g <- matrix(, 1, ng)
  for (j in 1:ng) {
    nobs_g[, j] <- as.numeric(table(group_var))[j]
  }
  case_index <- matrix(, ng, 2)
  kk <- 0
  for (j in 1:ng) {
    k <- kk + 1
    kk <- kk + nobs_g[j]
    case_index[j, 1] = k
    case_index[j, 2] = kk
  }
  nlv <- length(loadtype)
  ntv <- nvar + nlv
  for (j in 1:nlv) {
    if (loadtype[j] == 0) {
      C00[j, ] = matrix(0, 1, nvar)
    }
  }
  A00 <- cbind(C00, B00)
  V001 <- diag(1, nvar)
  V00 <- cbind(V001, W00)
  Wi <- W00
  Ai <- A00
  windex0 <- which(W00 == 99)
  aindex0 <- which(A00 == 99)
  W0 <- matrix(0, ng * nvar, ng * nlv)
  A0 <- matrix(0, ng * nlv, ntv)
  V0 <- matrix(0, ng * nvar, ng * ntv)
  W <- W0
  A <- A0
  V <- V0
  kk <- 0
  ss <- 0
  ll <- 0
  for (j in 1:ng) {
    k <- kk + 1
    kk <- kk + nvar
    s <- ss + 1
    ss <- ss + nlv
    l <- ll + 1
    ll <- ll + ntv
    W0[k:kk, s:ss] <- W00
    Wi[windex0] <- runif(length(windex0))
    W[k:kk, s:ss] <- Wi
    A0[s:ss, ] <- A00
    Ai[aindex0] <- runif(length(aindex0))
    A[s:ss, ] <- Ai
    V0[k:kk, l:ll] <- V00
    V[k:kk, l:ll] <- cbind(V001, Wi)
  }
  I <- matrix(0, ntv * ng, ntv)
  kk <- 0
  for (g in 1:ng) {
    k <- kk + 1
    kk <- kk + ntv
    I[k:kk, ] <- diag(1, ntv)
  }
  output.constmat <- gesca:::constmat(A0)
  PHT <- output.constmat$PHT
  num_nzct <- output.constmat$num_nzct
  num_const <- output.constmat$num_const
  num_nnz_W00 <- length(W00[!W00 == 0])
  num_nnz_C00 <- length(C00[!C00 == 0])
  num_nnz_B00 <- length(B00[!B00 == 0])
  vec_FIT <- matrix(0, nbt, 1)
  vec_FIT_m <- matrix(0, nbt, 1)
  vec_FIT_s <- matrix(0, nbt, 1)
  vec_AFIT <- matrix(0, nbt, 1)
  vec_GFI <- matrix(0, nbt, 1)
  vec_SRMR <- matrix(0, nbt, 1)
  MatW <- matrix(0, nbt, num_nnz_W00 * ng)
  Matload <- matrix(0, nbt, num_nnz_C00 * ng)
  Matbeta <- matrix(0, nbt, num_nnz_B00 * ng)
  Matsmc <- matrix(0, nbt, num_nnz_C00 * ng)
  MatcorF <- matrix(0, nbt, nlv^2 * ng)
  MatTE_S <- c()
  MatID_S <- c()
  MatTE_M <- c()
  MatID_M <- c()
  for (b in 0:nbt) {
    if (b == 0) {
      if (moption > 1) {
        output.bootsample.imp <- gesca:::bootsample.imp(z0, case_index, 
                                                nvar, nobs_g, ng, b, nobs_tot, moption, missingvalue)
        Z <- output.bootsample.imp$Z
        z0_meanimp <- output.bootsample.imp$z0_meanimp
        rawz0 <- output.bootsample.imp$rawz0
      }
      else {
        output.bootsample <- gesca:::bootsample(z0, case_index, 
                                        nvar, nobs_g, ng, b, nobs_tot)
        Z <- output.bootsample$Z
      }
    }
    else {
      output.bootsample <- gesca:::bootsample(z0, case_index, nvar, 
                                      nobs_g, ng, b, nobs_tot)
      Z <- output.bootsample$Z
    }
    if (b == 0) {
      if (moption == 3) {
        output.als.mg.imp <- gesca:::als.mg.imp(Z, rawz0, W0, 
                                        A0, W, A, V, I, PHT, nvar, nlv, ng, missingvalue, 
                                        itmax, ceps)
        W <- output.als.mg.imp$W
        A <- output.als.mg.imp$A
        ZZ <- output.als.mg.imp$Z
        Psi <- output.als.mg.imp$Psi
        Gamma <- output.als.mg.imp$Gamma
        f <- output.als.mg.imp$f
        it <- output.als.mg.imp$it
        imp <- output.als.mg.imp$imp
      }
      else {
        output.als.mg <- gesca:::als.mg(Z, W0, A0, W, A, V, I, 
                                PHT, nvar, nlv, ng, itmax, ceps)
        W <- output.als.mg$W
        A <- output.als.mg$A
        Psi <- output.als.mg$Psi
        Gamma <- output.als.mg$Gamma
        f <- output.als.mg$f
        it <- output.als.mg$it
        imp <- output.als.mg$imp
      }
    }
    else {
      output.als.mg <- gesca:::als.mg(Z, W0, A0, W, A, V, I, PHT, 
                              nvar, nlv, ng, itmax, ceps)
      W <- output.als.mg$W
      A <- output.als.mg$A
      Psi <- output.als.mg$Psi
      Gamma <- output.als.mg$Gamma
      f <- output.als.mg$f
      it <- output.als.mg$it
      imp <- output.als.mg$imp
    }
    corF <- t(Gamma) %*% Gamma
    CR <- t(A[, 1:nvar])
    BR <- t(A[, (nvar + 1):ntv])
    DF <- nobs_tot * nvar
    npw <- length(which(W0 == 99))
    dpht <- diag(PHT)
    if (num_nzct == 0) {
      cnzt <- length(which(dpht == 1))
    }
    else {
      cnzt <- num_const + length(which(dpht == 1))
    }
    NPAR <- cnzt + npw
    Fit <- 1 - f/sum(diag(t(Psi) %*% Psi))
    dif_m <- Psi[, 1:nvar] - Gamma %*% t(CR)
    dif_s <- Psi[, (nvar + 1):ntv] - Gamma %*% t(BR)
    Fit_m <- 1 - sum(diag(t(dif_m) %*% dif_m))/sum(diag(t(Z) %*% 
                                                          Z))
    Fit_s <- 1 - sum(diag(t(dif_s) %*% dif_s))/sum(diag(t(Gamma) %*% 
                                                          Gamma))
    Afit <- 1 - ((1 - Fit) * (DF)/(DF - NPAR))
    output.modelfit.mg <- gesca:::modelfit.mg(Z, W, A, nvar, nlv, 
                                      ng, case_index)
    Gfi <- output.modelfit.mg$GFI
    Srmr <- output.modelfit.mg$SRMR
    COR_RES <- output.modelfit.mg$COR_RES
    total_s <- matrix(0, ng * nlv, nlv)
    indirect_s <- matrix(0, ng * nlv, nlv)
    total_m <- matrix(0, ng * nlv, nvar)
    indirect_m <- matrix(0, ng * nlv, nvar)
    k <- kk <- 0
    for (g in 1:ng) {
      k = kk + 1
      kk = kk + nlv
      output.effects <- gesca:::effects(BR[, k:kk], CR[, k:kk])
      te_s <- output.effects$te_s
      ie_s <- output.effects$ie_s
      te_m <- output.effects$te_m
      ie_m <- output.effects$ie_m
      total_s[k:kk, ] <- te_s
      indirect_s[k:kk, ] <- ie_s
      total_m[k:kk, ] <- te_m
      indirect_m[k:kk, ] <- ie_m
    }
    nz_idx <- which(BR != 0, arr.ind = TRUE)
    nz_idx_C <- which(CR != 0, arr.ind = TRUE)
    nz_idx_W <- which(W != 0, arr.ind = TRUE)
    row_names_BR <- rownames(BR)  
    col_names_BR <- colnames(BR)
    relation_labels <- paste0(col_names_BR[nz_idx[,2]], " ~ ", row_names_BR[nz_idx[,1]])
    if (b == 0) {
      if (moption == 2) {
        z0 <- z0_meanimp
      }
      else if (moption == 3) {
        z0 <- matrix(0, nobs_tot, nvar)
        kk <- 0
        for (g in 1:ng) {
          k <- kk + 1
          kk <- kk + nvar
          z0[case_index[g, 1]:case_index[g, 2], ] = ZZ[case_index[g, 
                                                                  1]:case_index[g, 2], k:kk]
        }
      }
      if (it <= itmax) {
        if (imp <= ceps) {
          message(paste("The ALS algorithm converged in", 
                        it, "iterations (convergence criterion =", 
                        ceps, ")", "\n"))
        }
        else {
          message(paste("The ALS algorithm failed to converge in", 
                        it, "iterations (convergence criterion =", 
                        ceps, ")", "\n"))
        }
      }
      WR <- W
      Cr <- CR
      Br <- BR
      samplesizes <- nobs_g
      NPAR
      FIT <- Fit
      FIT_M <- Fit_m
      FIT_S <- Fit_s
      AFIT <- Afit
      GFI <- Gfi
      SRMR <- Srmr
      R2 <- matrix(0, ng, nlv)
      AVE <- matrix(0, ng, nlv)
      Alpha <- matrix(0, ng, nlv)
      rho <- matrix(0, ng, nlv)
      Dimension <- matrix(0, ng, nlv)
      lvmean <- matrix(0, ng, nlv)
      lvvar <- matrix(0, ng, nlv)
      corr_corres <- matrix(0, ng * nvar, nvar)
      ss <- 0
      kk <- 0
      for (g in 1:ng) {
        s <- ss + 1
        ss <- ss + nlv
        k <- kk + 1
        kk <- kk + nvar
        if (moption == 3) {
          z0_g <- z0_meanimp[case_index[g, 1]:case_index[g, 
                                                         2], ]
        }
        else {
          z0_g <- z0[case_index[g, 1]:case_index[g, 2], 
          ]
        }
        W_g <- W[k:kk, s:ss]
        CF_g <- corF[s:ss, s:ss]
        B <- t(BR[, s:ss])
        stdL <- CR[, s:ss]
        for (j in 1:nlv) {
          R2[g, j] <- t(B[, j, drop = FALSE]) %*% CF_g[, 
                                                       j, drop = FALSE]
          zind <- which(W00[, j] != 0)
          nnzload <- length(zind)
          if (nnzload > 0) {
            sumload <- sum(stdL[zind, j]^2)
            sumload_rho1 <- sum(stdL[zind, j])^2
            sumload_rho2 <- sum(1 - stdL[zind, j]^2)
            AVE[g, j] <- sumload/nnzload
            rho[g, j] <- sumload_rho1/(sumload_rho1 + 
                                         sumload_rho2)
          }
          nzj = length(zind)
          if (nzj > 1) {
            zsubset <- z0_g[, zind]
            Alpha[g, j] <- gesca:::cronbach.alpha(zsubset)
            eigval <- svd(cor(zsubset))$d
            kr <- length(which(eigval > 1))
            Dimension[g, j] <- kr
          }
          else {
            Alpha[g, j] <- 1
          }
        }
        lvscore_g <- gesca:::lvscore(z0_g, W_g)
        lvmean[g, ] <- apply(lvscore_g, 2, mean)
        lvvar[g, ] <- apply(lvscore_g, 2, var)
        sample_corr <- cor(z0_g)
        corr_corres[k:kk, ][upper.tri(corr_corres[k:kk, 
        ], diag = FALSE)] <- COR_RES[k:kk, ][upper.tri(COR_RES[k:kk, 
        ], diag = FALSE)]
        corr_corres[k:kk, ][lower.tri(corr_corres[k:kk, 
        ], diag = FALSE)] <- sample_corr[lower.tri(sample_corr, 
                                                   diag = FALSE)]
      }
      R2
      AVE
      Alpha
      rho
      LV_MEAN <- lvmean
      LV_VAR <- lvvar
      corr_corres
      Dimension
      mW <- data.frame(
        weight  = W[nz_idx_C],
        stringsAsFactors = FALSE
      )
      mC <- data.frame(
        Loading  = CR[nz_idx_C],
        stringsAsFactors = FALSE
      )
      mB <- data.frame(
        Estimate  = BR[nz_idx],
        stringsAsFactors = FALSE
      )
      mSMC <- mC^2
      mCF <- as.matrix(corF[which(!corF == 0)])
      latentcorr <- corF
      TE_S <- total_s
      ID_S <- indirect_s
      TE_M <- total_m
      ID_M <- indirect_m
    }
    else {
      vecw <- as.matrix(W[which(!W == 0)])
      vecload <- as.matrix(CR[which(!CR == 0)])
      vecbeta <- as.matrix(BR[which(!BR == 0)])
      veccorF <- as.matrix(corF[which(!corF == 0)])
      vec_FIT[b] <- Fit
      vec_FIT_m[b] <- Fit_m
      vec_FIT_s[b] <- Fit_s
      vec_AFIT[b] <- Afit
      vec_GFI[b] <- Gfi
      vec_SRMR[b] <- Srmr
      MatW[b, ] <- t(vecw)
      Matload[b, ] <- t(vecload)
      Matbeta[b, ] <- t(vecbeta)
      Matsmc[b, ] <- t(vecload^2)
      MatcorF[b, ] <- t(veccorF)
      MatTE_S <- rbind(MatTE_S, total_s[which(!total_s == 
                                                0)])
      MatID_S <- rbind(MatID_S, indirect_s[which(!indirect_s == 
                                                   0)])
      MatTE_M <- rbind(MatTE_M, total_m[which(!total_m == 
                                                0)])
      MatID_M <- rbind(MatID_M, indirect_m[which(!indirect_m == 
                                                   0)])
    }
  }
  if (nbt > 0) {
    se_beta <- apply(Matbeta, 2, sd)
    mB$Std_Error    <- se_beta
    lb <- ceiling(nbt * 0.025)
    ub <- ceiling(nbt * 0.975)
    sortFIT <- sort(vec_FIT)
    sortFIT_m <- sort(vec_FIT_m)
    sortFIT_s <- sort(vec_FIT_s)
    sortAFIT <- sort(vec_AFIT)
    sortGFI <- sort(vec_GFI)
    sortSRMR <- sort(vec_SRMR)
    sortw <- apply(MatW, 2, sort)
    se_w <- apply(MatW, 2, sd)
    mW$Std_Error  <- se_w
    mW$CI_95_LB <- sortw[lb, ]
    mW$CI_95_UB <- sortw[ub, ]
    sortload <- apply(Matload, 2, sort)
    se_load <- apply(Matload, 2, sd)
    mC$Std_Error  <- se_load
    mC$CI_95_LB <- sortload[lb, ]
    mC$CI_95_UB <- sortload[ub, ]
    sortbeta <- apply(Matbeta, 2, sort)
    mB$CI_95_LB <- sortbeta[lb, ]
    mB$CI_95_UB <- sortbeta[ub, ]
    sortsmc <- apply(Matsmc, 2, sort)
    sortcorF <- apply(MatcorF, 2, sort)
    sortte_s <- apply(MatTE_S, 2, sort)
    sortid_s <- apply(MatID_S, 2, sort)
    sortte_m <- apply(MatTE_M, 2, sort)
    sortid_m <- apply(MatID_M, 2, sort)
    output.gsca.mg <- list(WR = WR, CR = Cr, BR = Br, samplesizes = samplesizes, 
                           NPAR = NPAR, FIT = FIT, FIT_M = FIT_M, FIT_S = FIT_S, 
                           AFIT = AFIT, GFI = GFI, SRMR = SRMR, R2 = R2, AVE = AVE, 
                           Alpha = Alpha, rho = rho, LV_MEAN = LV_MEAN, LV_VAR = LV_VAR, 
                           corr_corres = corr_corres, Dimension = Dimension, 
                           latentcorr = latentcorr, TE_S = TE_S, ID_S = ID_S, 
                           TE_M = TE_M, ID_M = ID_M, mW = mW, mC = mC, mB = mB, 
                           mSMC = mSMC, mCF = mCF, lb = lb, ub = ub, vec_FIT = vec_FIT, 
                           vec_FIT_m = vec_FIT_m, vec_FIT_s = vec_FIT_s, vec_AFIT = vec_AFIT, 
                           vec_GFI = vec_GFI, vec_SRMR = vec_SRMR, MatW = MatW, 
                           Matload = Matload, Matbeta = Matbeta, Matsmc = Matsmc, 
                           MatcorF = MatcorF, MatTE_S = MatTE_S, MatID_S = MatID_S, 
                           MatTE_M = MatTE_M, MatID_M = MatID_M, sortFIT = sortFIT, 
                           sortFIT_m = sortFIT_m, sortFIT_s = sortFIT_s, sortAFIT = sortAFIT, 
                           sortGFI = sortGFI, sortSRMR = sortSRMR, sortw = sortw, 
                           sortload = sortload, sortbeta = sortbeta, sortsmc = sortsmc, 
                           sortcorF = sortcorF, sortte_s = sortte_s, sortid_s = sortid_s, 
                           sortte_m = sortte_m, sortid_m = sortid_m)
    output.gsca.mg
  }
  else {
    output.gsca.mg <- list(WR = WR, CR = Cr, BR = Br, samplesizes = samplesizes, 
                           NPAR = NPAR, FIT = FIT, FIT_M = FIT_M, FIT_S = FIT_S, 
                           AFIT = AFIT, GFI = GFI, SRMR = SRMR, R2 = R2, AVE = AVE, 
                           Alpha = Alpha, rho = rho, LV_MEAN = LV_MEAN, LV_VAR = LV_VAR, 
                           corr_corres = corr_corres, Dimension = Dimension, 
                           latentcorr = latentcorr, TE_S = TE_S, ID_S = ID_S, 
                           TE_M = TE_M, ID_M = ID_M, mW = mW, mC = mC, mB = mB, 
                           mSMC = mSMC, mCF = mCF)
    output.gsca.mg
  }
}


#-------------------------------------------------------------------------------

featmeasures2 <- function(object, print = TRUE) {
  
  lb <- object$lb
  ub <- object$ub
  vars <- c("FIT", "AFIT", "GFI", "SRMR", "FIT_M", "FIT_S")
  
  if (!is.null(object$vec_FIT)) {
    res <- data.frame(
      Measure   = vars,
      Estimate  = round(c(object$FIT, object$AFIT, object$GFI, 
                          object$SRMR, object$FIT_M, object$FIT_S), 4),
      Std.Error = round(c(
        sd(object$vec_FIT),
        sd(object$vec_AFIT),
        sd(object$vec_GFI),
        sd(object$vec_SRMR),
        sd(object$vec_FIT_m),
        sd(object$vec_FIT_s)
      ), 4),
      CI_LB = round(c(
        object$sortFIT[lb],
        object$sortAFIT[lb],
        object$sortGFI[lb],
        object$sortSRMR[lb],
        object$sortFIT_m[lb],
        object$sortFIT_s[lb]
      ), 4),
      CI_UB = round(c(
        object$sortFIT[ub],
        object$sortAFIT[ub],
        object$sortGFI[ub],
        object$sortSRMR[ub],
        object$sortFIT_m[ub],
        object$sortFIT_s[ub]
      ), 4),
      row.names = NULL
    )
  } else {
    res <- data.frame(
      Measure  = vars,
      Estimate = round(c(object$FIT, object$AFIT, object$GFI,
                         object$SRMR, object$FIT_M, object$FIT_S), 4),
      row.names = NULL
    )
  }
  
  if (print) {
    cat("  Number of parameters       ", object$NPAR, "\n")
    cat("  Number of bootstrap samples", object$nbt, "\n\n")
    print(res, row.names = FALSE)
  }
  
  invisible(res)
}

#------------------------------------------------------------------------------

summary_gesca2 <- function(object) {
  
  n.group <- ncol(object$samplesizes)
  lb <- object$lb
  ub <- object$ub
  gname <- object$gname
  grp  <- object$grp
  wname <- object$wname
  
  # ── helper : construit un data.frame à partir d'une matrice de résultats ──
  build_df <- function(res, row_names, group_idx, group_label) {
    has_se <- ncol(res) == 4
    df <- data.frame(
      group    = group_label,
      term     = row_names,
      estimate = res[, 1],
      stringsAsFactors = FALSE
    )
    if (has_se) {
      df$std.error <- res[, 2]
      df$ci_lb     <- res[, 3]
      df$ci_ub     <- res[, 4]
    }
    rownames(df) <- NULL
    df
  }
  
  # ── 1. WEIGHTS ────────────────────────────────────────────────────────────
  sub.wname <- wname[which(object$W00 != 0, arr.ind = TRUE)[, 1]]
  
  if (is.null(object$mW)) {
    # first-order weights
    if (is.null(object$MatW1))
      res_w <- round(cbind(object$mW1), 4)
    else
      res_w <- round(cbind(object$mW1,
                           apply(t(object$MatW1), 1, sd),
                           t(object$sortw1)[, lb],
                           t(object$sortw1)[, ub]), 4)
    
    pl <- 0
    w_list <- vector("list", n.group)
    for (j in 1:n.group) {
      idx  <- (pl + 1):(pl + length(sub.wname))
      lbl  <- if (n.group > 1) paste0(gname, "=", grp[j]) else "all"
      w_list[[j]] <- build_df(res_w[idx, , drop = FALSE], sub.wname, j, lbl)
      pl <- pl + length(sub.wname)
    }
    weights_df <- do.call(rbind, w_list)
    
    # second-order weights (si présents)
    weights2_df <- NULL
    if (!is.null(object$mW2)) {
      lname.a <- na.omit(object$lname[rowSums(object$W002) != 0][1:nrow(object$W002)])
      if (is.null(object$MatW2))
        res_w2 <- round(cbind(object$mW2), 4)
      else
        res_w2 <- round(cbind(object$mW2,
                              apply(t(object$MatW2), 1, sd),
                              t(object$sortw2)[, lb],
                              t(object$sortw2)[, ub]), 4)
      pl <- 0
      w2_list <- vector("list", n.group)
      for (j in 1:n.group) {
        idx  <- (pl + 1):(pl + length(lname.a))
        lbl  <- if (n.group > 1) paste0(gname, "=", grp[j]) else "all"
        w2_list[[j]] <- build_df(res_w2[idx, , drop = FALSE], lname.a, j, lbl)
        pl <- pl + length(lname.a)
      }
      weights2_df <- do.call(rbind, w2_list)
    }
    
  } else {
    # poids uniques (cas sans hiérarchie)
    if (is.null(object$MatW))
      res_w <- round(cbind(object$mW), 4)
    else
      res_w <- round(cbind(object$mW,
                           apply(t(object$MatW), 1, sd),
                           t(object$sortw)[, lb],
                           t(object$sortw)[, ub]), 4)
    pl <- 0
    w_list <- vector("list", n.group)
    for (j in 1:n.group) {
      idx  <- (pl + 1):(pl + length(sub.wname))
      lbl  <- if (n.group > 1) paste0(gname, "=", grp[j]) else "all"
      w_list[[j]] <- build_df(res_w[idx, , drop = FALSE], sub.wname, j, lbl)
      pl <- pl + length(sub.wname)
    }
    weights_df  <- do.call(rbind, w_list)
    weights2_df <- NULL
  }
  
  # ── 2. LOADINGS ───────────────────────────────────────────────────────────
  sub.wname_load <- wname[which(t(object$C00) != 0, arr.ind = TRUE)[, 1]]
  
  if (length(sub.wname_load) > 0) {
    if (is.null(object$Matload))
      res_l <- round(cbind(object$mC), 4)
    else
      res_l <- round(cbind(object$mC,
                           apply(t(object$Matload), 1, sd),
                           t(object$sortload)[, lb],
                           t(object$sortload)[, ub]), 4)
    pl <- 0
    l_list <- vector("list", n.group)
    for (j in 1:n.group) {
      idx  <- (pl + 1):(pl + length(sub.wname_load))
      lbl  <- if (n.group > 1) paste0(gname, "=", grp[j]) else "all"
      l_list[[j]] <- build_df(res_l[idx, , drop = FALSE], sub.wname_load, j, lbl)
      pl <- pl + length(sub.wname_load)
    }
    loadings_df <- do.call(rbind, l_list)
  } else {
    loadings_df <- NULL
  }
  
  # ── 3. PATH COEFFICIENTS ──────────────────────────────────────────────────
  idx_b <- which(t(object$B00) != 0, arr.ind = TRUE)
  
  if (nrow(idx_b)) {
    idx_b <- cbind(idx_b[, 2], idx_b[, 1])
    path_names <- paste0(object$lname[idx_b[, 2]], "~",
                         object$lname[idx_b[, 1]])
    
    if (is.null(object$Matbeta))
      res_b <- round(cbind(object$mB), 4)
    else
      res_b <- round(cbind(object$mB,
                           apply(t(object$Matbeta), 1, sd),
                           t(object$sortbeta)[, lb],
                           t(object$sortbeta)[, ub]), 4)
    pl <- 0
    b_list <- vector("list", n.group)
    for (j in 1:n.group) {
      idx_r <- (pl + 1):(pl + nrow(idx_b))
      lbl   <- if (n.group > 1) paste0(gname, "=", grp[j]) else "all"
      b_list[[j]] <- build_df(res_b[idx_r, , drop = FALSE], path_names, j, lbl)
      pl <- pl + nrow(idx_b)
    }
    paths_df <- do.call(rbind, b_list)
  } else {
    paths_df <- NULL
  }
  
  # ── sortie ────────────────────────────────────────────────────────────────
  out <- list(
    weights   = weights_df,
    weights2  = weights2_df,   # NULL si pas de 2nd order
    loadings  = loadings_df,
    paths     = paths_df
  )
  return(out)
}
