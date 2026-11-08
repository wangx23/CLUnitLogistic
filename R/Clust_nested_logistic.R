## not standardized ###
### nested error model clustered regression coefficient model

Clust_nested_logistic <- function(indexy, y, x, index_d, lam,
                                  initial0, vv0 = rep(0,length(unique(indexy))),
                                  var_v0 = max(0.001, var(vv0)), 
                                  nu = 1, gam = 3, lam0 = 0.01,
                                  maxiter= 500, inner_loop= 1,
                                  tol = 1e-3)
{
  n0 <- length(y)
  uindex <- unique(indexy)
  nr <- length(uindex)
  ns <- as.numeric(table(indexy)) # number of replicates
  nbar <- mean(ns)

  ## integer index map and precomputed per-group row indices
  idx_map   <- match(indexy, uindex)
  group_idx <- lapply(1:nr, function(i) which(idx_map == i))

  ## sparse function
  D <- matrix(0,nr*(nr-1)/2,nr)
  for(j in 1:(nr-1))
  {
    indexj <- (2*nr-j)*(j-1)/2
    indexvj <- indexj + (1:(nr-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nr)] <- -1
  }

  len_d <- length(index_d)

  x_d <- x[,index_d, drop = FALSE] ## subgroup part

  #### initial ###
  betam0   <- initial0$betam
  deltam   <- D %*% betam0
  um       <- matrix(0, nr*(nr-1)/2, len_d)
  beta_cur <- c(t(betam0))

  vv_cur <- vv0 ## random effects
  #var_v  <- max(0.0001, var(vv_cur))
  var_v <- var_v0

  ## Xm %*% b using block structure 
  xm_times_b <- function(b) {
    bm <- matrix(b, nr, len_d, byrow = TRUE)
    rowSums(x_d * bm[idx_map, , drop = FALSE])
  }
  ## t(Xm) %*% v 
  txm_times_v <- function(v) {
    c(t(rowsum(as.numeric(v) * x_d, idx_map)))
  }

  if(len_d == ncol(x))
  {
    ### initial
    lin <- xm_times_b(beta_cur) + vv_cur[idx_map]
    mu  <- 1/(1+exp(-lin))
    V   <- mu * (1-mu)
    #####

    for(mm in 1:maxiter)
    {
      ### inner loop
      for(r in 1:inner_loop)
      {
        Vv <- as.numeric(V)
        lin_xm <- xm_times_b(beta_cur)

        ## update beta
        temp1  <- txm_times_v(y - mu) / nbar
        temp2  <- c(t(rowsum(Vv * lin_xm * x_d, idx_map))) / nbar +
                  nu * c(t(crossprod(D, deltam - um/nu)))
        beta_r <- inverseR_apply(group_idx, x_d * sqrt(Vv)/sqrt(nbar), temp1 + temp2, nu)

        lin_xm <- xm_times_b(beta_r)
        lin <- lin_xm + vv_cur[idx_map]
        mu  <- 1/(1+exp(-lin))
        V <- mu * (1-mu)
        Vv<- as.numeric(V)

        ### update random effects
        tempv0 <- drop(rowsum(Vv, idx_map))
        tempv1 <- tempv0 + 1/var_v
        tempv2 <- drop(rowsum(as.numeric(y - mu), idx_map)) + tempv0 * vv_cur
        vv_r   <- tempv2 / tempv1

        diffinner <- sqrt(sum((beta_r - beta_cur)^2) + sum((vv_r - vv_cur)^2))

        ## updates
        beta_cur <- beta_r
        vv_cur   <- vv_r

        ## lin_xm 
        lin <- lin_xm + vv_cur[idx_map]
        mu <- 1/(1+exp(-lin))
        V <- mu * (1-mu)
        Vv <- as.numeric(V)

        ## update var_v
        var_v <- mean(1/(drop(rowsum(Vv, idx_map)) + 1/var_v) + vv_cur^2)

        if(diffinner < tol*sqrt(nr*len_d))
        {break}
      }

      betam  <- matrix(beta_r, nr, len_d, byrow = TRUE)

      ### update deltam
      Dbetam <- D %*% betam
      psim <- Dbetam + um/nu
      deltam_new <- mcp_mat(psim, lam, nu, gam)
      um <- um + nu * (Dbetam - deltam_new)

      sm <- sqrt(sum((deltam_new - deltam)^2))
      rm <- sqrt(sum((Dbetam - deltam_new)^2))

      deltam <- deltam_new

      if(rm <= tol*sqrt(nr*len_d))
      {break}
    }

    lin <- xm_times_b(beta_cur) + vv_cur[idx_map]
    muest <- 1/(1+exp(-lin))
    nllvalue  <- -sum((y * lin - log(1+exp(lin)))) + 0.5*sum(vv_cur^2)/var_v
    nllvalue1 <- -sum((y * lin - log(1+exp(lin))))

    cluster_est <- getgroup(t(deltam_new), nr)

    betac <- aggregate(betam, by = list(cluster_est), FUN = mean)
    output <- list(beta = beta_cur, betam = betam, betac = betac,
                   rand = vv_cur, var_v = var_v,
                   deltam = deltam_new, cluster = cluster_est,
                   lin = lin, mu = muest, nll = nllvalue, nll1 = nllvalue1,
                   sm = sm, rm = rm, niter = mm, betam0 = betam0)
  }

  if(len_d < ncol(x))
  {
    eta_cur <- initial0$eta
    ncz <- ncol(x) - len_d
    z   <- x[,-(index_d), drop = FALSE] ## common part

    ### initial
    lin <- drop(z %*% eta_cur) + xm_times_b(beta_cur) + vv_cur[idx_map]
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    #####

    for(mm in 1:maxiter)
    {
      ### inner loop
      for(r in 1:inner_loop)
      {
        Vv <- as.numeric(V)

        # update eta
        eta_r <- eta_cur + solve(crossprod(sqrt(Vv) * z)) %*% crossprod(z, y - mu)

        ## lin_xm for current beta_cur (reused in beta update below)
        lin_xm <- xm_times_b(beta_cur)
        lin <- drop(z %*% eta_r) + lin_xm + vv_cur[idx_map]
        mu  <- 1/(1+exp(-lin))
        V <- mu * (1-mu)
        Vv <- as.numeric(V)

    
        temp1  <- txm_times_v(y - mu) / nbar
        temp2  <- c(t(rowsum(Vv * lin_xm * x_d, idx_map))) / nbar +
                  nu * c(t(crossprod(D, deltam - um/nu)))
        beta_r <- inverseR_apply(group_idx, x_d * sqrt(Vv)/sqrt(nbar), temp1 + temp2, nu)

        lin_xm <- xm_times_b(beta_r)
        lin <- drop(z %*% eta_r) + lin_xm + vv_cur[idx_map]
        mu <- 1/(1+exp(-lin))
        V <- mu * (1-mu)
        Vv <- as.numeric(V)

        ### update random effects
        tempv0 <- drop(rowsum(Vv, idx_map))
        tempv1 <- tempv0 + 1/var_v
        tempv2 <- drop(rowsum(as.numeric(y - mu), idx_map)) + tempv0 * vv_cur
        vv_r   <- tempv2 / tempv1

        diffinner <- sqrt(sum((beta_r - beta_cur)^2) +
                          sum((vv_r   - vv_cur)^2)   +
                          sum((eta_r  - eta_cur)^2))

        ## updates
        eta_cur  <- eta_r
        beta_cur <- beta_r
        vv_cur   <- vv_r

        ## lin_xm 
        lin <- drop(z %*% eta_cur) + lin_xm + vv_cur[idx_map]
        mu <- 1/(1+exp(-lin))
        V <- mu * (1-mu)
        Vv <- as.numeric(V)

        if(diffinner < tol*sqrt(nr*len_d))
        {break}
      }

      var_v <- mean(1/(drop(rowsum(Vv, idx_map)) + 1/var_v) + vv_cur^2)

      betam <- matrix(beta_r, nr, len_d, byrow = TRUE)

      ### update deltam
      Dbetam    <- D %*% betam
      psim      <- Dbetam + um/nu
      deltam_new <- mcp_mat(psim, lam, nu, gam)
      um        <- um + nu * (Dbetam - deltam_new)

      sm <- sqrt(sum((deltam_new - deltam)^2))
      rm <- sqrt(sum((Dbetam - deltam_new)^2))

      deltam <- deltam_new

      if(rm <= tol*sqrt(nr*len_d))
      {break}
    }

    lin    <- drop(z %*% eta_cur) + xm_times_b(beta_cur) + vv_cur[idx_map]
    muest  <- 1/(1+exp(-lin))
    nllvalue  <- -sum((y * lin - log(1+exp(lin)))) + 0.5*sum(vv_cur^2)/var_v
    nllvalue1 <- -sum((y * lin - log(1+exp(lin))))

    cluster_est <- getgroup(t(deltam_new), nr)

    betac  <- aggregate(betam, by = list(cluster_est), FUN = mean)
    output <- list(eta = eta_cur, beta = beta_cur, betam = betam, betac = betac,
                   rand = vv_cur, var_v = var_v,
                   deltam = deltam_new, cluster = cluster_est,
                   lin = lin, mu = muest, nll = nllvalue, nll1 = nllvalue1,
                   sm = sm, rm = rm, niter = mm, betam0 = betam0)
  }

  return(output)
}




### warm start bic for vv0 ####
Clust_nested_logistic_bicw <- function(indexy, y, x, index_d, lamvec,
                                      initial0 = NULL, 
                                      vv0 = rep(0, length(unique(indexy))),
                                      var_v0 = max(0, var(vv0)),
                                      nu = 1, gam = 3, lam0 = 0.01, k0 = 10,
                                      maxiter= 500, inner_loop= 1,
                                      tol = 1e-3)
{
  nlam <- length(lamvec)
  K_est <- rep(0, nlam)
  nll_value <- rep(0, nlam)
  nll1_value <- rep(0, nlam)
  out <- list()

  ### initial values ###
  if(is.null(initial0))
  {
    k0 <- max(k0, round(sqrt(length(unique(indexy)))))
    initial0 <- initialEM(indexy, y, x, index_d, k0 = k0, lam0 = lam0)
    vv0 <- initial0$vv
  }

  initial00 <- initial0

  for(jj in 1:nlam)
  {
    resj <- Clust_nested_logistic(indexy, y, x, index_d, lam = lamvec[jj],
                                  initial0 = initial0, vv0 = vv0, var_v0 = var_v0,
                                  nu = nu, gam = gam, lam0 = lam0,
                                  maxiter = maxiter,
                                  inner_loop = inner_loop,
                                  tol = tol)
    clusterj <- resj$cluster
    betam0j <- resj$betac[clusterj,-1,drop=FALSE]
    rownames(betam0j) <- NULL

    initial0$eta  <- resj$eta
    initial0$betam <- as.matrix(betam0j)
    vv0 <- resj$rand
    var_v0 <- resj$var_v

    tryCatch({
      refit00 <- fit_clust(dom = indexy, y = y, x = x, index_d = index_d,
                             cluster = clusterj)

      vv0 <- refit00$rand
      var_v0 <- as.data.frame(VarCorr(refit00))$vcov
    }, error = function(e) {
    })

    out[[jj]] <- resj
    K_est[jj] <- length(unique(clusterj))
    nll_value[jj] <- resj$nll
    nll1_value[jj] <- resj$nll1
    if(K_est[jj]==1){break}
  }


  res_df <- data.frame(lambda= lamvec[1:jj], Khat = K_est[1:jj],
                       nll = nll_value[1:jj], nll1 = nll1_value[1:jj])
  outlist <- list(out = out, df = res_df, initial00 = initial00)
  return(outlist)
}


