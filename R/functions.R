# some functions used in the algorithm
#### MCP and SCAD ####
sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))
  thval <- 1 - th/xn
 
  if(xn ==0)
  {
    value <- xn
  }else{
    value <- thval*((thval) >0)*x
  }
}
mcp <- function(x,lam,nu,gam)
{
  temp <- gam*lam
  xn <- sqrt(sum(x^2))
  if(xn <= temp)
  {
    z <- sfun(x,lam/nu) / (1 - 1/(gam*nu))
  }else{
    z <- x
  }
  return(z)
}

scad <- function(x,lam,nu,gam)
{
  temp1 <- lam/nu
  temp2 <- gam * lam
  xn <- sqrt(sum(x^2))
  
  if(xn <= lam + temp1)
  {
    z <- sfun(x, temp1)
  }else if(xn <= temp2 & xn >= lam + temp1)
  {
    z <- sfun(x, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu))
  }else{
    z <- x
  }
  
  return(z)
  
}


getgroup = function(deltam, n, tol = 1e-2)
{
  p = nrow(deltam)
  b2value =sqrt(colMeans(deltam^2))
  b2value[b2value <= tol] = 0
  
  d2 = matrix(0, n, n)
  for(j in 1:(n-1))
  {
    indexj1 = (2*n -j)*(j-1)/2 + 1
    indexj2 = indexj1 + n - j - 1
    d2[(n - indexj2 + indexj1):n,j] = b2value[indexj1:indexj2]
  }
  d2 = t(d2) + d2
  
  
  ngj = 1:n
  groupest = rep(0,n)
  j = 0
  
  while (length(ngj) >0)
  {
    j = j + 1
    gj = (1:n)[d2[ngj[1],] ==0]
    indexj = ngj %in% gj
    gj = ngj[indexj]
    ngj = ngj[!indexj]
    groupest[gj] = j * rep(1,length(gj))
  }
  
  
  return(groupest)
  
}





### vectorized mcp applied to rows of a matrix
mcp_mat <- function(psim, lam, nu, gam)
{
  xn <- sqrt(rowSums(psim^2))
  temp <- gam * lam
  denom <- 1 - 1/(gam * nu)
  # default: identity (xn > temp)
  result <- psim
  mask <- xn <= temp
  if(any(mask)) {
    xn_m <- xn[mask]
    # sfun scale factor: max(0, 1 - (lam/nu)/xn) / denom; 0 when xn==0
    scale_m <- ifelse(xn_m == 0, 0,
                      pmax(0, 1 - (lam/nu) / xn_m) / denom)
    result[mask, ] <- psim[mask, , drop = FALSE] * scale_m
  }
  return(result)
}


### quick inverse — returns full (n*p)×(n*p) matrix (kept for compatibility)

inverseR <- function(indexy, x, nu = 1)
{
  uindex <- unique(indexy)
  n <- length(uindex)
  p <- ncol(x)
  Ip <- diag(1, p, p)

  DB <- matrix(0, p,   p)
  A0B <- matrix(0, n*p, p)
  matinv <- matrix(0, n*p, n*p)

  for(i in 1:n)
  {
    xmati <- x[indexy == uindex[i], , drop = FALSE]
    matj <- chol2inv(chol(crossprod(xmati) + nu*n*Ip))
    DB <- DB + matj
    indexi <- ((i-1)*p+1):(i*p)
    A0B[indexi, ]          <- matj
    matinv[indexi, indexi] <- matj
  }

  temp1 <- chol2inv(chol(1/nu*Ip - DB))
  A0tmp <- A0B %*% temp1
  matinv <- matinv + tcrossprod(A0tmp, A0B)
  return(matinv)
}

### fast version: directly computes M^{-1} %*% v without forming the full matrix.

inverseR_apply <- function(group_idx, x, v, nu = 1)
{
  n  <- length(group_idx)
  p  <- ncol(x)
  Ip <- diag(1, p, p)

  DB <- matrix(0, p, p)
  matj_list <- vector("list", n)

  for(i in 1:n)
  {
    xmati  <- x[group_idx[[i]], , drop = FALSE]
    matj <- chol2inv(chol(crossprod(xmati) + nu*n*Ip))
    DB <- DB + matj
    matj_list[[i]] <- matj
  }

  temp1 <- chol2inv(chol(1/nu*Ip - DB))

  # t(A0B) %*% v = sum_i matj_i %*% v_i  (matj_i is symmetric)
  sumterm <- numeric(p)
  for(i in 1:n)
  {
    indexi <- ((i-1)*p+1):(i*p)
    sumterm <- sumterm + drop(matj_list[[i]] %*% v[indexi])
  }

  corr <- drop(temp1 %*% sumterm)

  result <- numeric(n * p)
  for(i in 1:n)
  {
    indexi <- ((i-1)*p+1):(i*p)
    result[indexi] <- drop(matj_list[[i]] %*% (v[indexi] + corr))
  }

  return(result)
}



