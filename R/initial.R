
### ignore random effects
### no adjustment of the sample size
### and then adjust
initial <- function(indexy, y, x, index_d,
                    lam0  = 0.01, k0 = 10,
                    tol = 1e-5, 
                    maxiter = 500)
{

  ### first fit a logistic regression with all are
  #in one group to get the starting values
  x <- as.matrix(x)
  x_d <- x[,index_d, drop = FALSE]
 
  len_d <- length(index_d)
  n0 <- length(indexy)
  uindex <- unique(indexy)
  nr <- length(uindex)
  ns <- as.numeric(table(indexy))
  

  Ip <- diag(1,len_d,len_d)
  D <- matrix(0,nr*(nr-1)/2,nr)
  for(j in 1:(nr-1))
  {
    indexj <- (2*nr-j)*(j-1)/2
    indexvj <- indexj + (1:(nr-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nr)] <- -1
  }
  Am <- D %x% Ip
  AtA <- t(D)%*%D %x% Ip
  
  Xm <- matrix(0, n0, nr*len_d)
  for(i in 1:nr)
  {
    Xm[indexy == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[indexy == uindex[i],]
  }
  
  res0 <- glm(y ~ 0+x, family ="binomial")
 
  
  if(len_d == ncol(x))
  {
    beta00 <- rep(coef(res0),nr)
    lin <- Xm %*% beta00
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    beta_cur <- beta00
    for(r in 1:maxiter)
    {
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Xm) %*% Xm +  lam0 *AtA) %*% (t(Xm) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
      diffnorm <-  norm(beta_r - beta_cur,"F")
      beta_cur <- beta_r
      lin <-  Xm %*% beta_cur
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list( beta = beta_cur, betam = matrix(beta_cur, nr, byrow = TRUE))
  }
  
  if(len_d < ncol(x))
  {
    ncz <- ncol(x) - len_d
    
    z <- x[,-(index_d), drop = FALSE] ## common part

    eta00 <- coef(res0)[-(index_d)]
    beta00 <- rep(coef(res0)[index_d],nr)
    
    lin <- z %*% eta00 + Xm %*% beta00
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    eta_cur <- eta00
    beta_cur <- beta00
    for(r in 1:maxiter)
    {
      eta_r <- eta_cur + solve(t(as.numeric(V)*z) %*% z) %*% t(z) %*% ((y-mu))
      
      lin <- z %*% eta_r + Xm %*% beta_cur
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Xm) %*% Xm + lam0 *AtA) %*% (t(Xm) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
      diffnorm <- norm(eta_cur - eta_r,"F") + norm(beta_r - beta_cur,"F")
      
      eta_cur <- eta_r
      beta_cur <- beta_r
      lin <- z %*% eta_cur + Xm %*% beta_cur
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list(eta = eta_cur,beta = beta_cur,
                betam = matrix(beta_cur, nr, byrow = TRUE))
  }

  kk <- max(k0, round(sqrt(nr)))
  reskm <- pam(est$betam, kk, nstart = 50, stand = TRUE)
  cluster00 <- reskm$cluster
  
  ### refit this by ridge 
  initial0 <- fit_clust_reg_ridge(dom = indexy, y = y, x = x,
                                  index_d = index_d, cluster = cluster00,
                                  lam0 = lam0)
  
  # initial0 <- fit_clust_reg_ridge_pql(dom = indexy, y = y, x = x,
  #                                     index_d = index_d, cluster = cluster00,
  #                                     lam0 = lam0)
  # 

  return(initial0)
}


initialEM0 <- function(indexy, y, x, index_d, k0 = 10, lam0 = 0.01)
{
  df0 <- as.data.frame(cbind(as.factor(indexy), y, 1-y, x))
  colnames(df0) <- c("domain","y1","y2",paste("x",1:ncol(x),sep=""))
  form <- paste("cbind(y1,y2)~0+", paste(paste("x",1:ncol(x[,-1,drop = FALSE]),sep=""), collapse = "+"),"|domain",sep="")
  
  set.seed(123456)
  res_em <- flexmix(as.formula(form), data = df0, k = k0, 
                    model = FLXMRglm(family = "binomial"), 
                    control = list(
                      minprior = 0.01,
                      verbose = 0, 
                      nrep = 20))
  cluster00 <- as.numeric(as.factor(by(res_em@cluster, indexy,unique)))
  
  ### then ridge
  initial0 <- fit_clust_reg_ridge(dom = indexy, y = y, x = x,
                                  index_d = index_d, cluster = cluster00,
                                  lam0 = lam0)
  return(initial0)
}

### initial using EM and then PQL

initialEM <- function(indexy, y, x, index_d, k0 = 10, lam0 = 0.01)
{
  df0 <- as.data.frame(cbind(as.factor(indexy), y, 1-y, x))
  colnames(df0) <- c("domain","y1","y2",paste("x",1:ncol(x),sep=""))
  form <- paste("cbind(y1,y2)~0+", paste(paste("x",1:ncol(x[,-1,drop = FALSE]),sep=""), collapse = "+"),"|domain",sep="")
  
  set.seed(123456)
  res_em <- flexmix(as.formula(form), data = df0, k = k0, 
                           model = FLXMRglm(family = "binomial"), 
                    control = list(
                      minprior = 0.01,
                      verbose = 0, 
                      nrep = 20))
  cluster00 <- as.numeric(as.factor(by(res_em@cluster, indexy,unique)))
  
  ### then use pql 
  # initial0 <- fit_clust_reg_ridge(dom = indexy, y = y, x = x,
  #                                 index_d = index_d, cluster = cluster00, 
  #                                 lam0 = lam0)
  
  initial0 <- fit_clust_reg_ridge_pql(dom = indexy, y = y, x = x,
                                      index_d = index_d, cluster = cluster00,
                                      lam0 = lam0)
  
  
  return(initial0)
}



##### initial_pql 

initial_pql <- function(indexy, y, x, index_d,
                    lam0  = 0.01, k0 = 10,
                    tol = 1e-5, 
                    maxiter = 500)
{
  ### first fit a logistic regression with all are
  #in one group to get the starting values
  x <- as.matrix(x)
  x_d <- x[,index_d, drop = FALSE]
  
  len_d <- length(index_d)
  n0 <- length(indexy)
  uindex <- unique(indexy)
  nr <- length(uindex)
  ns <- as.numeric(table(indexy))
  idx_map <- match(indexy, uindex)
  
  Ip <- diag(1,len_d,len_d)
  D <- matrix(0,nr*(nr-1)/2,nr)
  for(j in 1:(nr-1))
  {
    indexj <- (2*nr-j)*(j-1)/2
    indexvj <- indexj + (1:(nr-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nr)] <- -1
  }
  Am <- D %x% Ip
  AtA <- t(D)%*%D %x% Ip
  
  Xm <- matrix(0, n0, nr*len_d)
  for(i in 1:nr)
  {
    Xm[indexy == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[indexy == uindex[i],]
  }
  
  res0 <- glmer(y ~ 0+x+(1|indexy), family ="binomial")
  
  #vv_cur <- rep(0, nr)
  re_df   <- ranef(res0)$indexy
  vv_cur  <- re_df[as.character(uindex), 1]
  
  tau_cur <- as.numeric(VarCorr(res0))

  
  if(len_d == ncol(x))
  {
    beta00 <- as.numeric(rep(fixef(res0),nr))
    lin <- Xm %*% beta00  + vv_cur[indexy]
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    beta_cur <- beta00
  
    for(r in 1:maxiter)
    {
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Xm) %*% Xm +  lam0 *AtA) %*% (t(Xm) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
      lin_xm <- Xm %*% beta_r
      lin <- lin_xm + vv_cur[idx_map]
      mu  <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      Vv<- as.numeric(V)
      
      ### update random effects
      tempv0 <- drop(rowsum(Vv, idx_map))
      tempv1 <- tempv0 + 1/tau_cur
      tempv2 <- drop(rowsum(as.numeric(y - mu), idx_map)) + tempv0 * vv_cur
      vv_r   <- tempv2 / tempv1
      
      ### update tau 
      lin <- lin_xm + vv_r[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      Vv <- as.numeric(V)
      
      ## update var_v
      tau_cur <- mean(1/(drop(rowsum(Vv, idx_map)) + 1/tau_cur) + vv_r^2)
      #diffnorm <-  norm(beta_r - beta_cur,"F")
      diffnorm <- norm(beta_r - beta_cur,"F") + sqrt(sum((vv_r - vv_cur)^2))
  
      beta_cur <- beta_r
      vv_cur <- vv_r
      
      lin <-  Xm %*% beta_cur + vv_cur[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list( beta = beta_cur, vv = vv_cur, betam = matrix(beta_cur, nr, byrow = TRUE))
  }
  
  if(len_d < ncol(x))
  {
    ncz <- ncol(x) - len_d
    
    z <- x[,-(index_d), drop = FALSE] ## common part
    
    eta00 <- fixef(res0)[-(index_d)]
    beta00 <- rep(fixef(res0)[index_d],nr)
    
    lin <- z %*% eta00 + Xm %*% beta00 + vv_cur[indexy]
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    eta_cur <- eta00
    beta_cur <- beta00
  
    
    for(r in 1:maxiter)
    {
      eta_r <- eta_cur + solve(t(as.numeric(V)*z) %*% z) %*% t(z) %*% ((y-mu))
      
      lin <- z %*% eta_r + Xm %*% beta_cur + vv_cur[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Xm) %*% Xm + lam0 *AtA) %*% (t(Xm) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
    
      lin_xm <- Xm %*% beta_r
      lin <- drop(z %*% eta_r) + lin_xm + vv_cur[idx_map]
      mu  <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      Vv <- as.numeric(V)
      
      ### update random effects
      tempv0 <- drop(rowsum(Vv, idx_map))
      tempv1 <- tempv0 + 1/tau_cur
      tempv2 <- drop(rowsum(as.numeric(y - mu), idx_map)) + tempv0 * vv_cur
      vv_r   <- tempv2 / tempv1
      
      lin <- lin_xm + vv_r[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      Vv <- as.numeric(V)
      
      tau_cur <- mean(1/(drop(rowsum(Vv, idx_map)) + 1/tau_cur) + vv_r^2)
      
      diffnorm <- norm(eta_cur - eta_r,"F") + norm(beta_r - beta_cur,"F")
      diffnorm <- diffnorm + sqrt(sum((vv_r - vv_cur)^2))
      
      eta_cur <- eta_r
      beta_cur <- beta_r
      vv_cur <- vv_r
      
      lin <- z %*% eta_cur + Xm %*% beta_cur + vv_cur[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list(eta = eta_cur,beta = beta_cur,  vv = vv_cur,
                betam = matrix(beta_cur, nr, byrow = TRUE))
  }
  
  set.seed(12345)
  kk <- max(k0, round(sqrt(nr)))
  reskm <- pam(est$betam, kk, nstart = 50, stand = TRUE)
  cluster00 <- reskm$cluster
  
  ### refit this by ridge 
  initial0 <- fit_clust_reg_ridge_pql(dom = indexy, y = y, x = x,
                                  index_d = index_d, cluster = cluster00,
                                  lam0 = lam0)
  
  
  return(initial0)
}


# ### adaptive using pam ###
# start_adaptive <- function(beta00, k0)
# {
# library(cluster)
#   nr <- nrow(beta00)
#   set.seed(123456)
#   #kk <- max(k0, round(sqrt(nr)))
#   kk <- min(k0, round(sqrt(nr)))
#   reskm <- pam(beta00, kk, nstart = 100, stand = TRUE)
#   #beta00 <- as.matrix(reskm$centers[reskm$cluster,])
#   cluster00 <- reskm$cluster
#   
#   #out <- list(beta00 = beta00, cluster00 = cluster00)
#   out <- list(cluster00 = cluster00)
#   return(out)
# }


##### fit the model without random effects and given cluster structure ### 

fit_clust_reg <- function(dom, y, x, index_d, cluster)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  len_d <- length(index_d)
  
  x_d <- x[,index_d, drop = FALSE]
  
  ng <- length(unique(cluster))
  Ux <- matrix(0, n0, ng*len_d)
  for(j in 1:ng)
  {
    uindexj <- uindex[cluster==j]
    indexj <- dom %in% uindexj
    Ux[indexj,(len_d*(j-1) + 1) : (len_d*j)] <- x_d[indexj,]
  }
  
  names_d <- paste(rep(paste("x",1:len_d,sep=""),ng), rep(paste("c",1:ng,sep=""), each = len_d),sep="_")
  
  ncx <- ncol(x)
  
  if(len_d == ncx)
  {
    Ux <- Ux
    colnames(Ux) <- names_d
  }else{
    Ux <- cbind(x[,-index_d, drop = FALSE], Ux)
    colnames(Ux) <- c(paste("z",1:(ncx-len_d),sep=""), names_d)
  }
  
  df_f <- as.data.frame(cbind(y, Ux))
  
  for1 <- as.formula(paste("y", "~0+", paste(colnames(Ux), collapse = "+")))
  res1 <- glm(for1, data = df_f, family = binomial)

  est_mat <- summary(res1)$coefficients
  
  if(len_d == ncx)
  {
    est_beta <- matrix(est_mat[,1], ng, byrow = TRUE)[cluster,,drop = FALSE]
    out <- list(betam = est_beta, est = est_mat, fit = res1)
  }else{
    
    est_beta <- matrix(est_mat[-(1:(ncx-len_d)),1], ng, byrow = TRUE)[cluster,,drop = FALSE]
    out <- list(eta = est_mat[1:(ncx-len_d),1],
                betam = est_beta, est = est_mat, fit = res1)
  }


  return(out)
}


##### fit the model without random effects plus a ridge regression form
####and given cluster structure ### 

fit_clust_reg_ridge <- function(dom, y, x, index_d, cluster, 
                                lam0 = 0.01, tol = 1e-5, maxiter = 500)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  len_d <- length(index_d)
  
  x_d <- x[,index_d, drop = FALSE]
  
  ng <- length(unique(cluster)) ### number of clusters
  Ux <- matrix(0, n0, ng*len_d)
  for(j in 1:ng)
  {
    uindexj <- uindex[cluster==j]
    indexj <- dom %in% uindexj
    Ux[indexj,(len_d*(j-1) + 1) : (len_d*j)] <- x_d[indexj,]
  }
  
  names_d <- paste(rep(paste("x",1:len_d,sep=""),ng), rep(paste("c",1:ng,sep=""), each = len_d),sep="_")
  

  ncx <- ncol(x)

  Ip <- diag(1,len_d,len_d)
  D <- matrix(0,ng*(ng-1)/2,ng)
  for(j in 1:(ng-1))
  {
    indexj <- (2*ng-j)*(j-1)/2
    indexvj <- indexj + (1:(ng-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):ng)] <- -1
  }
  Am <- D %x% Ip
  AtA <- t(D)%*%D %x% Ip
  
  res0 <- glm(y ~ 0+x, family ="binomial")
  
  if(len_d == ncol(x))
  {
    beta00 <- rep(coef(res0),ng)
    lin <- Ux %*% beta00
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    beta_cur <- beta00
    for(r in 1:maxiter)
    {
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Ux) %*% Ux +  lam0 *AtA) %*% (t(Ux) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
      diffnorm <-  norm(beta_r - beta_cur,"F")
      beta_cur <- beta_r
      lin <-  Ux %*% beta_cur
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list( beta = beta_cur)
  }

  
  if(len_d < ncol(x))
  {
    ncz <- ncol(x) - len_d
    z <- x[,-(index_d), drop = FALSE] ## common part
    
    eta00 <- coef(res0)[-(index_d)]
    beta00 <- rep(coef(res0)[index_d],ng)
    
    lin <- z %*% eta00 + Ux %*% beta00
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    eta_cur <- eta00
    beta_cur <- beta00
    for(r in 1:maxiter)
    {
      eta_r <- eta_cur + solve(t(as.numeric(V)*z) %*% z) %*% t(z) %*% ((y-mu))
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Ux) %*% Ux + lam0 *AtA) %*% (t(Ux) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
      diffnorm <- norm(eta_cur - eta_r,"F") + norm(beta_r - beta_cur,"F")
      eta_cur <- eta_r
      beta_cur <- beta_r
      lin <- z %*% eta_cur + Ux %*% beta_cur
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      
      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list(eta = eta_cur,beta = beta_cur)
  }
  
  betam <-  matrix(beta_cur, ng, byrow = TRUE)[cluster,,drop=FALSE]
  est$betam <- betam
  est$cluster <- cluster

  return(est)
}


### fit ridge with random effects 
fit_clust_reg_ridge_pql <- function(dom, y, x, index_d, cluster, 
                                lam0 = 0.01, tol = 1e-5, maxiter = 500)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  len_d <- length(index_d)
  idx_map <- match(dom, uindex)
  
  x_d <- x[,index_d, drop = FALSE]
  
  ng <- length(unique(cluster)) ### number of clusters
  Ux <- matrix(0, n0, ng*len_d)
  for(j in 1:ng)
  {
    uindexj <- uindex[cluster==j]
    indexj <- dom %in% uindexj
    Ux[indexj,(len_d*(j-1) + 1) : (len_d*j)] <- x_d[indexj,]
  }
  
  names_d <- paste(rep(paste("x",1:len_d,sep=""),ng), rep(paste("c",1:ng,sep=""), each = len_d),sep="_")
  
  
  ncx <- ncol(x)
  
  Ip <- diag(1,len_d,len_d)
  D <- matrix(0,ng*(ng-1)/2,ng)
  for(j in 1:(ng-1))
  {
    indexj <- (2*ng-j)*(j-1)/2
    indexvj <- indexj + (1:(ng-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):ng)] <- -1
  }
  Am <- D %x% Ip
  AtA <- t(D)%*%D %x% Ip
  
  res0 <- glmer(y ~ 0+x+ (1|dom), family ="binomial")
  vv_cur <- rep(0, nr)
  tau_cur <- as.numeric(VarCorr(res0))
  
  if(len_d == ncol(x))
  {
    beta00 <- rep(fixef(res0),ng)
    lin <- Ux %*% beta00
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    beta_cur <- beta00
    for(r in 1:maxiter)
    {
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Ux) %*% Ux +  lam0 *AtA) %*% (t(Ux) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
      lin_xm <- Ux %*% beta_r
      lin <- lin_xm + vv_cur[idx_map]
      mu  <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      Vv <- as.numeric(V)
      
      ### update random effects
      tempv0 <- drop(rowsum(Vv, idx_map))
      tempv1 <- tempv0 + 1/tau_cur
      tempv2 <- drop(rowsum(as.numeric(y - mu), idx_map)) + tempv0 * vv_cur
      vv_r   <- tempv2 / tempv1
      
      ### update tau 
      lin <- lin_xm + vv_r[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      Vv <- as.numeric(V)
      
      ## update var_v
      tau_cur <- mean(1/(drop(rowsum(Vv, idx_map)) + 1/tau_cur) + vv_r^2)
      diffnorm <-  norm(beta_r - beta_cur,"F")
      
      beta_cur <- beta_r
      vv_cur <- vv_r
      
      lin <-  Ux %*% beta_cur + vv_cur[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      

      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list(beta = beta_cur, vv = vv_cur)
  }
  
  
  if(len_d < ncol(x))
  {
    ncz <- ncol(x) - len_d
    z <- x[,-(index_d), drop = FALSE] ## common part
    
  
    eta00 <- fixef(res0)[-(index_d)]
    beta00 <- rep(fixef(res0)[index_d],ng)
    
    lin <- z %*% eta00 + Ux %*% beta00 + vv_cur[idx_map]
    mu <- 1/(1+exp(-lin))
    V <- mu * (1-mu)
    
    eta_cur <- eta00
    beta_cur <- beta00
    
    for(r in 1:maxiter)
    {
      eta_r <- eta_cur + solve(t(as.numeric(V)*z) %*% z) %*% t(z) %*% ((y-mu))
      beta_r <- beta_cur +
        solve(t(as.numeric(V)*Ux) %*% Ux + lam0 *AtA) %*% (t(Ux) %*% ((y-mu)) - lam0*AtA %*%beta_cur)
      
      lin_xm <- Ux %*% beta_r
      lin <- drop(z %*% eta_r) + lin_xm + vv_cur[idx_map]
      mu  <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
      Vv <- as.numeric(V)
      
      ### update random effects
      tempv0 <- drop(rowsum(Vv, idx_map))
      tempv1 <- tempv0 + 1/tau_cur
      tempv2 <- drop(rowsum(as.numeric(y - mu), idx_map)) + tempv0 * vv_cur
      vv_r   <- tempv2 / tempv1
      
      tau_cur <- mean(1/(drop(rowsum(Vv, idx_map)) + 1/tau_cur) + vv_r^2)
      
      diffnorm <- norm(eta_cur - eta_r,"F") + norm(beta_r - beta_cur,"F")
      
      eta_cur <- eta_r
      beta_cur <- beta_r
      vv_cur <- vv_r

      lin <- z %*% eta_cur + Ux %*% beta_cur + vv_cur[idx_map]
      mu <- 1/(1+exp(-lin))
      V <- mu * (1-mu)
    
      if(diffnorm <= tol)
      {
        break
      }
      
    }
    
    est <- list(eta = eta_cur,beta = beta_cur, vv = vv_cur)
  }
  
  betam <-  matrix(beta_cur, ng, byrow = TRUE)[cluster,,drop=FALSE]
  est$betam <- betam
  est$cluster <- cluster
  
  return(est)
}
