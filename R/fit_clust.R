##### fit the model with a given cluster structure ###

fit_clust <- function(dom, y, x, index_d, cluster)
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
  
  if(len_d == ncol(x))
  {
    Ux <- Ux
    colnames(Ux) <- names_d
  }else{
    Ux <- cbind(x[,-index_d, drop = FALSE], Ux)
    colnames(Ux) <- c(paste("z",1:(ncol(x)-len_d),sep=""), names_d)
  }
  
  df_f <- as.data.frame(cbind(y, Ux))
  
  for1 <- as.formula(paste("y", "~0+", paste(c(colnames(Ux),"(1|dom)"), collapse = "+")))
  # res1 <- glmer(for1, data = df_f, family = binomial, 
  #               control = glmerControl(optimizer = "nloptwrap",
  #                                      check.conv.singular = .makeCC(action = "ignore", tol = 1e-10)))
  # 
  res1 <- glmer(for1, data = df_f, family = binomial, 
                control = glmerControl(optimizer = "nlminbwrap",
                                       check.conv.grad = "ignore",
                                       check.conv.singular = "ignore",
                                       check.conv.hess = "ignore"))
  
  randv <- ranef(res1)
  beta_est <- fixef(res1)
  est_mat <- summary(res1)$coefficients
  
  out <- list(est = est_mat, rand = randv$dom[,1], fit = res1)
  return(out)
}


#### for a given cluster, give the estimate and give the area estimate 

fit_clust_area <- function(dom, y, x, xunit, index_d, cluster)
{
  m <- length(unique(dom))
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
  
  if(len_d == ncol(x))
  {
    Ux <- Ux
    colnames(Ux) <- names_d
  }else{
    Ux <- cbind(x[,-index_d, drop = FALSE], Ux)
    colnames(Ux) <- c(paste("z",1:(ncol(x)-len_d),sep=""), names_d)
  }
  
  df_f <- as.data.frame(cbind(y, Ux))
  
  for1 <- as.formula(paste("y", "~0+", paste(c(colnames(Ux),"(1|dom)"), collapse = "+")))
  #res1 <- glmer(for1, data = df_f, family = binomial, 
   #             control = glmerControl(optimizer = "nloptwrap",
    #                                   check.conv.singular = .makeCC(action = "ignore", tol = 1e-10)))
  
  res1 <- glmer(for1, data = df_f, family = binomial, 
                control = glmerControl(optimizer = "nlminbwrap",
                                       check.conv.grad = "ignore",
                                       check.conv.singular = "ignore",
                                       check.conv.hess = "ignore"))
  
  
  
  randv <- ranef(res1)$dom[,1]
  beta_est <- fixef(res1)
  est_mat <- summary(res1)$coefficients ### estimate 
  
  xall <-  xunit[,-1]
  len_d <- length(index_d)
  ncx <- ncol(xall)
  Ni <- as.numeric(table(xunit[,1]))
  
  if(len_d == ncx)
  {
    subfun1 <- function(betam_est, vi_est)
    {
      prob_est <- rep(0,m)
      for(i in 1:m)
      {
        etai_est <- xall[xunit[,1]==i,] %*% betam_est[i,] + vi_est[i]
        prob_est[i] <- sum(1/(1+exp(-etai_est)))
      }
      return(prob_est)
    }
    
    subfun2 <- function(betam_est, vi_est)
    {
      prob_est <- rep(0,m)
      for(i in 1:m)
      {
        etai_est <- x[dom==i,] %*% betam_est[i,] + vi_est[i]
        prob_est[i] <- sum(1/(1+exp(-etai_est)))
      }
      return(prob_est)
    }
    
 
    betam <- matrix(beta_est, ncol = len_d, byrow= TRUE)[cluster,,drop=FALSE]
      
    estimator1 <- subfun1(betam, randv)/Ni
    estimator2 <- (subfun1(betam, randv) +  as.numeric(by(y, dom, sum))-subfun2(betam, randv))/Ni
      
      
    estimator_mat <- cbind(estimator1, estimator2)
    
    
  }
  
  if(len_d < ncx)
  {
    x_d <- xall[,index_d, drop=FALSE]
    z <- xall[,-(index_d), drop = FALSE]
    
    xs_d <- x[,index_d, drop = FALSE]
    zs <- x[,-c(index_d),drop=FALSE]
    
    ind_all <- xunit[,1]
    subfun1 <- function(eta_est, betam_est, vi_est)
    {
      prob_est <- rep(0,m)
      for(i in 1:m)
      {
        etai_est <- z[ind_all==i,] %*% eta_est + x_d[ind_all==i,,drop = FALSE] %*% betam_est[i,] + vi_est[i]
        prob_est[i] <- sum(1/(1+exp(-etai_est)))
      }
      return(prob_est)
    }
    
    subfun2 <- function(eta_est, betam_est, vi_est)
    {
      prob_est <- rep(0,m)
      for(i in 1:m)
      {
        etai_est <- zs[dom==i,] %*% eta_est + xs_d[dom==i,,drop = FALSE] %*% betam_est[i,] + vi_est[i]
        prob_est[i] <- sum(1/(1+exp(-etai_est)))
      }
      return(prob_est)
    }
    
   
      betam <- matrix(beta_est[-(1:(ncx-len_d))], ncol = len_d, byrow= TRUE)[cluster,,drop=FALSE]
      eta <- beta_est[1:(ncx-len_d)]
    
      estimator1 <- subfun1(eta, betam, randv)/Ni
      estimator2 <- (subfun1(eta,betam, randv) +  
                       as.numeric(by(y, dom, sum))-subfun2(eta,betam, randv))/Ni
      
      estimator_mat <- cbind(estimator1, estimator2)
  
  }
  return(estimator_mat)
}


