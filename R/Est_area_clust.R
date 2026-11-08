######## unit level model based estimator in survey sampling ######
### area index 
## y x are sampled y and x, x should be include 1 if intercept is needed. 
# xunit: population level data with the first column as area index,  should include 1 if intercept is needed, 
# and the order should be same as x
## intercept should be included if intercept is included in the model
# bic based on l_pql

Est_area_clust <- function(area, y, x, xunit,
                      index_d, lamvec,
                      initial0 = NULL, 
                      vv0 = rep(0, length(unique(area))),
                      var_v0 = max(0.001, var(vv0)), 
                      nu = 1, gam = 3, lam0 = 0.01, 
                      maxiter= 300, inner_loop= 1, 
                      tol = 1e-3)
{
  #### fit model ###
  m <- length(unique(area))
  res <- Clust_nested_logistic_bicw(indexy = area, y, x,
                                    index_d = index_d , lamvec = lamvec,
                                    initial0 = initial0, vv0 = vv0, var_v0 = var_v0,
                                    nu = nu, gam = gam, lam0 = lam0, 
                                    maxiter= maxiter, inner_loop= inner_loop, 
                                    tol = tol)
  len_d <- length(index_d)
  n0 <- length(y)
  Ni <- as.numeric(table(xunit[,1]))
  
  
  bic_df <- res$df
  bic_df$bic <- 2*bic_df$nll + log(m*len_d)*log(m)*bic_df$Khat*len_d
  
  lamselect <- which.min(bic_df$bic)
  out_select <- res$out[[lamselect]] #### only one 
  cluster_est <- out_select$cluster ### estimated cluster
  
  betam <- out_select$betam
  rand <- out_select$rand
  
  
  ## use the result to construct area level estimators
 xall <-  xunit[,-1]
 len_d <- length(index_d)
 ncx <- ncol(x)
 
 if(len_d == ncx)
 {
   ### refit if can 
   res_refit <- NULL
   tryCatch({
     res_refit <- fit_clust(dom = area, y = y, x = x, index_d = index_d,
                            cluster = cluster_est)
     rand <- res_refit$rand
     betam <- matrix(res_refit$est[,1], ncol = len_d, byrow= TRUE)[cluster_est,,drop=FALSE]
   }, error = function(e) {
   })
   
   
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
       etai_est <- x[area==i,] %*% betam_est[i,] + vi_est[i]
       prob_est[i] <- sum(1/(1+exp(-etai_est)))
     }
     return(prob_est)
   }
   
   estimator1 <- subfun1(betam, rand)/Ni
   estimator2 <- (subfun1(betam, rand) +  as.numeric(by(y, area, sum))-subfun2(betam, rand))/Ni
   
   estimatorm <- cbind(estimator1, estimator2)
 }
 
 if(len_d < ncx)
 {
   
   eta <- out_select$eta
   res_refit <- NULL
   tryCatch({
     res_refit <- fit_clust(dom = area, y = y, x = x, index_d = index_d,
                            cluster = cluster_est)
     betam <- matrix(res_refit$est[-(1:(ncx-len_d)),1], ncol = len_d, byrow= TRUE)[cluster_est,,drop=FALSE]
     rand <- res_refit$rand
     eta <- res_refit$est[1:(ncx-len_d),1]
   }, error = function(e) {})
   
   
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
       etai_est <- z[ind_all==i,,drop=FALSE] %*% eta_est + x_d[ind_all==i,,drop=FALSE] %*% betam_est[i,] + vi_est[i]
       prob_est[i] <- sum(1/(1+exp(-etai_est)))
     }
     return(prob_est)
   }
   
   subfun2 <- function(eta_est, betam_est, vi_est)
   {
     prob_est <- rep(0,m)
     for(i in 1:m)
     {
       etai_est <- zs[area==i,,drop=FALSE] %*% eta_est + xs_d[area==i,,drop=FALSE] %*% betam_est[i,] + vi_est[i]
       prob_est[i] <- sum(1/(1+exp(-etai_est)))
     }
     return(prob_est)
   }
   

   estimator1 <- subfun1(eta, betam, rand)/Ni
   estimator2 <- (subfun1(eta,betam, rand) +  
                    as.numeric(by(y, area, sum))-subfun2(eta,betam, rand))/Ni
   
   estimatorm <- cbind(estimator1, estimator2)
 }

  out <- list(bic = bic_df, 
              refit = res_refit,
              out = out_select, 
              est = estimatorm,
              lam = lamvec,
              initial00 = res$initial00)
  return(out)
}

