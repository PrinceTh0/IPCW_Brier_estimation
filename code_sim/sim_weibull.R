
get.WeibullData <- function(n, n_cov_informative, n_cov_noninformative, s, s_c, rate, offset, admin_cens, structure, cens_distr, seed){
  
  ##packages
  library(evd)
  library(MASS)
  library(data.table)
  
  if(seed == T){
    set.seed(1234)
  }
  
  cutoff <- 0.9 # cutoff point for administrative censoring (quantile of observed times)
  
  n_cov <- n_cov_informative + n_cov_noninformative

  ##regression coefficients
  beta <- c(seq(-0.5,-0.1,0.1), seq(0.1,0.5,0.1))
  beta_c <- c(seq(-0.5,-0.1,0.1), seq(0.1,0.5,0.1))

  ##standard extreme value distribution
  e   <- - rgumbel(n, loc = 0, scale = 1)
  
  ##specify covariance structure
  if(structure == "compoundSymmetry"){
    # uniform covariance values (0.5)
    Sigma  <- diag(n_cov)
    n_triangle <- sum(1:(n_cov-1)) 
    Sigma[lower.tri(Sigma)] <- Sigma[upper.tri(Sigma)] <- rep(0.5, n_triangle)
    
  } else if(structure == "unstructured"){
    # non-uniform covariance values
    pos_def_corr_mat <- function(d,k){
      # returns a positive-definite correlation matrix of dimension dxd
      # k determines distribution of correlation values
      if(k>=d){
        stop("d must be greater than k")
      } else{
        W <- matrix(rnorm(d**2),nrow=d,ncol=k)
        S <- W%*%t(W) + diag(runif(d))
        S <- diag(1/sqrt(diag(S))) %*% S %*% diag(1/sqrt(diag(S)))
      }
      
      return(S)
    }
    
    Sigma <- pos_def_corr_mat(n_cov, 4)
  }
  
  ##matrix of predictor variables
  X <- matrix(mvrnorm(n, mu = rep(0,n_cov), Sigma = Sigma), ncol = n_cov)
  
  if(n_cov_noninformative == 0){
    X_inf <- X
  }else{
    X_inf <- X[,1:n_cov_informative]
  }
  
  d <- model.matrix(~ X)
  d <- data.table(d)
  
  ##generate event times with Weibull baseline hazard
  d$logT   <- X_inf%*%beta + e*s # on log scale
  d$lambda <- exp(-X_inf%*%beta/s)
  R2       <- c(var(X_inf%*%beta)/var(d$logT))
  d$Ttrue  <- exp(d$logT) # on arithmetic scale
  
  #administrative censoring
  d[, Cadmin := max(d[Ttrue <= quantile(d$Ttrue,cutoff), Ttrue])]
  
  #indicator variable for uncensored data
  d[, event_Ttrue := 1]

  ##generate censoring times
  if(cens_distr == "none"){
    if(admin_cens == T){
      d$event <- (d$Ttrue <= d$Cadmin) + 0
      d$Ttilde <- pmin(d$Ttrue, d$Cadmin)
    }else if(admin_cens == F){
      d[, event := 1]
      d[, Ttilde := Ttrue]
    }
    }else if(cens_distr == "exponential"){
      #censoring times, exponential baseline hazard (no covariates)
      d$Ctrue <- rexp(n, rate = rate)
      R2c     <- NA
      if(admin_cens == T){
        d$event <- (d$Ttrue <= d$Ctrue & d$Ttrue <= d$Cadmin) + 0
        d$Ttilde <- pmin(d$Ttrue, d$Ctrue, d$Cadmin) # observed times
      }else if(admin_cens == F){
        d$event <- (d$Ttrue <= d$Ctrue) + 0
        d$Ttilde <- pmin(d$Ttrue, d$Ctrue) # observed times
      }
    }else if(cens_distr == "Weibull"){
      #censoring times, Weibull baseline hazard (with covariates)
      d$Ctrue <- exp(offset + X_inf%*%beta_c + ((- rgumbel(n, loc = 0, scale = 1))*s_c))
      R2c     <- c(var(X_inf%*%beta_c)/var(log(d$Ctrue)))
      if(admin_cens == T){
        d$event <- (d$Ttrue <= d$Ctrue & d$Ttrue <= d$Cadmin) + 0
        d$Ttilde <- pmin(d$Ttrue, d$Ctrue, d$Cadmin) # observed times
      }else if(admin_cens == F){
        d$event <- (d$Ttrue <= d$Ctrue) + 0
        d$Ttilde <- pmin(d$Ttrue, d$Ctrue) # observed times
      }
    }
  
  return(list("data" = d, "R2" = R2, "R2c" = R2c))
  
}

##example
#validation <- get.WeibullData(n=100000,n_cov_informative=10,n_cov_noninformative=0,s=0.578,s_c=0.579,rate=0.864,offset=-0.002,admin_cens=F,structure="compoundSymmetry",cens_distr="Weibull",seed = F)$data
#setorder(validation, Ttilde)
