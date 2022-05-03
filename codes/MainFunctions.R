  
  
  
  estimating_psi = function(X, y){
    
    ##### Input: The independent data (X) and the target variable (y).
    ##### X-data is stored as a n x p matrix, y is kept as a vector ( n x 1).
    
    ##### Output: Getting \psi_{ij} from all p-marginal runs, only for -
    ##### the slope-based effects.
    
    ###### Objective: This function performs a MMM-based estimation.
    
    n = dim(X)[1] 
    
    p = dim(X)[2]
    Psi = c()
    
    betas <- rep(0, p)
    
    for(j in 1:p){
      
        g <- glm(y ~ X[, j], family = binomial(link = 'logit'))
        
        ######################################## 
        ## Gather all marginal estimates     ##
        ######################################## 
        
        betas[j] <- g$coefficients[2]
      
        ######################################## 
        ## Score contributions                ##
        ######################################## 
        
        scores =  (as.vector(y) - g$fitted.values) * cbind(1, X[, j])
        
        ######################################## 
        ## Extracting the relevant Fisher row ##
        ######################################## 
        
        scaling =  summary(g)$cov.unscaled
        
        ######################################## 
        ## Take out only the second row - this corresponds to the slope! ##
        ######################################## 
      
        iidcontribution  = scaling %*% t(scores) * sqrt(n)
        
        Psi = rbind(Psi, t(iidcontribution[2,]))
        
    }
    
    return(list(betas = betas, Psi = Psi, num_obs = n, num_feat = p))
  }
  ################################################################################
  marginal_covarinace_matrix = function(psi, n = num_obs, p = num_feat){
  
    ##### Input: All \psi_{ij} from the 'estimating_psi' function.
   
    ##### Output: The covariance matrix from among the test statistics. 
    
    ###### Objective: We aim at getting the null distribution.
    
    
    Sigma <- matrix(rep(0, p * p), nrow = p)
    tmp_psi = vector()
    
    for(i in 1:p){
      for(j in 1:p){
        
        if(i == j ){
          
          Sigma[i, i] = (t(psi[i,]) %*% psi[i,])/(n)
        }
        
        if( i < j ){
          
          Sigma[i, j] = Sigma[j, i] =  (t(psi[i,]) %*% psi[j,])/(n)
  
        }
      }
      
    }
    
    return(Sigma)
  }
  
  
  
  ###########################################################
  # Compute the quantity Ahat from Eq. (55) in Efron (2007) #
  ###########################################################
  
  Ahat <- function(x0, z){
    N  <- length(z)
    Y0 <- sum((z < x0) & (z > -x0));
    P0 <- 2 * pnorm(x0) - 1;
    Q0 <- sqrt(2) * x0 * dnorm(x0);
    P0hat <- Y0 / N;
    return((P0-P0hat) / Q0);
  }
  
  
  
