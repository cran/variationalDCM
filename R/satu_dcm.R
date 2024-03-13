

satu_dcm = function(X,
                    Q,
                    max_it  = 500,
                    epsilon = 1e-04,
                    verbose = TRUE,
                    # hyperparameter
                    delta_0 = NULL,
                    A_0 = NULL,
                    B_0 = NULL
){


  if(!inherits(X, "matrix")){
    X <- as.matrix(X)
  }
  if(!inherits(Q, "matrix")){
    Q <- as.matrix(Q)
  }

  if(!all(X %in% c(0,1)))
    stop("item response data should only contain 0/1. \n")

  if(!all(Q %in% c(0,1)))
    stop("Q-matrix should only contain 0/1. \n")

  I <- nrow(X);  J <- nrow(Q);  K <- ncol(Q);  L <- 2^K
  not_zero_q <- apply(Q, 1, function(x) which(x != 0))


  # Attribute pattern matrix

  A <- as.matrix(expand.grid(lapply(1:K, function(x)rep(0:1))))
  A_j <- lapply(not_zero_q, function(x) A[,x,drop=F])
  A_red <- lapply(A_j , function(x) apply(x,1,function(y) paste0(y,collapse = "")))
  A_red_uni <- lapply(A_j , function(x) unique(apply(x,1,function(y) paste0(y,collapse = ""))))

  # Make G-matrix

  G_j <- lapply(1:J, function(j) t(sapply(A_red_uni[[j]], function(x) x == A_red[[j]] ))*1 )
  att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
  for(j in 1:J) {
    colnames(G_j[[j]]) <- att_pat
  }

  # Hyper parameters

  if(is.null(delta_0)){
    delta_0 = rep(1,L) # For π
  }

  # Weekly informative prior

  number_of_attributes <- lapply(A_red_uni, function(x) sapply(strsplit(x, ""), function(y)sum(as.numeric(y))) )
  if(is.null(A_0)){
    A_0_hyperparam <- seq(from = 1+epsilon, to = 2, length.out = max(unlist( number_of_attributes))+1)
    A_0 <- lapply( number_of_attributes , function(x){A_0_hyperparam[x + 1] })
  }
  if(is.null(B_0)){
    B_0_hyperparam <- seq(from = 2, to = 1+epsilon, length.out = max(unlist( number_of_attributes))+1)
    B_0 <- lapply( number_of_attributes , function(x){B_0_hyperparam[x + 1] })
  }

  # Initialization

  r_il <- matrix(1/L ,ncol=L, nrow = I)
  one_vec  = matrix(1,nrow=I,ncol=1)

  # lower bound of log marginal likelihood

  llb_fun <- function(X,G_j,delta_ast,delta_0, A_ast, A_0, B_ast, B_0, r_il){
    tmp1 <- 0
    for(j in 1:length(G_j)){
      tmp1 <- tmp1+ sum(((X[,j]%*%(digamma(A_ast[[j]]) - digamma(A_ast[[j]]+B_ast[[j]])) + (1-X[,j])%*%(digamma(B_ast[[j]]) - digamma(A_ast[[j]]+B_ast[[j]]))) %*% G_j[[j]]) * r_il)
    }
    tmp2 <- sum(r_il * ( one_vec%*%(digamma(delta_ast) - digamma(sum(delta_ast))) - log(r_il)))
    tmp3 <- sum(lgamma(delta_ast))- lgamma(sum(delta_ast)) - (sum(lgamma(delta_0))- lgamma(sum(delta_0))) + sum((delta_0 - delta_ast)*(digamma(delta_ast) - digamma(sum(delta_ast))))

    A_ast_unlist <- unlist(A_ast)
    B_ast_unlist <- unlist(B_ast)
    A_0_unlist <- unlist(A_0)
    B_0_unlist <- unlist(B_0)
    tmp4 <- sum(lbeta(a =A_ast_unlist,  b=B_ast_unlist) - lbeta(a =A_0_unlist,  b=B_0_unlist) +  (A_0_unlist - A_ast_unlist)*(digamma(A_ast_unlist) - digamma(A_ast_unlist + B_ast_unlist ) ) + (B_0_unlist - B_ast_unlist)*(digamma(B_ast_unlist) - digamma(A_ast_unlist + B_ast_unlist ) ) )

    tmp1 + tmp2 + tmp3 + tmp4  }
  l_lb = rep(0, max_it+1)
  l_lb[1] = -Inf

  m = 1
  for(m in 1:max_it){

    # VM-step
    delta_ast <- colSums(r_il) + delta_0
    A_ast  <- lapply(1:J, function(j) t(G_j[[j]] %*% t(r_il) %*% X[,j]     + A_0[[j]] ))
    B_ast  <- lapply(1:J, function(j) t(G_j[[j]] %*% t(r_il) %*% (1-X[,j]) + B_0[[j]] ))

    E_log_pi      = digamma(delta_ast) - digamma(sum(delta_ast))
    E_log_theta   = lapply(1:J, function(j) digamma(A_ast[[j]]) - digamma(A_ast[[j]] + B_ast[[j]]))
    E_log_1_theta = lapply(1:J, function(j) digamma(B_ast[[j]]) - digamma(A_ast[[j]] + B_ast[[j]]))

    # VE-step: r_il
    temp <- matrix(0, ncol=I, nrow = L)
    for(j in 1:J) temp <- temp + t(G_j[[j]] ) %*% (t(E_log_theta[[j]]) %*% t(X[,j,drop=F]) +  t(E_log_1_theta[[j]]) %*% t(1 - X[,j,drop=F]))
    log_rho_il <- t(temp + E_log_pi)
    temp <- exp(log_rho_il)
    r_il =  temp / rowSums(temp)

    l_lb[m+1] <- llb_fun(X,G_j,delta_ast,delta_0, A_ast, A_0, B_ast,B_0, r_il)

    if(verbose){
      cat("\riteration = ", m+1,sprintf(",last change = %.05f", abs(l_lb[m] - l_lb[m+1])))
    }

    if(abs(l_lb[m] - l_lb[m+1]) < epsilon){
      if(verbose){
        cat("\nreached convergence.\n")
      }
      break()
    }
  }
  l_lb <- l_lb[-1]

  delta_sum <- sum(delta_ast)
  pi_est <-  delta_ast/delta_sum
  pi_sd <-sqrt(delta_ast*(delta_sum - delta_ast)/(delta_sum^2*(delta_sum+1)) )
  theta_est <- mapply(function(x,y) x/(x+y), A_ast, B_ast)
  theta_sd <- mapply(function(x,y) sqrt((x*y)/(((x+y)^2) *(x+y+1)) ),  A_ast, B_ast)

  model_params = list(
    theta_est = theta_est,
    theta_sd = theta_sd
  )

  res = list(model_params = model_params,
       pi_est = pi_est,
       pi_sd = pi_sd,
       #r_il  = r_il,
       A_ast = A_ast,
       B_ast = B_ast,
       delta_ast = delta_ast,
       A_0 = A_0,
       B_0 = B_0,
       delta_0 = delta_0,
       l_lb = l_lb[l_lb != 0],
       att_pat_est = A[apply(r_il, 1, which.max),],
       G_j = G_j,
       m = m)

  return(res)
}
