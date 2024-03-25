#
# Make G matrix funciton -----
#
make_G_mat <- function(Q){
  J <- max(Q[,1])
  K <- ncol(Q) -2
  L <- 2^K

  H <- lapply(1:J, function(x) Q[Q[,1] == x ,2])

  #
  # All attribute pattern matrix
  #
  A <- as.matrix(expand.grid(lapply(1:K, function(x)rep(0:1))))

  G_mat <- vector("list", J)
  for(j in 1:J){
    Q_j <- as.matrix(Q[Q[,1] == j,-c(1:2)])

    if( all(rowSums(Q_j) != 0)){
      Q_j <- rbind(Q_j, 0)
    }

    not_0_col <- colSums(Q_j) != 0
    Q_j    <- Q_j[,not_0_col, drop = F]
    A_temp <- A[,not_0_col,drop =F]

    G_j <- matrix(0,nrow=nrow(Q_j), ncol=L)

    q_patt <- apply(Q_j[rowSums(Q_j) != 0, ,drop = F], 1, function(x)paste0(x,collapse=""))
    alp_patt <- apply(A_temp, 1, function(x)paste0(x,collapse=""))
    alp_patt_not_in_q_patt <- !(alp_patt %in% q_patt)

    for(h in H[[j]]){
      for(l in 1:L){
        if(any(Q_j[h,] == 1)){
          G_j[h,l] <- all(Q_j[h,] == A_temp[l,])*1
        } else {
          G_j[h,l] <- alp_patt_not_in_q_patt[l]*1
        }
      }
    }
    G_mat[[j]] <- G_j
  }

  list(G_mat = G_mat, H = H,  A  = A)
}

#' @title Artificial data generating function for the multiple-choice DINA model based on the given Q-matrix
#'
#' @description \code{mc_dina_data_gen()} returns the artificially generated item response data for the MC-DINA model
#'
#' @param Q the \eqn{J \times K} binary matrix
#' @param I the number of assumed respondents
#' @param att_cor the true value of the correlation among attributes (default: 0.1)
#' @param seed the seed value used for random number generation (default: 17)
#' @return A list including:
#' \describe{
#'   \item{X}{the generated artificial item response data}
#'   \item{att_pat}{the generated true vale of the attribute mastery pattern}
#' }
#' @references Yamaguchi, K. (2020). Variational Bayesian inference for the
#'   multiple-choice DINA model. \emph{Behaviormetrika}, 47(1), 159-187.
#'   \doi{10.1007/s41237-020-00104-w}
#'
#' @examples
#' # load a simulated Q-matrix
#' mc_Q = mc_sim_Q
#' mc_sim_data = mc_dina_data_gen(Q=mc_Q,I=200)
#'
#' @export

#
# Data generation function------
#
mc_dina_data_gen <- function(I, Q, att_cor = 0.1, seed = 17){
  set.seed(seed)

  J <- max(Q[,1])
  K = ncol(Q) - 2
  L = 2^K
  match_p = 0.80

  att_threshold = "rand"
  tmp <- make_G_mat(Q)
  A <- tmp$A
  H <- tmp$H
  G_mat <- tmp$G_mat

  #
  # item parameters
  #
  i_par <- vector("list", J)

  not_match = "equal"
  if( not_match == "equal"){
    for(j in 1:J){
      Q_j <- as.matrix(Q[Q[,1] == j,-c(1:2)])

      not_0_col <- colSums(Q_j) != 0
      Q_j    <- Q_j[,not_0_col, drop = F]
      A_temp <- A[,not_0_col,drop =F]

      q_patt <- apply(Q_j[rowSums(Q_j) != 0, ,drop = F], 1, function(x)paste0(x,collapse=""))
      alp_patt <- apply(A_temp, 1, function(x)paste0(x,collapse=""))
      alp_patt_not_in_q_patt <- !(alp_patt %in% q_patt)


      dim_G_mat <- dim(G_mat[[j]])
      i_par_temp <- array(NA, dim_G_mat)

      H_j <- H[[j]]
      for(h in H_j){
        for(l in 1:L){
          #
          # not included
          #
          if(alp_patt_not_in_q_patt[l]){
            i_par_temp[h,l] <- 1/nrow(Q_j)
          }else if(!alp_patt_not_in_q_patt[l] & (G_mat[[j]][h,l] == 1)){
            #
            # match
            #
            i_par_temp[h,l] <- match_p
          }else{
            # others
            i_par_temp[h,l] <- (1-match_p)/(nrow(Q_j) -1)
          }
        }
      }

      i_par[[j]] <- i_par_temp
    }

  }else if (not_match == "rand"){

    for(j in 1:J){

      dim_G_mat <- dim(G_mat[[j]])
      i_par_temp <- array(NA, dim_G_mat)

      match_p_temp <- stats::runif(dim_G_mat[1], min=0.5, max=1.0)
      i_par_temp  <- diag(match_p_temp)
      i_par_temp[lower.tri(i_par_temp)] <- stats::runif(sum(lower.tri(i_par_temp)), min=0, max=0.5)
      i_par_temp[upper.tri(i_par_temp)] <- stats::runif(sum(upper.tri(i_par_temp)), min=0, max=0.5)
      i_par_temp[,ncol(i_par_temp)] <- stats::runif(nrow(i_par_temp),min=0.2,max=0.8)
      i_par_temp <- t(t(i_par_temp)/colSums(i_par_temp))

      i_par[[j]]  <-  i_par_temp %*% G_mat[[j]]
    }
  } else {
    stop("Error: You should specify not_match as \"equal\" or \"rand\" ")
  }


  sigma <- diag(K)
  sigma[sigma == 0] <- att_cor

  if(any(att_threshold == "rand")){
    att_threshold = stats::rnorm(K, 0, 1)
  } else if(is.null(att_threshold)) {
    att_threshold = rep(0,K)
  } else if(is.vector(att_threshold)){

  }
  att_value <- mvtnorm::rmvnorm(n = I,mean = rep(0,K),sigma = sigma)
  att_threshold
  att_pat <- t((t(att_value) > att_threshold)*1)
  att_all_pat <- apply(A, 1, function(x)paste0(x,collapse = ""))
  cluss_num <- apply(att_pat, 1, function(x) which(paste0(x,collapse = "") == att_all_pat))

  Z <- matrix(0, ncol=L,nrow=I)
  for(i in 1:I){
    Z[i,cluss_num[i]] <- 1
  }

  X <- matrix(NA, ncol=J, nrow = I)
  for(i in 1:I){
    for(j in 1:J){
      jh_prob <- H[[j]]
      jh_prob <- i_par[[j]] %*% Z[i, ]
      X[i,j] <- sample(x = H[[j]],size = 1, prob = jh_prob)
    }
  }
  list(X = X,att_pat = att_pat)
}


#
extend_X <- function(X){
  I <- nrow(X)
  J <- ncol(X)
  H_max <- max(X, na.rm = T)
  X_ijh <- array(NA, dim = c(I,J,H_max))
  for(i in 1:I){
    tmp <- matrix(0, nrow=J, ncol=H_max)
    for(j in 1:J){
      tmp[j,X[i,j]] <- 1
    }
    X_ijh[i,,] <- tmp
  }

  X_ijh
}



#
# VB script
#
mc_dina = function(
    X,
    Q,
    max_it  = 500,
    epsilon = 1e-04,
    verbose = TRUE,
    # hyperparameter
    delta_0 = NULL,
    a_0 = NULL
){


  if(!inherits(X, "matrix")){
      X <- as.matrix(X)
    }
    if(!inherits(Q, "matrix")){
      Q <- as.matrix(Q)
    }

  # Index
  I <- nrow(X)
  J <- ncol(X)
  K <- ncol(Q) - 2
  L <- 2^K

  #
  # G_matrix
  #
  tmp <- make_G_mat(Q)
  A <- tmp$A
  H <- tmp$H
  G_mat <- tmp$G_mat

  #
  # Data comvert
  #
  X_ijh <- extend_X(X)

  #
  # Hyper parameters
  #
  if(is.null(delta_0)){
    delta_0 = rep(1, L)
  }
  if(is.null(a_0)){
    a_0 = lapply(1:J, function(j) matrix(1, nrow = max(H[[j]]), ncol=nrow(G_mat[[j]]) ))
  }

  #
  # Objects for estimates
  #
  delta_ast  <- as.numeric(table(sample(1:L,replace = T,size = I)))
  a_ast      <- lapply(1:J, function(j) matrix(0, nrow = max(H[[j]]), ncol=nrow(G_mat[[j]]) ))
  for(j in 1:J){
    n_h <- ncol(a_ast[[j]])
    for(h in 1:n_h){
      if(h != n_h){
        prob <- rep(1,n_h)
        prob[h] <- 2*prob[h]
        prob <- prob/sum(prob)
        a_ast[[j]][,h] <- table(sample(1:n_h,replace = T,size = I, prob = prob))

      } else{
        a_ast[[j]][,h] <- table(sample(1:n_h,replace = T,size = I))
      }
    }

  }

  # r_il <- matrix(runif(I*L) ,ncol=L, nrow = I)
  # r_il <- diag(1/rowSums(r_il)) %*% r_il
  r_il <- matrix(1/L ,ncol=L, nrow = I)
  one_vec  = matrix(1,nrow=I,ncol=1)


  m = 1
  l_lb = rep(0, max_it+1)
  l_lb[1] = -Inf
  for(m in 1:max_it){

    #
    # Expectations
    #
    E_log_pi  = digamma(delta_ast) - digamma(sum(delta_ast))

    E_log_theta_jhhd = lapply(1:J, function(j) digamma(a_ast[[j]]) - matrix(1,ncol=1, nrow=nrow(a_ast[[j]])) %*% digamma(colSums(a_ast[[j]])))

    #
    # E-step: r_il
    #

    temp <- matrix(0, ncol=L, nrow = I)
    for(j in 1:J) temp <- temp +  X_ijh[,j,H[[j]]] %*% E_log_theta_jhhd[[j]] %*% G_mat[[j]]
    log_rho_il <- t(t(temp) + E_log_pi)
    temp <- exp(log_rho_il)
    r_il =  temp / rowSums(temp)


    #
    # M-step
    #
    delta_ast <- colSums(r_il) + delta_0
    a_ast <- lapply(1:J, function(j) t(X_ijh[,j,H[[j]]]) %*% t(G_mat[[j]] %*% t(r_il))+ a_0[[j]] )


    temp_order_l <- sample(1:L,L)
    for (l in temp_order_l) {

      f = function(x) {
        lgamma(delta_ast[l] - delta_0[l] + x) - lgamma(sum(delta_ast) - delta_0[l]+ x) - lgamma(x) + digamma(sum(delta_ast) - delta_0[l] + x)
      }
      temp_max <- stats::optimise(f, maximum = T,lower = 10E-5,upper = 10E5)$maximum
      delta_ast[l] <- delta_ast[l] - delta_0[l] + temp_max
      delta_0[l] <- temp_max
    }

    temp_order_j <- sample(1:J,J)
    for(j in temp_order_j){

      for(hd in H[[j]]){
        a_0_temp <- a_0[[j]][,hd]
        a_ast_temp <- a_ast[[j]][,hd]

        temp_order_h <- sample(H[[j]], max(H[[j]]))
        for(h in temp_order_h){
          f = function(x) {
            lgamma(a_ast_temp[h] - a_0_temp[h] + x) - lgamma(sum(a_ast_temp) - a_0_temp[h]+ x) - lgamma(x) + digamma(sum(a_ast_temp) - a_0_temp[h] + x)
          }
          temp_max <- stats::optimise(f, maximum = T,lower = 10E-5,upper = 10E5)$maximum
          a_ast_temp[h] <- a_ast_temp[h] - a_0_temp[h] + temp_max
          a_0_temp[h] <- temp_max

        }

        a_ast[[j]][,hd] <- a_ast_temp
        a_0[[j]][,hd] <- a_0_temp


      }


    }


    #
    # Evidense of Lower Bound
    #

    tmp1 <- (sum(lgamma(delta_ast)) - lgamma(sum(delta_ast)))  - (sum(lgamma(delta_0)) + lgamma(sum(delta_0)))
    tmp2 <- sum(unlist(lapply(1:J, function(j) (colSums(lgamma(a_ast[[j]])) - lgamma(colSums(a_ast[[j]]))) - (colSums(lgamma(a_0[[j]])) + lgamma(colSums(a_0[[j]]))) ) ) )
    tmp3 <- -sum(r_il*log(r_il+10E-100))

    l_lb[m+1] <- tmp1 + tmp2 + tmp3

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

  a_ast_sum <- lapply(1:J, function(j) colSums(a_ast[[j]]))
  theta_est <- lapply(1:J, function(j) a_ast[[j]] %*% diag(1/a_ast_sum[[j]]))
  theta_sd <- lapply(1:J, function(j) sqrt(a_ast[[j]]*( matrix(rep(1,max(H[[j]]),ncol=1 )) %*% a_ast_sum[[j]] - a_ast[[j]]) %*% diag(1/(a_ast_sum[[j]]^2*(a_ast_sum[[j]]+1)) )))

  model_params = list(
    theta_est = theta_est,
    theta_sd = theta_sd
  )

  res = list(model_params = model_params,
       pi_est = pi_est,
       pi_sd = pi_sd,
       #r_il  = r_il,
       a_ast = a_ast,
       delta_ast   = delta_ast,
       a_0 = a_0,
       delta_0 = delta_0,
       l_lb = l_lb[l_lb != 0],
       att_pat_est = A[apply(r_il, 1, which.max),],
       G_mat = G_mat,
       m = m)

  return(res)
}



