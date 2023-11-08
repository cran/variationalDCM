#' the artificial data generation for the hidden-Markov DCM based on the given Q-matrix
#' \code{hmdcm_data_gen()} returns the artificially generated item response data for the HM-DCM
#' @param Q the J by K binary matrix
#' @param I the number of assumed respondents
#' @param min_theta the minimum value of the item parameter theta
#' @param max_theta the maximum value of the item parameter theta
#' @param att_cor the true value of the correlation among attributes (default: 0.1)
#' @param seed the seed value used for random number generation (default: 17)
#' @return A list including:
#' \describe{
#'   \item{X}{the generated artificial item response data}
#'   \item{alpha_true}{the generated true vale of the attribute mastery pattern, matrix form}
#'   \item{alpha_patt_true}{the generated true vale of the attribute mastery pattern, string form}
#' }
#' @references Yamaguchi, K., & Martinez, A. J. (2023). Variational Bayes
#' inference for hidden Markov diagnostic classification models. \emph{British Journal
#' of Mathematical and Statistical Psychology}, 00, 1– 25. \doi{10.1111/bmsp.12308}
#'
#' @examples
#' indT = 3
#' Q = sim_Q_J30K3
#' hm_sim_Q = lapply(1:indT,function(time_point) Q)
#' hm_sim_data = hmdcm_data_gen(Q=hm_sim_Q,I=200)
#'
#' @export

hmdcm_data_gen <- function(I = 500,
                           Q,
                           min_theta = 0.2,
                           max_theta = 0.8,
                           att_cor = 0.7,
                           seed = 17
){
  indI = I
  indK <- ncol(Q[[1]])
  indT <- length(Q)
  indJt <- sapply(Q,nrow)
  indL <- 2^indK
  #  cut_offs <- stats::qnorm(c(1:indK)/(indK+1))
  cut_offs <- stats::qnorm((1:indK+0.5)/(indK+1))
  item_par = "random"

  mean_norm <- rep(0,indK)
  vcov_norm <- matrix(att_cor,indK,indK)
  diag(vcov_norm) <- 1

  alpha_cont_t_1 <- mvtnorm::rmvnorm(n = indI, mean = mean_norm,sigma = vcov_norm)

  alpha_true_t_1 <- t((t(alpha_cont_t_1) > cut_offs)*1)

  #
  # All attribute pattern matrix
  #
  K_jt <- sapply(Q,rowSums)
  H_jt <- 2^K_jt
  not_zero_q_t <- lapply(Q,function(y)apply(y, 1, function(x) which(x != 0)))
  A <- as.matrix(expand.grid(lapply(1:indK, function(x)rep(0:1))))
  A_jt <- lapply(not_zero_q_t, function(y)lapply(y, function(x) A[,x,drop=F]))

  A_red <- lapply(  A_jt, function(z)lapply(z , function(x) apply(x,1,function(y) paste0(y,collapse = ""))))

  #
  # Unique correct item response probability label for each time point and item.
  #
  A_red_uni <- lapply(  A_jt,function(z)lapply(z , function(x) unique(apply(x,1,function(y) paste0(y,collapse = "")))))

  #
  # Make G-matrix
  #
  G_jt <- lapply(1:indT,function(time_t)lapply(1:indJt[[time_t]], function(j) outer(A_red_uni[[time_t]][[j]], A_red[[time_t]][[j]], function(x,y) (x == y)*1  )))

  att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
  for(time_t in 1:indT){
    for(j in 1:indJt[time_t]) {
      colnames(G_jt[[time_t]][[j]]) <- att_pat
      row.names(G_jt[[time_t]][[j]]) <- A_red_uni[[time_t]][[j]]
    }
  }

  #
  # Transition matrix: Tau
  #
  Tau_true <- matrix(0, indL, indL)
  colnames(Tau_true) <- att_pat
  row.names(Tau_true) <- att_pat

  for(l in 1:indL){
    temp <- rep(0, indL)
    for(ld in 1:indL){
      temp <- A[ld,] - A[l,]
      temp[temp == 1]  <- 3/10
      temp[temp == -1] <- 1/10
      temp[temp == 0] <- 6/10
      Tau_true[l, ld] <- prod(temp)

    }

  }

  Tau_true <- Tau_true / rowSums(Tau_true)
  # Tau_true

  #
  # Alpha true, z_i
  #
  z_i_true <- vector("list",indT)
  alpha_patt_true <- vector("list",indT)
  alpha_true <- vector("list",indT)

  alpha_true[[1]] <- alpha_true_t_1
  alpha_patt_true[[1]] <- apply(alpha_true[[1]], 1, function(x) paste0(x,collapse = ""))
  z_i_true[[1]] <- outer(alpha_patt_true[[1]], att_pat, function(x,y) (x == y)*1)
  colnames(z_i_true[[1]]) <- att_pat

  for(time_t in 1:(indT-1)){
    tansition_prob <- z_i_true[[time_t]] %*% Tau_true
    z_i_true[[time_t+1]] <- t(apply(tansition_prob, 1, function(p)stats::rmultinom(1,1,p) ))
    alpha_true[[time_t+1]] <- z_i_true[[time_t+1]] %*% A
    alpha_patt_true[[time_t+1]] <- apply(alpha_true[[time_t+1]], 1, function(x) paste0(x,collapse = ""))
  }

  n_of_att <- lapply(A_red_uni, function(y)lapply(y,function(x) sapply(strsplit(x, ""), function(y)sum(y == "1") )))

  #
  # π parameter
  #
  cut_off_k_dim_list <- lapply(1:indL, function(l){
    cut_off_k_dim <- matrix(NA, ncol=2, nrow = indK)
    cut_off_k_dim[A[l,] == 1,1] <- cut_offs[A[l,] == 1]
    cut_off_k_dim[A[l,] == 1,2] <- Inf
    cut_off_k_dim[A[l,] == 0,1] <- -Inf
    cut_off_k_dim[A[l,] == 0,2] <- cut_offs[A[l,] == 0]
    cut_off_k_dim
  })

  names(cut_off_k_dim_list) <- att_pat

  pi_true <- sapply(cut_off_k_dim_list, function(x)mvtnorm::pmvnorm(lower = x[,1],upper = x[,2],mean = mean_norm,sigma = vcov_norm))
  pi_true <- pi_true/sum(pi_true)
  # pi_true

  #
  # True theta_jh.
  #
  theta_jht_true <- lapply(indJt, function(Jt)vector("list",Jt))

  for(time_t in 1:indT){
    for(j in 1:indJt[time_t]){

      if(item_par == "random"){
        true_par_temp <- stats::runif(min = min_theta, max = max_theta, max(n_of_att[[time_t]][[j]])+1)
        true_par_temp <- true_par_temp[order(true_par_temp)]
        theta_jht_true[[time_t]][[j]] <-  n_of_att[[time_t]][[j]]

      }else{
        true_par_temp <- seq(from = min_theta, to = max_theta, length.out = max(n_of_att[[time_t]][[j]])+1)
      }

      for(k in min(n_of_att[[time_t]][[j]]):max( n_of_att[[time_t]][[j]])){
        theta_jht_true[[time_t]][[j]][n_of_att[[time_t]][[j]] == k] <- true_par_temp[k+1]
      }
      names(theta_jht_true[[time_t]][[j]]) <- A_red_uni[[time_t]][[j]]
    }

  }
  # theta_jht_true

  #
  # Generate item response matrix
  #
  X <- rand_mat <- item_res_prob <- vector("list", indT)

  for(time_t in 1:indT){
    item_res_prob[[time_t]] <- sapply(1:indJt[[time_t]], function(j){z_i_ast <- z_i_true[[time_t]] %*% t(G_jt[[time_t]][[j]])
    z_i_ast %*%  theta_jht_true[[time_t]][[j]]
    })

    rand_mat[[time_t]] <- matrix(stats::runif(indI*indJt[[time_t]]),ncol=indJt[[time_t]],nrow=indI)

    X[[time_t]] <- (item_res_prob[[time_t]] > rand_mat[[time_t]])*1

  }

  list(X = X,
       alpha_patt_true = alpha_patt_true,
       alpha_true = alpha_true)
}

hmdcm_vb = function(
    X,
    Q,
    model = "General",
    A_0 = NULL,
    B_0 = NULL,
    delta_0 = NULL,
    ommega_0 = NULL,
    max_it  = 500,
    epsilon = 1e-04,
    random_start = FALSE,
    verbose = TRUE,
    random_block_design=FALSE,
    Test_versions=NULL,
    Test_order=NULL
){

  if((random_block_design)){
    stop("sorry, current version does not soppot the case where random_block_design is true when nondecreesing_atribute is false.\n")
  }

  indI <- sapply(X, nrow)[1] # Assume all individuals take all tests
  indK <- ncol(Q[[1]])
  indT <- length(Q)
  indJt <- sapply(Q,nrow)
  indL <- 2^indK

  #
  # All attribute pattern matrix
  #

  not_zero_q_t <- lapply(Q,function(y)apply(y, 1, function(x) which(x != 0)))
  A <- as.matrix(expand.grid(lapply(1:indK, function(x)rep(0:1))))
  A_jt <- lapply(not_zero_q_t, function(y)lapply(y, function(x) A[,x,drop=F]))

  A_red <- lapply(  A_jt, function(z)lapply(z , function(x) apply(x,1,function(y) paste0(y,collapse = ""))))

  #
  # Unique correct item response probability label for each time point and item.
  #
  A_red_uni <- lapply(  A_jt,function(z)lapply(z , function(x) unique(apply(x,1,function(y) paste0(y,collapse = "")))))

  #
  # Make G-matrix
  #
  if(model == "General"){

    G_jt <- lapply(1:indT,function(time_t)lapply(1:indJt[[time_t]], function(j) outer(A_red_uni[[time_t]][[j]], A_red[[time_t]][[j]], function(x,y) (x == y)*1  )))

    att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
    for(time_t in 1:indT){
      for(j in 1:indJt[time_t]) {
        colnames(G_jt[[time_t]][[j]]) <- att_pat
        row.names(G_jt[[time_t]][[j]]) <- A_red_uni[[time_t]][[j]]
      }
    }

  }else if(model == "DINA"){

    G_jt <- lapply(1:indT,function(time_t)lapply(1:indJt[time_t], function(j)matrix(0,ncol=indL,nrow=2)))

    att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
    for(time_t in 1:indT){
      for(j in 1:indJt[time_t]) {
        temp_eta <- apply(t(t(A) ^ Q[[time_t]][j,]),1, prod)
        G_jt[[time_t]][[j]][1,] <- 1 - temp_eta
        G_jt[[time_t]][[j]][2,] <- temp_eta

        colnames(G_jt[[time_t]][[j]]) <- att_pat
        row.names(G_jt[[time_t]][[j]]) <- c("0","1")
      }
    }

  } else {

    stop("Error: Specify model General or DINA.\n")

  }


  #
  # Hyper parameter
  #
  if(is.null(delta_0) ){
    delta_0 = rep(1,indL)# For π

  }

  if(is.null(ommega_0) ){
    ommega_0 = matrix(1,indL,indL)# For Tau matrix

  }

  #
  # Weekly Monotonicity constraint
  #
  number_of_attributes <- lapply(A_red_uni,function(y)lapply(y, function(x) sapply(strsplit(x, ""), function(y)sum(as.numeric(y))) ) )

  if(model == "DINA") {number_of_attributes <- lapply(1:indT,function(time_t){lapply(1:indJt[time_t],function(j)c(0,1))})}

  if(is.null(A_0)){
    A_0_hyperparam <- lapply(number_of_attributes, function(x)seq(from = 1+epsilon, to = 2, length.out = max(unlist( x))+1) )
    A_0 <- vector("list", length = indT)
    for(time_t in 1:indT){
      A_0[[time_t]] <-lapply( number_of_attributes[[time_t]] , function(x){A_0_hyperparam[[time_t]][x + 1] })
    }
  }

  if(is.null(B_0)){
    B_0_hyperparam <- lapply(number_of_attributes, function(x)seq(from = 2, to = 1+epsilon, length.out = max(unlist( x))+1) )
    B_0 <- vector("list", length = indT)
    for(time_t in 1:indT){
      B_0[[time_t]] <-lapply( number_of_attributes[[time_t]] , function(x){B_0_hyperparam[[time_t]][x + 1] })
    }
  }

  #
  # Initialization
  #
  if(random_start == TRUE){
    E_z_itl_temp <- lapply(1:indT, function(time_t)matrix(stats::runif(indI*indL) ,ncol=indL, nrow = indI))
    E_z_itl_temp <- lapply(E_z_itl_temp, function(x)diag(1/rowSums(x)) %*% x)

    E_z_itl <-  array(0, dim=c(indI, indL, indT))
    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- E_z_itl_temp[[time_t]]
    }


    E_z_itl_z_itm1l <- array(0, dim=c(indI, indL, indL, indT-1))

    for(i in 1:indI){
      for(time_t in 1:(indT-1)){

        E_z_itl_z_itm1l_temp <- matrix(stats::runif(indL*indL) ,ncol=indL, nrow = indL)
        E_z_itl_z_itm1l_temp <- E_z_itl_z_itm1l_temp/ sum(E_z_itl_z_itm1l_temp)

        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l_temp
      }
    }

  }else{

    E_z_itl <-  array(0, dim=c(indI, indL, indT))
    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- matrix(1/indL, nrow=indI, ncol=indL)
    }

    E_z_itl_z_itm1l <- array(0, dim=c(indI, indL, indL, indT-1))

    for(i in 1:indI){
      for(time_t in 1:(indT-1)){

        E_z_itl_z_itm1l[i,,,time_t] <- 1/(indL*indL)
      }
    }

  }

  #
  # Evidence of Lower Bound
  #


  llb_fun <- function(delta_ast,delta_0,ommega_ast,ommega_0, A_ast, A_0, B_ast, B_0, log_zeta_sum){

    A_0_unlist <- unlist(A_0)
    B_0_unlist <- unlist(B_0)
    A_ast_unlist <- unlist(A_ast)
    B_ast_unlist <- unlist(B_ast)

    tmp1 <- sum( lbeta(A_ast_unlist, B_ast_unlist) - lbeta(A_0_unlist, B_0_unlist) + (A_0_unlist - A_ast_unlist)*(digamma(A_ast_unlist) - digamma(A_ast_unlist+B_ast_unlist)) + (B_0_unlist - B_ast_unlist)*( digamma(B_ast_unlist)-digamma(A_ast_unlist+B_ast_unlist)) )

    tmp2 <- (sum(lgamma(delta_ast)) - lgamma(sum(delta_ast)))  - (sum(lgamma(delta_0)) - lgamma(sum(delta_0))) + sum((delta_0 - delta_ast)*(digamma(delta_ast) - digamma(sum(delta_ast))) )

    tmp3 <- 0
    for(l in 1:indL){
      tmp3 <- tmp3 + (sum(lgamma(ommega_ast[l,])) - lgamma(sum(ommega_ast[l,])))  - (sum(lgamma(ommega_0[l,])) - lgamma(sum(ommega_0[l,]))) + sum((ommega_0[l,] - ommega_ast[l,])*(digamma(ommega_ast[l,]) - digamma(sum(ommega_ast[l,]))) )
    }


    tmp1 + tmp2 + tmp3 + log_zeta_sum
  }


  #
  # Make objects for variational parameters
  #
  E_log_theta <- E_log_1_theta <- B_ast <- A_ast <- A_0
  delta_ast <- delta_0
  ommega_ast <- ommega_0


  b_z_it <- f_z_it <- array(0, dim=c(indI, indL ,indT) )
  gamma_t_x_it <- matrix(0, nrow=indI,ncol=indT)
  P_til_x_it_z_it <- array(0, dim=c(indI,indL,indT))


  m = 1

  l_lb = rep(NA_real_, max_it+1)
  l_lb[1] = 100

  for(m in 1:max_it){
    #
    # M-step and Calculation of Expectations
    #
    delta_ast <- colSums(E_z_itl[,,1]) + delta_0
    ommega_ast <- apply(E_z_itl_z_itm1l, c(2,3),sum) + ommega_0 #Check this point

    E_log_pi      = digamma(delta_ast) - digamma(sum(delta_ast))
    E_log_tau     = digamma(ommega_ast) - digamma(rowSums(ommega_ast))
    for(time_t in 1:indT){

      A_ast[[time_t]]  <- lapply(1:indJt[[time_t]], function(j) t(G_jt[[time_t]][[j]] %*% (t(E_z_itl[,,time_t]) %*% X[[time_t]][,j])     + A_0[[time_t]][[j]] ))
      B_ast[[time_t]]  <- lapply(1:indJt[[time_t]], function(j) t(G_jt[[time_t]][[j]] %*% (t(E_z_itl[,,time_t]) %*% (1-X[[time_t]][,j])) + B_0[[time_t]][[j]] ))

      E_log_theta[[time_t]]   = lapply(1:indJt[[time_t]], function(j) digamma(A_ast[[time_t]][[j]]) - digamma(A_ast[[time_t]][[j]] + B_ast[[time_t]][[j]]))
      E_log_1_theta[[time_t]] = lapply(1:indJt[[time_t]], function(j) digamma(B_ast[[time_t]][[j]]) - digamma(A_ast[[time_t]][[j]] + B_ast[[time_t]][[j]]))
    }

    #
    # E-step
    #
    # f_z_it
    # b_z_it

    for(time_t in 1:indT){
      temp <- matrix(0, nrow = indI, ncol=indL)
      for(j in 1:indJt[time_t]){
        temp <- temp +  ( X[[time_t]][,j]%*% E_log_theta[[time_t]][[j]] + (1-X[[time_t]][,j]) %*% E_log_1_theta[[time_t]][[j]]) %*% G_jt[[time_t]][[j]]
      }
      P_til_x_it_z_it[,,time_t] <-  exp(temp)
    }


    f_z_it[,,1] <- exp(t(log(t(P_til_x_it_z_it[,,1])) + E_log_pi))
    gamma_t_x_it[,1] <- rowSums(f_z_it[,,1])
    f_z_it[,,1] <- f_z_it[,,1]/gamma_t_x_it[,1] # Normarize

    b_z_it[,,indT] <- 1

    #
    # Recursive calculation
    #
    for(time_t in 2:indT){
      #
      # calc f
      #

      f_z_it[,,time_t] <-  P_til_x_it_z_it[,,time_t] * (f_z_it[,,time_t-1] %*% exp(E_log_tau))
      gamma_t_x_it[,time_t] <- rowSums(f_z_it[,,time_t])
      f_z_it[,,time_t] <- f_z_it[,,time_t]/gamma_t_x_it[,time_t] # Normarize

      #
      # calc b
      #
      b_z_it[,,indT - time_t + 1] <-  (P_til_x_it_z_it[,,indT - time_t + 2] * b_z_it[,,indT - time_t + 2]) %*% t(exp(E_log_tau))

      b_z_it[,,indT - time_t + 1] <- b_z_it[,,indT - time_t + 1] / rowSums(b_z_it[,,indT - time_t + 1])

    }



    E_z_itl <- f_z_it*b_z_it
    E_z_itl_temp <-  apply(E_z_itl, c(1,3),sum)

    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- E_z_itl[,,time_t]/E_z_itl_temp[,time_t]
    }

    for(l in 1:indL){
      for(time_t in 2:indT){

        E_z_itl_z_itm1l[,l,,time_t-1] <- t(t(P_til_x_it_z_it[,,time_t]*b_z_it[,,time_t]*f_z_it[,l,time_t-1]) * exp(E_log_tau[l,]))

      }

    }


    E_z_itl_z_itm1l_temp <- apply(E_z_itl_z_itm1l, c(1,4),sum)


    for(i in 1:indI){
      for(time_t in 1:(indT-1)){
        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l[i,,,time_t]/E_z_itl_z_itm1l_temp[i,time_t]
      }
    }

    log_zeta_sum <- sum(log(gamma_t_x_it))

    l_lb[m+1] <- llb_fun(delta_ast,delta_0,ommega_ast,ommega_0, A_ast, A_0, B_ast, B_0, log_zeta_sum)

    if(verbose){
      cat("\riteration = ", m+1, sprintf(",last change = %.05f", abs(l_lb[m] - l_lb[m+1])))
    }

    if(abs(l_lb[m] - l_lb[m+1]) < epsilon){
      break()
    }

  }
  l_lb <- l_lb[-1]
  #  plot(l_lb,type="l")

  #
  # Calculation of mean and sd of VB posteriors
  #
  delta_sum <- sum(delta_ast)
  pi_est <-  delta_ast/delta_sum
  pi_sd <-sqrt(delta_ast*(delta_sum - delta_ast)/(delta_sum^2*(delta_sum+1)) )

  names(pi_est) <- att_pat
  names(pi_sd)  <- att_pat

  ommega_sum <- rowSums(ommega_ast)
  Tau_est <-  ommega_ast/ommega_sum
  Tau_sd <- matrix(0, indL, indL)
  for(l in 1:indL) Tau_sd[,l] <- sqrt(ommega_ast[,l]*(ommega_sum - ommega_ast[,l])/(ommega_sum^2*(ommega_sum+1)) )

  colnames(Tau_est) <- att_pat
  row.names(Tau_est) <- att_pat
  colnames(Tau_sd) <- att_pat
  row.names(Tau_sd) <- att_pat

  theta_sd <- theta_est <- vector("list",indT)
  for(time_t in 1:indT){
    theta_est[[time_t]] <- mapply(function(x,y) x/(x+y), A_ast[[time_t]], B_ast[[time_t]])
    theta_sd[[time_t]] <- mapply(function(x,y) sqrt((x*y)/(((x+y)^2) *(x+y+1)) ),  A_ast[[time_t]], B_ast[[time_t]])

  }

  #
  # MAP and EAP of attribute mastery.
  #
  post_max_class <- matrix(0, nrow=indI, ncol=indT)
  EAP_att_pat  <- att_master_prob  <- MAP_att_pat  <- lapply(1:indT, function(time_t) matrix(0, nrow=indI, ncol=indK))
  for(time_t in 1:indT){
    post_max_class[,time_t] <- apply(E_z_itl[,,time_t], 1, function(x)which.max(x) )
    MAP_att_pat[[time_t]] <- A[post_max_class[,time_t],]

    att_master_prob[[time_t]] <- E_z_itl[,,time_t] %*% A
    EAP_att_pat[[time_t]] <- (att_master_prob[[time_t]] > 0.5)*1

  }




  list(theta_est = theta_est,
       theta_sd = theta_sd,
       pi_est = pi_est,
       pi_sd = pi_sd,
       Tau_est = Tau_est,
       Tau_sd = Tau_sd,
       post_max_class = post_max_class,
       MAP_att_pat = MAP_att_pat,
       att_master_prob = att_master_prob,
       EAP_att_pat = EAP_att_pat,
       A_ast = A_ast,
       B_ast = B_ast,
       delta_ast   = delta_ast,
       ommega_ast   = ommega_ast,
       E_z_itl = E_z_itl,
       E_z_itl_z_itm1l = E_z_itl_z_itm1l,
       A_0 = A_0,
       B_0 = B_0,
       delta_0 = delta_0,
       ommega_0 = ommega_0,
       l_lb = l_lb,
       gamma_t_x_it = gamma_t_x_it,
       log_zeta_sum = log_zeta_sum,
       E_z_itl = E_z_itl,
       E_z_itl_z_itm1l = E_z_itl_z_itm1l,
       A = A,
       Q = Q,
       X = X,
       G_jt = G_jt,
       m = m)
}

hmdcm_vb_nondec= function(
    X,
    Q,
    A_0 = NULL,
    B_0 = NULL,
    delta_0 = NULL,
    ommega_0 = NULL,
    max_it  = 500,
    epsilon = 10E-4,
    random_block_design=FALSE,
    Test_versions=NULL,
    Test_order=NULL,
    model="General",
    random_start = FALSE,
    verbose=TRUE
){

  if(!(random_block_design)){
    stop("sorry, current version does not soppot the case where random_block_design is false when nondecreesing_atribute is true.\n")
  }

  indI <- sapply(X, nrow)[1] # Assume all individuals take all tests.
  indK <- ncol(Q[[1]]) # Assume attributes are all same across time and individual.
  indT <- length(Q) # Assume time points is all same across individual.
  indJt <- sapply(Q,nrow) # Assume items presented at a time is just same across time.
  indL <- 2^indK

  #
  # All attribute pattern matrix
  #
  not_zero_q_t <- lapply(Q,function(y)apply(y, 1, function(x) which(x != 0)))
  A <- as.matrix(expand.grid(lapply(1:indK, function(x)rep(0:1))))
  A_jt <- lapply(not_zero_q_t, function(y)lapply(y, function(x) A[,x,drop=F]))

  A_red <- lapply(  A_jt, function(z)lapply(z , function(x) apply(x,1,function(y) paste0(y,collapse = ""))))

  #
  # Unique correct item response probability label for each time point and item.
  #
  A_red_uni <- lapply(  A_jt,function(z)lapply(z , function(x) unique(apply(x,1,function(y) paste0(y,collapse = "")))))

  #
  # Make G-matrix
  #
  # G_j <- lapply(1:J, function(j) t(sapply(A_red_uni[[j]], function(x) x == A_red[[j]] ))*1 )

  if(model == "General"){

    G_jt <- lapply(1:indT,function(time_t)lapply(1:indJt[[time_t]], function(j) outer(A_red_uni[[time_t]][[j]], A_red[[time_t]][[j]], function(x,y) (x == y)*1  )))


    att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
    for(time_t in 1:indT){
      for(j in 1:indJt[time_t]) {
        colnames(G_jt[[time_t]][[j]]) <- att_pat
        row.names(G_jt[[time_t]][[j]]) <- A_red_uni[[time_t]][[j]]
      }
    }

  }else if(model == "DINA"){

    G_jt <- lapply(1:indT,function(time_t)lapply(1:indJt[time_t], function(j)matrix(0,ncol=indL,nrow=2)))

    att_pat <- apply(A, 1, function(x) paste0(x, collapse=""))
    for(time_t in 1:indT){
      for(j in 1:indJt[time_t]) {
        temp_eta <- apply(t(t(A) ^ Q[[time_t]][j,]),1, prod)
        G_jt[[time_t]][[j]][1,] <- 1 - temp_eta
        G_jt[[time_t]][[j]][2,] <- temp_eta

        colnames(G_jt[[time_t]][[j]]) <- att_pat
        row.names(G_jt[[time_t]][[j]]) <- c("0","1")
      }
    }


  } else {
    stop("Error: Specify model General or DINA.\n")
  }


  #
  # Hyper parameter
  #
  if(is.null(delta_0) ){
    delta_0 = rep(1,indL)# For π
    #delta_0 = rep(1/indL,indL)# For π

  }

  if(is.null(ommega_0) ){
    ommega_0 = matrix(1,indL,indL)# For Tau matrix
    #ommega_0 = matrix(1/indL,indL,indL)# For Tau matrix

    for(l in 1:indL){
      for(ld in 1:indL){
        dif_pat <- A[l,] - A[ld,]
        ommega_0[l,ld] <- ifelse(any(dif_pat > 0), 0, 1)
      }
    }


  }

  #
  # Weekly Monotonicity constraint
  #
  number_of_attributes <- lapply(A_red_uni,function(y)lapply(y, function(x) sapply(strsplit(x, ""), function(y)sum(as.numeric(y))) ) )

  if(model == "DINA") {number_of_attributes <- lapply(1:indT,function(time_t){lapply(1:indJt[time_t],function(j)c(0,1))})}

  if(is.null(A_0)){

    A_0_hyperparam <- lapply(number_of_attributes, function(x)seq(from = 1+epsilon, to = 2, length.out = max(unlist( x))+1) )
    A_0 <- vector("list", length = indT)
    for(time_t in 1:indT){
      A_0[[time_t]] <-lapply( number_of_attributes[[time_t]] , function(x){A_0_hyperparam[[time_t]][x + 1] })
    }

  }

  if(is.null(B_0)){
    B_0_hyperparam <- lapply(number_of_attributes, function(x)seq(from = 2, to = 1+epsilon, length.out = max(unlist( x))+1) )
    B_0 <- vector("list", length = indT)
    for(time_t in 1:indT){
      B_0[[time_t]] <-lapply( number_of_attributes[[time_t]] , function(x){B_0_hyperparam[[time_t]][x + 1] })
    }

  }




  #
  # Initialization
  #
  if(random_start == TRUE){
    E_z_itl_temp <- lapply(1:indT, function(time_t)matrix(stats::runif(indI*indL) ,ncol=indL, nrow = indI))
    E_z_itl_temp <- lapply(E_z_itl_temp, function(x)diag(1/rowSums(x)) %*% x)

    E_z_itl <-  array(0, dim=c(indI, indL, indT))
    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- E_z_itl_temp[[time_t]]
    }


    E_z_itl_z_itm1l <- array(0, dim=c(indI, indL, indL, indT-1))

    for(i in 1:indI){
      for(time_t in 1:(indT-1)){

        E_z_itl_z_itm1l_temp <- matrix(stats::runif(indL*indL) ,ncol=indL, nrow = indL)
        E_z_itl_z_itm1l_temp[ommega_0==0] <-0
        E_z_itl_z_itm1l_temp <- E_z_itl_z_itm1l_temp/ sum(E_z_itl_z_itm1l_temp)

        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l_temp
      }
    }

  }else{

    E_z_itl <-  array(0, dim=c(indI, indL, indT))
    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- matrix(1/indL, nrow=indI, ncol=indL)

    }

    E_z_itl_z_itm1l <- array(0, dim=c(indI, indL, indL, indT-1))

    for(i in 1:indI){
      for(time_t in 1:(indT-1)){

        E_z_itl_z_itm1l_temp <- matrix(1/(indL*indL),indL,indL)
        E_z_itl_z_itm1l_temp[ommega_0==0] <-0
        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l_temp/ sum(E_z_itl_z_itm1l_temp)

      }
    }

  }

  #
  # Evidence of Lower Bound
  #
  llb_fun <- function(delta_ast,delta_0,ommega_ast,ommega_0, A_ast, A_0, B_ast, B_0, log_zeta_sum){

    A_0_unlist <- unlist(A_0)
    B_0_unlist <- unlist(B_0)
    A_ast_unlist <- unlist(A_ast)
    B_ast_unlist <- unlist(B_ast)

    tmp1 <- sum( lbeta(A_ast_unlist, B_ast_unlist) - lbeta(A_0_unlist, B_0_unlist) + (A_0_unlist - A_ast_unlist)*(digamma(A_ast_unlist) - digamma(A_ast_unlist+B_ast_unlist)) + (B_0_unlist - B_ast_unlist)*( digamma(B_ast_unlist)-digamma(A_ast_unlist+B_ast_unlist)) )

    tmp2 <- (sum(lgamma(delta_ast)) - lgamma(sum(delta_ast)))  - (sum(lgamma(delta_0)) - lgamma(sum(delta_0))) + sum((delta_0 - delta_ast)*(digamma(delta_ast) - digamma(sum(delta_ast))) )

    tmp3 <- 0
    for(l in 1:indL){
      ommega_not_0   <- ommega_0[l,]!=0
      tmp3 <- tmp3 + (sum(lgamma(ommega_ast[l,ommega_not_0])) - lgamma(sum(ommega_ast[l,ommega_not_0])))  - (sum(lgamma(ommega_0[l,ommega_not_0])) - lgamma(sum(ommega_0[l,ommega_not_0]))) + sum((ommega_0[l,ommega_not_0] - ommega_ast[l,ommega_not_0])*(digamma(ommega_ast[l,ommega_not_0]) - digamma(sum(ommega_ast[l,ommega_not_0]))) )
    }


    tmp1 + tmp2 + tmp3 + log_zeta_sum
  }


  #
  # Make objects for variational parameters
  #
  E_log_theta <- E_log_1_theta <- B_ast <- A_ast <- A_0
  delta_ast <- delta_0
  ommega_ast <- ommega_0
  ommega_zero_elem <- ommega_0 == 0


  b_z_it <- f_z_it <- array(0, dim=c(indI, indL ,indT) )
  gamma_t_x_it <- matrix(0, nrow=indI,ncol=indT)
  P_til_x_it_z_it <- array(0, dim=c(indI,indL,indT))



  X_reord <- X
  for(i in 1:indI){
    for(time_t in 1:indT){
      X_reord[[Test_order[Test_versions[i],time_t ]]][i,] <- X[[time_t]][i,]
    }
  }

  m = 1

  l_lb = rep(NA_real_, max_it+1)
  l_lb[1] = 100

  for(m in 1:max_it){
    #
    # M-step and Calculation of Expectations
    #
    delta_ast <- colSums(E_z_itl[,,1]) + delta_0
    ommega_ast <- apply(E_z_itl_z_itm1l, c(2,3),sum) + ommega_0 #Check this point


    E_log_pi      = digamma(delta_ast) - digamma(sum(delta_ast))
    #digamma_omega <- try(digamma(ommega_ast))
    #digamma_omega[ommega_zero_elem]  <- 0
    E_log_tau     = try(digamma(ommega_ast), silent = T)  - digamma(rowSums(ommega_ast))
    E_log_tau[ommega_zero_elem]  <- 0


    #
    # Reorder
    #
    E_z_itl_reord <- E_z_itl
    for(i in 1:indI){
      for(time_t in 1:indT){
        E_z_itl_reord[i,,Test_order[Test_versions[i],time_t]] <-  E_z_itl[i,,time_t]

      }
    }


    for(time_t in 1:indT){

      A_ast[[time_t]]  <- lapply(1:indJt[[time_t]], function(j) t(G_jt[[time_t]][[j]] %*% (t(E_z_itl_reord[,,time_t]) %*% X_reord[[time_t]][,j])     + A_0[[time_t]][[j]] ))
      B_ast[[time_t]]  <- lapply(1:indJt[[time_t]], function(j) t(G_jt[[time_t]][[j]] %*% (t(E_z_itl_reord[,,time_t]) %*% (1-X_reord[[time_t]][,j])) + B_0[[time_t]][[j]] ))

      E_log_theta[[time_t]]   = lapply(1:indJt[[time_t]], function(j) digamma(A_ast[[time_t]][[j]]) - digamma(A_ast[[time_t]][[j]] + B_ast[[time_t]][[j]]))
      E_log_1_theta[[time_t]] = lapply(1:indJt[[time_t]], function(j) digamma(B_ast[[time_t]][[j]]) - digamma(A_ast[[time_t]][[j]] + B_ast[[time_t]][[j]]))
    }

    #
    # E-step
    #
    for(time_t in 1:indT){
      temp <- matrix(0, nrow = indI, ncol=indL)
      for(i in 1:indI){

        for(j in 1:indJt[time_t]){
          #          temp <- temp +  ( X[[time_t]][,j]%*% E_log_theta[[time_t]][[j]] + (1-X[[time_t]][,j]) %*% E_log_1_theta[[time_t]][[j]]) %*% G_jt[[time_t]][[j]]
          temp[i,] <- temp[i,] +  ( X[[time_t]][i,j]* E_log_theta[[Test_order[Test_versions[i],time_t]]][[j]] + (1-X[[time_t]][i,j]) * E_log_1_theta[[Test_order[Test_versions[i],time_t]]][[j]]) %*% G_jt[[Test_order[Test_versions[i],time_t]]][[j]]
        }
      }
      P_til_x_it_z_it[,,time_t] <-  exp(temp)
    }


    f_z_it[,,1] <- exp(t(log(t(P_til_x_it_z_it[,,1])) + E_log_pi))
    gamma_t_x_it[,1] <- rowSums(f_z_it[,,1])
    f_z_it[,,1] <- f_z_it[,,1]/gamma_t_x_it[,1] # Normarize

    b_z_it[,,indT] <- 1

    #
    # Recursive calculation
    #

    exp_E_log_tau <- exp(E_log_tau)
    exp_E_log_tau[ommega_zero_elem] <- 0

    for(time_t in 2:indT){
      #
      # calc f
      #

      f_z_it[,,time_t] <-  P_til_x_it_z_it[,,time_t] * (f_z_it[,,time_t-1] %*% exp_E_log_tau)
      gamma_t_x_it[,time_t] <- rowSums(f_z_it[,,time_t])
      f_z_it[,,time_t] <- f_z_it[,,time_t]/gamma_t_x_it[,time_t] # Normarize

      #
      # calc b
      #
      b_z_it[,,indT - time_t + 1] <-  (P_til_x_it_z_it[,,indT - time_t + 2] * b_z_it[,,indT - time_t + 2]) %*% t(exp_E_log_tau)

      b_z_it[,,indT - time_t + 1] <- b_z_it[,,indT - time_t + 1] / rowSums(b_z_it[,,indT - time_t + 1])

    }



    E_z_itl <- f_z_it*b_z_it
    E_z_itl_temp <-  apply(E_z_itl, c(1,3),sum)

    for(time_t in 1:indT){
      E_z_itl[,,time_t] <- E_z_itl[,,time_t]/E_z_itl_temp[,time_t]
    }


    for(l in 1:indL){
      for(time_t in 2:indT){
        E_z_itl_z_itm1l[,l,,time_t-1] <- t(t(P_til_x_it_z_it[,,time_t]*b_z_it[,,time_t]*f_z_it[,l,time_t-1]) * exp_E_log_tau[l,])
      }
    }

    E_z_itl_z_itm1l_temp <- apply(E_z_itl_z_itm1l, c(1,4),sum)


    for(i in 1:indI){
      for(time_t in 1:(indT-1)){
        E_z_itl_z_itm1l[i,,,time_t] <- E_z_itl_z_itm1l[i,,,time_t]/E_z_itl_z_itm1l_temp[i,time_t]
      }
    }

    log_zeta_sum <- sum(log(gamma_t_x_it))

    l_lb[m+1] <- llb_fun(delta_ast,delta_0,ommega_ast,ommega_0, A_ast, A_0, B_ast, B_0, log_zeta_sum)

    if(verbose){
      cat("\riteration = ", m+1, sprintf(",last change = %.05f", abs(l_lb[m] - l_lb[m+1])))
    }

    if(abs(l_lb[m] - l_lb[m+1]) < epsilon){
      break()
    }

  }
  l_lb <- l_lb[-1]
  #  plot(l_lb,type="l")


  #
  # Calculation of mean and sd of VB posteriors
  #
  delta_sum <- sum(delta_ast)
  pi_est <-  delta_ast/delta_sum
  pi_sd <-sqrt(delta_ast*(delta_sum - delta_ast)/(delta_sum^2*(delta_sum+1)) )
  names(pi_est) <- att_pat
  names(pi_sd)  <- att_pat

  ommega_sum <- rowSums(ommega_ast)
  Tau_est <-  ommega_ast/ommega_sum
  Tau_sd <- matrix(0, indL, indL)
  for(l in 1:indL) Tau_sd[,l] <- sqrt(ommega_ast[,l]*(ommega_sum - ommega_ast[,l])/(ommega_sum^2*(ommega_sum+1)) )

  colnames(Tau_est) <- att_pat
  row.names(Tau_est) <- att_pat
  colnames(Tau_sd) <- att_pat
  row.names(Tau_sd) <- att_pat


  theta_sd <- theta_est <- vector("list",indT)
  for(time_t in 1:indT){
    theta_est[[time_t]] <- mapply(function(x,y) x/(x+y), A_ast[[time_t]], B_ast[[time_t]])
    theta_sd[[time_t]] <- mapply(function(x,y) sqrt((x*y)/(((x+y)^2) *(x+y+1)) ),  A_ast[[time_t]], B_ast[[time_t]])

  }

  #
  # MAP and EAP of attribute mastery.
  #
  post_max_class <- matrix(0, nrow=indI, ncol=indT)
  EAP_att_pat  <- att_master_prob  <- MAP_att_pat  <- lapply(1:indT, function(time_t) matrix(0, nrow=indI, ncol=indK))
  for(time_t in 1:indT){
    post_max_class[,time_t] <- apply(E_z_itl[,,time_t], 1, function(x)which.max(x) )
    MAP_att_pat[[time_t]] <- A[post_max_class[,time_t],]

    att_master_prob[[time_t]] <- E_z_itl[,,time_t] %*% A
    EAP_att_pat[[time_t]] <- (att_master_prob[[time_t]] > 0.5)*1

  }




  list(theta_est = theta_est,
       theta_sd = theta_sd,
       pi_est = pi_est,
       pi_sd = pi_sd,
       Tau_est = Tau_est,
       Tau_sd = Tau_sd,
       post_max_class = post_max_class,
       MAP_att_pat = MAP_att_pat,
       att_master_prob = att_master_prob,
       EAP_att_pat = EAP_att_pat,
       A_ast = A_ast,
       B_ast = B_ast,
       delta_ast   = delta_ast,
       ommega_ast   = ommega_ast,
       E_z_itl = E_z_itl,
       E_z_itl_z_itm1l = E_z_itl_z_itm1l,
       A_0 = A_0,
       B_0 = B_0,
       delta_0 = delta_0,
       ommega_0 = ommega_0,
       l_lb = l_lb,
       gamma_t_x_it = gamma_t_x_it,
       log_zeta_sum = log_zeta_sum,
       E_z_itl = E_z_itl,
       E_z_itl_z_itm1l = E_z_itl_z_itm1l,
       A = A,
       Q = Q,
       X = X,
       G_jt = G_jt,
       m = m)
}



#' for the hidden Markov DCM.
#'
#' \code{hm_dcm()} returns variational Bayesian estimates for the hidden
#' Markov DCM.
#'
#' @param X  T-length list or 3-dim array whose each element is I by J/T binary item response data matrix
#' @param Q  T-length list or 3-dim array whose each element is J/T by K Q-matrix
#' @param random_block_design logical value; whether the test design adopts different item ordering or not (default: FALSE)
#' @param A_0 the value of hyperparameter \eqn{A^0} (default: NULL)
#' @param B_0 the value of hyperparameter \eqn{B^0} (default: NULL)
#' @param delta_0 the value of hyperparameter delta_0 (default: NULL)
#' @param ommega_0 the value of hyperparameter ommega_0 (default: NULL)
#' @param Test_versions  indicates the test module each respondent is assigned (default: NULL)
#' @param Test_order the square matrix of order T that represents item order (default: NULL)
#' @param model "General" or "DINA" (default: "General"), specifies the measurement model
#' @param max_it Maximum number of iterations (default: 500)
#' @param epsilon convergence tolerance for iterations (default: 1e-4)
#' @param verbose logical, controls whether to print progress (default: TRUE)
#' @param random_start logical (default: FALSE)
#' @param nondecreasing_attribute logical; whether the assumption that mastered attributes are not forgotten is adopted or not (default: FALSE)
#'
#' @return A list including:
#' \describe{
#'   \item{theta_est}{the estimate of conditional response probability parameter \eqn{\Theta}}
#'   \item{theta_sd}{the posterior standard deviation of parameter \eqn{\Theta}}
#'   \item{pi_est}{the estimate of class mixing parameter \eqn{\pi}}
#'   \item{pi_sd}{the posterior standard deviation of parameter \eqn{\pi}}
#'   \item{Tau_est}{the estimate of class-transition probability parameter \eqn{\tau}}
#'   \item{Tau_sd}{the posterior standard deviation of parameter \eqn{\tau}}
#'   \item{post_max_class}{the result of class analysis}
#'   \item{MAP_att_pat}{the MAP estimate of attribute mastery patterns}
#'   \item{att_master_prob}{the estimated attribute mastery probabilities}
#'   \item{EAP_att_pat}{the EAP estimate of attribute mastery patterns}
#'   \item{A_ast}{the estimate of variational parameter \eqn{A^*}}
#'   \item{B_ast}{the estimate of variational parameter \eqn{B^*}}
#'   \item{delta_ast}{the estimate of variational posterior \eqn{\delta^*}}
#'   \item{ommega_ast}{the estimate of variational posterior \eqn{\Omega^*}}
#'   \item{E_z_itl}{the resulted expectations of  attribute mastery pattern}
#'   \item{E_z_itl_z_itm1l}{the resulted expectations of  attribute mastery pattern at the time point t and t-1}
#'   \item{A_0}{the value of hyperparameter \eqn{A^0}}
#'   \item{B_0}{the value of hyperparameter \eqn{B^0}}
#'   \item{delta_0}{the computed or entered value of delta_0}
#'   \item{ommega_0}{the computed or entered value of ommega_0}
#'   \item{l_lb}{the list of the values of evidence lower bound in each itertion}
#   \item{gamma_t_x_it}{the computed value of the normalizing constant for calculating zeta}
#   \item{log_zeta_sum}{the computed value of the normalizing constant for calculating variational posterior for class indicator}
#   \item{A}{all the possible attribute mastery patterns}
#   \item{Q}{the entered Q-matrix}
#   \item{X}{the entered data matrix}
#'   \item{G_jt}{the computed G-matrices}
#'   \item{m}{the number of performed iterations}
#' }
#'
#' @references Yamaguchi, K., & Martinez, A. J. (2023). Variational Bayes
#' inference for hidden Markov diagnostic classification models. \emph{British Journal
#' of Mathematical and Statistical Psychology}, 00, 1– 25. \doi{10.1111/bmsp.12308}
#'
#' @examples
#' indT = 3
#' Q = sim_Q_J30K3
#' hm_sim_Q = lapply(1:indT,function(time_point) Q)
#' hm_sim_data = hmdcm_data_gen(Q=hm_sim_Q,I=200)
#' res_hm = hm_dcm(X=hm_sim_data$X,Q=hm_sim_Q)
#'
#' @export

hm_dcm = function(
    X,
    Q,
    max_it  = 500,
    epsilon = 1e-04,
    nondecreasing_attribute=FALSE,
    model="General",
    verbose = TRUE,
    random_block_design=FALSE,
    Test_versions=NULL,
    Test_order=NULL,
    random_start = FALSE,
    # hyperparameters
    A_0 = NULL,
    B_0 = NULL,
    delta_0 = NULL,
    ommega_0 = NULL
){

  # convert X,Q to list
  if(length(dim(X)) == 3){
    X = lapply(1:dim(X)[3], function(i) X[,,i])
  }
  if(length(dim(Q)) == 3){
    Q = lapply(1:dim(Q)[3], function(i) Q[,,i])
  }

  if(random_block_design){
    if(is.null(Test_versions) || is.null(Test_order)){
      stop("if random_block_design is true, Test_versions and Test_order must be entered.\n")
    }
  }

  if(!nondecreasing_attribute){
    res = hmdcm_vb(X=X,Q=Q,max_it=max_it,epsilon=epsilon,model=model,
                   random_block_design=random_block_design,
                   Test_versions=Test_versions,Test_order=Test_order,
                   verbose=verbose,random_start=random_start,A_0=A_0,B_0=B_0,
                   delta_0=delta_0,ommega_0=ommega_0)
  }else{
    res = hmdcm_vb_nondec(X=X,Q=Q,max_it=max_it,epsilon=epsilon,model=model,
                          random_block_design=random_block_design,
                          Test_versions=Test_versions,Test_order=Test_order,
                          verbose=verbose,random_start=random_start,A_0=A_0,B_0=B_0,
                          delta_0=delta_0,ommega_0=ommega_0)
  }

  list(theta_est = res$theta_est,
       theta_sd = res$theta_sd,
       pi_est = res$pi_est,
       pi_sd = res$pi_sd,
       Tau_est = res$Tau_est,
       Tau_sd = res$Tau_sd,
       post_max_class = res$post_max_class,
       MAP_att_pat = res$MAP_att_pat,
       att_master_prob = res$att_master_prob,
       EAP_att_pat = res$EAP_att_pat,
       A_ast = res$A_ast,
       B_ast = res$B_ast,
       delta_ast = res$delta_ast,
       ommega_ast = res$ommega_ast,
       E_z_itl = res$E_z_itl,
       E_z_itl_z_itm1l = res$E_z_itl_z_itm1l,
       A_0 = res$A_0,
       B_0 = res$B_0,
       delta_0 = res$delta_0,
       ommega_0 = res$ommega_0,
       l_lb = res$l_lb,
       #gamma_t_x_it = res$gamma_t_x_it,
       #log_zeta_sum = res$log_zeta_sum,
       E_z_itl = res$E_z_itl,
       E_z_itl_z_itm1l = res$E_z_itl_z_itm1l,
       #A = res$A,
       #Q = res$Q,
       #X = res$X,
       G_jt = res$G_jt,
       m = res$m)
}
