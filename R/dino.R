#' for the deterministic input noisy OR gate (DINO) model.
#'
#' \code{dino()} returns variational Bayesian estimates for the DINO model.
#'
#' @param X I by J binary matrix, item response data
#' @param Q J by K binary matrix, Q-matrix
#' @param max_it the maximum number of iterations (default: 500)
#' @param epsilon the convergence tolerance for iterations (default: 1e-4)
#' @param verbose Logical, controls whether to print progress (default: TRUE)
#' @param delta_0 L by 1 vector, hyperparameter of prior dirichlet distribution
#'   for the class mixing parameter \eqn{\pi} (default: NULL).
#' @param alpha_s A positive scalar, hyperparameter that determines the
#'   shape of  prior beta distribution for slip parameter (default: NULL).
#' @param beta_s A positive scalar, hyperparameter that determines the
#'   shape of  prior beta distribution for slip parameter (default: NULL).
#' @param alpha_g A positive scalar, hyperparameter that determines the
#'   shape of  prior beta distribution for guessing parameter (default: NULL).
#' @param beta_g A positive scalar, hyperparameter that determines the
#'   shape of  prior beta distribution for guessing parameter (default: NULL).
#'
#' @return A list including:
#' \describe{
#'   \item{s_est}{the posterior mean of slip parameter.}
#'   \item{g_est}{the posterior mean of guessing parameter.}
#'   \item{s_sd}{the posterior standard diviation of slip parameter.}
#'   \item{g_sd}{the posterior standard diviation of guessing parameter.}
#   \item{r_il}{the estimates of parameter for categorical distribution that is optimal variational posterior for \emph{\strong{z}_i}.}
#'   \item{alpha_s_ast}{the estimates of variational parameter for slip parameter}
#'   \item{beta_s_ast}{the estimates of variational parameter for slip parameter}
#'   \item{alpha_g_ast}{the estimates of variational parameter for guessing parameter}
#'   \item{beta_g_ast}{the estimates of variational parameter for guessing parameter}
#'   \item{pi_est}{the estimates of class mixing parameter \eqn{\pi}}
#'   \item{delta_ast}{the estimates of variational parameter \eqn{\delta^*}}
#'   \item{delta_sd}{the standard diviation of variational parameter \eqn{\delta^*}}
#'   \item{l_lb}{the list of the values of evidence lower bound in each itertion}
#'   \item{att_pat_est}{the estimated attribute mastery patterns}
#   \item{A}{all of the possible attribute mastery patterns}
#   \item{Q}{the entered Q-matrix}
#   \item{X}{the entered data matrix}
#'   \item{eta_lj}{the computed ideal responce}
#'   \item{m}{the number of performed iterations}
#' }
#'
#' @examples
#' # load Q-matrix and create artificial item response data
#' Q = sim_Q_J80K5
#' sim_data = dina_data_gen(Q=Q,I=200)
#' # fit DINO model
#' res_dino = dino(X=sim_data$X, Q=Q)
#'
#' @export

#
# DINO VB
#

dino = function(
    X,
    Q,
    max_it  = 500,
    epsilon = 1e-04,
    verbose = TRUE,
    # Hyperparameters
    delta_0 = NULL, # For π
    alpha_s = NULL, # For s_j
    beta_s  = NULL, # For s_j
    alpha_g = NULL, # For g_j
    beta_g  = NULL # For g_j
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

  # Index
  I <- nrow(X)
  J <- ncol(X)
  K <- ncol(Q)
  L <- 2^K

  #
  # All attribute pattern matrix
  #
  A <- as.matrix(expand.grid(lapply(1:K, function(x)rep(0:1))))
  eta_lj <- (1-A) %*% t(Q)
  QQ <- diag(Q  %*% t(Q))

  #
  # hyperparameter
  #
  if(is.null(delta_0)){
    delta_0 = rep(1, L) # For π
  }
  if(is.null(alpha_s)){
    alpha_s = 1 # For s_j
  }
  if(is.null(beta_s)){
    beta_s = 1 # For s_j
  }
  if(is.null(alpha_g)){
    alpha_g = 1 # For g_j
  }
  if(is.null(beta_g)){
    beta_g = 1 # For g_j
  }

  #
  # Convert
  #
  for(j in 1:J)eta_lj[,j] <- 1*(eta_lj[,j] == QQ[j])
  eta_lj <- 1 - eta_lj
  #
  # Initialization
  #
  # r_il <- matrix(runif(I*L) ,ncol=L, nrow = I)
  # r_il <- diag(1/rowSums(r_il)) %*% r_il
  r_il <- matrix(1/L ,ncol=L, nrow = I)
  # r_il <- r_il + matrix(runif(I*L,-1/L,1/L) ,ncol=L, nrow = I)
  # r_il <- 1/rowSums(r_il) * r_il
  # rowSums(r_il)
  log_rho_il <- matrix(0,ncol=L, nrow = I)
  # z_il <- r_il



  delta_ast  <- rep(0,L)
  alpha_s_ast <- rep(0,J)
  beta_s_ast  <- rep(0,J)
  alpha_g_ast <- rep(0,J)
  beta_g_ast  <- rep(0,J)

  #
  one_vec  = matrix(1,nrow=I,ncol=1)


  m = 1
  #
  l_lb = rep(NA, max_it+1)
  l_lb[1] = 100
  for(m in 1:max_it){
    if(verbose){
      cat("\riteration = ", m, sprintf(": l_lb = %.05f", l_lb[m]))
    }

    #
    # M-step
    #
    delta_ast <- colSums(r_il) + delta_0

    alpha_s_ast <- colSums((r_il %*% eta_lj) * (1-X))    + alpha_s
    beta_s_ast  <- colSums((r_il %*% eta_lj) * X)        + beta_s
    alpha_g_ast <- colSums((1- r_il %*% eta_lj) * X)     + alpha_g
    beta_g_ast  <- colSums((1- r_il %*% eta_lj) * (1-X)) + beta_g

    #
    # Calculate Expectations
    #
    E_log_s   = digamma(alpha_s_ast) - digamma(alpha_s_ast + beta_s_ast)
    E_log_1_s = digamma(beta_s_ast)  - digamma(alpha_s_ast + beta_s_ast)
    E_log_g   = digamma(alpha_g_ast) - digamma(alpha_g_ast + beta_g_ast)
    E_log_1_g = digamma(beta_g_ast)  - digamma(alpha_g_ast + beta_g_ast)
    E_log_pi  = digamma(delta_ast)   - digamma(sum(delta_ast))


    #
    # E-step
    #

    #
    #
    # Fast
    #
    log_rho_il <- t(eta_lj %*% t(t((E_log_1_s - E_log_g) * t(X)) + t((E_log_s - E_log_1_g) * t(1-X)) ) ) + one_vec %*% E_log_pi
    # Slow
    # log_rho_il <- t(eta_lj %*% t(X %*% diag(E_log_1_s - E_log_g) + (1-X)%*% diag(E_log_s - E_log_1_g)) ) + one_vec %*% E_log_pi

    temp <- exp(log_rho_il)
    # Fast
    r_il <- 1/rowSums(temp) * temp
    # Slow
    # r_il = diag(1/rowSums(temp )) %*% temp



    #
    # Evidense of Lower Bound
    #
    #l_lb[m+1] = sum(lgamma(delta_ast) -lgamma(delta_0)) + lgamma(sum(delta_0)) - lgamma(sum(delta_ast))  + sum(lbeta(alpha_s_ast,beta_s_ast)+ lbeta(alpha_g_ast,beta_g_ast))  -length(alpha_s_ast)*(lbeta(alpha_s,beta_s) + lbeta(alpha_g,beta_g)) - sum(r_il*log(r_il))

    r_eta <- r_il %*% eta_lj
    one_m_r_eta <- 1 - r_eta
    tmp1 <- sum(X * t(t(r_eta)*E_log_1_s +t(one_m_r_eta)*E_log_g) + (1-X)*t(t(r_eta)*E_log_s +t(one_m_r_eta)*E_log_1_g)   )

    tmp2 <- sum(r_il*(one_vec%*%E_log_pi - log(r_il)))

    tmp3 <- sum(lgamma(delta_ast) -lgamma(delta_0)) + lgamma(sum(delta_0)) - lgamma(sum(delta_ast)) +sum(E_log_pi)

    tmp4 <- sum(lbeta(alpha_s_ast,beta_s_ast)-lbeta(alpha_s,beta_s) + (alpha_s - alpha_s_ast)*E_log_s + (beta_s - beta_s_ast)*E_log_1_s )
    tmp5 <- sum(lbeta(alpha_g_ast,beta_g_ast)-lbeta(alpha_g,beta_g) + (alpha_g - alpha_g_ast)*E_log_g + (beta_g - beta_g_ast)*E_log_1_g )

    l_lb[m+1] = tmp1 + tmp2 +tmp3+tmp4+tmp5

    if(verbose){
      cat("\riteration = ", m+1, sprintf(",last change = %.05f", abs(l_lb[m] - l_lb[m+1])))
    }

    if(abs(l_lb[m] -l_lb[m+1]) < epsilon){
      if(verbose){
        cat("\nreached convergence.")
      }
      break()
    }

  }
  l_lb <- l_lb[-1]
  # plot(l_lb,type="l")

  s_est <- alpha_s_ast /(alpha_s_ast + beta_s_ast )
  g_est <- alpha_g_ast /(alpha_g_ast + beta_g_ast )

  s_sd <-  sqrt(alpha_s_ast*beta_s_ast /(((alpha_s_ast + beta_s_ast )^2)*(alpha_s_ast + beta_s_ast +1)))
  g_sd <-  sqrt(alpha_g_ast*beta_g_ast /(((alpha_g_ast + beta_g_ast )^2)*(alpha_g_ast + beta_g_ast +1)))


  #
  # variance of estimates
  #
  delta_sum <- sum(delta_ast)
  pi_est <- delta_ast/delta_sum
  delta_sd <-sqrt(delta_ast*(delta_sum - delta_ast)/(delta_sum^2*(delta_sum+1)) )

  list(s_est = s_est,
       g_est = g_est,
       s_sd  = s_sd,
       g_sd  = g_sd,
       #r_il  = r_il,
       alpha_s_ast = alpha_s_ast,
       beta_s_ast  = beta_s_ast,
       alpha_g_ast = alpha_g_ast,
       beta_g_ast  = beta_g_ast,
       pi_est      = pi_est,
       delta_ast   = delta_ast,
       delta_sd    = delta_sd,
       l_lb = l_lb,
       att_pat_est = A[apply(r_il, 1, which.max),],
       #A = A,
       #Q = Q,
       #X = X,
       eta_lj = eta_lj,
       m = m)
}
