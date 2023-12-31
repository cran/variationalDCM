% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/satu_dcm_vb.R
\name{satu_dcm}
\alias{satu_dcm}
\title{for the saturated DCM.}
\usage{
satu_dcm(
  X,
  Q,
  max_it = 500,
  epsilon = 1e-04,
  seed = 123,
  verbose = TRUE,
  delta_0 = NULL,
  A_0 = NULL,
  B_0 = NULL
)
}
\arguments{
\item{X}{I by J binary matrix, item response data}

\item{Q}{J by K binary matrix, Q-matrix}

\item{max_it}{The maximum number of iterations (default: 500)}

\item{epsilon}{The convergence tolerance for iterations (default: 1e-4)}

\item{seed}{The seed value (default: 123)}

\item{verbose}{Logical, controls whether to print progress (default: TRUE)}

\item{delta_0}{the value of hyperparameter \eqn{\delta^0} (default: NULL)}

\item{A_0}{the value of hyperparameter \eqn{A^0} (default: NULL)}

\item{B_0}{the value of hyperparameter \eqn{B^0} (default: NULL)}
}
\value{
A list including:
\describe{
  \item{theta_est}{the estimate of the conditional response probability paramter \eqn{\Theta}}
  \item{theta_sd}{the posterior standard deviation of parameter \eqn{\Theta}}
  \item{pi_est}{the estimates of class mixing parameter \eqn{\pi}}
  \item{pi_sd}{the posterior standard deviations of class mixing parameter \eqn{\pi}}
  \item{A_ast}{the estimate of variational parameter \eqn{A^*}}
  \item{B_ast}{the estimate of variational parameter \eqn{B^*}}
  \item{delta_ast}{the estimate of variational parameter \eqn{\delta^*}}
  \item{A_0}{the value of hyperparameter \eqn{A^0}}
  \item{B_0}{the value of hyperparameter \eqn{B^0}}
  \item{delta_0}{the value of hyperparameter \eqn{\delta^0}}
  \item{l_lb}{the list of the values of evidence lower bound at each itertion}
  \item{att_pat_est}{the estimated attribute mastery patterns}
  \item{G_j}{the computed G-matrix}
  \item{m}{the number of the performed iterations}
  \item{seed}{the entered seed number}
}
}
\description{
\code{satu_dcm()} returns variational Bayesian estimates for the saturated DCM.
}
\examples{
# load Q-matrix and create artificial item response data
Q = sim_Q_J30K3
sim_data = dina_data_gen(Q=Q,I=200)
# fit saturated DCM
res_satu = satu_dcm(X=sim_data$X, Q=Q)

}
\references{
Yamaguchi, K., Okada, K. (2020). Variational Bayes Inference
  Algorithm for the Saturated Diagnostic Classification Model.
  \emph{Psychometrika}, 85(4), 973–995. \doi{10.1007/s11336-020-09739-w}
}
