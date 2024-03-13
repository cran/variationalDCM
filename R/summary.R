#' @include variationalDCM.R
#' @export
#' @describeIn variationalDCM print summary information

summary.variationalDCM = function(object,...){

  output = list(
    attr_mastery_pat = object$att_pat_est,
    ELBO = object$l_lb[length(object$l_lb)],
    time = object$time
    )
  output = c(object$model_params, output)
  class(output) = "summary.variationalDCM"
  output
}
