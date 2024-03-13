## ----eval=FALSE---------------------------------------------------------------
#  Q = sim_Q_J30K3
#  sim_data = dina_data_gen(Q=Q,I=200)
#  res = variationalDCM(X=sim_data$X, Q=Q, model="satu_dcm")
#  summary(res)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("variationalDCM")

## ---- eval=FALSE--------------------------------------------------------------
#  if(!require(devtools)){
#    install.packages("devtools")
#    }
#  devtools::install_github("khijikata/variationalDCM")

## ----eval=FALSE---------------------------------------------------------------
#  Q = sim_Q_J30K3
#  sim_data = dina_data_gen(Q=Q,I=200)

## ----eval=FALSE---------------------------------------------------------------
#  res = variationalDCM(X=sim_data$X, Q=Q, model="satu_dcm")
#  summary(res)

## ----eval=FALSE---------------------------------------------------------------
#  res = variationalDCM(X=sim_data$X, Q=Q, model="satu_dcm")
#  summary(res)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

