w_KS_stat_fast <- function(RCB_ctrl,RCB_exp, scale=0.5) {
  
  # create weighted eCDFs and calculated their difference
  ecdf_ctrl <- w_eCDF(RCB_ctrl, scale=scale)
  ecdf_exp <- w_eCDF(RCB_exp, scale=scale)
  ecdf_diff <- function(x){ecdf_exp(x) - ecdf_ctrl(x)}
  
  # calculate Treatment Efficacy score (TES) as supremum of ecdf_diff
  x <- sort(unique(c(0,RCB_ctrl,RCB_exp)))
  ecdf_diff_calc <- ecdf_diff(x)
  D_range <- range(ecdf_diff_calc)
  if (abs(D_range[1]) > abs(D_range[2])){
    TES <- D_range[1]
  } else {
    TES <- D_range[2]
  }
  
  # get AUC from ecdf_diff
  AUC <- sum(ecdf_diff_calc[-length(ecdf_diff_calc)] * (x[-1] - x[-length(x)]))/max(x)
  
  return(c(TES,AUC))
}