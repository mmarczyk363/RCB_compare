DensDiff <- function(RCB_ctrl, RCB_exp, n_perm=10000, sig="auto", lam=10^seq(-2,1,length.out=50),
                          ifplot=T, plots_folder=".", prj_title="Project", n_cores=1){
# Estimate Treatment Efficacy Score (TES) by comparing densities of RCB in experimental treatment cohort vs control cohort
#
# Input:
# RCB_ctrl - vector of RCB values in control cohort
# RCB_exp - vector of RCB values in experimental treatment cohort
# sig - kernel bandwidth (vector with numeric values > 0)
# lam - regularization parameter (vector with numeric values > 0)
# ifplot - if plot should be generated and saved (logical)
# plots_folder - plot saving destination
# prj_title - name of the project to include in plot file name
#
# Output: TES with 95% confidence intervals
# 
# Author: Michal Marczyk (michal.marczyk@yale.edu; michal.marczyk@polsl.edu)
# Based on Matlab implementation of Least-Squares Density-Difference Estimation by Masashi Sugiyama
# sugi@cs.titech.ac.jp, sugi@k.u-tokyo.ac.jp
  
  require(ggplot2)
  
  if (sig == "auto"){
    sig <- bw.nrd0(c(RCB_ctrl, RCB_exp)) 
  }
  
  # compute density difference
  res <- RCB_diff(RCB_ctrl, RCB_exp, sig=sig, lam=lam, ifplot=ifplot, plots_folder=plots_folder, prj_title=prj_title)
  
  # prepare data for permutation test 
  n1 <- length(RCB_exp)
  n2 <- length(RCB_ctrl)
  RCB_all <- c(RCB_ctrl, RCB_exp)
  perm_list <- lapply(1:n_perm,function(x){sample.int(n1+n2, n1)})
  
  # run permutation test in parallel with fixed density ratio parameters
  res_perm <- mclapply(perm_list, function(d){
    RCB_diff_fast(RCB_all[d], RCB_all[-d], sigma=res$Parameters["Sigma"], lambda=res$Parameters["Lambda"])
  },  mc.cores = n_cores)
  res_perm <- do.call(c,res_perm)
  TES_p <- sum(res_perm > res$TES)/n_perm
  TES_p <- max(1/n_perm, TES_p)
  
  if (ifplot){
    # plot distribution of TES
    data_plot <- data.frame(res_perm)
    p <- ggplot(data_plot, aes(x=res_perm)) + geom_histogram(bins=30, fill="lightblue", col="black") + 
      geom_vline(xintercept=res$TES, col="red") + annotate("text",x=res$TES, y=0, label="TES", col="red", hjust=0, vjust=1) +
      labs(x="TES null distribution", y="No. of permutations", title=paste0(prj_title,". TES = ",signif(res$TES,2), ". P = ",signif(TES_p,2))) + 
      theme_bw() + theme(plot.title=element_text(hjust=0.5))
    pdf(file=paste0(plots_folder,"/DensDiff_perm_",prj_title,".pdf"), width=6, height=4)
    print(p)
    dev.off()
  }

  # return results
  res <- c(res$TES, TES_p)
  names(res) <- c("TES", "TES_p")
  return(res)
}
