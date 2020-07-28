wKS <- function(RCB_ctrl,RCB_exp, n_perm=10000, scale=0.5, ifplot=T, plots_folder=".", prj_title="Project", n_cores=1) {
# Estimate Treatment Efficacy Score (TES) using two sample weighted Kolmogorov-Smirnov test
# by comparing RCB score from experimental treatment cohort vs control cohort and Area Beetwen CDF curves (ABC)
#
# Input:
# RCB_ctrl - vector of RCB values in control cohort
# RCB_exp - vector of RCB values in experimental treatment cohort
# n_perm - number of test permutations (positive integer)
# scale - weight function parameter. Higher value gives higher weights to low RCB scores (non-negative value)
# ifplot - if plot should be generated and saved (logical)
# plots_folder - plot saving destination
# prj_title - name of the project to include in plot file name
#
# Output: TES, ABC and p-values for both statistics
# 
# Author: Michal Marczyk (michal.marczyk@yale.edu; michal.marczyk@polsl.edu)
  
  # calculate Treatment Efficacy Score
  res <- w_KS_stat(RCB_ctrl,RCB_exp, scale=scale, ifplot=ifplot, plots_folder=plots_folder, prj_title=prj_title)
  TES <- res[1]
  ABC <- res[2]
  
  # prepare data for permutation test 
  n1 <- length(RCB_exp)
  n2 <- length(RCB_ctrl)
  RCB_all <- c(RCB_ctrl, RCB_exp)
  perm_list <- lapply(1:n_perm,function(x){sample.int(n1+n2, n1)})
  
  # run permutation test in parallel
  res_perm <- mclapply(perm_list, function(d){
    res_KS <- w_KS_stat_fast(RCB_all[d], RCB_all[-d], scale=scale)
  },  mc.cores = n_cores)
  res_perm <- do.call(rbind,res_perm)
  colnames(res_perm) <- c("TES","ABC")

  # calculate p-value of one sided test
  TES_p <- sum(res_perm[,1] > TES)/n_perm  
  TES_p <- max(1/n_perm, TES_p)
  ABC_p <- sum(res_perm[,2] > ABC)/n_perm  
  ABC_p <- max(1/n_perm, ABC_p)
  
  # plot null distribution of TES  
  if(ifplot){
    # plot null distribution of TES and ABC
    data_plot <- data.frame(res_perm)
    p1 <- ggplot(data_plot, aes(x=TES)) + geom_histogram(bins=30, fill="lightblue", col="black") + 
      geom_vline(xintercept=TES, col="red") + annotate("text",x=TES, y=0, label="TES", col="red", hjust=0, vjust=1) +
      labs(x="TES null distribution", y="No. of permutations", title=paste0(prj_title,". TES = ",signif(TES,2), ". P = ",signif(TES_p,2))) + 
      theme_bw() + theme(plot.title=element_text(hjust=0.5))
    p2 <- ggplot(data_plot, aes(x=ABC)) + geom_histogram(bins=30, fill="lightblue", col="black") + 
      geom_vline(xintercept=ABC, col="red") + annotate("text",x=ABC, y=0, label="ABC", col="red", hjust=0, vjust=1) +
      labs(x="ABC null distribution", y="No. of permutations", title=paste0(prj_title,". ABC = ",signif(ABC,2), ". P = ",signif(ABC_p,2))) + 
      theme_bw() + theme(plot.title=element_text(hjust=0.5))
    
    pdf(file=paste0(plots_folder,"/wKS_perm_",prj_title,".pdf"), width=6, height=6)
    grid.arrange(p1,p2,nrow=2)
    dev.off()
  }
  
  # return results
  res <- c(TES, ABC, TES_p, ABC_p)
  names(res) <- c("TES", "ABC" ,"TES_p","ABC_p")
  return(res)
}