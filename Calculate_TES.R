rm(list = ls())
gc()

# Please set the working directory first !!!
setwd("c:/Users/mm3525/Software/RCB_compare")
tmp <- sapply(paste0("Rscripts/",list.files(path="Rscripts")), source)

# Please set the number of permutations, method to calculate TES 
# and names for control and experimental treatment
# RCB scores should be saved as 1-column .txt file in data folder (one file per treatment)
method <- "wKS" #"wKS","DensRatio" or "DensDiff
n_perm <- 10000
n_cores <- 1
control_treatment <- "Ctrl_treat"
experimental_treatment <- "Exp_treat"

# Loading or installing necessary R libraries
libs <- c("ggplot2","densratio","parallel","reshape","pcg","gridExtra")
check_libs(libs)

# Creating folder to store the results
f_name <- paste0(control_treatment,"_vs_",experimental_treatment)
if(!dir.exists(paste0("res/",f_name))) dir.create(paste0("res/",f_name), recursive=T)

# Loading RCB score values
RCB_ctrl <- as.numeric(read.delim(file=paste0("data/",control_treatment,".txt"), header=F)$V1)
RCB_exp <- as.numeric(read.delim(file=paste0("data/",experimental_treatment,".txt"), header=F)$V1)
cat("Control treatment cohort size:",length(RCB_ctrl),"\n")
cat("Experimental treatment cohort size:",length(RCB_exp),"\n")

# Calculating Treatment Efficacy Score and its p-value
res_TES <- switch(method,
  wKS = wKS(RCB_ctrl, RCB_exp, n_perm=n_perm, plots_folder=paste0("res/",f_name), prj_title=f_name, n_cores=n_cores),
  DensRatio = DensRatio(RCB_ctrl, RCB_exp, n_perm=n_perm, plots_folder=paste0("res/",f_name), prj_title=f_name, n_cores=n_cores),
  DensDiff = DensDiff(RCB_ctrl, RCB_exp, n_perm=n_perm, sig="auto", plots_folder=paste0("res/",f_name), prj_title=f_name, n_cores=n_cores))
print(res_TES)
