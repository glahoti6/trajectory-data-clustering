library("ggplot2")
library("caret")
library("e1071")
library("backports")
library("seqHMM")
library("matrixStats")
library("proto")
library("gsubfn")
library("RSQLite")
library("sqldf")
library("gtools")
library("reshape2")
library("lattice")
library("ClusterR")
library("Biostrings")

source('func_clust_smm_cens.R')

#Import the data
full_data = read.table('data.csv', header = TRUE, sep = ",")          
#Import the data. If there are 4 transient states, then they will be denoted by 1,2,3,4. And if there are two absorbing states, then 
#they will be denoted by 5,6. Therefore, the censored state will be denoted by 7.

data_samples = sort(unique(full_data[,"id"]))
print(length(data_samples))

set.seed(0)
samp_frac = 1
if (samp_frac != 1) {
  data_samples_sel = sort(sample(data_samples,size = samp_frac*length(data_samples),replace = FALSE))
  print(length(data_samples_sel))
  data = full_data[which(full_data[,"id"] %in% data_samples_sel),]
} else {
  data = full_data
}
data[which(data[,"state"] == 4),"time"] = 0

#States
n_tran_st = 3           #Provide the total number of transient states 
n_abs_st = 0            #Provide the total number of absorbing states

#Clustering

allresult = list()

max_iter = 75 
extra_it = 10
Kp_init = 15 # Try different values
m_avg = ceiling(mean(unlist(sqldf("SELECT max(seq) FROM data GROUP BY id"))))         # 7
clust_iter = 25 #Number of consecutive iterations for which # of clusters should remain stable
eta = min(1,(0.5)^(floor((m_avg/2)-1)))

set.seed(0)
result1_1 = func_clust_smm_cens(data,n_tran_st,n_abs_st,max_iter,eta,Kp_init,clust_iter)
X1_1 = result1_1$X

result1_1$ntrans = n_tran_st
result1_1$nabs = n_abs_st
result1_1$maxiter = max_iter
result1_1$extra_it = extra_it

result1_1$m_avg = m_avg
result1_1$den = den
result1_1$eta = eta
result1_1$Kp_init = Kp_init
result1_1$clust_iter = clust_iter

allresult[[1]] = result1_1

saveRDS(allresult,file = "results.rds")