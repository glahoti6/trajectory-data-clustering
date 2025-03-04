func_clust_smm_cens = function(data,n_tran_st,n_abs_st,max_iter,eta,Kp,clust_iter) {
  
  s1 = n_tran_st
  s = s1 + n_abs_st
  cen_st = s + 1
  
  ##Clustering
  
  #Uncensored & Censored Sojourn Times
  
  if (n_abs_st>0) {
    uncen_time = sort(unique(data[-union(union(which(data$state %in% cen_st)-1,which(data$state %in% cen_st)),which(data$state %in% c((s1+1):s))),"time"]))
    Nu = uncen_time #Set of uncensored sojourn times
    M = length(Nu)
  } else if (n_abs_st==0) {
    uncen_time = sort(unique(data[-union(which(data$state %in% cen_st)-1,which(data$state %in% cen_st)),"time"]))
    Nu = uncen_time #Set of uncensored sojourn times
    M = length(Nu)
  }
  
  cen_time = sort(unique(data[(which(data$state %in% cen_st)-1),"time"]))       #Set of censored sojourn times
  
  #Initialization
  IDs = sort(unique(data[,"id"]))
  L = length(IDs)
  X = cbind(IDs, c(1:Kp, sample(1:Kp, size = L-Kp, replace = TRUE))) #Random cluster assignment to each ID ensuring that each cluster gets atleast one ID
  
  # print(table(X[,2]))
  
  cap_omega = matrix(1/Kp, L, Kp) #Membership function
  pi_old = rep(1/Kp,Kp)
  alpha = 1
  
  #M-Step
  
  smm_parameters = function(data, X, L, s1, s, M, Kp, Nu, cen_time, eta, cap_omega, pi_old, alpha) {
    cap_theta = list()
    Nu_ext = c(0,Nu)
    pi_em = colSums(cap_omega) / L
    pi = pi_em + (alpha*pi_old*(log(pi_old)-sum(pi_old*log(pi_old))))
    print(pi)
    alpha = min(sum(exp(-eta*L*abs(pi-pi_old)))/Kp,(1-max(pi_em))/(-max(pi_old)*(sum(pi_old*log(pi_old)))))
    print(alpha)
    discard_clust = which(pi<=(1/L))
    if (length(discard_clust)>0) {
      pi = pi[-discard_clust]
      pi = pi/sum(pi)
      Kp = Kp - length(discard_clust)
      cap_omega = cap_omega[,-discard_clust]
      for (i in c(1:dim(cap_omega)[1])) {
        if (sum(cap_omega[i,]) == 0) {
          cap_omega[i,sample(1:Kp,size=1)] = 1
        }
      }
      cap_omega = t(t(cap_omega) %*% diag(1/rowSums(cap_omega)))
      X = cbind(IDs,max.col(cap_omega))
    }
    rho = matrix(0, s1, Kp)
    m = array(0, c(s1, s, M, Kp)) 
    m_s_plus_1 = array(0, c(s1, 1, M+1, Kp))
    for (l in c(1:L)) {
      traj = data[which(data$id %in% X[l,1]),]
      rho[traj[1,"state"],X[l,2]] = rho[traj[1,"state"],X[l,2]] + cap_omega[l,X[l,2]]
      for (n_h in c(1:(nrow(traj)-1))) {
        if (traj[n_h+1, "state"] != s+1) {
          time_ind = which(Nu %in% traj[n_h,"time"])
          m[traj[n_h,"state"],traj[n_h+1,"state"],time_ind,X[l,2]] = m[traj[n_h,"state"],traj[n_h+1,"state"],time_ind,X[l,2]] + cap_omega[l,X[l,2]]
        } else {
          cen_time_lvals_ind = which(Nu_ext <= traj[n_h,"time"])
          cen_time_lb_ind = cen_time_lvals_ind[length(cen_time_lvals_ind)]
          m_s_plus_1[traj[n_h,"state"],1,cen_time_lb_ind,X[l,2]] = m_s_plus_1[traj[n_h,"state"],1,cen_time_lb_ind,X[l,2]] + cap_omega[l,X[l,2]]
        }
      }
    }
    
    rho = rho %*% diag(1/colSums(rho))
    
    n = array(0, c(s1, M, Kp))
    lambda_npmle = array(0, c(s1, s, M, Kp))
    q = array(0, c(s1, M, Kp))
    for (i in c(1:s1)) {
      for (k in c(1:M)) {
        for (kp in c(1:Kp)) {
          n[i,k,kp] = sum(colSums(m[i,,,kp])[k:M]) + sum(m_s_plus_1[i,1,,kp][(k+1):(M+1)])
          if (n[i,k,kp]>0) {
            q[i,k,kp] = 1 - (sum(m[i,,k,kp])/n[i,k,kp])
            lambda_npmle[i,1:s,k,kp] = m[i,1:s,k,kp]/n[i,k,kp]
          } else if (n[i,k,kp]==0) {
            q[i,k,kp] = 1
            lambda_npmle[i,1:s,k,kp] = 0
          }
        }
      }
    }
    
    S = array(0, c(M+1+1,s1,Kp))
    S[1,1:s1,1:Kp] = 1
    for (k in c(2:(M+1))) {
      for (i in c(1:s1)) {
        for (kp in c(1:Kp)) {
          S[k,i,kp] = prod(q[i,1:(k-1),kp])
        }
      }
    }
    S[M+1+1,1:s1,1:Kp] = 0
    
    X_temp = X
    for (kp in c(1:Kp)) {
      traj_id_grp_kp = X_temp[which(X_temp[,2]==kp),1]
      traj_grp_kp = data[data[,"id"] %in% traj_id_grp_kp,]
      entries_cen_tran_st = traj_grp_kp[which(traj_grp_kp$state == cen_st)-1,]
      if (nrow(entries_cen_tran_st) != 0) {            
        for (i in c(1:n_tran_st)) {
          entries_cen_tran_st_i = entries_cen_tran_st[which(entries_cen_tran_st$state==i),]
          if (nrow(entries_cen_tran_st_i) != 0) {
            entries_cen_tran_st_i_soj = sort(unique(entries_cen_tran_st_i[,"time"]))
            if (entries_cen_tran_st_i_soj[length(entries_cen_tran_st_i_soj)] > uncen_time[length(uncen_time)]) {
              S[M+1+1,i,kp] = 0.0001*S[M+1,i,kp]
            }
          }
        }
      }
    }
    
    cap_theta$pi = pi
    cap_theta$rho = rho
    cap_theta$lambda_npmle = lambda_npmle
    cap_theta$S = S
    cap_theta$Kp = Kp
    cap_theta$pi_old = pi
    cap_theta$alpha = alpha
    
    return(cap_theta)
  }
  
  #E-Step
  
  membership_prob = function(data, X, cap_theta, L, Kp, s, Nu) {
    cap_omega = matrix(0, L, Kp)
    pi = cap_theta$pi 
    rho = cap_theta$rho
    lambda_npmle = cap_theta$lambda_npmle
    S = cap_theta$S 
    Nu_ext = c(0,Nu,Nu[length(Nu)]+1e-05)
    for (kp in c(1:Kp)) {
      for (l in c(1:L)) {
        traj = data[which(data$id %in% X[l,1]),]
        cap_omega[l,kp] = pi[kp] * rho[traj[1,"state"],kp]
        for (n_h in c(1:(nrow(traj)-1))) {
          if (traj[n_h+1,"state"] != s+1) {
            uncen_time_val_ind = which(Nu %in% traj[n_h,"time"])
            uncen_time_lvals_ind = which(Nu_ext < traj[n_h,"time"])
            uncen_time_lb_ind = uncen_time_lvals_ind[length(uncen_time_lvals_ind)]
            cap_omega[l,kp] = cap_omega[l,kp] * lambda_npmle[traj[n_h,"state"],traj[n_h+1,"state"],uncen_time_val_ind,kp] * S[uncen_time_lb_ind,traj[n_h,"state"],kp]
          }
          else {
            cen_time_lvals_ind = which(Nu_ext <= traj[n_h,"time"])
            cen_time_lb_ind = cen_time_lvals_ind[length(cen_time_lvals_ind)]
            cap_omega[l,kp] = cap_omega[l,kp] * S[cen_time_lb_ind,traj[n_h,"state"],kp]
          }
        }
      }
    }
    cap_omega = t(t(cap_omega) %*% diag(1/rowSums(cap_omega))) #If all the values in a row are zero, this normalization will generate NAs
    cap_omega[is.nan(cap_omega)] = 0
    
    return(cap_omega)
  }
  
  #EM
  
  track_Kp = NULL
  
  cap_theta = smm_parameters(data, X, L, s1, s, M, Kp, Nu, cen_time, eta, cap_omega, pi_old, alpha) 
  Kp = cap_theta$Kp
  pi_old = cap_theta$pi_old
  alpha = cap_theta$alpha
  track_Kp[1] = Kp 
  
  for (i in c(2:max_iter)) {
    #print(i)
    cap_omega = membership_prob(data, X, cap_theta, L, Kp, s, Nu) 
    X = cbind(IDs,max.col(cap_omega))
    #print(table(X[,2]))
    cap_theta = smm_parameters(data, X, L, s1, s, M, Kp, Nu, cen_time, eta, cap_omega, pi_old, alpha) 
    Kp = cap_theta$Kp
    pi_old = cap_theta$pi_old
    alpha = cap_theta$alpha
    track_Kp[i] = Kp
    if (i>clust_iter) {
      if ((track_Kp[i]-track_Kp[i-clust_iter]) == 0) {
        print(i)
        cap_theta$clus_iter = i
        break
      }
    }
  }
  
  lambda_npmle = cap_theta$lambda_npmle
  S = cap_theta$S
  
  theta = array(0,c(s1,s,M,Kp))
  trans_prob = array(0,c(s1,s,Kp))
  for (i in c(1:s1)) {
    for (j in c(1:s)) {
      for (kp in c(1:Kp)) {
        for (k in c(1:M)) {
          theta[i,j,k,kp] = lambda_npmle[i,j,k,kp] * S[k,i,kp]
        }
        trans_prob[i,j,kp] = sum(theta[i,j,,kp])
      }
    }
  }
  
  H = array(0,c(M,s1,s,Kp))
  for (k in c(1:M)) {
    for (i in c(1:s1)) {
      for (j in c(1:s)) {
        for (kp in c(1:Kp)) {
          if (trans_prob[i,j,kp]>0) {
            H[k,i,j,kp] = theta[i,j,k,kp] / trans_prob[i,j,kp]
          } else if (trans_prob[i,j,kp]==0) {
            H[k,i,j,kp] = 0
          }
        }
      }
    }
  }
  
  #Normalization of transition probability matrices for the clusters
  for (i in c(1:dim(trans_prob)[length(dim(trans_prob))])) {
    trans_prob[,,i] = t(t(trans_prob[,,i]) %*% diag(1/rowSums(trans_prob[,,i])))          
  }
  trans_prob[is.nan(trans_prob)] = 0
  
  Q = array(0,c(M,s1,s,Kp))
  for (k in c(1:M)) {
    for (i in c(1:s1)) {
      for (j in c(1:s)) {
        for (kp in c(1:Kp)) {
          if (is.nan(trans_prob[i,j,kp])) {
            Q[k,i,j,kp] = 1
          } else if (trans_prob[i,j,kp]>0) {
            Q[k,i,j,kp] = (trans_prob[i,j,kp] - sum(theta[i,j,1:k,kp])) / (trans_prob[i,j,kp])
          } else if (trans_prob[i,j,kp]==0) {
            Q[k,i,j,kp] = 1
          }
        }
      }
    }
  }
  
  #Normalization of survival function
  for (j in c(1:s)) {
    for (kp in c(1:Kp)) {
      Q[,,j,kp] = prop.table(Q[,,j,kp],2)
    }
  }
  
  cap_theta$trans_prob = trans_prob
  cap_theta$H = H
  cap_theta$Q = Q
  
  result= list()
  result$cap_omega = cap_omega
  result$X = X
  result$cap_theta = cap_theta
  
  return(result)
  
}