# **************************************************************************
# Test reciprocity effect
# Data input: alpha - if linear part degenerates, sample n^alpha sub-network 
#                     for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, p-value, 
#             and degenerate test result of linear part (C=1)
# **************************************************************************

NET_re <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]

  # test linear part degenerate or not
  xi.matr <- adj_mat
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements)
  
  # Var(e_ij)_mmt_estimation
  zero_dig <- adj_mat
  diag(zero_dig) <- 0
  
  hat_g11 <- rep(NA,n)
  for (i in 1:n) {
    hat_g11[i] <- sum(zero_dig[i,] + zero_dig[,i])/(n-1)/2
  }
  hat_g11 <- hat_g11 - mean(xi.vec)
  
  # test using concentration result
  tmp_1 <- xi.matr
  tmp_2 <- xi.matr
  tmp_1[lower.tri(tmp_1)] <- NA
  tmp_2[upper.tri(tmp_2)] <- NA
  
  a_2_clt <- na.omit(as.vector(t(tmp_1)))
  b_2_clt <- na.omit(as.vector(tmp_2))
  numerator_2 <- mean(a_2_clt * b_2_clt) - mean(xi.vec) * mean(xi.vec) 
  
  hat_g21 <- rep(NA,n)
  for (i in 1:n) {
    hat_g21[i] <- sum(zero_dig[i,] * zero_dig[,i])/(n-1)
  }
  hat_g21 <- hat_g21 - mean(a_2_clt * b_2_clt)
  sigma_square_21 <- mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11)^2)  
  denominator_2 <- sqrt(sigma_square_21/n)
  
  if(sigma_square_21 > sqrt(log(n)/n)){
    
    result <- matrix(c(numerator_2, denominator_2, numerator_2/denominator_2, 2*(1-pnorm(abs(numerator_2/denominator_2))), 'non-degenerate'), ncol=5, byrow=TRUE)
    colnames(result) <- c('estimate','std.error','statistic','pvalue', 'linear part')
    rownames(result) <- c('reciprocity effect')
    result <- as.table(result)
    
    list(estimate = numerator_2,
         std.error = denominator_2,
         statistic = numerator_2/denominator_2,
         pvalue = 2*(1-pnorm(abs(numerator_2/denominator_2))),
         sumry = result)
    
  }else{

    a_n <- 1/n/(n-1)
    b_n <- (n-2)/n/(n-1)
    c_n <- (n-2)*(n-3)/n/(n-1)
    
    single_kernel_ord4 <- function(adj_mat_4t4){
      
      h_stars.save <- rep(NA,6)
      
      xi.matr <- adj_mat_4t4
      diag(xi.matr) <- NA
      elements <- as.vector(t(xi.matr))
      xi.vec <- na.omit(elements)
      
      # h1_star
      h_stars.save[1] <- mean(xi.vec^2)
      
      # h2_star
      tmp_1 <- xi.matr
      tmp_2 <- xi.matr
      tmp_1[lower.tri(tmp_1)] <- NA
      tmp_2[upper.tri(tmp_2)] <- NA
      
      a_2_clt <- na.omit(as.vector(t(tmp_1)))
      b_2_clt <- na.omit(as.vector(tmp_2))
      
      h_stars.save[2] <- mean(a_2_clt * b_2_clt)
      
      # h3_star
      a_3 <- rep(xi.vec, each = 2)
      B.temp <- matrix(xi.vec, 3, 4)
      combine.temp_3 <- B.temp[-1,]
      for(i in 2:3){
        combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
      }
      b_3 <- as.vector(combine.temp_3)
      
      h_stars.save[3] <- mean(a_3 * b_3)
      
      # h4_star
      xi.vec.transp <- na.omit(as.vector(xi.matr))
      a_4 <- rep(xi.vec.transp, each = 2)
      D.temp <- matrix(xi.vec.transp, 3, 4)
      combine.temp_4 <- D.temp[-1,]
      for(i in 2:3){
        combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
      }
      b_4 <- as.vector(combine.temp_4)
      
      h_stars.save[4] <- mean(a_4 * b_4)
      
      # h5_star
      h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
      
      # h6_star
      e <- adj_mat_4t4
      h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                             e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                             e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
      
      h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
        b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
      
      # combined to be estimator's kernel
      single_kernel_ord4_ouput <- rep(NA,5)
      single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
      single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
      single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
      single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
      
      return(single_kernel_ord4_ouput) 
      
    }
    
    replctn <- floor(n^alpha)
    incomplete_Us <- matrix(NA,replctn,5)
    
    for (i in 1:replctn){
      ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
      adj_mat_dim4 <- adj_mat[ijkl,ijkl]
      incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
    }
    
    eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
    eta_est_mean <- apply(incomplete_Us,2,mean)
    eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
    
    result <- matrix(c(eta_est_mean[2], eta_deno[2], eta_est[2], 2*(1-pnorm(abs(eta_est[2]))), 'degenerate'), ncol=5, byrow=TRUE)
    colnames(result) <- c('estimate','std.error','statistic','pvalue', 'linear part')
    rownames(result) <- c('reciprocity effect')
    result <- as.table(result)
    
    list(estimate = eta_est_mean[2],
         std.error = eta_deno[2],
         statistic = eta_est[2],
         pvalue = 2*(1-pnorm(abs(eta_est[2]))),
         sumry = result)
  }
}



# **************************************************************************
# Test reciprocity effect given that linear part degenerate
# Data input: alpha - sample n^alpha sub-network for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, and p-value
# **************************************************************************

NET_re_degen <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]
  
  a_n <- 1/n/(n-1)
  b_n <- (n-2)/n/(n-1)
  c_n <- (n-2)*(n-3)/n/(n-1)
  
  single_kernel_ord4 <- function(adj_mat_4t4){
    
    h_stars.save <- rep(NA,6)
    
    xi.matr <- adj_mat_4t4
    diag(xi.matr) <- NA
    elements <- as.vector(t(xi.matr))
    xi.vec <- na.omit(elements)
    
    # h1_star
    h_stars.save[1] <- mean(xi.vec^2)
    
    # h2_star
    tmp_1 <- xi.matr
    tmp_2 <- xi.matr
    tmp_1[lower.tri(tmp_1)] <- NA
    tmp_2[upper.tri(tmp_2)] <- NA
    
    a_2_clt <- na.omit(as.vector(t(tmp_1)))
    b_2_clt <- na.omit(as.vector(tmp_2))
    
    h_stars.save[2] <- mean(a_2_clt * b_2_clt)
    
    # h3_star
    a_3 <- rep(xi.vec, each = 2)
    B.temp <- matrix(xi.vec, 3, 4)
    combine.temp_3 <- B.temp[-1,]
    for(i in 2:3){
      combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
    }
    b_3 <- as.vector(combine.temp_3)
    
    h_stars.save[3] <- mean(a_3 * b_3)
    
    # h4_star
    xi.vec.transp <- na.omit(as.vector(xi.matr))
    a_4 <- rep(xi.vec.transp, each = 2)
    D.temp <- matrix(xi.vec.transp, 3, 4)
    combine.temp_4 <- D.temp[-1,]
    for(i in 2:3){
      combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
    }
    b_4 <- as.vector(combine.temp_4)
    
    h_stars.save[4] <- mean(a_4 * b_4)
    
    # h5_star
    h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
    
    # h6_star
    e <- adj_mat_4t4
    h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                           e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                           e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
    
    h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
      b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
    
    # combined to be estimator's kernel
    single_kernel_ord4_ouput <- rep(NA,5)
    single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
    single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
    single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
    single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
    
    return(single_kernel_ord4_ouput) 
    
  }
  
  replctn <- floor(n^alpha)
  incomplete_Us <- matrix(NA,replctn,5)
  
  for (i in 1:replctn){
    ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
    adj_mat_dim4 <- adj_mat[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
  }
  
  eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
  eta_est_mean <- apply(incomplete_Us,2,mean)
  eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
  
  result <- matrix(c(eta_est_mean[2], eta_deno[2], eta_est[2], 2*(1-pnorm(abs(eta_est[2])))), ncol=4, byrow=TRUE)
  colnames(result) <- c('estimate','std.error','statistic','pvalue')
  rownames(result) <- c('reciprocity effect')
  result <- as.table(result)
  
  list(estimate = eta_est_mean[2],
       std.error = eta_deno[2],
       statistic = eta_est[2],
       pvalue = 2*(1-pnorm(abs(eta_est[2]))),
       sumry = result)
}



# **************************************************************************
# Test same sender effect
# Data input: alpha - sample n^alpha sub-network for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, and p-value
# **************************************************************************

NET_ss <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]
  
  a_n <- 1/n/(n-1)
  b_n <- (n-2)/n/(n-1)
  c_n <- (n-2)*(n-3)/n/(n-1)
  
  single_kernel_ord4 <- function(adj_mat_4t4){
    
    h_stars.save <- rep(NA,6)
    
    xi.matr <- adj_mat_4t4
    diag(xi.matr) <- NA
    elements <- as.vector(t(xi.matr))
    xi.vec <- na.omit(elements)
    
    # h1_star
    h_stars.save[1] <- mean(xi.vec^2)
    
    # h2_star
    tmp_1 <- xi.matr
    tmp_2 <- xi.matr
    tmp_1[lower.tri(tmp_1)] <- NA
    tmp_2[upper.tri(tmp_2)] <- NA
    
    a_2_clt <- na.omit(as.vector(t(tmp_1)))
    b_2_clt <- na.omit(as.vector(tmp_2))
    
    h_stars.save[2] <- mean(a_2_clt * b_2_clt)
    
    # h3_star
    a_3 <- rep(xi.vec, each = 2)
    B.temp <- matrix(xi.vec, 3, 4)
    combine.temp_3 <- B.temp[-1,]
    for(i in 2:3){
      combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
    }
    b_3 <- as.vector(combine.temp_3)
    
    h_stars.save[3] <- mean(a_3 * b_3)
    
    # h4_star
    xi.vec.transp <- na.omit(as.vector(xi.matr))
    a_4 <- rep(xi.vec.transp, each = 2)
    D.temp <- matrix(xi.vec.transp, 3, 4)
    combine.temp_4 <- D.temp[-1,]
    for(i in 2:3){
      combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
    }
    b_4 <- as.vector(combine.temp_4)
    
    h_stars.save[4] <- mean(a_4 * b_4)
    
    # h5_star
    h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
    
    # h6_star
    e <- adj_mat_4t4
    h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                           e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                           e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
    
    h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
      b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
    
    # combined to be estimator's kernel
    single_kernel_ord4_ouput <- rep(NA,5)
    single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
    single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
    single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
    single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
    
    return(single_kernel_ord4_ouput) 
    
  }
  
  replctn <- floor(n^alpha)
  incomplete_Us <- matrix(NA,replctn,5)
  
  for (i in 1:replctn){
    ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
    adj_mat_dim4 <- adj_mat[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
  }
  
  eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
  eta_est_mean <- apply(incomplete_Us,2,mean)
  eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
  
  result <- matrix(c(eta_est_mean[3], eta_deno[3], eta_est[3], 2*(1-pnorm(abs(eta_est[3])))), ncol=4, byrow=TRUE)
  colnames(result) <- c('estimate','std.error','statistic','pvalue')
  rownames(result) <- c('same sender effect')
  result <- as.table(result)
  
  list(estimate = eta_est_mean[3],
       std.error = eta_deno[3],
       statistic = eta_est[3],
       pvalue = 2*(1-pnorm(abs(eta_est[3]))),
       sumry = result)
}

st_a_null_pois <- function(n){
  
  ###
  ### Generate_Epsilon_ij
  ###
  epsilon_ij.save <- matrix(rpois(n^2, 1), n, n)
  
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- epsilon_ij.save
  diag(e_ij.save) <- NA
  
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save)
}  

Error_term.example <- st_a_null_pois(n = 10^4)

result <- NET_ss(alpha = 1, adj_mat = Error_term.example$e_ij.save)

# **************************************************************************
# Test same receiver effect
# Data input: alpha - sample n^alpha sub-network for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, and p-value
# **************************************************************************

NET_sr <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]
  
  a_n <- 1/n/(n-1)
  b_n <- (n-2)/n/(n-1)
  c_n <- (n-2)*(n-3)/n/(n-1)
  
  single_kernel_ord4 <- function(adj_mat_4t4){
    
    h_stars.save <- rep(NA,6)
    
    xi.matr <- adj_mat_4t4
    diag(xi.matr) <- NA
    elements <- as.vector(t(xi.matr))
    xi.vec <- na.omit(elements)
    
    # h1_star
    h_stars.save[1] <- mean(xi.vec^2)
    
    # h2_star
    tmp_1 <- xi.matr
    tmp_2 <- xi.matr
    tmp_1[lower.tri(tmp_1)] <- NA
    tmp_2[upper.tri(tmp_2)] <- NA
    
    a_2_clt <- na.omit(as.vector(t(tmp_1)))
    b_2_clt <- na.omit(as.vector(tmp_2))
    
    h_stars.save[2] <- mean(a_2_clt * b_2_clt)
    
    # h3_star
    a_3 <- rep(xi.vec, each = 2)
    B.temp <- matrix(xi.vec, 3, 4)
    combine.temp_3 <- B.temp[-1,]
    for(i in 2:3){
      combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
    }
    b_3 <- as.vector(combine.temp_3)
    
    h_stars.save[3] <- mean(a_3 * b_3)
    
    # h4_star
    xi.vec.transp <- na.omit(as.vector(xi.matr))
    a_4 <- rep(xi.vec.transp, each = 2)
    D.temp <- matrix(xi.vec.transp, 3, 4)
    combine.temp_4 <- D.temp[-1,]
    for(i in 2:3){
      combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
    }
    b_4 <- as.vector(combine.temp_4)
    
    h_stars.save[4] <- mean(a_4 * b_4)
    
    # h5_star
    h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
    
    # h6_star
    e <- adj_mat_4t4
    h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                           e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                           e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
    
    h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
      b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
    
    # combined to be estimator's kernel
    single_kernel_ord4_ouput <- rep(NA,5)
    single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
    single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
    single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
    single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
    
    return(single_kernel_ord4_ouput) 
    
  }
  
  replctn <- floor(n^alpha)
  incomplete_Us <- matrix(NA,replctn,5)
  
  for (i in 1:replctn){
    ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
    adj_mat_dim4 <- adj_mat[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
  }
  
  eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
  eta_est_mean <- apply(incomplete_Us,2,mean)
  eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
  
  result <- matrix(c(eta_est_mean[4], eta_deno[4], eta_est[4], 2*(1-pnorm(abs(eta_est[4])))), ncol=4, byrow=TRUE)
  colnames(result) <- c('estimate','std.error','statistic','pvalue')
  rownames(result) <- c('same receiver effect')
  result <- as.table(result)
  
  list(estimate = eta_est_mean[4],
       std.error = eta_deno[4],
       statistic = eta_est[4],
       pvalue = 2*(1-pnorm(abs(eta_est[4]))),
       sumry = result)
}



# **************************************************************************
# Test sender-receiver effect
# Data input: alpha - if linear part degenerates, sample n^alpha sub-network 
#                     for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, p-value, 
#             and degenerate test result of linear part (C=1)
# **************************************************************************

NET_s.r <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]
  
  # test linear part degenerate or not
  xi.matr <- adj_mat
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements)
  
  # Var(e_ij)_mmt_estimation
  zero_dig <- adj_mat
  diag(zero_dig) <- 0
  
  hat_g11 <- rep(NA,n)
  for (i in 1:n) {
    hat_g11[i] <- sum(zero_dig[i,] + zero_dig[,i])/(n-1)/2
  }
  hat_g11 <- hat_g11 - mean(xi.vec)
  
  # test using concentration result
  a_3 <- rep(xi.vec, each = n-2)
  xi.vec.transp <- na.omit(as.vector(xi.matr))
  D.temp <- matrix(xi.vec.transp, (n-1), n)
  combine.temp_4 <- D.temp[-1,]
  for(i in 2:(n-1)){
    combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
  }
  b_4 <- as.vector(combine.temp_4)
  numerator_5 <- mean(c(a_3,b_4) * c(b_4,a_3)) - mean(xi.vec) * mean(xi.vec)
  
  hat_g51_1 <- rep(NA,n)
  hat_g51_2 <- rep(NA,n)
  hat_g51_3 <- rep(NA,n)
  for (i in 1:n) {
    hat_g51_1[i] <- sum(c(a_3,b_4)[((i-1)*(n-1)*(n-2)+1):(i*(n-1)*(n-2))] * c(b_4,a_3)[((i-1)*(n-1)*(n-2)+1):(i*(n-1)*(n-2))])
  }
  for (i in 1:n) {
    hat_g51_2[i] <- sum(rep(na.omit(xi.matr[i,]), each = n-2) * na.omit(as.vector(t(xi.matr[-i,-i]))))
  }
  for (i in 1:n) {
    hat_g51_3[i] <- sum(rep(na.omit(xi.matr[,i]), each = n-2) * na.omit(as.vector((xi.matr[-i,-i]))))
  }
  hat_g51 <- (hat_g51_1 + hat_g51_2 + hat_g51_3)/(n-1)/(n-2)/3
  hat_g51 <- hat_g51 - mean(c(a_3,b_4) * c(b_4,a_3))
  sigma_square_51 <- mean((3*hat_g51 - 4*mean(xi.vec)*hat_g11)^2)  
  denominator_5 <- sqrt(sigma_square_51/n)
  
  if(sigma_square_51 > sqrt(log(n)/n)){
    
    result <- matrix(c(numerator_5, denominator_5, numerator_5/denominator_5, 2*(1-pnorm(abs(numerator_5/denominator_5))), 'non-degenerate'), ncol=5, byrow=TRUE)
    colnames(result) <- c('estimate','std.error','statistic','pvalue', 'linear part')
    rownames(result) <- c('sender-receiver effect')
    result <- as.table(result)
    
    list(estimate = numerator_5,
         std.error = denominator_5,
         statistic = numerator_5/denominator_5,
         pvalue = 2*(1-pnorm(abs(numerator_5/denominator_5))),
         sumry = result)
    
  }else{
    
    a_n <- 1/n/(n-1)
    b_n <- (n-2)/n/(n-1)
    c_n <- (n-2)*(n-3)/n/(n-1)
    
    single_kernel_ord4 <- function(adj_mat_4t4){
      
      h_stars.save <- rep(NA,6)
      
      xi.matr <- adj_mat_4t4
      diag(xi.matr) <- NA
      elements <- as.vector(t(xi.matr))
      xi.vec <- na.omit(elements)
      
      # h1_star
      h_stars.save[1] <- mean(xi.vec^2)
      
      # h2_star
      tmp_1 <- xi.matr
      tmp_2 <- xi.matr
      tmp_1[lower.tri(tmp_1)] <- NA
      tmp_2[upper.tri(tmp_2)] <- NA
      
      a_2_clt <- na.omit(as.vector(t(tmp_1)))
      b_2_clt <- na.omit(as.vector(tmp_2))
      
      h_stars.save[2] <- mean(a_2_clt * b_2_clt)
      
      # h3_star
      a_3 <- rep(xi.vec, each = 2)
      B.temp <- matrix(xi.vec, 3, 4)
      combine.temp_3 <- B.temp[-1,]
      for(i in 2:3){
        combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
      }
      b_3 <- as.vector(combine.temp_3)
      
      h_stars.save[3] <- mean(a_3 * b_3)
      
      # h4_star
      xi.vec.transp <- na.omit(as.vector(xi.matr))
      a_4 <- rep(xi.vec.transp, each = 2)
      D.temp <- matrix(xi.vec.transp, 3, 4)
      combine.temp_4 <- D.temp[-1,]
      for(i in 2:3){
        combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
      }
      b_4 <- as.vector(combine.temp_4)
      
      h_stars.save[4] <- mean(a_4 * b_4)
      
      # h5_star
      h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
      
      # h6_star
      e <- adj_mat_4t4
      h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                             e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                             e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
      
      h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
        b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
      
      # combined to be estimator's kernel
      single_kernel_ord4_ouput <- rep(NA,5)
      single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
      single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
      single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
      single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
      
      return(single_kernel_ord4_ouput) 
      
    }
    
    replctn <- floor(n^alpha)
    incomplete_Us <- matrix(NA,replctn,5)
    
    for (i in 1:replctn){
      ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
      adj_mat_dim4 <- adj_mat[ijkl,ijkl]
      incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
    }
    
    eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
    eta_est_mean <- apply(incomplete_Us,2,mean)
    eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
    
    result <- matrix(c(eta_est_mean[5], eta_deno[5], eta_est[5], 2*(1-pnorm(abs(eta_est[5]))), 'degenerate'), ncol=5, byrow=TRUE)
    colnames(result) <- c('estimate','std.error','statistic','pvalue', 'linear part')
    rownames(result) <- c('reciprocity effect')
    result <- as.table(result)
    
    list(estimate = eta_est_mean[5],
         std.error = eta_deno[5],
         statistic = eta_est[5],
         pvalue = 2*(1-pnorm(abs(eta_est[5]))),
         sumry = result)
  }
}



# **************************************************************************
# Test sender-receiver effect given that linear part degenerate
# Data input: alpha - sample n^alpha sub-network for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, and p-value
# **************************************************************************

NET_s.r_degen <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]
  
  a_n <- 1/n/(n-1)
  b_n <- (n-2)/n/(n-1)
  c_n <- (n-2)*(n-3)/n/(n-1)
  
  single_kernel_ord4 <- function(adj_mat_4t4){
    
    h_stars.save <- rep(NA,6)
    
    xi.matr <- adj_mat_4t4
    diag(xi.matr) <- NA
    elements <- as.vector(t(xi.matr))
    xi.vec <- na.omit(elements)
    
    # h1_star
    h_stars.save[1] <- mean(xi.vec^2)
    
    # h2_star
    tmp_1 <- xi.matr
    tmp_2 <- xi.matr
    tmp_1[lower.tri(tmp_1)] <- NA
    tmp_2[upper.tri(tmp_2)] <- NA
    
    a_2_clt <- na.omit(as.vector(t(tmp_1)))
    b_2_clt <- na.omit(as.vector(tmp_2))
    
    h_stars.save[2] <- mean(a_2_clt * b_2_clt)
    
    # h3_star
    a_3 <- rep(xi.vec, each = 2)
    B.temp <- matrix(xi.vec, 3, 4)
    combine.temp_3 <- B.temp[-1,]
    for(i in 2:3){
      combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
    }
    b_3 <- as.vector(combine.temp_3)
    
    h_stars.save[3] <- mean(a_3 * b_3)
    
    # h4_star
    xi.vec.transp <- na.omit(as.vector(xi.matr))
    a_4 <- rep(xi.vec.transp, each = 2)
    D.temp <- matrix(xi.vec.transp, 3, 4)
    combine.temp_4 <- D.temp[-1,]
    for(i in 2:3){
      combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
    }
    b_4 <- as.vector(combine.temp_4)
    
    h_stars.save[4] <- mean(a_4 * b_4)
    
    # h5_star
    h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
    
    # h6_star
    e <- adj_mat_4t4
    h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                           e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                           e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
    
    h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
      b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
    
    # combined to be estimator's kernel
    single_kernel_ord4_ouput <- rep(NA,5)
    single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
    single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
    single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
    single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
    
    return(single_kernel_ord4_ouput) 
    
  }
  
  replctn <- floor(n^alpha)
  incomplete_Us <- matrix(NA,replctn,5)
  
  for (i in 1:replctn){
    ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
    adj_mat_dim4 <- adj_mat[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
  }
  
  eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
  eta_est_mean <- apply(incomplete_Us,2,mean)
  eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
  
  result <- matrix(c(eta_est_mean[5], eta_deno[5], eta_est[5], 2*(1-pnorm(abs(eta_est[5])))), ncol=4, byrow=TRUE)
  colnames(result) <- c('estimate','std.error','statistic','pvalue')
  rownames(result) <- c('sender-receiver effect')
  result <- as.table(result)
  
  list(estimate = eta_est_mean[5],
       std.error = eta_deno[5],
       statistic = eta_est[5],
       pvalue = 2*(1-pnorm(abs(eta_est[5]))),
       sumry = result)
}



# **************************************************************************
# Test four effects
# Data input: alpha - for reciprocity and sender-receiver effect, 
#                     if linear part degenerates, sample n^alpha
#                     sub-network for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, p-value, 
#             and degenerate test result of linear part (C=1)
# **************************************************************************

NET_four_effects <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]
  
  a_n <- 1/n/(n-1)
  b_n <- (n-2)/n/(n-1)
  c_n <- (n-2)*(n-3)/n/(n-1)
  
  single_kernel_ord4 <- function(adj_mat_4t4){
    
    h_stars.save <- rep(NA,6)
    
    xi.matr <- adj_mat_4t4
    diag(xi.matr) <- NA
    elements <- as.vector(t(xi.matr))
    xi.vec <- na.omit(elements)
    
    # h1_star
    h_stars.save[1] <- mean(xi.vec^2)
    
    # h2_star
    tmp_1 <- xi.matr
    tmp_2 <- xi.matr
    tmp_1[lower.tri(tmp_1)] <- NA
    tmp_2[upper.tri(tmp_2)] <- NA
    
    a_2_clt <- na.omit(as.vector(t(tmp_1)))
    b_2_clt <- na.omit(as.vector(tmp_2))
    
    h_stars.save[2] <- mean(a_2_clt * b_2_clt)
    
    # h3_star
    a_3 <- rep(xi.vec, each = 2)
    B.temp <- matrix(xi.vec, 3, 4)
    combine.temp_3 <- B.temp[-1,]
    for(i in 2:3){
      combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
    }
    b_3 <- as.vector(combine.temp_3)
    
    h_stars.save[3] <- mean(a_3 * b_3)
    
    # h4_star
    xi.vec.transp <- na.omit(as.vector(xi.matr))
    a_4 <- rep(xi.vec.transp, each = 2)
    D.temp <- matrix(xi.vec.transp, 3, 4)
    combine.temp_4 <- D.temp[-1,]
    for(i in 2:3){
      combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
    }
    b_4 <- as.vector(combine.temp_4)
    
    h_stars.save[4] <- mean(a_4 * b_4)
    
    # h5_star
    h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
    
    # h6_star
    e <- adj_mat_4t4
    h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                           e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                           e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
    
    h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
      b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
    
    # combined to be estimator's kernel
    single_kernel_ord4_ouput <- rep(NA,5)
    single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
    single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
    single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
    single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
    
    return(single_kernel_ord4_ouput) 
    
  }
  
  replctn <- floor(n^alpha)
  incomplete_Us <- matrix(NA,replctn,5)
  
  for (i in 1:replctn){
    ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
    adj_mat_dim4 <- adj_mat[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
  }
  
  eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
  eta_est_mean <- apply(incomplete_Us,2,mean)
  eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
  
  result <- matrix(c(eta_est_mean[2], eta_deno[2], eta_est[2], 2*(1-pnorm(abs(eta_est[2]))), 'degenerate',
                     eta_est_mean[3], eta_deno[3], eta_est[3], 2*(1-pnorm(abs(eta_est[3]))), '',
                     eta_est_mean[4], eta_deno[4], eta_est[4], 2*(1-pnorm(abs(eta_est[4]))), '',
                     eta_est_mean[5], eta_deno[5], eta_est[5], 2*(1-pnorm(abs(eta_est[5]))), 'degenerate'), ncol=5, byrow=TRUE)
  colnames(result) <- c('estimate','std.error','statistic','pvalue', 'linear part')
  rownames(result) <- c('reciprocity effect','same sender effect','same receiver effect','sender-receiver effect')

  
  # reciprocity effect: test linear part degenerate or not
  xi.matr <- adj_mat
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements)
  
  # Var(e_ij)_mmt_estimation
  zero_dig <- adj_mat
  diag(zero_dig) <- 0
  
  hat_g11 <- rep(NA,n)
  for (i in 1:n) {
    hat_g11[i] <- sum(zero_dig[i,] + zero_dig[,i])/(n-1)/2
  }
  hat_g11 <- hat_g11 - mean(xi.vec)
  
  # test using concentration result
  tmp_1 <- xi.matr
  tmp_2 <- xi.matr
  tmp_1[lower.tri(tmp_1)] <- NA
  tmp_2[upper.tri(tmp_2)] <- NA
  
  a_2_clt <- na.omit(as.vector(t(tmp_1)))
  b_2_clt <- na.omit(as.vector(tmp_2))
  numerator_2 <- mean(a_2_clt * b_2_clt) - mean(xi.vec) * mean(xi.vec) 
  
  hat_g21 <- rep(NA,n)
  for (i in 1:n) {
    hat_g21[i] <- sum(zero_dig[i,] * zero_dig[,i])/(n-1)
  }
  hat_g21 <- hat_g21 - mean(a_2_clt * b_2_clt)
  sigma_square_21 <- mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11)^2)  
  denominator_2 <- sqrt(sigma_square_21/n)
  

  # sender-receiver effect: test linear part degenerate or not
  xi.matr <- adj_mat
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements)
  
  # Var(e_ij)_mmt_estimation
  zero_dig <- adj_mat
  diag(zero_dig) <- 0
  
  hat_g11 <- rep(NA,n)
  for (i in 1:n) {
    hat_g11[i] <- sum(zero_dig[i,] + zero_dig[,i])/(n-1)/2
  }
  hat_g11 <- hat_g11 - mean(xi.vec)
  
  # test using concentration result
  a_3 <- rep(xi.vec, each = n-2)
  xi.vec.transp <- na.omit(as.vector(xi.matr))
  D.temp <- matrix(xi.vec.transp, (n-1), n)
  combine.temp_4 <- D.temp[-1,]
  for(i in 2:(n-1)){
    combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
  }
  b_4 <- as.vector(combine.temp_4)
  numerator_5 <- mean(c(a_3,b_4) * c(b_4,a_3)) - mean(xi.vec) * mean(xi.vec)
  
  hat_g51_1 <- rep(NA,n)
  hat_g51_2 <- rep(NA,n)
  hat_g51_3 <- rep(NA,n)
  for (i in 1:n) {
    hat_g51_1[i] <- sum(c(a_3,b_4)[((i-1)*(n-1)*(n-2)+1):(i*(n-1)*(n-2))] * c(b_4,a_3)[((i-1)*(n-1)*(n-2)+1):(i*(n-1)*(n-2))])
  }
  for (i in 1:n) {
    hat_g51_2[i] <- sum(rep(na.omit(xi.matr[i,]), each = n-2) * na.omit(as.vector(t(xi.matr[-i,-i]))))
  }
  for (i in 1:n) {
    hat_g51_3[i] <- sum(rep(na.omit(xi.matr[,i]), each = n-2) * na.omit(as.vector((xi.matr[-i,-i]))))
  }
  hat_g51 <- (hat_g51_1 + hat_g51_2 + hat_g51_3)/(n-1)/(n-2)/3
  hat_g51 <- hat_g51 - mean(c(a_3,b_4) * c(b_4,a_3))
  sigma_square_51 <- mean((3*hat_g51 - 4*mean(xi.vec)*hat_g11)^2)  
  denominator_5 <- sqrt(sigma_square_51/n)
  

  if(sigma_square_51 > sqrt(log(n)/n)){
    result[4,] <- c(numerator_5, denominator_5, numerator_5/denominator_5, 2*(1-pnorm(abs(numerator_5/denominator_5))), 'non-degenerate')
  }
  if(sigma_square_21 > sqrt(log(n)/n)){
    result[1,] <- c(numerator_2, denominator_2, numerator_2/denominator_2, 2*(1-pnorm(abs(numerator_2/denominator_2))), 'non-degenerate')
  }
  
  list(estimate = result[,1],
         std.error = result[,2],
         statistic = result[,3],
         pvalue = result[,4],
         sumry = as.table(result))
}



# **************************************************************************
# Test four effects given that linear part degenerate for reciprocity 
# and sender-receiver effect
# Data input: alpha - sample n^alpha sub-network for reduced network moments
#             adj_mat - adjacency matrix of a network
# Output:     point estimate, std.error,test statistic, and p-value
# **************************************************************************

NET_four_effects_degen <-function(alpha = 1, adj_mat){
  
  diag(adj_mat) <- NA
  n <- dim(adj_mat)[1]
  
  a_n <- 1/n/(n-1)
  b_n <- (n-2)/n/(n-1)
  c_n <- (n-2)*(n-3)/n/(n-1)
  
  single_kernel_ord4 <- function(adj_mat_4t4){
    
    h_stars.save <- rep(NA,6)
    
    xi.matr <- adj_mat_4t4
    diag(xi.matr) <- NA
    elements <- as.vector(t(xi.matr))
    xi.vec <- na.omit(elements)
    
    # h1_star
    h_stars.save[1] <- mean(xi.vec^2)
    
    # h2_star
    tmp_1 <- xi.matr
    tmp_2 <- xi.matr
    tmp_1[lower.tri(tmp_1)] <- NA
    tmp_2[upper.tri(tmp_2)] <- NA
    
    a_2_clt <- na.omit(as.vector(t(tmp_1)))
    b_2_clt <- na.omit(as.vector(tmp_2))
    
    h_stars.save[2] <- mean(a_2_clt * b_2_clt)
    
    # h3_star
    a_3 <- rep(xi.vec, each = 2)
    B.temp <- matrix(xi.vec, 3, 4)
    combine.temp_3 <- B.temp[-1,]
    for(i in 2:3){
      combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
    }
    b_3 <- as.vector(combine.temp_3)
    
    h_stars.save[3] <- mean(a_3 * b_3)
    
    # h4_star
    xi.vec.transp <- na.omit(as.vector(xi.matr))
    a_4 <- rep(xi.vec.transp, each = 2)
    D.temp <- matrix(xi.vec.transp, 3, 4)
    combine.temp_4 <- D.temp[-1,]
    for(i in 2:3){
      combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
    }
    b_4 <- as.vector(combine.temp_4)
    
    h_stars.save[4] <- mean(a_4 * b_4)
    
    # h5_star
    h_stars.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3))
    
    # h6_star
    e <- adj_mat_4t4
    h_stars.save[6] <- ((e[1,2]*e[3,4]+e[1,2]*e[4,3]+e[2,1]*e[3,4]+e[2,1]*e[4,3]+
                           e[1,3]*e[2,4]+e[1,3]*e[4,2]+e[3,1]*e[2,4]+e[3,1]*e[4,2]+
                           e[1,4]*e[2,3]+e[1,4]*e[3,2]+e[4,1]*e[2,3]+e[4,1]*e[3,2])/12)
    
    h_star <- a_n*(h_stars.save[1] + h_stars.save[2]) + c_n*h_stars.save[6] +
      b_n*(h_stars.save[3] + h_stars.save[4] + 2*h_stars.save[5]) 
    
    # combined to be estimator's kernel
    single_kernel_ord4_ouput <- rep(NA,5)
    single_kernel_ord4_ouput[5] <- h_stars.save[5]-h_star
    single_kernel_ord4_ouput[4] <- h_stars.save[4]-h_star
    single_kernel_ord4_ouput[3] <- h_stars.save[3]-h_star
    single_kernel_ord4_ouput[2] <- h_stars.save[2]-h_star
    
    return(single_kernel_ord4_ouput) 
    
  }
  
  replctn <- floor(n^alpha)
  incomplete_Us <- matrix(NA,replctn,5)
  
  for (i in 1:replctn){
    ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
    adj_mat_dim4 <- adj_mat[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(adj_mat_dim4)
  }
  
  eta_est <- (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2))
  eta_est_mean <- apply(incomplete_Us,2,mean)
  eta_deno <- sqrt(apply(incomplete_Us,2,var))/n^(alpha/2)
  
  result <- matrix(c(eta_est_mean[2], eta_deno[2], eta_est[2], 2*(1-pnorm(abs(eta_est[2]))), 
                     eta_est_mean[3], eta_deno[3], eta_est[3], 2*(1-pnorm(abs(eta_est[3]))), 
                     eta_est_mean[4], eta_deno[4], eta_est[4], 2*(1-pnorm(abs(eta_est[4]))), 
                     eta_est_mean[5], eta_deno[5], eta_est[5], 2*(1-pnorm(abs(eta_est[5])))), ncol=4, byrow=TRUE)
  colnames(result) <- c('estimate','std.error','statistic','pvalue')
  rownames(result) <- c('reciprocity effect','same sender effect','same receiver effect','sender-receiver effect')
  
  list(estimate = result[,1],
       std.error = result[,2],
       statistic = result[,3],
       pvalue = result[,4],
       sumry = as.table(result))
}



# *******
# example
# *******

mat_example <- matrix(rnorm(10000,0,1),100,100) + matrix(rep(rnorm(100,0,1),100),100,100)
set.seed(68)
fit1 <- NET_four_effects_degen(1,mat_example)
fit1$sumry

#                           estimate   std.error   statistic      pvalue
# reciprocity effect     -0.05429095  0.06556567 -0.82803934  0.40764822
# same sender effect      0.94020604  0.09681378  9.71149028  0.00000000
# same receiver effect   -0.03954222  0.03334653 -1.18579712  0.23570243
# sender-receiver effect -0.03424310  0.04359755 -0.78543644  0.43219772


