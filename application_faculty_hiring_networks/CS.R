library(readr)
library(dplyr)

estimate_eta_use_eij <- function(n, error_mari){
  
  ###
  ### Data Manipulation
  ###
  
  xi.matr <- error_mari
  diag(xi.matr) <- NA
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements)
  
  lam_est.save <- rep(NA,10)
  
  ###
  ###'@Var(e_ij)_mmt_estimation
  ###
  lam_est.save[1] <- mean(xi.vec^2) - mean(xi.vec)*mean(xi.vec)  
  
  zero_dig <- error_mari
  diag(zero_dig) <- 0
  
  hat_g11 <- rep(NA,n)
  for (i in 1:n) {
    hat_g11[i] <- sum(zero_dig[i,] + zero_dig[,i])/(n-1)/2
  }
  hat_g11 <- hat_g11 - mean(xi.vec)
  
  ###
  ###'@Cov(e_ij,e_ji)=eta_2_with_count=(n^2-n)/2_(for_independent_CLT)
  ###
  tmp_1 <- xi.matr
  tmp_2 <- xi.matr
  tmp_1[lower.tri(tmp_1)] <- NA
  tmp_2[upper.tri(tmp_2)] <- NA
  
  a_2_clt <- na.omit(as.vector(t(tmp_1)))
  b_2_clt <- na.omit(as.vector(tmp_2))
  lam_est.save[2] <- mean(a_2_clt * b_2_clt) - mean(xi.vec) * mean(xi.vec) 
  
  
  hat_g21 <- rep(NA,n)
  for (i in 1:n) {
    hat_g21[i] <- sum(zero_dig[i,] * zero_dig[,i])/(n-1)
  }
  hat_g21 <- hat_g21 - mean(a_2_clt * b_2_clt)
  sigma_square_21 <- mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11)^2)  
  lam_est.save[7] <- sqrt(sigma_square_21/n)
  
  ###
  ###'@Cov(e_ij,e_il)=eta_3_with_count=(n^2-n)(n-2)
  ###
  a_3 <- rep(xi.vec, each = n-2)
  B.temp <- matrix(xi.vec, (n-1), n)
  combine.temp_3 <- B.temp[-1,]
  for(i in 2:(n-1)){
    combine.temp_3 <- rbind(combine.temp_3, B.temp[-i,])
  }
  b_3 <- as.vector(combine.temp_3)
  lam_est.save[3] <- mean(a_3 * b_3) - mean(xi.vec) * mean(xi.vec)  
  
  ###
  ###'@Cov(e_ij,e_kj)=eta_4_with_count=(n^2-n)(n-2)
  ###
  xi.vec.transp <- na.omit(as.vector(xi.matr))
  a_4 <- rep(xi.vec.transp, each = n-2)
  D.temp <- matrix(xi.vec.transp, (n-1), n)
  combine.temp_4 <- D.temp[-1,]
  for(i in 2:(n-1)){
    combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
  }
  b_4 <- as.vector(combine.temp_4)
  lam_est.save[4] <- mean(a_4 * b_4) - mean(xi.vec) * mean(xi.vec)  
  
  
  ###
  ###'@Cov(e_ij,e_ki)((e_ij,e_jk))
  ###
  lam_est.save[5] <- mean(c(a_3,b_4) * c(b_4,a_3)) - mean(xi.vec) * mean(xi.vec)  
  
  
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
  lam_est.save[10] <- sqrt(sigma_square_51/n)
  
  
  ###
  ### Write Output 
  ###  
  list(lam_est.save=lam_est.save, sigma_square_21=sigma_square_21, sigma_square_51=sigma_square_51) 
  
}

incomplete_Us_standz <- function(n, alpha, error_mari){
  
  diag(error_mari) <- NA
  
  a_n <- 1/n/(n-1)
  b_n <- (n-2)/n/(n-1)
  c_n <- (n-2)*(n-3)/n/(n-1)
  
  single_kernel_ord4 <- function(error_mari_4t4){
    
    h_stars.save <- rep(NA,6)
    
    xi.matr <- error_mari_4t4
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
    e <- error_mari_4t4
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
    error_mari_dim4 <- error_mari[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(error_mari_dim4)
  }
  
  list(eta_est = (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2)),
       eta_est_mean = apply(incomplete_Us,2,mean), eta_deno = sqrt(apply(incomplete_Us,2,var))/n^(alpha/2))
  
}


mat2 <- as.matrix(read_csv("ready_to_use_data_matrices/CS General Adjacency Matrix.csv",
                           col_names = FALSE))

colnames(mat2) <- NULL
mat2 <- unname(mat2)
Bs112_mat <- mat2
diag(Bs112_mat) <- NA


#'@_Test_eta_2,5_when_g2,5=\=0,_respectively

city_size <- 205 

etas_reslt <- estimate_eta_use_eij(city_size, Bs112_mat)

sqrt(log(city_size)/city_size)  
etas_reslt$sigma_square_21  # FTR linear part degenerate when testing reciprocity effect
etas_reslt$sigma_square_51  # FTR linear part degenerate when testing sender-receiver effect

etas <- etas_reslt$lam_est.save
etas[2:5]  # estimates of network effects


#'@_Test_eta_2,3,4,5,_respectively

replication <- 10000
check_2 <- rep(NA,replication)
check_3 <- rep(NA,replication)
check_4 <- rep(NA,replication)
check_5 <- rep(NA,replication)

set.seed(6)
for (i in 1:replication){
  
  print(paste(city_size,"-", i))
  
  #'@_e_ij
  etas_estimation <- incomplete_Us_standz(city_size, alpha=1.2, Bs112_mat)
  
  #'@-----------------eta3_estimates
  check_2[i] <- etas_estimation$eta_est[2]
  check_3[i] <- etas_estimation$eta_est[3]
  check_4[i] <- etas_estimation$eta_est[4]
  check_5[i] <- etas_estimation$eta_est[5]
}

# Z-average
2*(1-pnorm(abs(1*mean(check_2))))  # p-value of testing reciprocity effect
2*(1-pnorm(abs(1*mean(check_3))))  # p-value of testing the same sender effect
2*(1-pnorm(abs(1*mean(check_4))))  # p-value of testing the same receiver effect
2*(1-pnorm(abs(1*mean(check_5))))  # p-value of testing the sender receiver effect


