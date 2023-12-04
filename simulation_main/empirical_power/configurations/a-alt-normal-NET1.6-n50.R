library(doParallel)

eta2_nondegen_concentration_test <- function(n, error_mari){
  
  ###
  ### Data Manipulation
  ###
  xi.matr <- error_mari
  diag(xi.matr) <- NA
  elements <- as.vector(t(xi.matr))
  xi.vec <- na.omit(elements)
  
  ###
  ###'@Var(e_ij)_mmt_estimation
  ###
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
  numerator_2 <- mean(a_2_clt * b_2_clt) - mean(xi.vec) * mean(xi.vec) 
  
  hat_g21 <- rep(NA,n)
  for (i in 1:n) {
    hat_g21[i] <- sum(zero_dig[i,] * zero_dig[,i])/(n-1)
  }
  hat_g21 <- hat_g21 - mean(a_2_clt * b_2_clt)
  sigma_square_21 <- mean((2*hat_g21 - 4*mean(xi.vec)*hat_g11)^2)  
  denominator_2 <- sqrt(sigma_square_21/n)
  
  ###
  ### Write Output 
  ###  
  return(list(t=numerator_2/denominator_2,sigma_square_21=sigma_square_21))
  
}


incomplete_Us_NO_debias_standz <- function(n, alpha, error_mari){
  
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
    # set.seed(i)  # fixed design/samplying
    ijkl <- sort(sample(1:n, 4, replace = FALSE, prob = NULL))                
    error_mari_dim4 <- error_mari[ijkl,ijkl]
    incomplete_Us[i,] <- single_kernel_ord4(error_mari_dim4)
  }
  
  list(test_stat = (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2)))
  
}


WE_alt_normal <- function(n,C){
  
  ###
  ###'@Generate_Gamma_ij__Paras:a1,b1
  ###
  gamma_ij.save <- matrix(rnorm(n^2, 1, 1), n, n)
  gamma_ij.save[upper.tri(gamma_ij.save)] <- t(gamma_ij.save)[upper.tri(gamma_ij.save)]
  diag(gamma_ij.save) <- 0
  
  a_i.save <- rnorm(n, 1, 1)
  b_j.save <- rnorm(n, 1, 1)
  
  ###
  ###'@Generate_Epsilon_ij__Paras:a2,b2
  ###
  epsilon_ij.save <- matrix(rnorm(n^2, 0, 1), n, n)  
  diag(epsilon_ij.save) <- 0
  
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- matrix(rep(a_i.save, n), n, n) + 
    matrix(rep(b_j.save, n), n, n, byrow = TRUE) + C * gamma_ij.save + epsilon_ij.save
  diag(e_ij.save) <- NA
  
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save)
} 



num_cores <- 20 #  detectCores() - 4  
registerDoParallel(cores=num_cores)  
# cl <- makeCluster(num_cores, type="FORK")  


#'@_Simulation
city_size <- 50
C_list <- sqrt(c(0.05,0.2,0.5,1,5)) 
replication <- 1000
test_stat <- matrix(NA,replication,5)


coln <- 0
set.seed(68)
for (c_val in C_list){
  coln <- coln + 1

  r2 <- foreach(i=1:replication, .combine=rbind) %dopar% {
    
    #'@_e_ij
    Error_term.example <- WE_alt_normal(n = city_size, C = c_val)

    #'@_test_stat_of_eta2
    nondegentest <- eta2_nondegen_concentration_test(city_size, Error_term.example$e_ij.save)

    if(nondegentest$sigma_square_21 > sqrt(log(city_size)/city_size) ){
      nondegentest$t
    }else if(nondegentest$sigma_square_21 <= sqrt(log(city_size)/city_size)){
      studentuzed_t <- incomplete_Us_NO_debias_standz(city_size, alpha=1.6, Error_term.example$e_ij.save)
      studentuzed_t$test_stat[2]
    }
  }
  
  test_stat[,coln] <- r2
}

write.csv(test_stat,file='a-alt-normal-Net1.6-n50.csv', row.names=FALSE)
# stopCluster(cl)  

