library(TripleR)
library(srm)
library(doParallel)

eta5_nondegen_concentration_test <- function(n, error_mari){
  
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
  ###'@Cov(e_ij,e_il)=eta_3_with_count=(n^2-n)(n-2)
  ###
  a_3 <- rep(xi.vec, each = n-2)
  
  
  ###
  ###'@Cov(e_ij,e_kj)=eta_4_with_count=(n^2-n)(n-2)
  ###
  xi.vec.transp <- na.omit(as.vector(xi.matr))
  D.temp <- matrix(xi.vec.transp, (n-1), n)
  combine.temp_4 <- D.temp[-1,]
  for(i in 2:(n-1)){
    combine.temp_4 <- rbind(combine.temp_4, D.temp[-i,])
  }
  b_4 <- as.vector(combine.temp_4)
  
  
  ###
  ###'@Cov(e_ij,e_ki)((e_ij,e_jk))
  ###
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
  
  
  ###
  ### Write Output 
  ###  
  return(list(t=numerator_5/denominator_5,sigma_square_51=sigma_square_51))
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
  
  list(test_stat = (apply(incomplete_Us,2,mean)/sqrt(apply(incomplete_Us,2,var))*n^(alpha/2)))
  
}


TripleR_test <- function(city_size,error_mari){
  
  net <- matrix(NA,city_size^2-city_size,3)
  
  net[,1] <- rep(c(1:city_size), each = city_size - 1)
  indx <- matrix(rep(c(1:city_size), each = city_size),city_size,city_size)
  diag(indx) <- NA
  net[,2] <- na.omit(as.vector(t(indx)))
  net[,3] <- na.omit(as.vector(t(error_mari)))
  
  dataframe_net <- as.data.frame(net)
  
  fit_TripleR <- RR(V3 ~ V1*V2, data = dataframe_net)
  return(fit_TripleR$varComp$p.value)
}


srm_test <- function(city_size,error_mari){
  
  net <- matrix(NA,city_size^2-city_size,3)
  
  net[,1] <- rep(c(1:city_size), each = city_size - 1)
  indx <- matrix(rep(c(1:city_size), each = city_size),city_size,city_size)
  diag(indx) <- NA
  net[,2] <- na.omit(as.vector(t(indx)))
  net[,3] <- na.omit(as.vector(t(error_mari)))
  
  dataframe_net <- as.data.frame(net)
  colnames(dataframe_net) <- c('Actor','Partner','values')
  
  mf <- '
  %Person
  F1@A =~ 1*values@A
  F1@P =~ 1*values@P
  values@A ~~ 0*values@A + 0*values@P
  values@P ~~ 0*values@P 
  
  %Dyad
  F1@AP =~ 1*values@AP
  F1@PA =~ 1*values@PA
  values@AP ~~ 0*values@AP + 0*values@PA
  values@PA ~~ 0*values@PA 
  '
  
  mod1 <- srm::srm(mf, data = dataframe_net, conv_par=1e-4, maxiter=20)
  return(2*(1-pnorm(abs(mod1$coef/mod1$se))))
}


st_a_alt_normal <- function(n,C){

  a_i.save <- rnorm(n, 1, 1)
  
  ###
  ### Generate_Epsilon_ij
  ###
  epsilon_ij.save <- matrix(rnorm(n^2, 0, 1), n, n)  
  
  
  ###
  ### Generate Error Terms
  ###
  e_ij.save <- C * matrix(rep(a_i.save, n), n, n) +  epsilon_ij.save
  diag(e_ij.save) <- NA
  
  ###
  ### Write Output 
  ###  
  list(e_ij.save = e_ij.save)
}  


num_cores <- 20 
registerDoParallel(cores=num_cores)  


#'@_Simulation
city_size <- 50
C_list <- sqrt(c(0.05,0.2,0.5,1,5)) 
replication <- 1000
test_stat <- matrix(NA,replication,5)


coln <- 0
set.seed(68)
for (c_val in C_list){
  coln <- coln + 1

  r3 <- foreach(i=1:replication, .combine=rbind) %dopar% {
    
    #'@_e_ij
    Error_term.example <- st_a_alt_normal(n = city_size, C = c_val)
    
    #'@_Our_method_eta3
    studentuzed_t <- incomplete_Us_standz(city_size, alpha=1, Error_term.example$e_ij.save)
    
    studentuzed_t$test_stat[3]
  }

  test_stat[,coln] <- r3
}

write.csv(test_stat,file='b-alt-normal-Net1-n50.csv', row.names=FALSE)
 

