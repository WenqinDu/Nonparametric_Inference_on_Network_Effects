library(readr)

# **********************************************
# empirical type-I error rates for SRMA and SRML
# **********************************************

test_stat <- read_csv("outputs/b-null-normal-SRMA.csv")
# please iterate over all following configurations
# b-null-normal-SRMA.csv
# b-null-normal-SRML.csv
# b-null-pois-SRMA.csv
# b-null-pois-SRML.csv
# c-null-normal-SRMA.csv
# c-null-normal-SRML.csv
# c-null-pois-SRMA.csv
# c-null-pois-SRML.csv

test_stat <- data.matrix(test_stat)
round(apply(test_stat < 0.05, 2, sum) / apply(test_stat != 666, 2, sum),3)
# p-values for testing \eta_3 and \eta_5 for n=25,50,100,200,400(SRMA) 
# p-values for testing \eta_3 and \eta_5 for n=25,50(SRML) 

apply(test_stat == 666, 2, sum)
# number of unsuccessfully compiled tests applying SRMA and SRML


# ************************************
# empirical type-I error rates for NET
# ************************************

test_stat <- read_csv("outputs/a-null-pois-Net1.csv")
# please iterate over all following configurations (for 1, 1.2, 1.6)
# a-null-normal-Net1.csv
# a-null-pois-Net1.csv
# b-null-normal-Net1.csv
# b-null-pois-Net1.csv
# c-null-normal-Net1.csv
# c-null-pois-Net1.csv

# for setting (a)
test_stat <- data.matrix(test_stat)
apply(abs(test_stat) > qnorm(0.975,0,1), 2, sum) / 1000 
# p-values for testing \eta_2 for n=25,50,100,200,400(NET) 

# for setting (b) and (c)
test_stat <- data.matrix(test_stat)
apply(abs(test_stat) > qnorm(0.975,0,1), 2, sum) / 1000 
# p-values for testing \eta_3 and \eta_5 for n=25,50,100,200,400(NET) 


# ****************
# QQ-plots for NET
# ****************

# use *-null-normal-Net*.csv for normal configurations
# use *-null-pois-Net*.csv for poisson configurations

test_stat <- read_csv("outputs/a-null-normal-Net1.csv")
test_stat <- data.matrix(test_stat)

# NET test eta_2 qq-plot
qqnorm(test_stat[,1], pch = 16, col = "pink", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), 
       main = expression(paste("(a) ", lambda, "=1.2")),
       cex.lab=1.4, cex.axis=1.3, cex.main=1.6)
par(new=TRUE)
qqnorm(test_stat[,2], pch = 16, col = "red", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "", xlab="", ylab="", xaxt="n", yaxt="n")  # xaxt="n", yaxt="n"
par(new=TRUE)
qqnorm(test_stat[,3], pch = 16, col = "orange", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "", xlab="", ylab="", xaxt="n", yaxt="n")
par(new=TRUE)
qqnorm(test_stat[,4], pch = 16, col = "blue", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "", xlab="", ylab="", xaxt="n", yaxt="n")
par(new=TRUE)
qqnorm(test_stat[,5], pch = 16, col = "purple", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "", xlab="", ylab="", xaxt="n", yaxt="n")
abline(0,1, lty = 2, lwd = 1.5)
legend(1, 0.2, legend=c("n=25", "n=50", "n=100", "n=200", "n=400"),
       col=c("pink", "red", "orange", "blue", "purple"), pch = 16, cex=0.86)

# NET test eta_3 qq-plot
qqnorm(test_stat[,1], pch = 16, col = "pink", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "(a)")
par(new=TRUE)
qqnorm(test_stat[,3], pch = 16, col = "red", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "") 
par(new=TRUE)
qqnorm(test_stat[,5], pch = 16, col = "orange", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "")
par(new=TRUE)
qqnorm(test_stat[,7], pch = 16, col = "blue", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "")
par(new=TRUE)
qqnorm(test_stat[,9], pch = 16, col = "purple", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "")
abline(0,1, lty = 2, lwd = 1.5)
legend(1, 0.2, legend=c("n=25", "n=50", "n=100", "n=200", "n=400"),
       col=c("pink", "red", "orange", "blue", "purple"), pch = 16, cex=0.86)

# NET test eta_5 qq-plot
qqnorm(test_stat[,2], pch = 16, col = "pink", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "(b)")
par(new=TRUE)
qqnorm(test_stat[,4], pch = 16, col = "red", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "") 
par(new=TRUE)
qqnorm(test_stat[,6], pch = 16, col = "orange", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "")
par(new=TRUE)
qqnorm(test_stat[,8], pch = 16, col = "blue", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "")
par(new=TRUE)
qqnorm(test_stat[,10], pch = 16, col = "purple", frame = FALSE, cex=0.4, xlim = c(-3,3), ylim = c(-3,3), main = "")
abline(0,1, lty = 2, lwd = 1.5)






