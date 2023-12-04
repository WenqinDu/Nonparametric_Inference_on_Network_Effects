library(readr)
test_stat <- read_csv("b-c0.2-normal-Net1.8.csv")
# b-c0-normal-Net1
# b-c0-normal-Net1.2
# b-c0-normal-Net1.4
# b-c0-normal-Net1.6
# b-c0-normal-Net1.8

# b-c0.2-normal-Net1
# b-c0.2-normal-Net1.2
# b-c0.2-normal-Net1.4
# b-c0.2-normal-Net1.6
# b-c0.2-normal-Net1.8

apply(abs(test_stat) > qnorm(0.975,0,1), 2, sum) / 1000 

