library(readr)
test_stat <- read_csv("outputs/a-alt-normal-Net1-n50.csv")
# please iterate over all following configurations 
# (and iterate over Net1, Net1.2, Net1.6)

# a-alt-normal-Net1-n50.csv
# a-alt-normal-Net1-n100.csv
# a-alt-pois-Net1-n50.csv
# a-alt-pois-Net1-n100.csv

# b-alt-normal-Net1-n50.csv
# b-alt-normal-Net1-n100.csv
# b-alt-normal-SRMA-n50.csv
# b-alt-normal-SRMA-n100.csv
# b-alt-normal-SRML-n50.csv

# b-alt-pois-Net1-n50.csv
# b-alt-pois-Net1-n100.csv
# b-alt-pois-SRMA-n50.csv
# b-alt-pois-SRMA-n100.csv
# b-alt-pois-SRML-n50.csv

# c-alt-normal-Net1-n50.csv
# c-alt-normal-Net1-n100.csv
# c-alt-normal-SRMA-n50.csv
# c-alt-normal-SRMA-n100.csv
# c-alt-normal-SRML-n50.csv

# c-alt-pois-Net1-n50.csv
# c-alt-pois-Net1-n100.csv
# c-alt-pois-SRMA-n50.csv
# c-alt-pois-SRMA-n100.csv
# c-alt-pois-SRML-n50.csv


# ***********************
# empirical powers of NET
# ***********************

test_stat <- data.matrix(test_stat)
apply(abs(test_stat) > qnorm(0.975,0,1), 2, sum) / 1000 
# p-values for c=sqrt(c(0.05,0.2,0.5,1,5)) 


# *********************************
# empirical powers of SRMA and SRML
# *********************************

test_stat <- data.matrix(test_stat)
round(apply(test_stat < 0.05, 2, sum) / apply(test_stat != 666, 2, sum),3)
# p-values for c=sqrt(c(0.05,0.2,0.5,1,5)) 

apply(test_stat == 666, 2, sum)
# number of unsuccessfully compiled tests applying SRMA and SRML




