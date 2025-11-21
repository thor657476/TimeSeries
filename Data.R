# Timo Stuij 612337
# Pien Küthe 648127
# Job van Onna 638473
# Thor Hogerbrugge 657476

# Data inladen 
library(readr)
df <- read_csv("growth_data.csv")

T <- nrow(df)
num_holdout <- 12

# Alle data, behalve laatste twaalf voor estimation
estimation_data <- df[1:(T - num_holdout), ]
holdout_data <- df[(T - num_holdout + 1):T, ]

