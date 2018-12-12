#Example of how to calculate a Z-score

setwd('#####') #set your working directory

#DATA LOAD
data <- read.csv('Height_data.csv')
height <- data$Height
hist(height) #histogram

#population parameter calculations
pop_sd <- sd(height)*sqrt((length(height)-1)/(length(height)))
pop_mean <- mean(height)

#z-score calculation
z <- (72 - pop_mean)/pop_sd

p_yellow1 <- pnorm(72, pop_mean, pop_sd)    #using x, mu, and sigma
p_yellow2 <- pnorm(z)                       #using z-score of 2.107

p_blue1 <- 1 - p_yellow1
p_blue2 <- 1 - p_yellow2

p_blue1
p_blue2

#if x is a vector with raw scores then scale(x) is a vector with standardized scores.
(x-mean(x))/sd(x)

x <- c(1,3,4,5,7)
scale(x)
