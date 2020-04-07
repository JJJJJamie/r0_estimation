alpha_1 <- (1/0.86)^2
beta_1 <- alpha_1/5.1 # distribution 1
alpha_2 <- (1/0.45)^2
beta_2 <- alpha_2/18.8 # distribution 2
gamma_data <- rgamma(10000000,alpha_2,beta_2) + rgamma(10000000,alpha_1,beta_1)
h <- hist(gamma_data,breaks = seq(0,round(max(gamma_data))+1))
p <- h$density[1:t_total]
save(p,file = 'p')