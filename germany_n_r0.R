# starting from 25th Feb
t_total <- 53
alpha_1 <- (1/0.86)^2
beta_1 <- alpha_1/5.1

# c_result from germany_c.r
germany_0.0025   <- c(1.72, 1.71, 1.62, 1.04, 1.04, 1.04) # c = 1.608, 1.384, 1.038, 1.009
germany_0.01   <- c(1.66, 1.56, 1.56, 1.04, 1.04, 1.04) # c = 1.661, 1.47, 1.04, 1.01
germany_0.03 <- c(1.56, 1.52, 1.51, 1.04, 1.04, 1.04) # c = 1.722, 1.605, 1.043, 1.016

c0123 <- germany_0.0025
#c0123 <- germany_0.01
#c0123 <- germany_0.03

# compute c(t)
cut_off_0 <- 10
cut_off_1 <- 18 # 3.6 policy 1
cut_off_2 <- 26 # 3.14 policy 2
cut_off_3 <- 34 # 3.22 policy 3
cut_off_4 <- 42

c_t <- rep(NA,(t_total-1))
c_t[1:cut_off_0] <- c0123[1]
c_t[(cut_off_0+1):cut_off_1] <- c0123[2]
c_t[(cut_off_1+1):cut_off_2] <- c0123[3]
c_t[(cut_off_2+1):cut_off_3] <- c0123[4]
c_t[(cut_off_3+1):cut_off_4] <- c0123[5]
c_t[(cut_off_4+1):(t_total-1)] <- c0123[6]

# calculate N_cum and N
start_n <- 8
N_total <- rep(NA,t_total)
N_total[1] <- 1
for (i in 2:t_total) {
  N_total[i] <- N_total[i-1]*c_t[i-1]
}

N <- rep(1, t_total)
for (i in 2:t_total) {
  N[i] <- N_total[i] - N_total[i-1]
}

N_cum_germany_0.0025 <- N_total[start_n:t_total]
N_germany_0.0025 <- N[start_n:t_total]

#N_cum_germany_0.01 <- N_total[start_n:t_total]
#N_germany_0.01 <- N[start_n:t_total]

#N_cum_germany_0.03 <- N_total[start_n:t_total]
#N_germany_0.03 <- N[start_n:t_total]

germany_N <- data.frame(N_cum_germany_0.0025,N_cum_germany_0.01,N_cum_germany_0.03
                        ,N_germany_0.0025,N_germany_0.01,N_germany_0.03)
write.csv(germany_N, file = "germany_N.csv")

# interpolating c(t)
c_t[7:14] <- seq(c_t[7],c_t[14],length.out = 8)
c_t[17:20] <- seq(c_t[17],c_t[20],length.out = 4)
c_t[23:32] <- seq(c_t[23],c_t[32],length.out = 10)
c_t[33:36] <- seq(c_t[33],c_t[36],length.out = 4)
c_t[41:44] <- seq(c_t[41],c_t[44],length.out = 4)
plot(c_t)

# N_cum and N with interpolated c(t)
N_total <- rep(NA,t_total)
N_total[1] <- 1
for (i in 2:t_total) {
  N_total[i] <- N_total[i-1]*c_t[i-1]
}

N <- rep(1, t_total)
for (i in 2:t_total) {
  N[i] <- N_total[i] - N_total[i-1]
}

N <- N[start_n:t_total] # only keeping N > 10
n <- length(N)

# R0 estimation from N(t) # repeat and find average
n_repeat <- 100 # of repeats

array_I <- array(NA,dim = c(n_repeat,n,n))
for (nn in 1:n_repeat) {
  for (i in 1:n){
    array_I[nn,,i]<-rep(N[i],n)
  }
}

for (nn in 1:n_repeat) {
  for (i in 1:n){
    d_iso <- rep(N[i],n)
    iso <- ceiling(rgamma(N[i],alpha_1,beta_1))
    period_iso <- min(max(iso)+1,n)
    d_iso[1:period_iso] <- d_iso[1:period_iso] - N[i]
    for (k in 1:period_iso){
      d_iso[k:period_iso] <- d_iso[k:period_iso] + sum(iso == (k-1))
    }
    array_I[nn,,i] <- array_I[nn,,i] - d_iso
  }
}

I <- colSums(array_I)/n_repeat

cul_inf <- matrix(rep( 0, len=(n-1)*(n-1)), nrow = (n-1)) # I matrix
for (i in 1:(n-1)){
  cul_inf[i:(n-1),i] <- I[1:(n-i),i]
}

w <- cul_inf/rowSums(cul_inf) # weight matrix
R0 <- N[2:n]*w 
R0 <-colSums(R0)/N[2:n]

# save R0
germany_0.0025 <- R0
#germany_0.01 <- R0
#germany_0.03 <- R0

# smooth curve
smooth_par <- 15
germany_0.01_smooth <- rep(NA,n-smooth_par)
germany_0.0025_smooth <- rep(NA,n-smooth_par)
germany_0.03_smooth <- rep(NA,n-smooth_par)
for (i in 1:(n-smooth_par)) {
  germany_0.01_smooth[i] <- mean(germany_0.01[i:(i+smooth_par-1)])
  germany_0.0025_smooth[i] <- mean(germany_0.0025[i:(i+smooth_par-1)])
  germany_0.03_smooth[i] <- mean(germany_0.03[i:(i+smooth_par-1)])
}
plot(germany_0.01_smooth,type = 'l')
lines(germany_0.0025_smooth)
lines(germany_0.03_smooth)

# save
germany_r0 <- data.frame(germany_0.0025_smooth, germany_0.01_smooth, germany_0.03_smooth)
write.csv(germany_r0, file = "germany_r0.csv")