# starting from 5th Feb
t_total <- 64
alpha_1 <- (1/0.86)^2
beta_1 <- alpha_1/5.1

# c_result from korea_c.r
korea_0.03   <- c(2.19, 1.08, 1.05, 1.05, 1.01, 1.002) # c = 1.608, 1.384, 1.038, 1.009
korea_0.01   <- c(2.39, 1.12, 1.05, 1.05, 1.01, 1.005) # c = 1.661, 1.47, 1.04, 1.01
korea_0.0025 <- c(2.59, 1.20, 1.05, 1.05, 1.01, 1.01) # c = 1.722, 1.605, 1.043, 1.016

c0123 <- korea_0.03
#c0123 <- korea_0.01
#c0123 <- korea_0.03

# compute c(t)
cut_off_0 <- 9
cut_off_1 <- 18 
cut_off_2 <- 27 # 2.27 policy 1
cut_off_3 <- 39 # 3.10 policy 2
cut_off_4 <- 51 # 3.22 policy 3

c_t <- rep(NA,(t_total-1))
c_t[1:cut_off_0] <- c0123[1]
c_t[(cut_off_0+1):cut_off_1] <- c0123[2]
c_t[(cut_off_1+1):cut_off_2] <- c0123[3]
c_t[(cut_off_2+1):cut_off_3] <- c0123[4]
c_t[(cut_off_3+1):cut_off_4] <- c0123[5]
c_t[(cut_off_4+1):(t_total-1)] <- c0123[6]

# calculate N_cum and N
start_n <- 4
N_total <- rep(NA,t_total)
N_total[1] <- 1
for (i in 2:t_total) {
  N_total[i] <- N_total[i-1]*c_t[i-1]
}

N <- rep(1, t_total)
for (i in 2:t_total) {
  N[i] <- N_total[i] - N_total[i-1]
}

N_cum_korea_0.0025 <- N_total[start_n:t_total]
N_korea_0.0025 <- N[start_n:t_total]

#N_cum_korea_0.01 <- N_total[start_n:t_total]
#N_korea_0.01 <- N[start_n:t_total]

#N_cum_korea_0.03 <- N_total[start_n:t_total]
#N_korea_0.03 <- N[start_n:t_total]

korea_N <- data.frame(N_cum_korea_0.0025,N_cum_korea_0.01,N_cum_korea_0.03
                      ,N_korea_0.0025,N_korea_0.01,N_korea_0.03)
write.csv(korea_N, file = "korea_N.csv")

# interpolating c(t)
c_t[5:14] <- seq(c_t[5],c_t[14],length.out = 10)
c_t[17:20] <- seq(c_t[17],c_t[20],length.out = 4)
c_t[26:29] <- seq(c_t[26],c_t[29],length.out = 4)
c_t[38:41] <- seq(c_t[38],c_t[41],length.out = 4)
c_t[50:53] <- seq(c_t[50],c_t[53],length.out = 4)
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
korea_0.0025 <- R0
#korea_0.01 <- R0
#korea_0.03 <- R0

# smooth curve
smooth_par <- 8
korea_0.01_smooth <- rep(NA,n-smooth_par)
korea_0.0025_smooth <- rep(NA,n-smooth_par)
korea_0.03_smooth <- rep(NA,n-smooth_par)
for (i in 1:(n-smooth_par)) {
  korea_0.01_smooth[i] <- mean(korea_0.01[i:(i+smooth_par-1)])
  korea_0.0025_smooth[i] <- mean(korea_0.0025[i:(i+smooth_par-1)])
  korea_0.03_smooth[i] <- mean(korea_0.03[i:(i+smooth_par-1)])
}
plot(korea_0.01_smooth,type = 'l')
lines(korea_0.0025_smooth)
lines(korea_0.03_smooth)

# save
korea_r0 <- data.frame(korea_0.0025_smooth, korea_0.01_smooth, korea_0.03_smooth)
write.csv(korea_r0, file = "korea_r0.csv")