# starting from 8th Feb
t_total <- 70
alpha_1 <- (1/0.86)^2
beta_1 <- alpha_1/5.1

# c_result from italy_c.r
italy_0.03 <-   c(1.67, 1.28, 1.05, 1.05, 1.05, 1.006) # average increase rate: c = 1.373, 1.04, 1.028
italy_0.01 <-   c(1.81, 1.26, 1.05, 1.05, 1.05, 1.01) # average increase rate: c = 1.424, 1.08, 1.035
italy_0.0025 <- c(1.99, 1.24, 1.06, 1.05, 1.05, 1.01) # average increase rate: c = 1.522, 1.1 , 1.047

c0123 <- italy_0.0025
#c0123 <- italy_0.01
#c0123 <- italy_0.03

# compute c(t)
cut_off_0 <- 17
cut_off_1 <- 34 # 3.5 policy 1
cut_off_2 <- 40 # 3.11 policy 2
cut_off_3 <- 50
cut_off_4 <- 60

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

N_cum_italy_0.0025 <- N_total[start_n:t_total]
N_italy_0.0025 <- N[start_n:t_total]

N_cum_italy_0.01 <- N_total[start_n:t_total]
N_italy_0.01 <- N[start_n:t_total]

N_cum_italy_0.03 <- N_total[start_n:t_total]
N_italy_0.03 <- N[start_n:t_total]

italy_N <- data.frame(N_cum_italy_0.0025,N_cum_italy_0.01,N_cum_italy_0.03
                      ,N_italy_0.0025,N_italy_0.01,N_italy_0.03)
write.csv(italy_N, file = "italy_N.csv")

# interpolating c(t)
c_t[14:21] <- seq(c_t[14],c_t[21],length.out = 8)
c_t[31:38] <- seq(c_t[31],c_t[38],length.out = 8)
c_t[39:42] <- seq(c_t[39],c_t[42],length.out = 4)
c_t[49:52] <- seq(c_t[49],c_t[52],length.out = 4)
c_t[59:62] <- seq(c_t[59],c_t[62],length.out = 4)
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
italy_0.0025 <- R0
#italy_0.01 <- R0
#italy_0.03 <- R0

# smooth curve
smooth_par <- 15
italy_0.01_smooth <- rep(NA,n-smooth_par)
italy_0.0025_smooth <- rep(NA,n-smooth_par)
italy_0.03_smooth <- rep(NA,n-smooth_par)
for (i in 1:(n-smooth_par)) {
  italy_0.01_smooth[i] <- mean(italy_0.01[i:(i+smooth_par-1)])
  italy_0.0025_smooth[i] <- mean(italy_0.0025[i:(i+smooth_par-1)])
  italy_0.03_smooth[i] <- mean(italy_0.03[i:(i+smooth_par-1)])
}
plot(italy_0.01_smooth,type = 'l')
lines(italy_0.0025_smooth)
lines(italy_0.03_smooth)

# save
italy_r0 <- data.frame(italy_0.0025_smooth, italy_0.01_smooth, italy_0.03_smooth)
write.csv(italy_r0, file = "italy_r0.csv")