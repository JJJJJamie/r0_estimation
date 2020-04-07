# starting from 20th Feb
t_total <- 62
alpha_1 <- (1/0.86)^2
beta_1 <- alpha_1/5.1

# c_result from uk_c.r
uk_0.0025 <- c(1.73, 1.65, 1.35, 1.27, 1.03, 1.03)
uk_0.01   <- c(1.60, 1.60, 1.32, 1.26, 1.03, 1.03)
uk_0.03   <- c(1.54, 1.51, 1.30, 1.26, 1.03, 1.03)

c0123 <- uk_0.0025
#c0123 <- uk_0.01
#c0123 <- uk_0.03

# compute c(t)
cut_off_0 <- 10
cut_off_1 <- 20
cut_off_2 <- 29 # 3.12 policy 1
cut_off_3 <- 41 # 3.24 policy 2
cut_off_4 <- 51

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

N_cum_uk_0.0025 <- N_total[start_n:t_total]
N_uk_0.0025 <- N[start_n:t_total]

#N_cum_uk_0.01 <- N_total[start_n:t_total]
#N_uk_0.01 <- N[start_n:t_total]

#N_cum_uk_0.03 <- N_total[start_n:t_total]
#N_uk_0.03 <- N[start_n:t_total]

uk_N <- data.frame(N_cum_uk_0.0025,N_cum_uk_0.01,N_cum_uk_0.03
                      ,N_uk_0.0025,N_uk_0.01,N_uk_0.03)
write.csv(uk_N, file = "uk_N.csv")

# interpolating c(t)
c_t[9:12] <- seq(c_t[9],c_t[12],length.out = 4)
c_t[18:23] <- seq(c_t[18],c_t[23],length.out = 6)
c_t[28:31] <- seq(c_t[28],c_t[31],length.out = 4)
c_t[39:44] <- seq(c_t[39],c_t[44],length.out = 6)
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
uk_0.0025 <- R0
#uk_0.01 <- R0
#uk_0.03 <- R0

# smooth curve
smooth_par <- 10
uk_0.01_smooth <- rep(NA,n-smooth_par)
uk_0.0025_smooth <- rep(NA,n-smooth_par)
uk_0.03_smooth <- rep(NA,n-smooth_par)
for (i in 1:(n-smooth_par)) {
  uk_0.01_smooth[i] <- mean(uk_0.01[i:(i+smooth_par-1)])
  uk_0.0025_smooth[i] <- mean(uk_0.0025[i:(i+smooth_par-1)])
  uk_0.03_smooth[i] <- mean(uk_0.03[i:(i+smooth_par-1)])
}
plot(uk_0.01_smooth,type = 'l')
lines(uk_0.0025_smooth)
lines(uk_0.03_smooth)

# save
uk_r0 <- data.frame(uk_0.0025_smooth, uk_0.01_smooth, uk_0.03_smooth)
write.csv(uk_r0, file = "uk_r0.csv")