library(readxl)

# prepare data
Predictive_death <- read_excel("Downloads/Mean Predictive  death per day 020420.xlsx")
deaths<- Predictive_death$`Korea, South`
deaths <- deaths[32:length(deaths)] # only keeping daily new deaths > 10 # 19th March - 30th March
deaths[38:42] <- c(4,3,4,5,3) # data correction
deaths <- deaths[1:42]
t_deaths <- length(deaths)
t_before_fd <- 22 # 1st feb
t_total <- t_deaths + t_before_fd

n_c <- 6

cut_off_0 <- 9
cut_off_1 <- 18 
cut_off_2 <- 27 # 2.27 policy 1
cut_off_3 <- 39 # 3.10 policy 2
cut_off_4 <- 51 # 3.22 policy 3

alpha_1 <- (1/0.86)^2
beta_1 <- alpha_1/5.1 # distribution infected to onset/isolation

load('p') # death delay pdf


# death rates
dr <- 0.01# dr <- seq(dr_min,dr_max, length.out = t_total)


# manually tune the grid by changing range and density
# better method to be proposed
nc0 <- 11
nc1 <- 11
nc2 <- 11
nc3 <- 11
nc4 <- 11
nc5 <- 11

c0 <- seq(1,2,length.out = nc0)
c1 <- seq(1,2,length.out = nc1)
c2 <- seq(1,2,length.out = nc2)
c3 <- seq(1,2,length.out = nc3)
c4 <- seq(1,2,length.out = nc4)
c5 <- seq(1,2,length.out = nc5)


# build a grid of c's
c_total <- nc0*nc1*nc2*nc3*nc4*nc5
N <- 0
a <- array(NA, dim = c(c_total,n_c))
for (n in 1:nc0) {
  for (i in 1:nc1) {
    for (j in 1:nc2) {
      for (k in 1:nc3) {
        for (m in 1:nc4) {
          for (q in 1:nc5) {
            N <- N+1
            a[N,] <- c(c0[n],c1[i],c2[j],c3[k],c4[m],c5[q])
          }
        }
      }
    }
  }
}

a <- a[which(a[,2] - a[,1] <= 0),] # remove c1 > c0
a <- a[which(a[,3] - a[,2] <= 0),] # remove c2 > c1
a <- a[which(a[,4] - a[,3] <= 0),] # remove c3 > c2
a <- a[which(a[,5] - a[,4] <= 0),] # remove c4 > c3
a <- a[which(a[,6] - a[,5] <= 0),] # remove c5 > c4
c_total <- dim(a)[1]


# find c1,c2,c3 with MSE
pb <- txtProgressBar(max = c_total, style = 3)
se <- rep(NA,t_total)
for (nn in 1:c_total) {
  c0123 <- a[nn,]
  c_t <- rep(NA,(t_total-1))
  c_t[1:cut_off_0] <- c0123[1]
  c_t[(cut_off_0+1):cut_off_1] <- c0123[2]
  c_t[(cut_off_1+1):cut_off_2] <- c0123[3]
  c_t[(cut_off_2+1):cut_off_3] <- c0123[4]
  c_t[(cut_off_3+1):cut_off_4] <- c0123[5]
  c_t[(cut_off_4+1):(t_total-1)] <- c0123[6]
  
  # Generate N_total
  N_total <- rep(NA,t_total)
  N_total[1] <- 1
  for (i in 2:t_total) {
    N_total[i] <- N_total[i-1]*c_t[i-1]
  }
  
  # Calculate N(t)
  N <- rep(1, t_total)
  for (i in 2:t_total) {
    N[i] <- N_total[i] - N_total[i-1]
  }
  
  # Cauculate d(t)
  d <- rep(NA,t_deaths)
  for (i in 1:t_deaths) {
    n_terms <- t_before_fd + i
    d[i] <- sum(dr*p[n_terms:1]*N[1:n_terms]) # change to dr while using varying dr
  }
  
  # Calculate square error
  d_error <- deaths - d
  se[nn] <- sum(d_error^2)
  setTxtProgressBar(pb, nn)
}
close(pb)

# daily death error
sqrt(min(se)/(t_deaths - bond_n + 1))/bond_n

# c values
a[which(se == min(se)),]
