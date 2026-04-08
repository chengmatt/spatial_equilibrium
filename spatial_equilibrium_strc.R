# Age Structured Extension ----------------------------------------------------------
# Dimensions
n_regions = 2
n_ages = 15

# Natural Mortality
natmort = array(c(0.2, 0.25), dim = c(n_regions, n_ages))

# Mean recruitment
mean_rec_r = c(0.6, 0.4)

# Movement
T_mat = {
  rate = 0.3
  arr = array(0, dim = c(n_regions, n_regions, n_ages)) # structured as from, to, n_ages, n_sexes
  arr[1,1,] = rate
  arr[1,2,] = 1 - rate
  arr[2,2,] = 1 - rate
  arr[2,1,] = rate
  arr
}

# Container
n = array(0, dim = c(n_regions, n_ages))

# projection initial abundance forward
n = array(0, dim = c(n_regions, n_ages))
for(i in 1:n_ages) {
  n[,1] = mean_rec_r
  for(a in 1:n_ages) n[,a] = t(T_mat[,,a]) %*% n[,a]
  for(r in 1:n_regions) {
    n[r, 2:n_ages] = n[r, 1:(n_ages-1)] * exp(-natmort[r, 1:(n_ages-1)])
  }
}

# now apply analytical plus group on top of penultimate age
Move_penult  = t(T_mat[,,n_ages-1])
Move_plus    = t(T_mat[,,n_ages])
s_penult     = exp(-natmort[, n_ages-1])
s_plus       = exp(-natmort[, n_ages])
init_penult  = n[, n_ages-1]
source       = (Move_penult %*% init_penult) * s_penult
O_mat        = diag(s_plus, n_regions) %*% Move_plus
I_mat        = diag(n_regions)
n[, n_ages] = solve(I_mat - O_mat, source)

# look at eigens
max(abs(eigen(O_mat)$values))

# iterative approach
n_iter = array(0, dim = c(n_regions, n_ages))

# project forward iteratively until convergence
for(i in 1:500) {
  n_iter[,1] = mean_rec_r
  # movement
  for(a in 1:n_ages) n_iter[,a] = t(T_mat[,,a]) %*% n_iter[,a]
  for(r in 1:n_regions) {
    # save plus group survivors
    plus_survivors = n_iter[r, n_ages] * exp(-natmort[r, n_ages])
    # ageing and mortality - this overwrites n_ages with inflow from penultimate
    n_iter[r, 2:n_ages] = n_iter[r, 1:(n_ages-1)] * exp(-natmort[r, 1:(n_ages-1)])
    # accumulate plus group
    n_iter[r, n_ages] = n_iter[r, n_ages] + plus_survivors
  }
}

# plot to compare
par(mfrow=c(1,2))
plot(n[1,], type="l", xlab="Age bin", ylab="N", main="Region 1", ylim = c(0, 1), lwd = 5)
lines(n_iter[1,], col="blue", lty=3, lwd = 5)
legend("topright", c("Analytical","Iterative"), col=c("black","blue"), lty=c(1,3), lwd = 5)

plot(n[2,], type="l", xlab="Age bin", ylab="N", main="Region 2", ylim = c(0, 1), lwd = 5)
lines(n_iter[2,], col="red", lty=3, lwd = 5)
legend("topright", c("Analytical","Iterative"), col=c("black","red"), lty=c(1,3), lwd = 5)


# Size Structured Extension ---------------------------------------------------------
n_region = 2
n_size   = 15

# growth transition
build_Xr = function(n, mean_growth, cv_growth=0.5) {
  Xr = matrix(0, n, n)
  for (i in 1:(n-1)) {
    mean_inc = mean_growth * (1 - 0.7*(i-1)/(n-1)) # assuming linear growth for simplicity
    sd_inc   = mean_inc * cv_growth
    shape    = (mean_inc/sd_inc)^2
    rate     = mean_inc/sd_inc^2
    probs    = rep(0, n)
    for (j in i:n) {
      probs[j] = pgamma(j-i+1, shape=shape, rate=rate) -
        pgamma(j-i,   shape=shape, rate=rate)
    }
    probs[n]  = probs[n] + (1 - sum(probs))  # terminal bin
    Xr[, i]   = probs
  }
  Xr[n, n] = 1.0  # terminal bin
  return(Xr)
}

# region 1: faster growth
# region 2: slower growth
X_list = list(
  build_Xr(n_size, mean_growth=1, cv_growth=0.5),  # X^1
  build_Xr(n_size, mean_growth=3, cv_growth=0.5)   # X^2
)

# X_full, block diagonal of X^r (n_size*n_region x n_size*n_region)
X_full = matrix(0, n_size*n_region, n_size*n_region)
for (r in 1:n_region) {
  idx = ((r-1)*n_size + 1):(r*n_size)
  X_full[idx, idx] = X_list[[r]]
}

# get survival
natmort_1    = rep(0.2, n_size)
natmort_2   = rep(0.25, n_size)
S        = exp(-cbind(natmort_1, natmort_2))  # n_size x n_region

# vec(S) stacked by region: [S[,1], S[,2]]
S_full = diag(as.vector(S))  # n_size*n_region x n_size*n_region

# get movement rates
m_rate <- 0.3
psi_list <- lapply(1:n_size, function(l) {
  matrix(c(1 - m_rate, m_rate,
           m_rate,     1 - m_rate), nrow = 2, byrow = TRUE)
})

# build T_full n_region x n_region blocks ... each n_size x n_size diagonal
T_full = matrix(0, n_size*n_region, n_size*n_region)
for (r in 1:n_region) {
  for (rp in 1:n_region) {
    # row block r, col block r'
    row_idx = ((r-1)*n_size  + 1):(r*n_size)
    col_idx = ((rp-1)*n_size + 1):(rp*n_size)
    # diagonal -  [psi^l]_{r,r'} for l = 1,...,n_size
    diag_vals = sapply(1:n_size, function(l) psi_list[[l]][r, rp])
    T_full[row_idx, col_idx] = diag(diag_vals)
  }
}

# recruitment vector
mean_rec_r = c(0.6, 0.4)
r_tilde = rep(0, n_size*n_region)
r_tilde[1]          = mean_rec_r[1]  # size bin 1, region 1
r_tilde[n_size + 1] = mean_rec_r[2]  # size bin 1, region 2

# get transition matrix
O = X_full %*% S_full %*% T_full

# look at eigens
max(abs(eigen(O)$values))

# analytical soln
I_mat   = diag(n_size * n_region)
n_equil = solve(I_mat - O) %*% r_tilde

# reshape to n_size x n_region (stacked by region)
n_equil = matrix(n_equil, n_size, n_region)

# iterative comparison
n_tilde = rep(0, n_size*n_region)
n_tilde[1]          = mean_rec_r[1]
n_tilde[n_size + 1] = mean_rec_r[2]

N_hist = list()
for (t in 1:500) {
  n_tilde        = O %*% n_tilde + r_tilde
  N_hist[[t]]    = matrix(n_tilde, n_size, n_region)
}

n_tilde_iter = matrix(n_tilde, n_size, n_region)

# plot to compare
par(mfrow=c(1,2))
plot(n_equil[,1], type="l", xlab="Size bin", ylab="N", main="Region 1", ylim = c(0, 2), lwd = 5)
lines(n_tilde_iter[,1], col="blue", lty=3, lwd=5)
legend("topright", c("Analytical","Iterative"), col=c("black","blue"), lty=c(1,3), lwd = 5)

plot(n_equil[,2], type="l", xlab="Size bin", ylab="N", main="Region 2", ylim = c(0, 2), lwd = 5)
lines(n_tilde_iter[,2], col="red", lty=3, lwd = 5)
legend("topright", c("Analytical","Iterative"), col=c("black","red"), lty=c(1,3), lwd = 5)

