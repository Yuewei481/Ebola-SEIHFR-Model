next_step <- function(state, beta_I, beta_F, alpha, delta, gamma, eta, sigma) {
  S <- state$S
  E <- state$E
  I <- state$I
  H <- state$H
  F <- state$F
  R <- state$R
  
  N <- S + E + I + H + F + R

  new_I <- beta_I * S * I / N
  new_F <- beta_F * S * F / N

  S1 <- S - new_I - new_F 
  E1 <- E + new_I + new_F - alpha * E
  I1 <- I - (gamma + delta) * I + alpha * E
  H1 <- H + delta * I - eta * H 
  F1 <- F + gamma * I - sigma * F 
  R1 <- R + eta * H + sigma * F 
  
  S <- S1
  E <- E1
  I <- I1
  H <- H1
  F <- F1
  R <- R1

  return(list(S = S, E = E, I = I, H = H, F = F, R = R))
}


simulate_model <- function(t, S0, E0, I0, H0, F0, R0, beta_I, beta_F, alpha, 
                           delta, gamma, eta, sigma) {

  time <- 0:t

  S <- numeric(t + 1)
  E <- numeric(t + 1)
  I <- numeric(t + 1)
  H <- numeric(t + 1)
  F <- numeric(t + 1)
  R <- numeric(t + 1)

  S[1] <- S0
  E[1] <- E0
  I[1] <- I0
  H[1] <- H0
  F[1] <- F0
  R[1] <- R0

  state <- list(S = S0, E = E0, I = I0, H = H0, F = F0, R = R0)

  for (i in 1:t) {
    state <- next_step(state, beta_I, beta_F, alpha, delta, gamma, eta, sigma)
    
    S[i + 1] <- state$S
    E[i + 1] <- state$E
    I[i + 1] <- state$I
    H[i + 1] <- state$H
    F[i + 1] <- state$F
    R[i + 1] <- state$R
  }

  return(list(time = time, S = S, E = E, I = I, H = H, F = F, R = R))
}

S0 <- 990
E0 <- 10
I0 <- 0
H0 <- 0
F0 <- 0
R0 <- 0

beta_I <- 0.4
beta_F <- 0.5

alpha <- 0.3
delta <- 0.1
gamma <- 0.1
eta   <- 0.12
sigma <- 0.2

t <- 100

result <- simulate_model(t, S0, E0, I0, H0, F0, R0, beta_I, beta_F, alpha, 
  delta, gamma, eta, sigma)

calc_R <- function(S, N, beta_I, beta_F, gamma, delta, sigma) {
  (S / N) * (beta_I * sigma + beta_F * gamma) / (sigma * (delta + gamma))
}

N2 <- result$S + result$E + result$I + result$H + result$F + result$R

D <- result$E + result$I + result$F

# find peak of D
D_peak_index <- which.max(D)
D_peak_time <- result$time[D_peak_index]
D_peak_value <- D[D_peak_index]

# calculate R
R_values <- calc_R(result$S, N2, beta_I, beta_F, gamma, delta, sigma)

# R at peak of D
R_at_D_peak <- R_values[D_peak_index]

# first time when R < 1
R_below_1_index <- which(R_values < 1)[1]

# find the point just before R < 1
R_before_1_index <- R_below_1_index - 1

# corresponding values
R_below_1_time <- result$time[R_below_1_index]
R_below_1_value <- R_values[R_below_1_index]

R_before_1_time <- result$time[R_before_1_index]
R_before_1_value <- R_values[R_before_1_index]

plot(result$time, result$S, type = "l", lwd = 2,
     ylim = range(c(result$S, result$E, result$I, result$H, result$F, result$R)),
     xlab = "Time(days)", ylab = "Population",
     main = "SEIHFR Model Simulation")

lines(result$time, result$S, lwd = 2, col = "black")
lines(result$time, result$E, lwd = 2, col = "yellow")
lines(result$time, result$I, lwd = 2, col = "red")
lines(result$time, result$H, lwd = 2, col = "blue")
lines(result$time, result$F, lwd = 2, col = "purple")
lines(result$time, result$R, lwd = 2, col = "darkgreen")

legend("left",
       legend = c("Susceptible", "Exposed", "Infected", "Hospitalization", 
                  "Funeral", "Removed"),
       col = c("black","yellow", "red", "blue", "purple", "darkgreen"),
       lwd = 2)

# save current graphic settings
old_par <- par(no.readonly = TRUE)

# make space for right y-axis
par(mar = c(5, 4, 4, 5))

# plot D on left axis
plot(result$time, D, type = "l", lwd = 2, col = "darkred",
     xlab = "Time(days)", ylab = "D",
     main = "D and R over Time")

# mark peak of D
points(D_peak_time, D_peak_value, pch = 19, col = "darkred")
text(D_peak_time, D_peak_value,
     labels = paste0("(", D_peak_time, ", ", round(D_peak_value, 4), ")"),
     pos = 3, col = "darkred")

# add R curve on right axis
par(new = TRUE)
plot(result$time, R_values, type = "l", lwd = 2, col = "blue",
     axes = FALSE, xlab = "", ylab = "")

axis(side = 4)
mtext("R", side = 4, line = 3)

# add horizontal line R = 1
abline(h = 1, lty = 2, col = "gray40")

# mark first point where R < 1
points(R_below_1_time, R_below_1_value, pch = 19, col = "forestgreen")
text(R_below_1_time, R_below_1_value,
     labels = paste0("First R < 1\n(",
                     R_below_1_time, ", ",
                     round(R_below_1_value, 4), ")"),
     pos = 1, col = "forestgreen")

# mark the point just before R < 1
points(R_before_1_time, R_before_1_value, pch = 19, col = "orange")
text(R_before_1_time, R_before_1_value,
     labels = paste0("Last R >= 1\n(",
                     R_before_1_time, ", ",
                     round(R_before_1_value, 4), ")"),
     pos = 3, col = "orange")

legend("topright",
       legend = c("D", "Reproduction(R)", "R = 1", "Peak of D", 
                  "First R < 1", "Last R >= 1"),
       col = c("darkred", "blue", "gray40", "darkred", "forestgreen", "orange"),
       lty = c(1, 1, 2, NA, NA, NA),
       pch = c(NA, NA, NA, 19, 19, 19),
       lwd = c(2, 2, 1, NA, NA, NA),
       bty = "n")

# restore graphic settings
par(old_par)

cat("Peak time of E + I + F =", D_peak_time, "\n")
cat("Maximum value of E + I + F =", D_peak_value, "\n")
cat("R at the peak of E + I + F =", R_at_D_peak, "\n")
cat(result$time[R_below_1_index], "\n")





