library(fitdistrplus)

n <- 152 #total number of index individuals
#vector of secondary cases of individuals for which secondary cases /= 0)
c1 <- c(1,2,2,5,14,1,4,4,1,3,3,8,2,1,1,4,9,9,1,1,17,2,1,1,1,4,3,3,4,2,5,1,2,2,1,9,1,3,1,2,1,1,2)
#concatenation of secondary cases from 43, non-zero individuals held in c1 with the remainder of individuals (n = 152 - 43) who produced no secondary cases
c0 <- c(c1, rep(0, n-length(c1)))

#c0 is a vector which contains the secondary cases produced by each of the total 152 index cases.
#Use plotdistr(c0) to visualise the distribution of cases
plotdistr(c0)

#use descdist to produce classical descriptive statistics
descdist(c0)
mean_cases <- descdist(c0)$mean

print(paste0("The mean number of secondary cases is : ", mean_cases))

#fit negative binomial distribution as used to describe offspring of overdispersed parameter
fit.cases <- fitdist(c0, "nbinom")
summary(fit.cases)
dispersion <- summary(fit.cases)$estimate[1]
mean_R <- summary(fit.cases)$estimate[2]

#maximum likelihood estimate of the mean(R0) = 0.95 and of the dispersion parameter (k) = 0.18

plot(fit.cases)
#this produces Empirial and theoretical density distribution and the cumulative density function (CDF) of the
#number of secondary cases



# Range of reported serial intervals for EBOLA cases
days <- 0:43
# Observed intervals for each day
frequency <- c(0,1,3,1,4,1,6,1,2,2,11,6,0,1,10,3,5,8,4,3,3,1,
               0,2,0,2,0,3,1,1,1,0,0,0,0,0,2,0,1,0,1,1,0,1)
d <- rep(days,frequency) #vector containing distribution of 92 cases and their serial interval
mean(d)

fit.serial <- fitdist(d, "gamma")
summary(fit.serial)


#####################################

### OUTBREAK SIMULATION

#random number generator
set.seed(645)
runs <- 100
#number of simulations

seed <- 1 #initial case

# Initialize plot
plot(NA,xlim=c(0,100),ylim=c(0,100),xlab="Time (days)",
     ylab="Cumulative number of EVD cases",frame=FALSE)
# Set color scheme for different trajectories
cols <- sample(terrain.colors(runs))

sims <- matrix(nrow = runs, ncol = 3)
colnames(sims) <- c("outbreak_no", "duration", "cases")

for (i in 1:runs) {
  cases <- seed # 1 if initial case
  t <- rep(0, seed) #
  times <- t
  while (cases > 0) {
    secondary <-
      rnbinom(cases,
              size = fit.cases$estimate[1],
              mu = fit.cases$estimate[2])
    t.new <- numeric()
    for (j in 1:length(secondary)) {
      t.new <-
        c(t.new,
          t[j] + rgamma(
            secondary[j],
            shape = fit.serial$estimate[1],
            rate = fit.serial$estimate[2]
          )
        )
    }
    cases <- length(t.new)
    t <- t.new
    times <- c(times, t.new)
  }
  lines(sort(times),
        1:length(times),
        col = cols[i],
        lwd = 1)
  points(max(times), length(times), pch = 16)
  
  outbreak_duration <- times[length(times)]
  outbreak_cases <- length (times)
  
  sims[i, ] <- c(i, outbreak_duration, outbreak_cases) #contains data for each simulated case import
}







