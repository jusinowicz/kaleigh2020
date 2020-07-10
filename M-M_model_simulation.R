#the goal of this script is to practice setting up michaelis-menten metabolic dynamics for diazotrophic C-fixers and algae and thne project those dynamics through time in order to predict competitive outcomes

#background: lotka-volterra competition, with temperature dependence incoroporated in the competition coefficients, is one way to model dynamics between algae and phototrophic diazotrophs across a temperature gradient. Another option is to use michaelis-menten kinetics, where the rate of metabolism depends on substrate (N) concentration. Here, temperature influences the maximum rate of reaction (Vmax).

library(drc) #drc = dose-response curve
library(tidyverse)
library(deSolve)
library(FME)


#model 1 ---------------------------------------------------
#here d_x (dx/dt)  is the standard MM equation, dependent on its kinetics parameters and the substrate concentration. d_s (ds/dt), where s = substrate concentration, is dependent on the density of x, a consumption rate parameter (C), and substrate concentration

#putting all of the above in one function
# parameters
parms <- c(Vmax = 10, Km = 3, C = 0.1)

#the function
mm_model <- function(parms, times = seq(0, 50, by = 1)) {
  # initial state 
  state <-c(x = 20, substrate = 1)
  # derivative
   deriv <- function(t, state, parms) {
    with(as.list(c(state, parms)), {
      d_x <- Vmax*(substrate/(Km + substrate))
      d_s <- x*C*substrate
      return(list(c(x = d_s, y = d_x)))
    })
  }
  # solve
  ode(y = state, times = times, func = deriv, parms = parms)
} 

#output results
mm_results <- mm_model(parms = parms, times = seq(0, 50, by = 1))

#studying the lotka-volterra example straight out of ODE documentation
## Example2: Substrate-Producer-Consumer Lotka-Volterra model

## Note:
## Function sigimp passed as an argument (input) to model
## (see also lsoda and rk examples)
SPCmod <- function(t, x, parms, input) {
  with(as.list(c(parms, x)), {
    import <- input(t)
    dS <- import - b*S*P + g*C # substrate -- import = some supply rate?
    dP <- c*S*P - d*C*P # producer
    dC <- e*P*C - f*C # consumer
    res <- c(dS, dP, dC)
    list(res)
  })
}

## The parameters
parms <- c(b = 0.001, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)
## vector of timesteps
times <- seq(0, 200, length = 101)
## external signal with rectangle impulse
signal <- data.frame(times = times,
                     import = rep(0, length(times)))
signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2
sigimp <- approxfun(signal$times, signal$import, rule = 2)
## Start values for steady state
xstart <- c(S = 1, P = 1, C = 1)
## Solve model
out <- ode(y = xstart, times = times,
           func = SPCmod, parms = parms, input = sigimp)

## Default plot method
plot(out)

## User specified plotting
# mf <- par(mfrow = c(1, 2))
# matplot(out[,1], out[,2:4], type = "l", xlab = "time", ylab = "state")
# legend("topright", col = 1:3, lty = 1:3, legend = c("S", "P", "C"))
# plot(out[,"P"], out[,"C"], type = "l", lwd = 2, xlab = "producer",
#      ylab = "consumer")
# par(mfrow = mf)

