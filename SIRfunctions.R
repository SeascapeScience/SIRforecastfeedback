library(deSolve) # using the "ode" function
#library(fields) # for image.plot

SIRalpha <- function(bet = 0.03, 
                     r = 1, # human response rate
                     alphaS = 0.001, # fastest infections contact rate (/person/day)
                     alpha0 = 0.0005,
                     S = 999,
                     I = 1,
                     R = 0,
                     alpha = 0.001,
                     time_values = seq(0, 50, by = 0.1))
{
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -alpha * I * S
      dI <-  alpha * I * S - bet * I
      dR <-  bet * I
      dalpha <- r * alpha * 
        (1 - (alpha / ((  (alpha0-alphaS)/((S+I+R)^2)*(I^2)+alphaS   ))))
      return(list(c(dS, dI, dR, dalpha)))
    })
  }
  parameters_values <- c(bet=bet, r=r, alphaS=alphaS, alpha0=alpha0)
  initial_values <- c(S=S,I=I,R=R,alpha=alpha)
  sir_values_1 <- ode(
    y = initial_values,
    times = time_values,
    func = sir_equations,
    parms = parameters_values 
  )
  sir_values_1 <- as.data.frame(sir_values_1)
  return(sir_values_1)
}

SIRalphaS <- function(bet = 0.03, 
                     r = 1, # human response rate
                     alphaS = 0.001, # fastest infections contact rate (/person/day)
                     alpha0 = 0.0005,
                     S = 999,
                     I = 1,
                     R = 0,
                     alpha = 0.001,
                     time_values = seq(0, 50, by = 0.1))
{
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -alpha * I * S
      dI <-  alpha * I * S - bet * I
      dR <-  bet * I
      dalpha <- r * alpha * 
        (1 - (alpha / (( (((S+I+R)*(1-(bet/alpha/(S+I+R))-(bet/alpha/(S+I+R))*log(alpha*(S+I+R)/bet)))^2)*(alpha0-alphaS)/((S+I+R)^2)+alphaS ))))
      return(list(c(dS, dI, dR, dalpha)))
    })
  }
  parameters_values <- c(bet=bet, r=r, alphaS=alphaS, alpha0=alpha0)
  initial_values <- c(S=S,I=I,R=R,alpha=alpha)
  sir_values_1 <- ode(
    y = initial_values,
    times = time_values,
    func = sir_equations,
    parms = parameters_values 
  )
  sir_values_1 <- as.data.frame(sir_values_1)
  return(sir_values_1)
}

plotSIR <- function(sir_values_1)
{
  with(sir_values_1, {
    # plotting the time series of susceptibles:
    plot(time, S, type = "l", col = "blue",
         xlab = "time (days)", ylab = "number of people",
         ylim = c(0,max(S)))
    # adding the time series of infectious:
    lines(time, I, col = "red")
    # adding the time series of recovered:
    lines(time, R, col = "green")
    lines(time, alpha, col = "black")
  })
  legend("right", c("susceptibles", "infectious", "recovered"),
         col = c("blue", "red", "green"), lty = 1, bty = "n")
}

#

SIRalphatheta <- function(bet = 0.03, 
                     r = 1, # human response rate
                     alphaS = 0.001, # fastest infections contact rate (/person/day)
                     alpha0 = 0.0005,
                     theta = 2,
                     S = 999,
                     I = 1,
                     R = 0,
                     alpha = 0.001,
                     time_values = seq(0, 50, by = 0.1))
{
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -alpha * I * S
      dI <-  alpha * I * S - bet * I
      dR <-  bet * I
      dalpha <- r * alpha * 
        (1 - (alpha / ((  (alpha0-alphaS)/((S+I+R)^2)*(I^theta)+alphaS   ))))
      return(list(c(dS, dI, dR, dalpha)))
    })
  }
  parameters_values <- c(bet=bet, r=r, alphaS=alphaS, alpha0=alpha0)
  initial_values <- c(S=S,I=I,R=R,alpha=alpha)
  sir_values_1 <- ode(
    y = initial_values,
    times = time_values,
    func = sir_equations,
    parms = parameters_values 
  )
  sir_values_1 <- as.data.frame(sir_values_1)
  return(sir_values_1)
}

SIRalphaE <- function(bet = 0.03, 
                     r = 1, # human response rate
                     alphaS = 0.001, # fastest infections contact rate (/person/day)
                     alpha0 = 0.0005,
                     S = 999,
                     I = 1,
                     R = 0,
                     alpha = 0.001,
                     time_values = seq(0, 50, by = 0.1))
{
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -alpha * I * S
      dI <-  alpha * I * S - bet * I
      dR <-  bet * I
      dalpha <- r * alpha * (1-(alpha / (alpha0+(alphaS-alpha0)/exp(1000*I)) ))
      return(list(c(dS, dI, dR, dalpha)))
    })
  }
  parameters_values <- c(bet=bet, r=r, alphaS=alphaS, alpha0=alpha0)
  initial_values <- c(S=S,I=I,R=R,alpha=alpha)
  sir_values_1 <- ode(
    y = initial_values,
    times = time_values,
    func = sir_equations,
    parms = parameters_values 
  )
  sir_values_1 <- as.data.frame(sir_values_1)
  return(sir_values_1)
}


