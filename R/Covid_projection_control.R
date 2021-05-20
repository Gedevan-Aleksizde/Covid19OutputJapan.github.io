Covid_projection_control <- function(
  InitialValues, # S, I, R, D [1, 4] vector
  alpha_on,  # status of emergency
  alpha_off,  # beta, gamma, delta, V; time-varying [T, 1] vector
  th_on,
  th_off,
  beta,
  gamma,
  delta,
  V,
  h,
  k,
  POP0,
  hconstant
){
  Time <- length(beta)
  SimData <- zeros(Time + 1, NCOL(InitialValues))
  SimData[1, ] <- InitialValues
  SimN <- zeros(Time, 1)
  alphapath <- zeros(Time, 1)
  alpha <- alpha_on
  for (i in 1:Time){
    alphapath[i] <- alpha
    if (hconstant == 0){
      SimN[i] <- ((1 + h * alpha) ^ k) * beta[i] * SimData[i, 1] * SimData[i, 2] * (1 / POP0)
    } else if (hconstant == 1){
      SimN[i] = ((1 + (h[2] / h[1]) * alpha) ^ k) * beta[i] * SimData[i, 1] * SimData[i, 2] * (1 / POP0)
    }
    SimData[i + 1, 1] <- SimData[i, 1] - SimN[i] - V[i]
    SimData[i + 1, 2] <- SimData[i, 2] + SimN[i] - gamma[i] * SimData[i, 2] - delta[i] * SimData[i, 2]
    SimData[i + 1, 3] <- SimData[i, 3] + gamma[i] * SimData[i, 2] + V[i];
    SimData[i + 1, 4] <- SimData[i, 4] + delta[i] * SimData[i, 2]
    if (hconstant == 0){
      if (th_on <= ((1 + h * alpha) ^ k) * beta[i] * SimData[i, 1] * SimData[i, 2] * (1 / POP0)){
        alpha <- alpha_on
      } else if (th_off >= ((1 + h * alpha) ^ k) * beta[i] * SimData[i,1] * SimData[i, 2] * (1 / POP0)){
        alpha <- alpha_off
      }
    } else if (hconstant == 1){
      if (th_on <= ((1 + (h[2]/h[1]) * alpha) ^ k) * beta[i] * SimData[i, 1] * SimData[i, 2] * (1 / POP0)){
        alpha <- alpha_on
      } else if (th_off >= ((1 + (h[2]/h[1]) * alpha) ^ k) * beta[i] * SimData[i,1] * SimData[i, 2] * (1 / POP0)){
        alpha <- alpha_off
      }
    }
  }
  CumD <- SimData[NROW(SimData), 4]  # Cumulative deaths during the simulation period
  GDPLoss <- mean(alphapath)         # Average output loss during the simulation period
  return(list(
    CumD,
    GDPLoss,
    alphapath,
    SimData,
    SimN
  ))
}