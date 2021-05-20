Covid_projection <- function(
  InitialValues,
  alpha,
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
  for(i in 1:Time){
    if(hconstant == 0){
      SimN[i] <- ((1 + h * alpha)^k) * beta[i] * SimData[i, 1] * SimData[i, 2] * (1 / POP0)
    } else if(hconstant == 1){
      # SimN[i] = ((1 + (h[2]/h[1]) * alpha[i])^k)*beta[i] * SimData[i, 1] * SimData[i, 2] * (1 / POP0)
      SimN[i] = ((1 + (h[2]/h[1]) * alpha[i])^k) * beta[i] * SimData[i, 1] * SimData[i, 2] * (1/POP0)
    }
    SimData[i + 1, 1] <- SimData[i, 1] - SimN[i] - V[i]
    SimData[i + 1, 2] <- SimData[i, 2] + SimN[i] - gamma[i] * SimData[i, 2] - delta[i] * SimData[i, 2]
    SimData[i + 1, 3] <- SimData[i, 3] + gamma[i] * SimData[i, 2] + V[i];
    SimData[i + 1, 4] <- SimData[i, 4] + delta[i] * SimData[i, 2]
  }
  CumD <- tail(SimData[, 4], 1)  # Cumulative deaths during the simulation period
  GDPLoss <- mean(alpha)         # Average output loss during the simulation period
  return(list(
    CumD,
    GDPLoss,
    SimData,
    SimN
  ))
}