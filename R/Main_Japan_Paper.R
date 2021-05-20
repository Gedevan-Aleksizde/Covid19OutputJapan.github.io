require(readr)
require(zeallot)
source("R/wrappers.R")
source("R/Covid_projection.R")
source("R/Covid_projection_control.R")

# --- Graphical parameters ---

figure_save <- T
dev.args <- list()
fontfamily <- "Noto Sans CJK JP"
savedir <- "R/img/Japan"
xlim_tradeoff <- c(1.5, 3)
iDrawUB <- T  # F = create quasi-UB without simulations, T = create UB with simulations
Nsim <- 40000  # if iDrawUB=1, this is the number of draws you use.

# ---- Model parameter values ----

POP0 <- 125710000  # initial population
gamma <- 7 / 5  # recovery rate from Covid
k <- 2  # exponent of (1-h*alpha)
hconstant <- 1  # 0 = without intercept, 1 = with intercept for h regression
TargetAlpha <- seq(1, 3, by = .1)  # values of alpha we simulate results
AlphaVals <- c(1.5, 0.85 * (1.5 + 3) / 2, 3)  # benchmark alpha we plot time-series figures
SimPeriod <- 52  # simulation period in weeks
VacStart <- 11  # time until the start of vaccination process
VacPace <- 1.33 * 0.8 * 1500000 # number of vaccinated persons per week (steady pace)
VacDuration <- 12  # time until vaccination pace reaches its steady pace
RetroPeriod <- 20  # retroactive periods used to estimate gamma and delta
wl <- 1:2  # Results as of these weeks ago

# ---- Import data ----
# Covid data are recorded at weekly frequencies (Mon-Sun)
# The first week start on January 20 (Mon), 2020
#covid = importdata([home 'Data/Covid_weekly2.csv']);  % Import weekly Covid data

covid <- read_csv("Covid19_Output/Data/Covid_weekly2.csv")  # Import weekly covid data

N <- covid$positive
dD <- covid$death
M <- covid$mobility
GDP <- covid$gdp
Tdata <- NROW(covid)  # Data period in weeks
TdataGDP <- Tdata - sum(is.na(GDP))
POP <- covid$manufacturing
Month <- covid$month[-1]
xtick1 <- c(seq(1, Tdata, by = 13), Tdata)
M <- 1 + 0.01 * M

RetroH <- TdataGDP - 4


# ---- Constructing the reference level of output ----

potentialGDP <- zeros(52 * 3, 1)                            # potential GDP for the next 3 years
potentialGDP[1] <- (548182 / (1.0122)) * (1.0063 ^ (1 / 12))
for (i in 2:length(potentialGDP)) {
  if (i <= 13) {
    potentialGDP[i] <- potentialGDP[i - 1] * (1.0063 ^ (1 / 52))
  } else if (i <= 52) {
    potentialGDP[i] <- potentialGDP[i - 1] * (1.0021 ^ (1 / 52))
  } else if (i <= 104) {
    potentialGDP[i] = potentialGDP[i - 1] * (1.0024 ^ (1 / 52))
  } else if (i <= 156) {
    potentialGDP[i] = potentialGDP[i - 1] * (1.0021 ^ (1 / 52))
  }
}

## TODO: 元コードでもここはコメントアウトされていた. 用途は何?

# gapGDP = zeros(length(potentialGDP),1);
# gapGDP(1) = 0.0166;
# for i = 2:length(gapGDP)
#     gapGDP(i) = gapGDP(i-1)*(0.975^(12/52));
# end

referenceGDP <- potentialGDP * (1 + 0.0166)
referenceGDP <- referenceGDP[-c(1, 2)]


# ---- Impute alpha (regress alpha on M) ----
## TODO: この辺のインデックスの定数は何を意味しているのか
Malt <- M
Malt[50] <- 0.5 * (Malt[49] + Malt[51])
alpha <- matrix((1 - GDP[1:TdataGDP] / referenceGDP[1:TdataGDP]), ncol = 1)  # output loss in percentage
X <- M[(TdataGDP - 17):TdataGDP]
Y <- matrix(alpha[(TdataGDP - 17):TdataGDP], ncol = 1)
XC <- cbind(ones(length(X), 1), X)
## matlab ではどうなのか知らないが R ではたぶん直接逆行列を求めるべきではない
s <- solve(t(XC) %*% XC, t(XC) %*% Y)  # OLS estimate of h with constant
reg <- XC %*% s
r <- Y - reg
SSE <- sum(r^2)
eps_p <- rep(0, Tdata - TdataGDP)
eps_p[1] <- tail(r, 1)
for(i in 1:(Tdata - TdataGDP - 1)){
  eps_p[i + 1] <- 1 * eps_p[i]
}
alpha_pred <- s[1] + s[2] * Malt[(TdataGDP+1):Tdata] + eps_p
alpha <- c(alpha, alpha_pred)

# ---- Plot mobility data ----
## TokyoMobilityGDPLine.png
par(family = fontfamily, mar = c(5, 4, 4, 4), xpd = T)
plot(M, type = "l", ylab = "Mobility", xaxt="n", xlab = "")
par(new = T)
plot(GDP, type = "l", lty = "dotdash" , col = "blue", axes = F, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(GDP, na.rm = T)), col = "blue", col.axis="blue")
axis(1, at = xtick1, labels = covid$month[xtick1], las=2)
mtext("GDP", side = 4, line = 2.5, col = "blue")
legend("top", inset=c(0,0),
       legend = c("Mobility (left axis)", "GDP (right axis)"), col = c("black", "blue"), lty = c(1, 4))
if (figure_save) saveas("MobilityGDPLine.png", savedir, dev.args = dev.args)
par(xpd = F)
# ---- Regress mobility on alpha to estimate the elasticity h ----

Y <- M[4:TdataGDP]
X <- matrix(alpha[4:TdataGDP], ncol = 1)
if (hconstant == 0) {
  Y <- Y - 1
  h_all <- solve(t(X) %*% X, t(X) %*% Y) # OLS estimate of h
  reg <- X %*% h_all
  r <- Y - reg                           # r is the residuals, which is the observed minus fitted values
  SSE <- sum(r ^ 2)                       # SSE is the sum of squared errors
  MSE <- SSE / (length(Y) - 1)            # mean squared error
  h_all_se <- sqrt(MSE / sum(X ^ 2))      # standard error of h
} else if (hconstant == 1) {
  XC <- cbind(ones(length(X), 1), X)
  h_all <- solve(t(XC) %*% XC, t(XC) %*% Y) # OLS estimate of h with constant
  reg <- XC %*% h_all
  r <- Y - reg
  SSE <- sum(r ^ 2)
  MSE <- SSE / (length(Y) - 1)            # mean squared error
  h_all_se <- zeros(2, 1)
  h_all_se[1] <- sqrt(MSE / sum(XC[, 1] ^ 2))   # standard error of h
  h_all_se[2] <- sqrt(MSE / sum(XC[, 2] ^ 2))
}

Y <- M[(TdataGDP - RetroH):TdataGDP]
X <- matrix(alpha[(TdataGDP - RetroH):TdataGDP], ncol = 1)
if (hconstant == 0) {
  Y <- Y - 1
  h <- solve(t(X) %*% X, t(X) %*% Y)  # OLS estimate of h
  reg = X %*% h
  r <- Y - reg                            # r is the residuals, which is the observed minus fitted values
  SSE <- sum(r ^ 2)                       # SSE is the sum of squared errors
  MSE <- SSE / (length(Y) - 1)            # mean squared error
  h_se <- sqrt(MSE / sum(X ^ 2))          # standard error of h
} else if (hconstant == 1) {
  XC <- cbind(ones(length(X), 1), X)
  h <- solve(t(XC) %*% XC, t(XC) %*% Y)   # OLS estimate of h with constant
  reg <- XC %*% h
  r <- Y - reg
  SSE <- sum(r ^ 2)
  MSE = SSE / (length(Y) - 1)                   # mean squared error
  h_se = zeros(2, 1)
  h_se[1] <- sqrt(MSE / sum(XC[, 1] ^ 2))     # standard error of h
  h_se[2] <- sqrt(MSE / sum(XC[, 2] ^ 2))
}

## TODO: グラフ完成品が保存されてない?
par(family = fontfamily, ann = F)
if (hconstant == 0) {
  plot(X, Y, bg = "blue", pch = 21)
  par(new = T)
  plot(X, reg, axes = F, type = "l")
} else if (hconstant == 1) {
  plot(XC[, 2], Y, bg = "blue", pch = 21)
  par(new = T)
  plot(XC[, 2], reg, type = "l", axes = F)
}
title(xlab = "Output loss (%)", ylab = "Mobility")
if(figure_save) saveas("MobilityGDPScatter.png", savedir, dev.args = dev.args)



# ---- Compute the history of S, I, R, D in the data period ----

S <- zeros(Tdata + 1, 1)
I <- zeros(Tdata + 1, 1)
R <- zeros(Tdata + 1, 1)
D <- zeros(Tdata + 1, 1)
S[1] <- POP0
for (i in 1:Tdata) {
  S[i + 1] <- S[i] - N[i]
  I[i + 1] <- I[i] + N[i] - gamma * I[i] - dD[i]
  R[i + 1] <- R[i] + gamma * I[i]
  D[i + 1] <- D[i] + dD[i]
  if (i > TdataGDP) {
    GDP[i] <- referenceGDP[i] * (1 - alpha[i])
  }
}

# ---- Compute the history of time-varying parameters ----

delta <- (D[2:(Tdata + 1)] - D[1:Tdata]) / I[1:Tdata]       # death rate
beta_tilde <- -POP0 * ((S[2:(Tdata + 1)] - S[1:Tdata]) / (S[1:Tdata] * I[1:Tdata]))  # overall infection rate
ERN = (S[1:(NROW(S) - 1)] / POP0) * beta_tilde / (gamma + delta)                        # effective reproduction number
if (hconstant == 0) {
  beta <- beta_tilde / (1 + h_all * alpha) ^ k              # raw infection rate
} else if (hconstant == 1) {
  ## TODO: このコメントアウトの意味
  # beta = beta_tilde./(h(1)+h(2)*alpha).^k
  beta <- beta_tilde / (1 + (h_all[2] / h_all[1]) * alpha) ^ k
}
beta <- matrix(beta, ncol = 1)
AverageAlpha2020 <- 100 * mean(c(alpha[1], alpha[2], alpha[1:(NROW(alpha) - 2)]))

# ---- Plot varabiles in the data period ----
## Variables.png
par(family = fontfamily)
plot(N, type = "l", xaxt = "n", yaxt = "n", xlab = "")
axis(1, at = xtick1, labels = covid$month[xtick1], las = 2)
axis(2, at = 0:4 * 10000, labels = formatC(0:4 * 10000, "%,6.0f"))
if(figure_save) saveas("NewCases.png", savedir, dev.args = dev.args)
par(new=F)

VarList <- c("S","I","R","D","N","GDP")
VarName <- c("S","I","R","D","N","Y")
par(family = fontfamily, mfrow = c(2, 3))
for(i in 1:length(VarList)){
  y_ <- get(VarList[i])
  if (i == length(VarList)) {
    ylim <- c(min(y_, referenceGDP), max(y_, referenceGDP))
  } else {
    ylim = c(min(y_), max(y_))
  }
  plot(y_, type = "l", lwd = 2, ann = F, xaxt = "n", yaxt = "n", xlim = c(1, Tdata), ylim = ylim)
  title(main = VarName[i], ylab = "")
  axis(1, at = xtick1, labels = covid$month[xtick1])
  if(i == 5){
    ## 指数表記はめんどうなので再現しない
    # ax.YAxis.Exponent = 0;
    axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
  } else {
    axis(2, at = axTicks(2))
  }
  if(i == length(VarList)){
    par(new = T)
    plot(referenceGDP[1:Tdata], lty = 2, lwd = 2, type = "l", axes = F, ann = F,
         xlim = c(1, Tdata), ylim = ylim)
  }
}
if(figure_save) saveas("Variables.png", savedir, dev.args = dev.args)
par(mfrow = c(1, 1), new = F)

# %--- Plot parameters in the data period ---%
## parameters.png
ParamList <- c("alpha", "ERN", "beta", "delta")
ParamName <- c(
  "$\\alpha$ (decline in Y)",
  "Effective reproduction number",
  "$\\beta$ (raw infection rate)",
  "$\\delta$ (death rate)"
)
par(family = fontfamily, mfrow = c(2, 2))
for (i in 1:4) {
  ## TODO: 手動で調整しているのはなぜ?
  if (i == 2) ylim <- c(0, 2.5)
  else if (i == 3) ylim <- c(0, 5)
  else if (i == 4) ylim <- c(0, .1)
  else ylim <- NULL
  plot(get(ParamList[i]), lty = 2, type = "l", xaxt = "n", ann = F, ylim = ylim)
  axis(1, at = xtick1, labels = covid$month[xtick1])
  title(main = ParamName[i], ylab = "")
}
if(figure_save) saveas("Parameters.png", savedir, dev.args = dev.args)
par(mfrow = c(1, 1))

################### Projection starts here ###################

# ---- Construct time series of parameters ----

InitialValues <- cbind(tail(S, 1), tail(I, 1), tail(R, 1), tail(D, 1))
beta_sample <- beta[(NROW(beta) - RetroPeriod + 1):NROW(beta)]
delta_sample <- delta[(NROW(beta) - RetroPeriod + 1):NROW(beta)]
beta_average <- mean(beta_sample)
delta_average <- mean(delta_sample)
betaT <- beta_average * ones(SimPeriod, 1)
deltaT <- delta_average * ones(SimPeriod, 1)
gammaT <- gamma * ones(SimPeriod, 1)
## TODO: このコメントアウトの意図
# V = zeros(SimPeriod,1);
# V(VacStart:end) = VacPace;
V <- zeros(SimPeriod, 1)
V[(VacStart - 1):(VacStart + VacDuration - 1)] <- seq(0, VacPace, by = VacPace / VacDuration) ## 0:VacPace/VacDuration:VacPace ## TODO
V[(VacStart + VacDuration):NROW(V)] <- VacPace
CumD <- zeros(1, length(TargetAlpha))
AverageAlpha <- zeros(1, length(TargetAlpha))
LagResults <- zeros(2, length(TargetAlpha), length(wl)[1])
ParamList2 <- c("alpha2", "ERN2", "beta2", "delta2")

# ---- Sensitivity parameters ----

BetaVals <- beta_average * c(1.05, 1, 0.95)
KVals <- c(1.5, 2, 2.5)
HVals <- c(0.8, 1, 1.2)

VacStartVals <- c(16, 12, 20)
VacPaceVals <- c(0.5 * VacPace, VacPace, 2 * VacPace)
FrontLoadingVals <- c(0.5, 1, 1.5)
SA <- zeros(3, length(TargetAlpha), 4)  #matrix to contain the results of sensitivity analysis
AlphaIndexVariables <- zeros(SimPeriod, 6, length(AlphaVals))    # matrix to contain the projected variables for certain alpha values

# ---- Simulation for different alpha values ----
for(i in 1:length(TargetAlpha)){
  # %%% Baseline projection %%%
  # TODO: このコメントアウトの意図
  # %     alphaT = flip(0:0.01*TargetAlpha(i)*2/(SimPeriod-1):(0.01*TargetAlpha(i)*2))';
  alphaT <- zeros(SimPeriod, 1)
  ## alphaT[1:26] <- t(flip(0:0.01*TargetAlpha[i]*4/(26-1):(0.01*TargetAlpha[i]*4)))
  alphaT[1:26] <- flip(seq(0, 0.01 * TargetAlpha[i] * 4, by = 0.01 * TargetAlpha[i] * 4 / (26 - 1)))  ## ここでは R は転置しなくても代入できる
  c(CumD[i], AverageAlpha[i], SimData, SimN) %<-% Covid_projection(InitialValues, alphaT, betaT, gammaT, deltaT, V, h, k, POP0, hconstant)  # this syntax requires zeallot package
  
  # %%% Results from previous weeks %%%
  for (j in 1:length(wl)){
    wlag <- wl[j]
    VacStartL <- VacStart + wlag
    V <- zeros(SimPeriod, 1)
    V[(VacStartL - 1):(VacStartL + VacDuration - 1)] <- seq(0, VacPace, by = VacPace/VacDuration) ## 0:VacPace/VacDuration:VacPace
    V[(VacStartL + VacDuration):NROW(V)] <- VacPace
    InitialValuesL <- cbind(S[NROW(S) - wlag], I[NROW(I) - wlag], R[NROW(R) - wlag], D[NROW(D) -wlag])
    betaL <- mean(beta[(NROW(beta) - RetroPeriod + 1 - wlag):(NROW(beta)-wlag)]) * ones(SimPeriod, 1)
    deltaL <- mean(delta[(NROW(delta) - RetroPeriod + 1 - wlag):(NROW(delta) - wlag)]) * ones(SimPeriod, 1)
    c(LagResults[1,i,j], LagResults[2,i,j]) %<-% Covid_projection(InitialValuesL, alphaT, betaL, gammaT, deltaL, V, h, k, POP0, hconstant)[1:2]
  }
  
  # %%% Beta sensitivity %%%
  for(j in 1:length(BetaVals)){
    betaS <- BetaVals[j] * ones(SimPeriod, 1)
    SA[j, i, 1] <- Covid_projection(InitialValues, alphaT, betaS, gammaT, deltaT, V, h, k, POP0, hconstant)[[1]]
  }
  # %%% h sensitivity %%%
  for(j in 1:length(HVals)){
    if(hconstant == 0){
      HS = HVals[j]
    } else if(hconstant == 1){
      HS <- h
      HS[2] <- h[2] * HVals[j]
    }
    SA[j, i, 2] <- Covid_projection(InitialValues, alphaT, betaT, gammaT, deltaT, V, HS, k, POP0, hconstant)[[1]]
  }
  # %%% k sensitivity %%%
  for(j in 1:length(KVals)){
    kS <- KVals[j]
    SA[j, i, 3] <- Covid_projection(InitialValues, alphaT, betaT, gammaT, deltaT, V, h, kS, POP0, hconstant)[[1]]
  }
  # %%% VacPace sensitivity %%%
  for(j in 1:length(VacPaceVals)){
    VacPaceS <- VacPaceVals[j]
    VS <- zeros(SimPeriod, 1)
    VS[(VacStart - 1):(VacStart + VacDuration - 1)] <- seq(0, VacPaceS, by = VacPaceS/VacDuration) ## 0:VacPaceS/VacDuration:VacPaceS
    VS[(VacStart + VacDuration):NROW(VS)] <- VacPaceS
    SA[j, i, 4] <- Covid_projection(InitialValues, alphaT, betaT, gammaT, deltaT, VS, h, k, POP0, hconstant)[[1]]
  }
}

# %%% Time series figures (variables and parameters) %%%

for(i in 1:length(AlphaVals)){
  alphaT <- zeros(SimPeriod, 1);    
  alphaT[1:26] <- t(flip(seq(0, 0.01 * AlphaVals[i] * 4, by = 0.01 * AlphaVals[i] * 4 / (26 - 1))))
  SimGDP <- referenceGDP[(Tdata + 1):(Tdata + NROW(alphaT))] * (1 - alphaT)
  c(CD, AA, SimData, SimN) %<-% Covid_projection(InitialValues, alphaT, betaT, gammaT, deltaT, V, h, k, POP0, hconstant)
  AlphaIndexVariables[, , i] <- cbind(SimData[1:(NROW(SimData) - 1), ], SimN, SimGDP)
}



Past <- cbind(N, GDP, ERN, zeros(Tdata, 1))
VarList2 <- c("N", "GDP", "ERN", "V")
VarName2 <- c("N", "Y", "Effective reproduction number", "Newly vaccinated persons")
par(family = fontfamily, mfrow = c(2,2), ann = F)
xlim <- c(1, Tdata+SimPeriod)
for(j in 1:4){
  for (i in 1:length(AlphaVals)){
    ## TODO: このコメントアウトの意図
    # %     alphaT = flip(0:0.01*AlphaVals(i)*2/(SimPeriod-1):(0.01*AlphaVals(i)*2))';
    alphaT <- zeros(SimPeriod, 1)
    alphaT[1:26] = flip(seq(0, 0.01 * AlphaVals[i] * 4, by = 0.01 * AlphaVals[i] * 4 / (26 - 1)))
    # %     alphaT = 0.01*AlphaVals(i)*ones(SimPeriod,1);
    
    if (hconstant == 0){
      ERNT <- (AlphaIndexVariables[, 1, i] / POP0) * (((1 + h * alphaT) ^ k) * betaT) / (gammaT + deltaT)
    } else if(hconstant == 1){
      ERNT <- (AlphaIndexVariables[, 1, i] / POP0) * (((1 + (h[2] / h[1]) * alphaT) ^ k) * betaT) / (gammaT + deltaT)
    }
    ProjectedData <- cbind(AlphaIndexVariables[, 5, i], AlphaIndexVariables[, 6, i], ERNT, V)
    CombinedData <- rbind(Past, ProjectedData)
    
    if (j == 2){
      ## base R では事後的にスケール変更できないのでmatlabに対応させるのが難しい
      ylim <- c(-12, 2)
      temp <- 100 * (CombinedData[, j] - potentialGDP[1 * NROW(CombinedData[, j])]) /  potentialGDP[length(CombinedData[, j])]
      if (i == 1){
        plot(temp, col = "red", lwd = 2, type = "l", xaxt = "n", xlim = xlim, ylim = ylim)
      } else if (i == 2){
        par(new = T)
        plot(temp,  type = "l", lwd = 2, xaxt = "n", xlim = xlim, ylim = ylim)
      } else if (i == 3){
        par(new = T)
        plot(temp, type = "l", col = "blue", lwd = 2, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
      }
      par(new = T)
      plot(temp[1:Tdata], type = "l", lwd = 2, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
    } else if (j == 4){
      ylim <- c(0, 15e5 * 1.1)
      if(i != 1) par(new = T)
      plot(CombinedData[, j], type = "l", lwd = 2, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim) 
    } else {
      ## base R では事後的にスケール変更できないので同上
      if(j == 1) ylim <- c(0, 140000)
      if(j == 3) ylim <- c(0.5, 2)
      if (i == 1){
        plot(CombinedData[, j], type = "l", col = 'red', lwd = 2, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
      } else if (i == 2){
        par(new = T)
        plot(CombinedData[, j], type = "l", lwd = 2, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
      } else if (i == 3){
        par(new = T)
        plot(CombinedData[, j], type = "l", col = 'blue', lwd = 2, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
      }
      par(new = T)
      plot(Past[, j], type = "l", lwd = 2, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
    }
  }
  if (j == 2){
    title(main = "output (deviation from potenial)")
  } else {
    title(main = VarName2[j])
  }
  abline(v = Tdata)
  axis(1, at = c(1, 27, 53, 79, 98),
       labels = c('Jan-20', 'Jul-20', 'Jan-21', 'Jul-21', 'Dec-21'))
  if (j == 1){
    axis(2, at = c(2e4, 4e4, 6e4, 8e4, 10e4, 12e4),
         labels = formatC(c(2e4, 4e4, 6e4, 8e4, 10e4, 12e4), "%6.0f", big.mark = ","))
  }
  if (j == 3){
    abline(h = 1)
  }
  if (j == 4){
    axis(2, at = c(5e5, 10e5, 15e5, 20e5, 25e5),
         labels = c('0.5M','1M','1.5M','2M','2.5M'))
  }
}
if(figure_save) saveas("VariablesProjection.png", savedir, dev.args = dev.args)
par(new = F, mfrow = c(1, 1))

## TODO:
#   %%%%%for i = 1:length(AlphaVals)
#   %%%%%%     alphaT = flip(0:0.01*AlphaVals(i)*2/(SimPeriod-1):(0.01*AlphaVals(i)*2))';
# %%%%%    alphaT = zeros(SimPeriod,1);    
# %%%%%alphaT(1:26) = flip(0:0.01*AlphaVals(i)*4/(26-1):(0.01*AlphaVals(i)*4))';
#   %%%%%    %     alphaT = 0.01*AlphaVals(i)*ones(SimPeriod,1);
#   %%%%%    if hconstant == 0
#   %%%%%        ERNT = (AlphaIndexVariables(:,1,i)/POP0).*(((1+h*alphaT).^k).*betaT)./(gammaT+deltaT);
#   %%%%%    elseif hconstant == 1
#   %%%%%        ERNT = (AlphaIndexVariables(:,1,i)/POP0).*(((1+(h(2)/h(1))*alphaT).^k).*betaT)./(gammaT+deltaT);
#   %%%%%    end
#   %%%%%    ProjectedData = [AlphaIndexVariables(:,5,i),AlphaIndexVariables(:,6,i),ERNT,V];
#   %%%%%    CombinedData = [Past;ProjectedData];
#   %%%%%    figure(5)
#   %%%%%    for j = 1:length(VarList2)
#   %%%%%        subplot(2,2,j)
#   %%%%%        plot(CombinedData(:,j),'k','LineWidth',2)
#   %%%%%        xlim([1 Tdata+SimPeriod])
#   %%%%%        title(VarName2(j),'FontWeight','normal')
#   %%%%%        xline(Tdata);
#   %%%%%        xticks([1 27 53 79 98])
#   %%%%%        xticklabels( {'Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'} )
#   %%%%%        xtickangle(45)
#   %%%%%        if j == 1
#   %%%%%            yticks([2e4 4e4 6e4 8e4 10e4 12e4]);
#   %%%%%            ax = gca;
#   %%%%%            ax.YAxis.Exponent = 0;
#   %%%%%            ytickformat('%,6.0f')
# %%%%%        end
# %%%%%        if j == 3
# %%%%%            ylim([0 3]);
# %%%%%            yline(1);
# %%%%%        end
# %%%%%        if j == 4
# %%%%%            yticks([5e5 10e5 15e5 20e5 25e5]);
# %%%%%            yticklabels({'0.5M','1M','1.5M','2M','2.5M'});
# %%%%%        end
# %%%%%        hold on
# %%%%%    end
# %%%%%    
# %%%%%    alpha2 = [alpha;alphaT];
# %%%%%    beta2 = [beta;betaT];
# %%%%%    if hconstant == 0
# %%%%%    beta_tildeT = betaT.*(1+h*alphaT).^k;                                      % raw infection rate
# %%%%%    elseif hconstant == 1
# %%%%%%     beta_tildeT = betaT.*(h(1)+h(2)*alpha).^k;
# %%%%%    beta_tildeT = betaT.*(1+(h(2)/h(1))*alphaT).^k;    
# %%%%%    end
# %%%%%    beta_tilde2 = [beta_tilde;beta_tildeT];
# %%%%%    ERN2 = [ERN;ERNT];
# %%%%%    delta2 = [delta;deltaT];
# %%%%%    V2 = [zeros(Tdata,1);V];
# %%%%%    figure(101)
# %%%%%    %set(gcf,'Position',[100,100,800,500])
# %%%%%    for j = 1:length(ParamList2)
# %%%%%        subplot(2,2,j)
# %%%%%        plot(eval(ParamList2(j)),'k','LineWidth',2)
# %%%%%        xlim([1 Tdata+SimPeriod])
# %%%%%        title(ParamName(j),'Interpreter','latex','FontSize',11,'FontWeight','normal')
# %%%%%        xline(Tdata);
# %%%%%        xticks([1 27 53 79 98])
# %%%%%        xticklabels( {'Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'} )
# %%%%%        xtickangle(45)
# %%%%%        if j == 2
# %%%%%            ylim([0,2])
# %%%%%        end
# %%%%%        if j == 3
# %%%%%            ylim([0,5])
# %%%%%        end
# %%%%%        if j == 4
# %%%%%            ylim([0,0.1])
# %%%%%        end
# %%%%%        hold on
# %%%%%    end
# %%%%%end


# %--- Trade-off figure (baseline) ---%
par(family = fontfamily)
plot(100 * AverageAlpha, CumD, type = "l",
     xlim = xlim_tradeoff, ann = F, yaxt = "n")
abline(h = D[NROW(D) - 2])
title(xlab = "Output Loss (%)", ylab = "Cumulative Deaths")
axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
if(figure_save) saveas("BaselineTradeoff.png", savedir, dev.args = dev.args)

# ---- Trade-off figure (lag) ----
par(family = fontfamily)
cols <- c("black", "red", "blue")
plot(100 * AverageAlpha, CumD,
     type = "l", lwd = 2,
     ann = F, xlim = xlim_tradeoff, ylim = c(0, 50000), yaxt = "n")
for(i in 1:length(wl)){
  par(new = T)
  plot(100 * LagResults[2, , i], LagResults[1, , i],
       type = "l", col = cols[i + 1], lwd = 2,
       xlim = xlim_tradeoff, ylim = c(0, 50000),
       ann = F, axes = F, yaxt = "n")
}
abline(v = AverageAlpha2020)
legend("top",
       legend = c("baseline", sprintf('%.0f week%s ago', 1:length(wl), ifelse(1:length(wl) == 1, "", "s"))),
       col = cols[c(1:(length(wl) + 1))], lty = 1, cex = .9)
title(xlab = "Output Loss (%)", ylab = "Cumlative Deaths")
## TODO: グラフには表示されてない??
# text(3.5,100000,{'Output loss';'      in 2020'},'FontSize',16);
axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
if(figure_save) saveas("LaggedTradeoff.png", savedir, dev.args = dev.args)

# ---- Trade-off figures with sensitivity analysis ----
SensitivityTitle <- c("beta","h","k","Vaccine pace")
SensitivitySaveTitle <- c("beta", "h", "k", "VacPace")
SensitivityLegend = list(
  c("5% higher", "baseline", "5% lower"),
  c("20% lower", "baseline", "20% higher"),
  c("1.5", "2 (baseline)", "2.5"),
  c("0.5*baseline", "baseline", "2*baseline")
)
col <- c("blue", "black", "red")
lty <- c(2, 1, 4)
par(family = fontfamily, mfrow = c(2, 2))
for (i in 1:length(SensitivityTitle)){
  ## そのまま for 文に置き換えるのはつらいので matplot 使用
  matplot(matrix(rep(TargetAlpha, 3), ncol = 3), t(SA[, , i]), type = "l",
          lwd = 1.5,  col = col, lty = lty,
          xlim = xlim_tradeoff,
          xlab = "Output Loss (%)",
          ylab = "Cumulative Deaths"
          )
  title(main = SensitivityTitle[i])
  legend("top", legend = SensitivityLegend[[i]], lwd = 2, col = col, lty = lty,
         cex = .9, x.intersp = 0, y.intersp=0.5, xpd = T, xjust=0, yjust=0, text.width = 1)
}
if(figure_save) saveas("Sentitivity.png", savedir, dev.args = dev.args)
par(mfrow = c(1, 1))


# --- Counterfactual analysis for 2020 (multiplicative deviation of alpha) ---
sw <- 2
IV_CF <- cbind(S[sw], I[sw], R[sw], D[sw])
AlphaDeviation2 <- seq(0.5, 1.3, by = .01)
CF2 <- zeros(2, length(AlphaDeviation2))
for (i in 1:length(AlphaDeviation2)){
  alpha_CF2 <- alpha[sw:(NROW(alpha)-2)] * AlphaDeviation2[i]
  c(CF2[1, i], CF2[2, i]) %<-% Covid_projection(IV_CF, alpha_CF2, beta[sw:(NROW(beta) - 2)], gamma * ones(Tdata - sw - 1, 1), delta[sw:(NROW(delta)-2)], zeros(Tdata - sw - 1, 1), h, k, POP0, hconstant)[1:2]
}

wedge <- 100 * mean(alpha[2:(NROW(alpha) - 2)]) - AverageAlpha2020
par(family = fontfamily)
x_ <- 100 * CF2[2, ] - wedge
y_ <- CF2[1, ]
xlim <- c(3, 5)
ylim <- c(0, 40000 * 1.1)
plot(
  x_, y_,
  type = "l", lwd = 2,
  xlim = xlim, ylim = ylim,
  xlab = 'Output Loss (%)',
  ylab = 'Cumulative Deaths'
)
par(new = T)
plot(x_[AlphaDeviation2 == 1], y_[AlphaDeviation2 == 1],
  lwd = 1.5, col = "red", pch = 19,
  axes = F, ann = F, xlim = xlim, ylim = ylim
)
if(figure_save) saveas("CF_multiplicative.png", savedir, dev.args = dev.args)

slope <-
  mean(cbind(CF2[1, 51] - CF2[1, 50]) / (100 * (CF2[2, 51] - CF2[2, 50])), (CF2[1, 52] -
                                                                              CF2[1, 51]) / (100 * (CF2[2, 52] - CF2[2, 51])))

SVL <- (POP0 / slope) / 100


# %%%%%%%%%%%%%%%%% Forecasting accuracy analysis %%%%%%%%%%%%%%%%%
StartDateF  <- 33
HorizonValsF <- matrix(c(1, 4), nrow = 1)
dNForecast <- zeros(Tdata, length(HorizonValsF))
dDForecast <- zeros(Tdata, length(HorizonValsF))
dNActual <- zeros(Tdata, length(HorizonValsF))
dDActual <- zeros(Tdata, length(HorizonValsF))


# --- Regress mobility on alpha to estimate the elasticity h ---


for (iH in 1:length(HorizonValsF)) {
  HorizonF  <- HorizonValsF[iH]
  EndtDateF <- Tdata + 1 - HorizonF
  
  for (iF in StartDateF:EndtDateF) {
    GDP_F <- GDP[1:(iF - 6)]
    
    referenceGDP_F <- referenceGDP[1:(iF - 1)]
    
    # --- Impute alpha (regress alpha on M)--- #
    Malt_F <- M[1:(iF - 1)]
    alpha_F <- matrix(1 - GDP[1:(iF - 6)] / referenceGDP[1:(iF - 6)])   # output loss in percentage
    X_F <- Malt_F[(iF - 6 - 17):(iF - 6)]
    Y_F <- alpha_F[(iF - 6 - 17):(iF - 6)]
    XC_F <- cbind(ones(length(X_F), 1), X_F)
    s_F <- solve(t(XC_F) %*% XC_F, t(XC_F) %*% Y_F) # OLS estimate of h with constant
    reg_F <- XC_F %*% s_F
    r_F <- Y_F - reg_F
    eps_p_F <- zeros(iF - 1 - (iF - 6), 1)
    eps_p_F[1] <- tail(r_F, 1)
    for (i in 1:(iF - 1 - (iF - 6) - 1)) {
      eps_p_F[i + 1] <- 1 * eps_p_F[i]
    }
    alpha_pred_F <- s_F[1] + s_F[2] * Malt_F[((iF - 6) + 1):(iF - 1)] + eps_p_F
    
    alpha_F <- rbind(alpha_F, alpha_pred_F)
    
    # --- Regress mobility on alpha to estimate the elasticity h ---#
    Y_F <- M[4:(iF - 6)]
    X_F <- alpha_F[4:(iF - 6)]
    if (hconstant == 0) {
      Y_F <- Y_F - 1
      h_all_F <- solve(t(X_F) %*% X_F, t(X_F) %*% Y_F) # OLS estimate of h
      reg_F <- X_F %*% h_all_F
      r_F <- Y_F - reg_F  # r is the residuals, which is the observed minus fitted values
      SSE_F <- sum(r_F ^ 2)            #SSE is the sum of squared errors
      MSE_F <- SSE_F / (length(Y_F) - 1)  # mean squared error
      h_all_se_F <- sqrt(MSE_F / sum(X_F ^ 2))  # standard error of h
    } else if (hconstant == 1) {
      XC_F <- cbind(ones(length(X_F), 1), X_F)
      h_all_F <- solve(t(XC_F) %*% XC_F, t(XC_F) %*% Y_F) # OLS estimate of h with constant
      reg_F <- XC_F %*% h_all_F
      r_F <- Y_F - reg_F
      SSE_F <- sum(r_F ^ 2)
    }
    h_F <- h_all_F
    
    # --- Compute the history of S, I, R, D in the data period --- #
    S_F <- zeros(iF, 1)
    I_F <- zeros(iF, 1)
    R_F <- zeros(iF, 1)
    D_F <- zeros(iF, 1)
    S_F[1] <- POP0
    for (i in 1:(iF - 1)) {
      S_F[i + 1] <- S_F[i] - N[i]
      I_F[i + 1] <- I_F[i] + N[i] - gamma * I_F[i] - dD[i]
      R_F[i + 1] <- R_F[i] + gamma * I_F[i]
      D_F[i + 1] <- D_F[i] + dD[i]
    }
    
    # --- Compute the history of time-varying parameters --- #
    
    delta_F <- (D_F[2:iF] - D_F[1:(iF - 1)]) / I_F[1:(iF - 1)]
    beta_tilde_F <- -POP0 * ((S_F[2:iF] - S_F[1:(iF - 1)]) / (S_F[1:(iF - 1)] * I_F[1:(iF - 1)]))   # overall infection rate
    
    if (hconstant == 0) {
      beta_F <- beta_tilde_F / (1 + h_all_F * alpha_F) ^ k                                      # raw infection rate
    }
    else if (hconstant == 1) {
      # %     beta = beta_tilde./(h(1)+h(2)*alpha).^k;
      beta_F <- beta_tilde_F / (1 + (h_all_F[2] / h_all_F[1]) * alpha_F) ^ k
    }
    
    IV_F <- cbind(S_F[iF], I_F[iF], R_F[iF], D_F[iF])
    alphaT_F <- alpha[iF:(iF + HorizonF - 1)] # Use actual values because we are interested in conditional forecasts.
    beta_average_F <- mean(beta_F[(iF - RetroPeriod):(iF - 1)])
    delta_average_F <- mean(delta_F[(iF - 1 - RetroPeriod):(iF - 1)])
    betaT_F <- beta_average_F * ones(HorizonF, 1)
    deltaT_F <- delta_average_F * ones(HorizonF, 1)
    V_F <- zeros(HorizonF, 1)
    gammaT_F <- gamma * ones(HorizonF, 1)
    
    c(CumD_F, AverageAlpha_F, SimData_F, SimN_F) %<-%
      Covid_projection(IV_F, alphaT_F, betaT_F, gammaT_F, deltaT_F, V_F, h_F, k, POP0, hconstant)
    
    dNForecast[iF, iH] <- sum(SimN_F[1:HorizonF])
    dDForecast[iF, iH] <- CumD_F - IV_F[4]
    dNActual[iF, iH] <- sum(N[iF:(iF + HorizonF - 1)])
    dDActual[iF, iH] <- D[iF + HorizonF - 1] - D[iF - 1]
    
  }
}
xtick1_F <- c(seq(StartDateF + 4, Tdata, by = 8), Tdata)
xtick4_F <- c(seq(StartDateF + 4, Tdata, by = 8), Tdata)
fs_F <- 10
fs_legend_F <- 10
fs_title <- 12

## ForecastErrorsN.png
par(family = fontfamily, mfrow = c(2, 2))
titles <- c(
  "Conditional Forecast vs. Actual\n (one-week horizon)",
  "Conditional Forecast vs. Actual\n (four-week horizon)",
  "Conditional Forecast Errors\n (one-week horizon)",
  "Conditional Forecast Errors\n (four-week horizon)"
  )
for (i in 1:4){
  if (i %in% 1:2){
    matplot(matrix(rep(StartDateF:Tdata, 2), ncol = 2),
            cbind(dNForecast[, i], dNActual[, i])[StartDateF:Tdata, ], type = "l",
            lwd = 1.5, col = c("red", "black"), lty = 1,
            xlim = c(StartDateF, Tdata), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    legend("topleft", legend = c('Forecast','Actual'), col = c("red", "black"), lty = 1,
           xjust=0, yjust=0, x.intersp=0, text.width=0.1, xpd = T)
  } else if (i == 3){
    plot(StartDateF:Tdata, dNForecast[StartDateF:Tdata,1]-dNActual[StartDateF:Tdata,1],
         type = "l", lwd = 1.5, xlim = c(StartDateF, Tdata), xaxt = "n", yaxt = "n", xlab = "", ylab = "") 
    abline(h = 0)
  } else if (i == 4){
    plot((StartDateF + 3):Tdata, dNForecast[StartDateF:(Tdata - 3),2]-dNActual[StartDateF:(Tdata - 3),2],
         type = "l", lwd = 1.5, xlim = c(StartDateF + 3, Tdata), xaxt = "n", yaxt = "n", xlab = "", ylab = "") 
    abline(h = 0)
  }
  title(main = titles[i])
  at_ <- if(i == 4) xtick4_F else xtick1_F
  axis(1, at = at_, labels = Month[at_])
  axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
}
if(figure_save) saveas("ForecastErrorN.png", savedir, dev.args = dev.args)

## ForecastErrorsD.png
par(family = fontfamily, mfrow = c(2, 2))
titles <- c(
  "Conditional Forecast vs. Actual\n (one-week horizon)",
  "Conditional Forecast vs. Actual\n (four-week horizon)",
  "Conditional Forecast Errors\n (one-week horizon)",
  "Conditional Forecast Errors\n (four-week horizon)"
)
for (i in 1:4){
  if (i %in% 1:2){
    matplot(matrix(rep(StartDateF:Tdata, 2), ncol = 2),
            cbind(dDForecast[, i], dDActual[, i])[StartDateF:Tdata, ], type = "l",
            lwd = 1.5, col = c("red", "black"), lty = 1,
            xlim = c(StartDateF, Tdata), xaxt = "n",  yaxt = "n", xlab = "", ylab = "")
    legend("topleft", legend = c('Forecast','Actual'), col = c("red", "black"), lty = 1,
           xjust=0, yjust=0, x.intersp=0, text.width=0.1, xpd = T)
  } else if (i == 3){
    plot(StartDateF:Tdata, dDForecast[StartDateF:Tdata,1]-dDActual[StartDateF:Tdata,1],
         type = "l", lwd = 1.5, xlim = c(StartDateF, Tdata), xaxt = "n", yaxt = "n", xlab = "", ylab = "") 
    abline(h = 0)
  } else if (i == 4){
    plot((StartDateF + 3):Tdata, dDForecast[StartDateF:(Tdata - 3),2]-dDActual[StartDateF:(Tdata - 3),2],
         type = "l", lwd = 1.5, xlim = c(StartDateF + 3, Tdata), xaxt = "n", yaxt = "n", xlab = "", ylab = "") 
    abline(h = 0)
  }
  title(main = titles[i])
  at_ <- if(i == 4) xtick4_F else xtick1_F
  axis(1, at = at_, labels = Month[at_])
  axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
}
if(figure_save) saveas("ForecastErrorD.png", savedir, dev.args = dev.args)
par(mfrow = c(1, 1))

# %%%%%%%%%%%%%%%%% Uncertainty band for the tradeoff curve %%%%%%%%%%%%%%%%%

# --- Uncertainty-band parameters --- #
beta_se  <- sd(beta_sample) / sqrt(length(beta_sample))
delta_se <- sd(delta_sample) / sqrt(length(delta_sample))
# "h_se" is computed above.
if (!iDrawUB) {
  ## TODO: 画像すらないのでこのブロック確認ができない
  BetaValsUB <-
    beta_average + beta_se * c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
  DeltaValsUB <-
    delta_average + delta_se * c(-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)
  HValsUB <- h_se %*% t(c(1, 0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75, -1)) + rep(h, 9)  ## たぶん sweep() のほうがいい
  UB <- zeros(length(BetaValsUB), length(TargetAlpha))
  # matrix to contain the results of uncertainty-band analysis
  
  for (i in 1:length(TargetAlpha)) {
    # alphaT = flip(0:0.01*TargetAlpha(i)*2/(SimPeriod-1):(0.01*TargetAlpha(i)*2))';
    alphaT <- zeros(SimPeriod, 1)
    alphaT[1:26] <-
      t(flip(seq(
        0, (0.01 * TargetAlpha[i] * 4), by = 0.01 * TargetAlpha[i] * 4 / (26 - 1)
      )))
    for (j in 1:length(BetaValsUB)) {
      betaS  <- BetaValsUB[j] * ones(SimPeriod, 1)
      deltaS <-  DeltaValsUB[j] * ones(SimPeriod, 1)
      HS <- HValsUB[, j]
      UB[j, i] <- Covid_projection(InitialValues, alphaT, betaS, gammaT, deltaS, V, HS, k, POP0, hconstant)[[1]]
    }
  }
  Y1UB <- zeros(length(AverageAlpha), 9)
  X1UB <- 100 * AverageAlpha
  Y1UB[, 1] <- UB[1, ]
  for (iUB in 2:9) {
    Y1UB[, iUB] <- UB[iUB, ] - UB[iUB - 1,]
  }
  # %--- Trade-off figure with UB (baseline) ---%
  # AreaUB <- area(X1UB,Y1UB,0)
  par(family = fontfamily)
  plot(X1UB[1, ], CumD, type = "l", lwd = 1.5,
       xlim = xlim_tradeoff, ylim = c(0, 50000), xlab = "", ylab = "", yaxt = "n")
  for(k in 1:9){
    grad <- abs((5 - k) / 10)
    polygon(x = c(X1UB, rev(X1UB)),
            y = c(if(k == 1) rep(0, NCOL(X1UB)) else  UB[k - 1, ], rev(UB[k, ])),
            xlim = xlim_tradeoff, ylim = c(0, 50000), lty = 0,
            col = rgb(grad, grad, grad, alpha = .3))
  }
  title(xlab = 'Output Loss (%)', ylab = 'Cumulative Deaths')
  abline(h = D[NROW(D) - 2], lwd = 1.5, col = "red", lty = 2)
  axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
  if(figure_save) saveas("BaselineTradeoffUB.png", savedir, dev.args = dev.args)
} else if (iDrawUB){
  # %%%%%%%%%%%%%%%%% Uncertainty band for the tradeoff curve (with simulation) %%%%%%%%%%%%%%%%%
  ## TODO: インデックスがおかしい?
  ## UB2 <- zeros(9, Nsim)
  UB2 <- zeros(Nsim, length(TargetAlpha))
  #  matrix to contain the results of uncertainty-band analysis
  
  h1 <- h[1]
  h2 <- h[2]
  h_se1 <- h_se[1]
  h_se2 <- h_se[2]
  
  
  for (i in 1:length(TargetAlpha)){
    # alphaT = flip(0:0.01*TargetAlpha(i)*2/(SimPeriod-1):(0.01*TargetAlpha(i)*2))';        
    alphaT <- zeros(SimPeriod, 1)
    alphaT[1:26] = flip(seq(0, (0.01 * TargetAlpha[i] * 4), by = 0.01 * TargetAlpha[i] * 4 / (26 - 1)))
    ## base R には parfor に対応するものがない (ないわけではないが) ので一旦 for にする 
    for(j in 1:Nsim){
      # for j = 1:Nsim
      betaS  <- rnorm(1, beta_average, beta_se) * ones(SimPeriod, 1)
      deltaS <- rnorm(1, delta_average, delta_se) * ones(SimPeriod, 1)
      ## matalab では動かないらしい
      # HS(1) = normrnd(h(1),h_se(1)); % It does not work well with parfor
      # HS(2) = normrnd(h(2),h_se(2)); % It does not work well with parfor
      HS1 <- rnorm(1, h1, h_se1)
      HS2 <- rnorm(1, h2, h_se2)
      ## TODO: なぜ h と同じ次元ではないのか
      ## HS <- cbind(HS1, HS2)
      HS <- rbind(HS1, HS2)
      UB2[j, i] <- Covid_projection(InitialValues, alphaT, betaS, gammaT, deltaS, V, HS, k, POP0, hconstant)[[1]]
    }
  }

  UBp <- zeros(9, length(TargetAlpha)) # matrix to contain the results of uncertainty-band analysis
  Y1UB <- zeros(length(AverageAlpha), 9)
  X1UB <- 100 * AverageAlpha
  for (i in 1:length(TargetAlpha)){
    UBp[, i] <- quantile(UB2[, i], 1:9 / 10, names = F, type = 2)  # TODO: https://stackoverflow.com/questions/24764966/numpy-percentile-function-different-from-matlabs-percentile-function ここによるとこの方法らしいが, 結果が合わない. Wikipedia には https://mathworks.com/help/stats/prctile.html?lang=en に補間を使っていると書いているが, わかりづらい...
  }
  Y1UBp <- matrix(NA, nrow = NCOL(UBp), ncol = NROW(UBp))
  Y1UBp[,1] <- UBp[1, ]
  for (iUB in 2:9){
    ## R の描画は積み上げではないので絶対値を使用
    Y1UBp[, iUB] <- UBp[iUB, ]
  }
  # figure(66)
  par(family = fontfamily)
  plot(X1UB[1, ], CumD, type = "l", lwd = 1.5,
       xlim = xlim_tradeoff, ylim = c(0, 50000), xlab = "", ylab = "", yaxt="n")
  for(k in 1:9){
    grad <- abs((5 - k) / 5)
    polygon(x = c(X1UB, rev(X1UB)),
            y = c(if(k == 1) rep(0, NCOL(X1UB)) else  Y1UBp[, k - 1], rev(Y1UBp[, k])),
            xlim = xlim_tradeoff, ylim = c(0, 50000), lty = 0,
            col = rgb(grad, grad, grad, alpha = .3))
  }
  axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
  title(xlab = 'Output Loss (%)', ylab = 'Cumulative Deaths')
  abline(h = D[NROW(D) - 2], lwd = 1.5, col = "red", lty = 2)
  
  # end of if-else statement for iDrawUB
}
if(figure_save) saveas("BaselineTradeoffUBp.png", savedir, dev.args = dev.args)

## 一括 save は不可能.
