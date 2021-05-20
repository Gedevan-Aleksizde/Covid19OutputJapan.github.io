require(readr)
require(zeallot)
source("R/wrappers.R")
source("R/Covid_projection.R")
source("R/Covid_projection_control.R")

# --- Graphical parameters ---

figure_save <- T
dev.args <- list()
fontfamily <- "Noto Sans CJK JP"
savedir <- "R/img/Tokyo"
xlim_tradeoff <- c(3, 4)

# ---- Model parameter values ----

pref <- 'Tokyo'  #prefecture to be analyzed
gamma <- 7 / 5  # recovery rate from Covid
k <- 2  # exponent of (1-h*alpha)
hconstant <- 1  # 0 = without intercept, 1 = with intercept for h regression
SimPeriod <- 52  # simulation period in weeks
TargetAlpha <- seq(.1, 8, by = .1) # values of alpha we simulate results
AlphaVals <- c(1.5, 0.85 * (1.5 + 3) / 2, 3)  # benchmark alpha we plot time-series figures
VacStart <- 11  # time until the start of vaccination process
VacPace <- 0.1 * 1.33 * 0.8 * 1500000  # number of vaccinations per week
VacDuration <- 12  # time until vaccination pace reaches its steady pace
RetroPeriod <- 17  # retroactive periods used to estimate gamma and delta
alpha_on <- 4.2  # 3.75 = 12 weeks, 4.2 = 8 weeks, 5.5 = 4 weeks
th_on <- 14000  # threshold to place the state of emergency (daily new infected persons in Tokyo = 2000)
th_off <- 3500  # threshold to remove the state of emergency (daily new infected persons in Tokyo = 500)
target_duration <- 8
wl <- c(1, 2)  # Results as of these weeks ago

# ---- import data ----

# Covid data are recorded at weekly frequencies (Mon-Sun)
# The first week start on January 20 (Mon), 2020

## covid = importdata([home '\Data\Covid_weekly_prefecture.csv']);  % Import weekly Covid data by prefecture
covid <- read_csv("Covid19_Output/Data/Covid_weekly_prefecture.csv")
Data <- subset(covid, prefecture == pref)
# Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility, 5 = GDP, 6 = population
date <- Data$date + 21916
N <- Data$positive
dD <- Data$death
M <- Data$mobility
GDP <- Data$gdp
POP <- Data$population
Tdata <- NROW(Data)  # Data period in weeks
POP0 <- POP[1]  # initial population
xtick1 <- c(seq(1, Tdata, by = 13), Tdata)
date <- as.Date(date, origin = "1899-12-30") ## TODO: Exel日付値でよい?

## 以下は使わない
# dateJP = string(datetime(date,'ConvertFrom','excel','Format','M? dd?'));
# dateEN = string(datetime(date,'ConvertFrom','excel','Format','MMM-dd'));
# Month = string(datetime(date,'ConvertFrom','excel','Format','MMM-yy'));

M <- 1 + 0.01 * M
TdataGDP <- Tdata - sum(is.na(GDP))
RetroH <- TdataGDP - 4

# ---- Constructing the reference level of output ----

potentialGDP <- zeros(52 * 3, 1)
potentialGDP[1] <- (100 / (1.0122)) * (1.0063 ^ (1 / 12))
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
## TODO: これは何
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
axis(1, at = xtick1, labels = format(date[xtick1], "%m-%Y"), las=2)
mtext("GDP", side = 4, line = 2.5, col = "blue")
legend("top", inset=c(0,0),
       legend = c("Mobility (left axis)", "GDP (right axis)"), col = c("black", "blue"), lty = c(1, 4))
if(figure_save) saveas("TokyoMobilityGDPLine.png", savedir, dev.args = dev.args)
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
AverageAlpha2020 <- 100 * mean(c(alpha[1], alpha[2], alpha))


# ---- Plot varabiles in the data period ----
## TODO: これも何のグラフ?
par(family = fontfamily)
plot(N, type = "l", xaxt = "n", yaxt = "n", xlab = "")
axis(1, at = xtick1, labels = format(date[xtick1], "%m-%Y"), las = 2)
axis(2, at = 0:4 * 10000, labels = formatC(0:4 * 10000, "%,6.0f", big.mark = ","))

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
  axis(1, at = xtick1, labels = format(date[xtick1], "%Y-%m"))
  if(i == 5){
    # ax.YAxis.Exponent = 0;
    axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
  } else {
    axis(2, at = axTicks(2))
  }
  if(i == length(VarList)){
    par(new = T)
    plot(referenceGDP, lty = 2, lwd = 2, type = "l", axes = F, ann = F,
         xlim = c(1, Tdata), ylim = ylim)
  }
}
par(mfrow = c(1, 1))

# %--- Plot parameters in the data period ---%

## TODO: これも確認手段がない?
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
  axis(1, at = xtick1, labels = format(date[xtick1], "%m-%Y"))
  title(main = ParamName[i], ylab = "")
}
if(figure_save) saveas("parameters.png", savedir, dev.args = dev.args)
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
VacPaceVals <- c(0.5*VacPace, VacPace, 2*VacPace)
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
for (i in 1:length(AlphaVals)){
  ## TODO: このコメントアウトの意図
  # %     alphaT = flip(0:0.01*AlphaVals(i)*2/(SimPeriod-1):(0.01*AlphaVals(i)*2))';
  alphaT <- zeros(SimPeriod, 1)
  alphaT[1:26] = flip(seq(0, 0.01 * AlphaVals[i] * 4, by = 0.01 * AlphaVals[i] * 4 / (26 - 1)))
  # %     alphaT = 0.01*AlphaVals(i)*ones(SimPeriod,1);
  
  if(hconstant == 0){
    ERNT <- (AlphaIndexVariables[, 1, i] / POP0) * (((1 + h * alphaT) ^ k) * betaT) / (gammaT + deltaT)
  } else if(hconstant == 1){
    ERNT <- (AlphaIndexVariables[, 1, i] / POP0) * (((1 + (h[2] / h[1]) * alphaT) ^ k) * betaT) / (gammaT + deltaT)
  }
  ProjectedData <- cbind(AlphaIndexVariables[, 5, i], AlphaIndexVariables[, 6, i], ERNT, V)
  CombinedData <- rbind(Past, ProjectedData)
  par(family = fontfamily, mfrow = c(2, 2))
  for (j in 1:length(VarList2)){
    plot(CombinedData[, j], type = "l", xaxt = "n", xlim = c(1, Tdata + SimPeriod), ylim = if(j == 3) c(0, 3))
    title(main = VarName2[j], xlab = "", ylab = "")
    abline(v = Tdata)
    axis(1, at = c(1, 27, 53, 79, 98), labels = c('Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'))
    if (j == 1){
      axis(2, at = c(2e4, 4e4, 6e4, 8e4, 10e4, 12e4))
      axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
    }
    if (j == 3){
      abline(h = 1)
    }
    if (j == 4){
      axis(2, at = c(5e5, 10e5, 15e5, 20e5, 25e5), labels = c('0.5M','1M','1.5M','2M','2.5M'))
    }
  }
  
  alpha2 <- rbind(alpha, alphaT)
  beta2 <- rbind(beta, betaT)
  if (hconstant == 0){
    beta_tildeT <- betaT * (1 + h * alphaT) ^ k             # raw infection rate 
  } else if (hconstant == 1){
    ## TODO: このコメントアウトの意図
    # % beta_tildeT = betaT.*(h(1)+h(2)*alpha).^k;
    beta_tildeT <- betaT * (1 + (h[2] / h[1]) * alphaT) ^ k    
  }
  beta_tilde2 <- rbind(beta_tilde, beta_tildeT)
  ERN2 <- rbind(ERN, ERNT)
  delta2 <- rbind(delta, deltaT)
  V2 <- rbind(zeros(Tdata, 1), V)
  # %set(gcf,'Position',[100,100,800,500])
  par(mfrow = c(2, 2))
  ## TODO:
  # %set(gcf,'Position',[100,100,800,500])
  for (j in 1:length(ParamList2)){
    if (j == 2) ylim <- c(0, 2)
    else if (j == 3) ylim <- c(0, 5)
    else if (j == 4) ylim <- c(0, .1)
    else ylim <- NULL
    plot(get(ParamList2[j]), type = "l", xlab = "",  ylab = "", xlim = c(1, Tdata + SimPeriod), xaxt = "n", ylim = ylim)
    title(main = ParamName[j], xlab = "", ylab = "")
    abline(v = Tdata)
    axis(1, at = c(1, 27, 53, 79, 98), labels = c('Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'))
  }
}
par(mfrow = c(1, 1))

# %--- Trade-off figure (baseline) ---%
par(family = fontfamily)
plot(100 * AverageAlpha, CumD, type = "l",
     xlim = xlim_tradeoff, ann = F, yaxt = "n")
abline(h = D[NROW(D) - 2])
title(xlab = "Output Loss (%)", ylab = "Cumulative Deaths")
axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))

 
# ---- Trade-off figure (lag) ----
par(family = fontfamily)
plot(100 * AverageAlpha, CumD, type = "l", ann = F,
     xlim = xlim_tradeoff, ylim = c(0, 60000), yaxt = "n")
for(i in 1:length(wl)){
  par(new = T)
  plot(100 * LagResults[2, , i], LagResults[1, , i], col = i + 1, type = "l",
       ann = F, axes = F, yaxt = "n")
}
abline(v = AverageAlpha2020)
legend("top", legend = c("baseline",
                         sprintf('%.0f week%s ago', 1:length(wl), ifelse(1:length(wl) == 1, "", "s"))),
       col = c(1:(length(wl) + 1)), lty = 1, cex = .9)
title(xlab = "Output Loss (%)", ylab = "Cumulative Deaths")
## TODO: グラフには表示されてない?
# text(3.5,100000,{'Output loss';'      in 2020'},'FontSize',16);
axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))



# %%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%
alpha_off <- mean(alpha[format(date, "%Y-%m") %in% c("2020-09", "2020-10", "2020-11")]) # output loss without the state of emergency
InitialValues <- cbind(tail(S, 1), tail(I, 1), tail(R, 1), tail(D, 1))
dateP <- seq(tail(date, 1) + 7, tail(date, 1) + 7 * (SimPeriod + 1), by = 7)
## R では不要
# dateP=string(datetime(dateP,'ConvertFrom','excel','Format','MMM'))

# ---- Construct time series of parameters ----
beta_sample <- beta[(length(beta) - RetroPeriod + 1):length(beta)]
betaT <- mean(beta_sample) * ones(SimPeriod, 1)
delta_sample <- delta[(length(delta) - RetroPeriod + 1):length(beta)]
deltaT <- mean(delta_sample) * ones(SimPeriod, 1)
gammaT <- gamma*ones(SimPeriod, 1)
V <- zeros(SimPeriod, 1)
V[(VacStart - 1):(VacStart + 3)] <- seq(0, VacPace, by = VacPace / 4)
V[(VacStart + 4):NROW(V)] <- VacPace

# ---- Construct time series of parameters ----
beta_sample <- beta[(length(beta) - RetroPeriod + 1):length(beta)]
betaT <- mean(beta_sample) * ones(SimPeriod, 1)
delta_sample <- delta[(length(delta) - RetroPeriod + 1):length(delta)]
deltaT <- mean(delta_sample) * ones(SimPeriod, 1)
gammaT <- gamma * ones(SimPeriod, 1)
V <- zeros(SimPeriod, 1)
V[(VacStart - 1):(VacStart + 3)] <- seq(0, VacPace, by = VacPace / 4) ## 0:VacPace/4:VacPace;
V[(VacStart + 4):NROW(V)] <- VacPace

# --- Projection for different th_off values ----
TH <- seq(100, 1000, by = 50) * 7
TH_index <- 500 * 7
DM <- zeros(1, length(TH))
AlphaM <- zeros(1, length(TH))
AlphaPath <- zeros(SimPeriod, length(TH))
NPath <- zeros(SimPeriod, length(TH))

# ---- Projection for different th_off values ----
TH <- seq(100, 1000, by = 50) * 7
TH_index <- 500 * 7
DM <- zeros(1, length(TH))
AlphaM <- zeros(1, length(TH))
AlphaPath <- zeros(SimPeriod, length(TH))
NPath <- zeros(SimPeriod, length(TH))


# ---- Plot three scenarios ----
altA <- c(0.085, 0.074, 0.0695)

figurenames <- list(
  "12" = "BaselineDecline.png",
  "13" = "RapidDecline.png",
  "14" = "GradualDecline.png",
  "15" = "ThreeScenariosDecline.png",
  "16" = "ThreeVaccinationsDecline.png"
)


for(y in 1:length(altA)) {
  alpha_on = altA[y]
  
  ## TODO:
  # set(gcf,'Position',[100,100,1200,500])
  par(family = fontfamily, mfrow = c(1, 2), cex = .8)
  for (i in 1:length(TH)) {
    c(DM[i], AlphaM[i], AlphaPath[,], SimData, NPath[, i]) %<-% Covid_projection_control(InitialValues, alpha_on, alpha_off, th_on, TH[i], betaT, gammaT, deltaT, V, h, k, POP0, hconstant)
    if (sum(TH[i] == TH_index) == 1) {
      plot(NPath[, i], lwd = 1.5, col = "red", lty = 3, type = "l", xaxt = "n", yaxt = "n", ann = F, xlim = c(1, SimPeriod))
    } else {
      plot(NPath[, i], lwd = .4, col = "blue", lty = 4, type = "l", xaxt = "n", yaxt = "n", ann = F, xlim = c(1, SimPeriod))
    }
    par(new = T)
  }
  legend("top", legend = sprintf('%.0f', TH / 7), lty = c(1, 3:4), cex = .4, ncol = 2)
  title(main = "Projected path of new cases", ylab = "Number of new cases per week")
  axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
  if(y == 1){
    ## TODO:
    # lgd = legend;
    # lgd.NumColumns = 2;
  }
  axis(
    1,
    at = c(1, 9, 17, 25, 33, 41, 49),
    labels = dateP[c(1, 9, 17, 25, 33, 41, 49)]
    )
  ## TODO:
  # % xticks([1 14 27 40 52])
  # % xticklabels(dateP([1 14 27 40 52]))
  
  ## --- Record how many times on and off are triggered ---
  
  waves <- zeros(1, length(TH))
  for (i in 1:length(TH)){
    svec <- zeros(SimPeriod - 1, 1)
    for (t in 1:(SimPeriod - 1)){
      svec[t] <- AlphaPath[t + 1, i] - AlphaPath[t, i]
    }
    waves <- nnz(svec)
  }
  ## --- Number of cumulative deaths ---
  AlphaM <- AlphaM * 100
  par(new = F)
  plot(AlphaM[waves == 1], DM[waves == 1], lty = 2, col = "blue", type = "b", pch = 21, axes = F)
  par(new = T)
  plot(AlphaM[waves == 3], DM[waves == 3], lty = 2, col = "gree", type = "b", pch = 21, axes = F)
  par(new = T)
  plot(AlphaM[waves == 5], DM[waves == 5], lty = 2, col = "magenta", type = "b", pch = 21, axes = F)
  text(AlphaM, DM, TH / 7, offset = 1, cex = .5)
  plot(AlphaM[TH == 500 * 7], DM[TH == 500 * 7], pch = 21, col = "red")
  title(
    xlab = "Output Loss (%)",
    ylab = "Cumlative Deaths",
    main = "Relationsihp between Covid-19 and output"
    )
  if (figure_save){
    if (10 + y %in% 13:16) {
      saveas(paste0(figurenames[[as.character(y)]], ".svg"), savedir, dev.args = dev.args)
    }
  }
  
  ## TODO: このコメントアウトの意図 
  # %     %--- Number of new infections ---%
  # %     NM = sum(NPath,1)+sum(N);
  # %     subplot(1,2,2)
  # %     plot(AlphaM(waves==1),NM(waves==1),'-bo','LineWidth',1.5,'MarkerSize',10)
  # %     hold on
  # %     plot(AlphaM(waves==3),NM(waves==3),'-go','LineWidth',1.5,'MarkerSize',10)
  # %     hold on
  # %     plot(AlphaM(waves==5),NM(waves==5),'-mo','LineWidth',1.5,'MarkerSize',10)
  # %     hold on
  # %     text(AlphaM,NM,string(TH/7),'VerticalAlignment','bottom','HorizontalAlignment','left')
  # %     hold on
  # %     scatter(AlphaM(TH==500*7),NM(TH==500*7),80,'red','filled');
  # %     xlabel('Output Loss (%)','FontSize',20)
  #     %     ylabel('Cumlative number of infected people','FontSize',20)
  #     %     title('Relationsihp between Covid-19 and output','FontSize',20,'FontWeight','normal')
  #     %     grid on
  #     %     ax = gca;
  #     %     ax.YAxis.FontSize = 20;
  #     %     ax.XAxis.FontSize = 20;
  #     %     ax.YAxis.Exponent = 0;
  #     %     ytickformat('%,6.0f')
}

# ---- Combine three scenarios ----

for(y in 1:length(altA)) {
  alpha_on <- altA[y]
  par(family = fontfamily, mfrow = c(1, 2), ann = F)
  for (i in 1:length(TH)) {
    c(DM[i], AlphaM[i], AlphaPath[, i], SimData, NPath[, i]) %<-% Covid_projection_control(InitialValues, alpha_on, alpha_off, th_on, TH[i], betaT, gammaT, deltaT, V, h, k, POP0, hconstant)
  }
  ## for で書くのはさすがにつらい
  matplot(NPath, type = "l",
          col = ifelse(TH == TH_index, "red", "blue"),
          lwd = ifelse(TH == TH_index, 2, .3),
          xaxt = "n"
          )
  legend("topright", legend = sprintf('%.0f', TH / 7),
         col = ifelse(TH == TH_index, "red", "blue"),
         lty = 2, ncol = 2,
         cex = .7, x.intersp = 0, y.intersp=0.5, xpd = T, xjust=0, yjust=0, text.width = 1.2)
  title(main = "Projected path of new cases",
        ylab = "Number of new cases per week")
  axis(
    1,
    at = c(1, 9, 17, 25, 33, 41, 49),
    labels = dateP[c(1, 9, 17, 25, 33, 41, 49)]
    )
  
  # ---- Record how many times on and off are triggered ----
  
  waves <- zeros(1, length(TH))
  for (i in 1:length(TH)) {
    svec <- zeros(SimPeriod - 1, 1)
    for (t in 1:(SimPeriod - 1)) {
      svec[t] <- AlphaPath[t + 1, i] - AlphaPath[t, i]
    }
    waves[i] <- nnz(svec)
  }
  # ---- Number of cumulative deaths ----
  AlphaM <- AlphaM * 100
  xlim <- c(min(AlphaM), max(AlphaM))
  ylim <- c(min(DM), max(DM))
  plot(AlphaM[waves == 1], DM[waves == 1], lty = 2, col = "blue", type = "b", pch = 21,
       xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
  par(new = T)
  plot(AlphaM[waves == 3], DM[waves == 3], lty = 2, col = "green", type = "b", pch = 21,
       xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
  par(new = T)
  plot(AlphaM[waves == 5], DM[waves == 5], lty = 2, col = "magenta", type = "b", pch = 21,
       xlim = xlim, ylim = ylim, yaxt = "n")
  text(AlphaM, DM, TH / 7, offset = 1, cex = .5)
  par(new = T)
  plot(AlphaM[TH == 500 * 7], DM[TH == 500 * 7], pch = 21, col = "red", 
       xlim = xlim, ylim = ylim, yaxt = "n")
  title(
    xlab = "Output Loss (%)",
    ylab = "Cumlative Deaths",
    main = "Relationsihp between Covid-19 and output"
    )
  axis(2, at = axTicks(2), labels = formatC(axTicks(2), "%6.0f", big.mark = ","))
}
par(mfrow = c(1, 1), new = F)

# ---- Three vaccination scenarios ----
altV <- c(0.5, 1, 2)
alpha_on <- 0.074
par(family = fontfamily, mfrow = c(1, 2))
for(y in 1:length(altV)){
  VacPaceAlt <- altV[y]*VacPace
  V[(VacStart-1):(VacStart+3)] <- seq(0, VacPaceAlt, by = VacPaceAlt/4)
  V[(VacStart+4):NROW(V)] <- VacPaceAlt
  for(i in 1:length(TH)){
    c(DM[i],AlphaM[i],AlphaPath[,i],SimData,NPath[,i]) %<-% Covid_projection_control(InitialValues,alpha_on,alpha_off,th_on,TH[i],betaT,gammaT,deltaT,V,h,k,POP0,hconstant)
  }
  matplot(NPath, type = "l",
          col = ifelse(TH == TH_index, "red", "blue"),
          lwd = ifelse(TH == TH_index, 2, .3),
          xaxt = "n"
  )
  legend("topright", legend = sprintf('%.0f', TH / 7),
         col = ifelse(TH == TH_index, "red", "blue"),
         lty = 2, ncol = 2,
         cex = .7, x.intersp = 0, y.intersp=0.5, xpd = T, xjust=0, yjust=0, text.width = 1.2)
  title(main = "Projected path of new cases",
        ylab = "Number of new cases per week")
  axis(1, at = c(1, 9, 17, 25, 33, 41, 49), labels = dateP[c(1, 9, 17, 25, 33, 41, 49)])

  # ---- Record how many times on and off are triggered ----
  waves <- zeros(1, length(TH))
  for(i in 1:length(TH)){
    svec <- zeros(SimPeriod-1,1)
    for(t in 1:(SimPeriod -1)){
      svec[t] <- AlphaPath[t+1,i]-AlphaPath[t,i]
    }
    waves[i] <- nnz(svec)
  }
  # ---- Number of cumulative deaths ----
  AlphaM <- AlphaM*100
  xlim <- c(min(AlphaM), max(AlphaM))
  ylim <- c(min(DM), max(DM))
  par(family = fontfamily, ann = F)
  plot(AlphaM[waves==1], DM[waves==1], lty = 2, col = "blue", type = "b", pch = 21, xaxt = "n", yaxt = "n", ann = F, xlim = xlim, ylim = ylim)
  par(new = T)
  plot(AlphaM[waves==3], DM[waves==3], lty = 2, col = "green", type = "b", pch = 21, xaxt = "n", yaxt = "n", ann = F, xlim = xlim, ylim = ylim)
  par(new = T)
  plot(AlphaM[waves==5], DM[waves==5], lty = 2, col = "magenta", type = "b", pch = 21, xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim)
  text(AlphaM, DM, TH/7, offset = 1, cex = .5)
  par(new = T)
  plot(AlphaM[TH==500*7], DM[TH==500*7], pch = 21, col = "red", xlim = xlim, ylim = ylim)
  title(xlab = "Output Loss (%)", ylab = "Cumlative Deaths",
        main = "Relationsihp between Covid-19 and output")
}
par(mfrow = c(1, 1), new = F)




