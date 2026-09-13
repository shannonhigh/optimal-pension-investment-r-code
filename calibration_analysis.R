# Initial calibration implementation supplied by the project supervisor. The accompanying anonymised P+ observations and market inputs relate to the 2023 project data. The 2023 market assumptions are from Radet for Afkastforventninger. Subsequent sample diagnostics, restrictions and extensions used in this thesis are documented below.

# Settings
library(deSolve)
library(expm)
library(scales)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(latex2exp)
library(limSolve)
#library(lamW)

# Options
options(scipen = 999)
options(digits = 16)
scale = 100
linetypes = c("solid", "dashed", "dotted", "dotdash")

data_file <- "Samlet information til JP.csv"
figure_dir <- "figures"
results_dir <- "results"

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(data_file)) {stop(paste0("Input file not found: '", data_file, "'. ", "Set the working directory to the project root and place the anonymised CSV there before running this script."))}

# Aggregate function to bonds and stocks
pf_fordelingsfunktion = function(A)
  rbind(sum(A[1:4,1]), sum(A[5:6,1]), sum(A[7:10,1]))

# Correlation matrix from Radet for Afkastforventninger
correlation_matrix = matrix(c(1, 0.6, 0.1, 0.3, -0.1, -0.1, -0.2, -0.1, -0.1, -0.1,
                              0.6, 1, 0.6, 0.6, 0.2, 0.2, 0.2, 0.1, 0.1, 0.3,
                              0.1, 0.6, 1, 0.7, 0.7, 0.6, 0.6, 0.4, 0.3, 0.7,
                              0.3, 0.6, 0.7, 1, 0.5, 0.6, 0.4, 0.2, 0.2, 0.5, 
                             -0.1, 0.2, 0.7, 0.5, 1, 0.7, 0.8, 0.4, 0.4, 0.8,
                             -0.1, 0.2, 0.6, 0.6, 0.7, 1, 0.7, 0.4, 0.4, 0.7,
                             -0.2, 0.2, 0.6, 0.4, 0.8, 0.7, 1, 0.4, 0.4, 0.7,
                             -0.1, 0.1, 0.4, 0.2, 0.4, 0.4, 0.4, 1, 0.3, 0.4,
                             -0.1, 0.1, 0.3, 0.2, 0.4, 0.4, 0.4, 0.3, 1, 0.4,
                             -0.1, 0.3, 0.7, 0.5, 0.8, 0.7, 0.7, 0.4, 0.4, 1), ncol = 10, nrow = 10)

# Volatilities from Radet for Afkastforventninger
std_dev = matrix(nrow = 10, ncol = 1, data = c(3.7, 5.5, 11.9, 10.7, 15.3, 21.7, 20.4, 14.0, 10.8, 9.4) / scale)

# Variance matrix
variance_matrix = function(t = 0) {if (t < 11) {variance_matrix = matrix(rep(0, 100), nrow = 10, ncol = 10)
for (i in 1:10) {for (j in 1:10) {variance_matrix[i,j] = correlation_matrix[i,j] * std_dev[i] * std_dev[j]}}
return(variance_matrix)}
  
if (t >= 11) {return(matrix(data = c(0.08^2, 0, 0, 0, 0.18^2, 0, 0, 0, 0.12^2), nrow = 3, ncol = 3))}}
returns = function(t = 0) {
if (t < 6) {return(matrix(nrow = 10, ncol = 1, data = c(1.9, 2.2, 4.9, 4.3, 6.1, 8.3, 10.2, 5.6, 4.1, 3.8) / scale))}
if (t >= 6 && t < 11) {return(matrix(nrow = 10, ncol = 1, data = c(2.1, 3.1, 5.0, 4.7, 6.3, 8.8, 10.5, 5.9, 5.1, 4.3) / scale))}
if (t >= 11) {return(matrix(nrow = 3, ncol = 1, data = c(3.5, 6.5, 6.5) / scale))}}

# Load data
data = read.csv(data_file, header = T, sep = ";", dec = ",", colClasses = c("character", "character", "character", "character", "character", "character", "numeric", "numeric", "numeric", "character"))

# Internal P+ profile coding used throughout the implementation:
# "hoj" -> low risk aversion / higher-risk investment profile
# "mellem" -> moderate risk aversion
# "lav" -> high risk aversion / lower-risk investment profile

data$investeringsProfil[data$investeringsProfil != "P+ Livscyklus mellem" & data$investeringsProfil != "P+ Livscyklus lav"] = "hoj"
data$investeringsProfil[data$investeringsProfil == "P+ Livscyklus mellem"] = "mellem"
data$investeringsProfil[data$investeringsProfil == "P+ Livscyklus lav"] = "lav"

# Plot data
# Active P211R policyholders with positive pension wealth.
test_data = subset(data, data$gruppeNavn == "P211R" & data$policeStatus == "Betalende" & data$depot > 0)

# Age plot
(d1 = ggplot(as.data.frame(test_data$alder), aes(x = test_data$alder)) + geom_histogram(bins = 30, col = "grey") + ylab("Frequency") + xlab("Age") + theme_bw())
ggsave(filename = file.path(figure_dir, "Age Distribution of P211R Policyholders.png"), plot = d1, width = 7, height = 5, dpi = 300)

# Risk profile plot
(d2 = ggplot(test_data, aes(x = fct_relevel(investeringsProfil, "hoj", "mellem", "lav"), group = investeringsProfil)) + geom_bar(col = "grey") + ylab("Frequency") + xlab("Risk profile") + scale_x_discrete(labels = c(hoj = "Low risk aversion", mellem = "Moderate risk aversion", lav = "High risk aversion")) + theme_bw() + theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "none"))
ggsave(filename = file.path(figure_dir, "Distribution of Risk Profiles, P+ Livscyklus.png"), plot = d2, width = 7, height = 5, dpi = 300)

# Savings accounts plot
(d3 = ggplot(test_data, aes(x = depot / 1000000)) + geom_histogram(binwidth = 0.25, boundary = 0, colour = "grey40", fill = "grey75") + coord_cartesian(xlim = c(0, 5)) + labs(x = "Pension savings (DKK millions)", y = "Frequency") + theme_bw())
ggsave(filename = file.path(figure_dir, "Distribution of Pension Savings (Depot).png"), plot = d3, width = 7, height = 5, dpi = 300)

# The recorded portfolio weights do not always sum exactly to one, so each supplied allocation is rescaled to sum to one.
rescale_function = function(A) {
B = matrix(nrow = 10, ncol = 1, data = c(A[1,1] / sum(A), A[2,1] / sum(A), A[3,1] / sum(A), A[4,1] / sum(A), A[5,1] / sum(A), A[6,1] / sum(A), A[7,1] / sum(A), A[8,1] / sum(A), A[9,1] / sum(A), A[10,1] / sum(A)))
return(as.matrix(B))}

pension_om_30_hoj_original = matrix(nrow = 10, ncol = 1, data = c(5.0, 14.0, 25.0, 6.0, 69.0, 11.0, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_30_hoj = rescale_function(pension_om_30_hoj_original)
pension_om_30_mellem_original = matrix(nrow = 10, ncol = 1, data = c(10.8, 15.8, 23.2, 6.0, 57.0, 9.0, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_30_mellem = rescale_function(pension_om_30_mellem_original)
pension_om_30_lav_original = matrix(nrow = 10, ncol = 1, data = c(23.3, 19.7, 19.3, 6.0, 31.2, 4.7, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_30_lav = rescale_function(pension_om_30_lav_original)
pension_om_15_hoj_original = matrix(nrow = 10, ncol = 1, data = c(5.0, 14.0, 25.0, 6.0, 69.0, 11.0, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_15_hoj = rescale_function(pension_om_15_hoj_original)
pension_om_15_mellem_original = matrix(nrow = 10, ncol = 1, data = c(10.8, 15.8, 23.2, 6.0, 57.0, 9.0, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_15_mellem = rescale_function(pension_om_15_mellem_original)
pension_om_15_lav_original = matrix(nrow = 10, ncol = 1, data = c(23.3, 19.7, 19.3, 6.0, 31.2, 4.7, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_15_lav = rescale_function(pension_om_15_lav_original)
pension_om_5_hoj_original = matrix(nrow = 10, ncol = 1, data = c(14.0, 16.8, 22.2, 6.0, 50.4, 7.9, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_5_hoj = rescale_function(pension_om_5_hoj_original)
pension_om_5_mellem_original = matrix(nrow = 10, ncol = 1, data = c(21.8, 19.2, 19.8, 6.0, 34.3, 5.2, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_5_mellem = rescale_function(pension_om_5_mellem_original)
pension_om_5_lav_original = matrix(nrow = 10, ncol = 1, data = c(30.4, 21.9, 17.1, 6.0, 16.5, 2.3, 4.0, 6.0, 10.0, 6.0) / scale)
pension_om_5_lav = rescale_function(pension_om_5_lav_original)

pension_efter_5_hoj_original = matrix(nrow = 10, ncol = 1, data = c(18.3, 18.1, 20.9, 6.0, 41.4, 6.4, 4.0, 6.0, 10.0, 6.0) / scale)
pension_efter_5_hoj = rescale_function(pension_efter_5_hoj_original)
pension_efter_5_mellem_original = matrix(nrow = 10, ncol = 1, data = c(25.2, 20.3, 18.7, 6.0, 27.1, 4.0, 4.0, 6.0, 10.0, 6.0) / scale)
pension_efter_5_mellem = rescale_function(pension_efter_5_mellem_original)
pension_efter_5_lav_original = matrix(nrow = 10, ncol = 1, data = c(32.9, 22.7, 16.3, 6.0, 11.2, 1.4, 4.0, 6.0, 10.0, 6.0) / scale)
pension_efter_5_lav = rescale_function(pension_efter_5_lav_original)

# Defines the observed P+ life-cycle portfolio as a function of age and risk profile using linear interpolation.
pf = function(alder, risiko) {prod = matrix(data = list(pension_om_30_hoj, pension_om_15_hoj, pension_om_5_hoj, pension_efter_5_hoj, pension_om_30_mellem, pension_om_15_mellem, pension_om_5_mellem, pension_efter_5_mellem, pension_om_30_lav, pension_om_15_lav, pension_om_5_lav, pension_efter_5_lav), nrow = 4, ncol = 3)

if (alder < 50 & risiko == "hoj")
  x = c(1, 1, 0, 50, 0)
  
if (alder < 50 & risiko == "mellem")
  x = c(1, 2, 0, 50, 0)
  
if (alder < 50 & risiko == "lav")
  x = c(1, 3, 0, 50, 0)
  
if (alder < 60 & alder >= 50 & risiko == "hoj")
  x = c(2, 1, 50, 60, 1)
  
if (alder < 60 & alder >= 50 & risiko == "mellem")
  x = c(2, 2, 50, 60, 1)
  
if (alder < 60 & alder >= 50 & risiko == "lav")
  x = c(2, 3, 50, 60, 1)
  
if (alder < 70 & alder >= 60 & risiko == "hoj")
  x = c(3, 1, 60, 70, 1)
  
if (alder < 70 & alder >= 60 & risiko == "mellem")
  x = c(3, 2, 60, 70, 1)
  
if (alder < 70 & alder >= 60 & risiko == "lav")
  x = c(3, 3, 60, 70, 1)
  
if (alder >= 70 & risiko == "hoj")
  x = c(4, 1, 70, 100, 0)
  
if (alder >= 70 & risiko == "mellem")
  x = c(4, 2, 70, 100, 0)
  
if (alder >= 70 & risiko == "lav")
  x = c(4, 3, 70, 100, 0)
  
test = function(t) {
  class_1 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][1,1], prod[x[1]+x[5],x[2]][[1]][1,1]), rule = 2)
  class_2 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][2,1], prod[x[1]+x[5],x[2]][[1]][2,1]), rule = 2)
  class_3 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][3,1], prod[x[1]+x[5],x[2]][[1]][3,1]), rule = 2)
  class_4 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][4,1], prod[x[1]+x[5],x[2]][[1]][4,1]), rule = 2)
  class_5 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][5,1], prod[x[1]+x[5],x[2]][[1]][5,1]), rule = 2)
  class_6 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][6,1], prod[x[1]+x[5],x[2]][[1]][6,1]), rule = 2)
  class_7 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][7,1], prod[x[1]+x[5],x[2]][[1]][7,1]), rule = 2)
  class_8 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][8,1], prod[x[1]+x[5],x[2]][[1]][8,1]), rule = 2)
  class_9 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][9,1], prod[x[1]+x[5],x[2]][[1]][9,1]), rule = 2)
  class_10 = approxfun(c(x[3], x[4]), c(prod[x[1],x[2]][[1]][10,1], prod[x[1]+x[5],x[2]][[1]][10,1]), rule = 2)
  return(matrix(data = c(class_1(t), class_2(t), class_3(t), class_4(t), class_5(t), class_6(t), class_7(t), class_8(t), class_9(t), class_10(t)), nrow = 10, ncol = 1))}
return(test(alder))}

# Plot strategy
stocks = data.frame(age = seq(25, 80, 1))

stocks$Low = sapply(stocks$age, function(x) sum(pf(x, "hoj")[5:10]))
stocks$Moderate = sapply(stocks$age, function(x) sum(pf(x, "mellem")[5:10]))
stocks$High = sapply(stocks$age, function(x) sum(pf(x, "lav")[5:10]))
scaleFUN <- function(x)
  sprintf("%.2f", x)
data_long <- melt(stocks, id = "age")
gfg_plot <- ggplot(data_long, aes(x = age, y = round(value, digits = 4), group = variable)) + geom_line(aes(linetype = variable)) + ylab(TeX("$w_{5-10}$")) + xlab("Age") + scale_linetype_manual(values = linetypes, name = "Risk aversion profile") + theme_bw() + scale_y_continuous(breaks = scales::pretty_breaks(n = 5), labels = scaleFUN ) + scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
gfg_plot
ggsave(filename = file.path(figure_dir, "P+ Livscyklus Growth-Asset Weight by Age and Risk Profile.png"), plot = gfg_plot, width = 8, height = 5, dpi = 300)

# CALCULATING GAMMA

# Data of people under 30
# test_unge = subset(data, data$gruppeNavn == "P211R" & data$policeStatus == "Betalende" & data$depot > 0 & data$alder < 30)

# Data on people between 50 and 70 
# mellemaldersgruppe <- subset(data, data$gruppeNavn == "P211R" & data$depot > 0 & data$alder >= 50 & data$alder < 70)

# Alternative restriction:
# mellemaldersgruppe <- subset(data, data$gruppeNavn == "P211R" & data$depot > 0 & data$alder >= 50 & data$alder < 70 & data$depot > 5 * (data$maanedligPraemie * 12))

# Check older members
# test_pensioneret = subset(data, data$gruppeNavn == "P211R" & data$depot > 0 & data$alder >= 70)

# Calibration sample retained from the supplied implementation. Unlike test_data, it does not restrict policy status, it also uses age >= 50 rather than a separate upper-age cutoff.
thepensiondata <- subset(data, data$gruppeNavn == "P211R" & data$depot > 0 & data$alder >= 50)

# Mortality intensity retained from supplied calibration
mu = function(x)
  0.000500 + 10^(5.88 + 0.038 * x - 10)

mu_bar = function(t1, t2)
  integrate(mu, t1, t2)$value

loop2 <- function(y, r = 0.015, t1 = 65) {for (i in y$policeNr) {print(i)
  alder = y$alder[y$policeNr == i]
  pension = t1 - alder
  maanedligPraemie = y$maanedligPraemie[y$policeNr == i]
  praemie = 12 * maanedligPraemie
  depot = y$depot[y$policeNr == i]
  risiko = y$investeringsProfil[y$policeNr == i]
    
# Contribution reserve
components <- function(t, state, parameters) {with(as.list(c(state, parameters)), {
  dXc <- (r + mu(alder + t)) * Xc - praemie
  list(c(dXc))})}
stepsize <- 1 / 12
parameters <- c(praemie = praemie, r = r)
state <- c(Xc = 0)

if (pension >= 2) {times <- seq(pension, 1, by = -stepsize)
ReserveSolution <- ode(y = state, times = times, func = components, parms = parameters, method = "rk4")
annuity <- approxfun(ReserveSolution[,"time"] - stepsize, ReserveSolution[,"Xc"], rule = 2)} else {
annuity = function(t)
  praemie}

account_c = function(t) {if (t < pension - 1 / 12) {annuity(t)} else {return(0)}}

# Inverse-calibration function
gamma_fct = function(alder, c, S0, c0, t, alpha, sigma) {
  if (alder < 50) {return(NA_real_)}
  if (alder >= 50 & alder < 70) {gamma = ((alpha - r) / sigma^2) * ((S0 + c0) / S0)
  return(gamma)}
  if (alder >= 70) {gamma = (alpha - r) / sigma^2
  return(gamma)}}

account = matrix(ncol = 1, nrow = 65)
account[1,1] <- depot
    
# One-year-ahead timing retained from the supplied calibration: portfolio and contribution reserve are evaluated at t = 1, while depot is the currently recorded balance.
# Mutual fund
x = as.matrix(pf(alder + 1, risiko))
alpha = as.numeric(t(x) %*% returns(1))
sigma = as.numeric(sqrt(t(x) %*% variance_matrix(1) %*% x))
    
if (alder >= 50) {y$gamma[y$policeNr == i] = gamma_fct(alder + 1, account_c(1), depot, account_c(1), 1, alpha, sigma)}

gamma = gamma_fct(alder + 1, account_c(1), depot, account_c(1), 1, alpha, sigma)
  
year = 1
    
while (year < 50 - alder) {year = year + 1
print(year)
account[year,1] <- (account[year - 1,1] + account_c(year - 1)) * exp((r + 1 / gamma * (alpha - r)^2 / sigma^2)) - account_c(year)
x = as.matrix(pf(alder + year, risiko))
alpha = as.numeric(t(x) %*% returns(1))
sigma = as.numeric(sqrt(t(x) %*% variance_matrix(1) %*% x))
gamma = gamma_fct(alder + 1, account_c(year), account[year - 1,1], account_c(year - 1), year, alpha, sigma)
y$gamma[y$policeNr == i] = gamma_fct(alder + 1, account_c(year), account[year - 1,1], account_c(year - 1), year, alpha, sigma)}}
  return(y)}

thepensiondatainclgamma = loop2(thepensiondata, r = 0.02, t1 = 65)

# Annual pension premium and wealth-to-premium ratio
thepensiondatainclgamma$annual_premium = 12 * thepensiondatainclgamma$maanedligPraemie
thepensiondatainclgamma$depot_to_annual_premium = thepensiondatainclgamma$depot / thepensiondatainclgamma$annual_premium

# Some observations have zero recorded premium.
thepensiondatainclgamma$depot_to_annual_premium[thepensiondatainclgamma$annual_premium <= 0] = NA
zero_premium_count = sum(thepensiondatainclgamma$annual_premium <= 0, na.rm = TRUE)
zero_premium_count

# Unrestricted gamma distribution
res_unrestricted = as.data.frame(thepensiondatainclgamma)
res_unrestricted$investeringsProfil <- factor(res_unrestricted$investeringsProfil, levels = c("mellem", "hoj", "lav"))
(d_gamma_unrestricted = ggplot(as.data.frame(res_unrestricted$gamma), aes(x = res_unrestricted$gamma, fill = res_unrestricted$investeringsProfil)) + geom_histogram(bins = 30, position = "identity", alpha = 0.6, col = "white") + ylab("Frequency") + xlab(TeX("$\\gamma$")) + scale_fill_manual(name = "Risk aversion profile", breaks = c("hoj", "mellem", "lav"), labels = list("Low", "Moderate", "High"), values = c("#858484", "gray", "black")) + theme_bw() + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)))
ggsave(filename = file.path(figure_dir, "Calibrated Risk Aversion - No Depot Restriction.png"), plot = d_gamma_unrestricted, width = 7, height = 5, dpi = 300)

# Distribution of pension wealth relative to annual premiums
summary(thepensiondatainclgamma$depot_to_annual_premium)
quantile(thepensiondatainclgamma$depot_to_annual_premium, probs = c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1), na.rm = TRUE)

# Alternative depot restrictions
gamma_no_restriction = thepensiondatainclgamma
gamma_1yr = subset(thepensiondatainclgamma, depot > 1 * annual_premium)
gamma_2yr = subset(thepensiondatainclgamma, depot > 2 * annual_premium)
gamma_5yr = subset(thepensiondatainclgamma, depot > 5 * annual_premium)
N_unrestricted = nrow(gamma_no_restriction)
restriction_summary = data.frame(Restriction = c("No restriction", "Depot > 1 annual premium", "Depot > 2 annual premiums", "Depot > 5 annual premiums"), N = c(nrow(gamma_no_restriction), nrow(gamma_1yr), nrow(gamma_2yr), nrow(gamma_5yr)))
restriction_summary$Removed = N_unrestricted - restriction_summary$N
restriction_summary$Removed_Percent = 100 * restriction_summary$Removed / N_unrestricted
restriction_summary

# Gamma statistics under each restriction
gamma_summary_function = function(dataset, restriction_name) {data.frame(Restriction = restriction_name, N = nrow(dataset), Mean_Gamma = mean(dataset$gamma, na.rm = TRUE), Median_Gamma = median(dataset$gamma, na.rm = TRUE), SD_Gamma = sd(dataset$gamma, na.rm = TRUE), Min_Gamma = min(dataset$gamma, na.rm = TRUE), Q25_Gamma = quantile(dataset$gamma, 0.25, na.rm = TRUE), Q75_Gamma = quantile(dataset$gamma, 0.75, na.rm = TRUE), Max_Gamma = max(dataset$gamma, na.rm = TRUE))}
gamma_restriction_summary = rbind(gamma_summary_function(gamma_no_restriction, "No restriction"), gamma_summary_function(gamma_1yr, "Depot > 1 annual premium"), gamma_summary_function(gamma_2yr, "Depot > 2 annual premiums"), gamma_summary_function(gamma_5yr, "Depot > 5 annual premiums"))
gamma_restriction_summary

# Compare unrestricted and final restricted samples
gamma_no_restriction_plot = gamma_no_restriction
gamma_5yr_plot = gamma_5yr
gamma_no_restriction_plot$Sample = "No depot restriction"
gamma_5yr_plot$Sample = "Depot > 5 annual premiums"
gamma_compare = rbind(gamma_no_restriction_plot, gamma_5yr_plot)
gamma_compare$Sample = factor(gamma_compare$Sample, levels = c("No depot restriction", "Depot > 5 annual premiums"))
gamma_compare$investeringsProfil = factor(gamma_compare$investeringsProfil, levels = c("mellem", "hoj", "lav"))
(d_gamma_compare = ggplot(gamma_compare, aes(x = gamma, fill = investeringsProfil)) + geom_histogram(bins = 30, position = "identity", alpha = 0.6, col = "white") + facet_wrap( ~ Sample, scales = "free", ncol = 1) + ylab("Frequency") + xlab(TeX("$\\gamma$")) + scale_fill_manual(name = "Risk aversion profile", breaks = c("hoj", "mellem", "lav"), labels = list("Low", "Moderate", "High"), values = c("#858484", "gray", "black")) + theme_bw())
ggsave(filename = file.path(figure_dir, "Calibrated Risk Aversion - No Restriction vs 5 Annual Premiums.png"), plot = d_gamma_compare, width = 7, height = 7, dpi = 300)

# Gamma statistics by risk profile: final 5-year restriction
gamma_profile_summary = gamma_5yr %>% group_by(investeringsProfil) %>% summarise(N = n(), Mean_Gamma = mean(gamma, na.rm = TRUE), Median_Gamma = median(gamma, na.rm = TRUE), SD_Gamma = sd(gamma, na.rm = TRUE), Q25_Gamma = quantile(gamma, 0.25, na.rm = TRUE), Q75_Gamma = quantile(gamma, 0.75, na.rm = TRUE), Min_Gamma = min(gamma, na.rm = TRUE), Max_Gamma = max(gamma, na.rm = TRUE))
gamma_profile_summary
gamma_profile_medians = gamma_5yr %>% group_by(investeringsProfil) %>% summarise(N = n(), Median_Gamma = median(gamma, na.rm = TRUE))
gamma_profile_medians

# Final restriction and plots
thepensiondatainclgammayeardeposit = subset(thepensiondatainclgamma, thepensiondatainclgamma$depot > 5 * (thepensiondatainclgamma$maanedligPraemie * 12))
res = as.data.frame(rbind(thepensiondatainclgammayeardeposit))
res$investeringsProfil <- factor(res$investeringsProfil, levels = c("mellem", "hoj", "lav"))
(d2 = ggplot(as.data.frame(res$gamma), aes(x = res$gamma, fill = res$investeringsProfil)) + geom_histogram(bins = 30, position = "identity", alpha = 0.6, col = "white") + ylab("Frequency") + xlab(TeX("$\\gamma$")) + scale_fill_manual(name = "Risk aversion profile", breaks = c("hoj", "mellem", "lav"), labels = list("Low", "Moderate", "High"), values = c("#858484", "gray", "black")) + theme_bw() + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)))

(d3 = ggplot(as.data.frame(res$gamma), aes(x = res$gamma, fill = res$investeringsProfil)) + geom_histogram(bins = 30, position = "identity", col = "white") + ylab("Frequency") + xlab(TeX("$\\gamma$")) + scale_fill_manual(name = "Risk aversion profile", breaks = c("hoj", "mellem", "lav"), labels = list("Low", "Moderate", "High"), values = c("#858484", "gray", "black")) + theme_bw() + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)))

# Save final calibration figure and numerical outputs
ggsave(filename = file.path(figure_dir, "Calibrated Risk Aversion by Risk Profile - 5 Annual Premium Restriction.png"), plot = d2, width = 7, height = 5, dpi = 300)
write.csv(restriction_summary, file.path(results_dir, "calibration_restriction_summary.csv"), row.names = FALSE)
write.csv(gamma_restriction_summary, file.path(results_dir, "calibration_gamma_summary_by_restriction.csv"), row.names = FALSE)
write.csv(as.data.frame(gamma_profile_summary), file.path(results_dir, "calibration_gamma_summary_by_profile.csv"), row.names = FALSE)
write.csv(gamma_5yr, file.path(results_dir, "calibration_final_sample.csv"), row.names = FALSE)

if (nrow(gamma_no_restriction) != 569) {warning("The unrestricted calibration sample is not 569. Check the input data and filters.")}
if (nrow(gamma_5yr) != 65) {warning("The final five-premium sample is not 65. Review the regenerated thesis values.")}

capture.output({
  cat("CALIBRATION RESULTS\n")
  cat("===================\n\n")
  cat("Zero-premium observations:", zero_premium_count, "\n\n")
  cat("Sample restrictions:\n")
  print(restriction_summary, row.names = FALSE)
  cat("\nGamma statistics by restriction:\n")
  print(gamma_restriction_summary, row.names = FALSE)
  cat("\nFinal gamma statistics by profile:\n")
  print(as.data.frame(gamma_profile_summary), row.names = FALSE)}, file = file.path(results_dir, "calibration_key_results.txt"))
cat("\nCalibration complete. Updated figures are in '", figure_dir, "' and numerical outputs are in '", results_dir, "'.\n", sep = "")