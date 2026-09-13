# RQ4: QUANTITATIVE COMPARISON WITH THE THEORETICAL OPTIMUM

# Uses gamma_5yr, the final restricted calibration sample. Member-level calculations use recorded age, depot, premium, calibrated gamma and risk profile.
required_objects <- c("gamma_5yr", "returns", "variance_matrix", "pf")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects) > 0) {stop(paste0("Run source('code/calibration_analysis.R') first. Missing objects: ", paste(missing_objects, collapse = ", ")))}

figure_dir <- "figures"
results_dir <- "results"

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Market objects reused from calibration
r_riskfree <- 0.02
alpha_vec <- returns(1)
Sigma_mat <- variance_matrix(1)
Sigma_inv <- solve(Sigma_mat)
excess_ret <- as.numeric(alpha_vec - r_riskfree)
L <- t(chol(Sigma_mat))
Theta2 <- as.numeric(t(excess_ret) %*% Sigma_inv %*% excess_ret)
retirement_age <- 65

# Deterministic present value of remaining contributions from Chapter 3
h_of_t <- function(age, ell, T_ret = retirement_age, r = r_riskfree) {tau <- pmax(T_ret - age, 0)
ifelse(tau <= 0, 0, (ell / r) * (1 - exp(-r * tau)))}

# Certainty-equivalent balance evaluated in log space for numerical stability
ceb <- function(gamma, X) {logX <- log(X)
a <- (1 - gamma) * logX
m <- max(a)
exp((m + log(mean(exp(a - m)))) / (1 - gamma))}
growth_weight <- function(theta)
  sum(theta[5:10])

# Pre-retirement observations
rq4_sample <- subset(gamma_5yr, alder < retirement_age & !is.na(gamma))
rq4_sample$ell_i <- 12 * rq4_sample$maanedligPraemie
rq4_sample$h_i <- h_of_t(rq4_sample$alder, rq4_sample$ell_i)
cat("RQ4 sample size (age 50-64, final restricted sample):", nrow(rq4_sample), "\n")

# Cross-sectional glide-path comparison
rq4_sample$growth_observed <- NA_real_
rq4_sample$growth_theoretical <- NA_real_
rq4_sample$riskfree_theoretical <- NA_real_
for (k in seq_len(nrow(rq4_sample))) {
  age_k <- rq4_sample$alder[k]
  x_k <- rq4_sample$depot[k]
  h_k <- rq4_sample$h_i[k]
  gamma_k <- rq4_sample$gamma[k]
  risiko_k <- rq4_sample$investeringsProfil[k]
  w_obs <- as.numeric(pf(age_k, risiko_k))
  w_theo <- as.numeric((1 / gamma_k) * (1 + h_k / x_k) * (Sigma_inv %*% excess_ret))
  rq4_sample$growth_observed[k] <- growth_weight(w_obs)
  rq4_sample$growth_theoretical[k] <- growth_weight(w_theo)
  rq4_sample$riskfree_theoretical[k] <- 1 - sum(w_theo)}
rq4_sample$age_band <- cut(rq4_sample$alder, breaks = c(50, 55, 60, 65), include.lowest = TRUE, right = FALSE)
rq4_glidepath_summary <- rq4_sample %>% group_by(investeringsProfil, age_band) %>% summarise(N = n(), mean_growth_observed = mean(growth_observed), mean_growth_theoretical = mean(growth_theoretical), mean_riskfree_theoretical = mean(riskfree_theoretical),.groups = "drop")
print(rq4_glidepath_summary)
write.csv(rq4_sample, file.path(results_dir, "rq4_crosssectional_glidepath.csv"), row.names = FALSE)
write.csv(rq4_glidepath_summary, file.path(results_dir, "rq4_glidepath_summary_by_profile_age.csv"), row.names = FALSE)
rq4_plot_long <- rbind(data.frame(age = rq4_sample$alder, profile = rq4_sample$investeringsProfil, weight = rq4_sample$growth_observed, series = "Observed P+ strategy"), data.frame(age = rq4_sample$alder, profile = rq4_sample$investeringsProfil, weight = rq4_sample$growth_theoretical, series = "Theoretical optimal strategy"))
rq4_plot_long$profile_lab <- recode(rq4_plot_long$profile, hoj = "Low risk aversion", mellem = "Moderate risk aversion", lav = "High risk aversion")
(rq4_glidepath_plot <- ggplot(rq4_plot_long, aes(x = age, y = weight, color = series)) + geom_point(alpha = 0.5, size = 1.6) + geom_smooth(se = FALSE, method = "loess", span = 1.2) + facet_wrap(~ profile_lab) + ylab(TeX("Growth-asset weight $w_{5-10}$")) + xlab("Age") + scale_color_manual(name = NULL, values = c("Observed P+ strategy" = "black", "Theoretical optimal strategy" = "#858484")) + scale_y_continuous(breaks = seq(0.2, 0.7, by = 0.1), labels = function(x)sprintf("%.2f", x)) + theme_bw())
ggsave(file.path(figure_dir, "RQ4 Theoretical vs Observed Growth Weight (Real Members, Age 50+).png"), rq4_glidepath_plot, width = 9, height = 5, dpi = 300)

# PRE-RETIREMENT CERTAINTY-EQUIVALENT WELFARE COMPARISON
set.seed(2026)
Nsim <- 3000
dt <- 1 / 12
welfare_results <- data.frame(policeNr = rq4_sample$policeNr, alder = rq4_sample$alder, investeringsProfil = rq4_sample$investeringsProfil, gamma = rq4_sample$gamma, depot = rq4_sample$depot, ell = rq4_sample$ell_i, CEB_observed = NA_real_, CEB_theoretical = NA_real_, extra_balance_required = NA_real_, extra_balance_pct = NA_real_)
for (k in seq_len(nrow(rq4_sample))) {
  age_k <- rq4_sample$alder[k]
  x0_k <- rq4_sample$depot[k]
  ell_k <- rq4_sample$ell_i[k]
  gamma_k <- rq4_sample$gamma[k]
  risiko_k <- rq4_sample$investeringsProfil[k]
  n_steps <- max(round((retirement_age - age_k) / dt), 1)
  ages <- age_k + (0:(n_steps - 1)) * dt
  w_dir <- as.numeric((1 / gamma_k) * (Sigma_inv %*% excess_ret))
  dir_vol_vec <- as.numeric(t(w_dir) %*% L)
  M <- rep(1, Nsim)
  C <- rep(0, Nsim)
  cumshock <- rep(0, Nsim)
  for (s in seq_len(n_steps)) {
    age_now <- ages[s]
    w_t <- as.numeric(pf(age_now, risiko_k))
    drift_t <- r_riskfree + sum(w_t * excess_ret)
    vol_t <- as.numeric(t(w_t) %*% L)
    var_t <- sum(vol_t^2)
    Z <- matrix(rnorm(Nsim * 10), nrow = Nsim, ncol = 10)
    shock_t <- as.numeric(Z %*% vol_t) * sqrt(dt)
    logret_t <- (drift_t - 0.5 * var_t) * dt + shock_t
    prem_t <-
      if (age_now < retirement_age) {ell_k * dt}
    else {0}
    M <- M * exp(logret_t)
    C <- C * exp(logret_t) + prem_t
    bshock_t <- as.numeric(Z %*% dir_vol_vec) * sqrt(dt)
    cumshock <- cumshock + bshock_t}
  Tt <- n_steps * dt
  X_obs65 <- function(x)
    x * M + C
  mu_Y <- r_riskfree + Theta2 / gamma_k
  bnorm2 <- Theta2 / gamma_k^2
  h0 <- h_of_t(age_k, ell_k)
  X_theo65 <- (x0_k + h0) * exp((mu_Y - 0.5 * bnorm2) * Tt + cumshock)
  CEB_obs <- ceb(gamma_k, X_obs65(x0_k))
  CEB_theo <- ceb(gamma_k, X_theo65)
  f <- function(delta)
    ceb(gamma_k, X_obs65(x0_k + delta)) - CEB_theo
  lo <- -0.99 * x0_k
  hi <- 50 * max(x0_k, 1)
  root <- tryCatch(uniroot(f, lower = lo, upper = hi, extendInt = "yes", tol = 1), error = function(e)
    NULL)
  delta_star <- if (is.null(root)) {NA_real_} else {root$root}
  
  welfare_results$CEB_observed[k] <- CEB_obs
  welfare_results$CEB_theoretical[k] <- CEB_theo
  welfare_results$extra_balance_required[k] <- delta_star
  welfare_results$extra_balance_pct[k] <- 100 * delta_star / x0_k
  
  cat(sprintf(paste0("[%d/%d] policeNr=%s ", "age=%.1f gamma=%.2f ", "CEB_obs=%.0f ", "CEB_theo=%.0f ", "extra=%.1f%%\n"), k, nrow(rq4_sample), welfare_results$policeNr[k], age_k, gamma_k, CEB_obs, CEB_theo, welfare_results$extra_balance_pct[k]))}

write.csv(welfare_results, file.path(results_dir, "rq4_welfare_comparison_per_person.csv"), row.names = FALSE)

welfare_summary_by_profile <- welfare_results %>% group_by(investeringsProfil) %>% summarise(N = n(), mean_gamma = mean(gamma), mean_CEB_ratio = mean(CEB_theoretical / CEB_observed, na.rm = TRUE), median_extra_balance_pct = median(extra_balance_pct, na.rm = TRUE), mean_extra_balance_pct = mean(extra_balance_pct, na.rm = TRUE),.groups = "drop")
print(welfare_summary_by_profile)

write.csv(welfare_summary_by_profile, file.path(results_dir, "rq4_welfare_summary_by_profile.csv"), row.names = FALSE)

(rq4_welfare_plot <- ggplot(welfare_results, aes(x = investeringsProfil, y = extra_balance_pct, fill = investeringsProfil)) + geom_boxplot(alpha = 0.6) + ylab("Extra starting balance required (%)") + xlab("Risk profile") + scale_x_discrete(labels = c(hoj = "Low risk aversion", mellem = "Moderate risk aversion", lav = "High risk aversion")) + scale_fill_manual(name = "Risk aversion profile", breaks = c("hoj", "mellem", "lav"), labels = list("Low", "Moderate", "High"), values = c("#858484","gray","black")) + theme_bw() + theme(legend.position = "none"))

ggsave(file.path(figure_dir, "RQ4 Extra Balance Required by Profile.png"), rq4_welfare_plot, width = 7, height = 5, dpi = 300)

cat("\nDone.\n")

# POST-RETIREMENT WELFARE COMPARISON
# Uses Chapter 4's post-retirement optimal investment and consumption solution as the theoretical benchmark. Both the theoretical and observed-investment comparisons use the same theoretical consumption feedback rule: c*(t,x) = x / abar(t). Only the investment policy differs.
N_terminal <- 100
rho_discount <- 0
nu_fn <- function(gamma, rho = rho_discount) {(rho - (1 - gamma) * (r_riskfree + Theta2 / (2 * gamma))) / gamma}
abar_fn <- function(t, gamma, N = N_terminal) {
  nu <- nu_fn(gamma)
  if (abs(nu) < 1e-8) {return(N - t)}
  (1 - exp(-nu * (N - t))) / nu}
V_fn <- function(t, x, gamma, N = N_terminal, rho = rho_discount) {exp(-rho * t) * abar_fn(t, gamma, N)^gamma * x^(1 - gamma) / (1 - gamma)}

# Expected discounted utility when the observed P+ investment weights are continued forward and paired with the same theoretical optimal consumption rule.
J_hybrid <- function(t, gamma, profile, N = N_terminal, rho = rho_discount, dt = 1 / 12) {
  n_steps <- max(round((N - t) / dt), 1)
  u_grid <- t + (0:n_steps) * dt
  u_grid[length(u_grid)] <- N
  mu_obs <- numeric(length(u_grid) - 1)
  var_obs <- numeric(length(u_grid) - 1)
  nu <- nu_fn(gamma)
  for (k in seq_len(length(u_grid) - 1)) {u_mid <- 0.5 * (u_grid[k] + u_grid[k + 1])
  w_u <- as.numeric(pf(u_mid, profile))
  mu_obs[k] <- r_riskfree + sum(w_u * excess_ret)
  var_obs[k] <- as.numeric(t(w_u) %*% Sigma_mat %*% w_u)}
  step_dt <- diff(u_grid)
  incr_mean <- (mu_obs - nu - 0.5 * var_obs) * step_dt
  incr_var <- var_obs * step_dt
  cum_mean <- c(0, cumsum(incr_mean))
  cum_var <- c(0, cumsum(incr_var))
  abar_t <- abar_fn(t, gamma, N)
  logY0 <- -log(abar_t)
  mean_lnY <- logY0 + cum_mean
  var_lnY <- cum_var
  integrand <- exp((1 - gamma) * mean_lnY + 0.5 * (1 - gamma)^2 * var_lnY) / (1 - gamma)
  integrand <- exp(-rho * u_grid) * integrand
  sum(0.5 * (integrand[-1] + integrand[-length(integrand)]) * step_dt)}
post_retirement_sample <- subset(gamma_5yr, alder >= retirement_age & !is.na(gamma))
cat("Post-retirement (age 65+) sample size:", nrow(post_retirement_sample), "\n")
post_results <- data.frame(policeNr = post_retirement_sample$policeNr, alder = post_retirement_sample$alder, investeringsProfil = post_retirement_sample$investeringsProfil, gamma = post_retirement_sample$gamma, depot = post_retirement_sample$depot, V_theoretical = NA_real_, J_hybrid_observed = NA_real_, extra_balance_pct = NA_real_)
for (k in seq_len(nrow(post_retirement_sample))) {
  t_k <- post_retirement_sample$alder[k]
  gamma_k <- post_retirement_sample$gamma[k]
  risiko_k <- post_retirement_sample$investeringsProfil[k]
  V_theo <- V_fn(t_k, 1, gamma_k)
  J_hyb <- J_hybrid(t_k, gamma_k, risiko_k)
  extra_pct <- 100 * ((V_theo / J_hyb)^(1 / (1 - gamma_k)) - 1)
  post_results$V_theoretical[k] <- V_theo
  post_results$J_hybrid_observed[k] <- J_hyb
  post_results$extra_balance_pct[k] <- extra_pct
  cat(sprintf(paste0("[%d/%d] policeNr=%s ", "age=%.1f gamma=%.2f ", "extra_balance=%.1f%%\n"), k, nrow(post_retirement_sample), post_results$policeNr[k], t_k, gamma_k, extra_pct))}
print(post_results)
write.csv(post_results, file.path(results_dir, "rq4_post_retirement_welfare.csv"), row.names = FALSE)
cat("\nMean extra balance required (post-retirement, age 65+):", round(mean(post_results$extra_balance_pct, na.rm = TRUE), 1), "%\n")
cat("Median extra balance required (post-retirement, age 65+):", round(median(post_results$extra_balance_pct, na.rm = TRUE), 1), "%\n")
as.data.frame(welfare_summary_by_profile)

# POST-RETIREMENT OPTIMAL CONSUMPTION VISUALISATION
# This section illustrates the Chapter 4 optimal consumption rule: c*(t,x) = x / abar(t) for one representative retired member from the final restricted sample.
# The representative member is the age-65+ individual whose calibrated gamma is closest to the median gamma of the retired subsample.
# The main thesis figure shows: 1. one simulated realisation of c*(t,X*(t)) and 2. expected optimal consumption.
# A second figure shows the deterministic optimal annual consumption-to-wealth ratio 100/abar(t).
# Consumption is defined on [t0,N), so t=N is excluded from the consumption plot.

# Check post-retirement sample
if (!exists("post_retirement_sample")) {post_retirement_sample <- subset(gamma_5yr, alder >= retirement_age & !is.na(gamma))}
cat("\nNumber of age-65+ observations available for consumption illustration:", nrow(post_retirement_sample), "\n")
if (nrow(post_retirement_sample) == 0) {stop(paste("No age-65+ observations are available", "in the final restricted calibration sample."))}

# Choose representative retired member
median_gamma_post <- median(post_retirement_sample$gamma, na.rm = TRUE)
representative_index <- which.min(abs(post_retirement_sample$gamma - median_gamma_post))
representative_person <- post_retirement_sample[representative_index,, drop = FALSE]
age0_cons <- as.numeric(representative_person$alder[1])
x0_cons <- as.numeric(representative_person$depot[1])
gamma_cons <- as.numeric(representative_person$gamma[1])
profile_cons <- as.character(representative_person$investeringsProfil[1])
police_cons <- representative_person$policeNr[1]
cat("\n--------------------------------------------------\n")
cat("REPRESENTATIVE MEMBER FOR CONSUMPTION FIGURE\n")
cat("--------------------------------------------------\n")
cat("Police number:", police_cons,"\n")
cat("Current age:", round(age0_cons, 2),"\n")
cat("Initial pension wealth:", format(round(x0_cons, 0), big.mark = ",", scientific = FALSE),"\n")
cat("Calibrated gamma:", round(gamma_cons, 4),"\n")
cat("Risk profile:", profile_cons,"\n")
cat("Median retired-sample gamma:", round(median_gamma_post, 4),"\n")
cat("--------------------------------------------------\n\n")

# Post-retirement parameters

if (!exists("N_terminal")) {N_terminal <- 100}
if (!exists("rho_discount")) {rho_discount <- 0}

# Effective rate from Chapter 4
nu_cons <- (rho_discount - (1 - gamma_cons) * (r_riskfree + Theta2 / (2 * gamma_cons))) / gamma_cons

# Annuity factor specialised to this representative member
abar_cons_fn <- function(age) {
  remaining_horizon <- N_terminal - age
  if (remaining_horizon <=0) {return(0)}
  if (abs(nu_cons) < 1e-8) {return(remaining_horizon)}
  (1 - exp(-nu_cons * remaining_horizon)) / nu_cons}
abar0_cons <- abar_cons_fn(age0_cons)
if (!is.finite(abar0_cons) || abar0_cons <= 0) {stop("The initial annuity factor is non-positive or undefined.")}
cat("Effective rate nu:", round(nu_cons, 6),"\n")
cat("Initial annuity factor abar(t0):", round(abar0_cons, 4),"\n")
cat("Initial optimal annual consumption:", format(round(x0_cons / abar0_cons, 0), big.mark = ",", scientific = FALSE),"\n\n")

# Chapter 4 Merton allocation
w_merton_cons <- as.numeric((1 / gamma_cons) * (Sigma_inv %*% excess_ret))
merton_excess_cons <- as.numeric(sum(w_merton_cons * excess_ret))
merton_variance_cons <- as.numeric(t(w_merton_cons) %*% Sigma_mat %*% w_merton_cons)
cat("Merton expected excess return:", round(merton_excess_cons, 6),"\n")
cat("Merton portfolio volatility:", round(sqrt(merton_variance_cons), 6),"\n\n")

# Monthly simulation grid
set.seed(2026)
dt_cons <- 1 / 12
age_grid_cons <- seq(from = age0_cons, to = N_terminal, by = dt_cons)
if (tail(age_grid_cons, 1) < N_terminal) {age_grid_cons <- c(age_grid_cons, N_terminal)}
age_grid_cons[length(age_grid_cons)] <- N_terminal
n_cons <- length(age_grid_cons)
tau_cons <- age_grid_cons - age0_cons

# Annuity factor over the horizon
abar_grid_cons <- sapply(age_grid_cons, abar_cons_fn)

# Closed-loop stochastic process M(t)
bnorm2_cons <- Theta2 / gamma_cons^2
bnorm_cons <- sqrt(bnorm2_cons)
M_drift_cons <- r_riskfree + Theta2 / gamma_cons - Theta2 / (2 * gamma_cons^2) - nu_cons
M_cons <- rep(1, n_cons)
if (n_cons > 1) {
  for (j in 2:n_cons) {
    dt_j <- age_grid_cons[j] - age_grid_cons[j - 1]
    Z_j <- rnorm(1)
    M_cons[j] <- M_cons[j - 1] * exp(M_drift_cons * dt_j + bnorm_cons * sqrt(dt_j) * Z_j)}}

# Optimal closed-loop wealth
X_opt_cons <- x0_cons * (abar_grid_cons / abar0_cons) * M_cons
X_opt_cons[length(X_opt_cons)] <- 0

# Realised optimal consumption
c_realised_cons <- (x0_cons / abar0_cons) * M_cons

# Expected optimal consumption
expected_M_cons <- exp((r_riskfree + Theta2 / gamma_cons - nu_cons) * tau_cons)
c_expected_cons <- (x0_cons / abar0_cons) * expected_M_cons

# Consumption-to-wealth ratio
consumption_ratio_cons <- rep(NA_real_, n_cons)
non_terminal_cons <- age_grid_cons < N_terminal
consumption_ratio_cons[non_terminal_cons] <- 1 / abar_grid_cons[non_terminal_cons]

# Store results
consumption_path_data <- data.frame(Age = age_grid_cons, Annuity_Factor = abar_grid_cons, Optimal_Wealth = X_opt_cons, Realised_Optimal_Consumption = c_realised_cons, Expected_Optimal_Consumption = c_expected_cons, Consumption_Wealth_Ratio = consumption_ratio_cons)
head(consumption_path_data)
tail(consumption_path_data)
write.csv(consumption_path_data, file.path(results_dir, "optimal_post_retirement_consumption_representative_member.csv"), row.names = FALSE)

# MAIN THESIS FIGURE: REALISED AND EXPECTED OPTIMAL CONSUMPTION
# Exclude t=N because consumption is defined on [t0,N).
consumption_plot_base <- subset(consumption_path_data, Age < N_terminal)
consumption_plot_data <- rbind(data.frame(Age = consumption_plot_base$Age, Consumption = consumption_plot_base$Realised_Optimal_Consumption, Series = "One simulated optimal path"), data.frame(Age = consumption_plot_base$Age, Consumption = consumption_plot_base$Expected_Optimal_Consumption, Series = "Expected optimal consumption"))
consumption_plot_data$Series <- factor(consumption_plot_data$Series, levels = c("One simulated optimal path", "Expected optimal consumption"))
optimal_consumption_plot <- ggplot(consumption_plot_data, aes(x = Age, y = Consumption, linetype = Series, group = Series)) + geom_line(linewidth = 0.85) + scale_linetype_manual(name = NULL, values = c("One simulated optimal path" = "solid", "Expected optimal consumption" = "dashed")) + scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 6)) + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + labs(x = "Age", y = "Optimal annual consumption (DKK)") + theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())
optimal_consumption_plot
ggsave(filename = file.path(figure_dir, "Optimal Post-Retirement Consumption - Representative Member.png"), plot = optimal_consumption_plot, width = 8, height = 5, dpi = 300)

# OPTIMAL CONSUMPTION-TO-WEALTH RATIO
# The ratio diverges as t -> N. The figure therefore ends one year before N so that the final divergence does not dominate the displayed scale.
ratio_plot_end_age <- N_terminal - 1
consumption_ratio_plot_data <- subset(consumption_path_data, Age <= ratio_plot_end_age & is.finite(Consumption_Wealth_Ratio))
consumption_ratio_plot_data$Consumption_Wealth_Percent <- 100 * consumption_ratio_plot_data$Consumption_Wealth_Ratio
optimal_consumption_ratio_plot <- ggplot(consumption_ratio_plot_data, aes(x = Age, y = Consumption_Wealth_Percent)) + geom_line(linewidth = 0.85) + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 6), labels = scales::label_number(accuracy = 1)) + labs(x = "Age", y = "Optimal annual consumption-to-wealth ratio (%)") + theme_bw()
optimal_consumption_ratio_plot
ggsave(filename = file.path(figure_dir, "Optimal Post-Retirement Consumption-to-Wealth Ratio.png"), plot = optimal_consumption_ratio_plot, width = 8, height = 5, dpi = 300)

# OPTIMAL WEALTH FIGURE
optimal_wealth_plot <- ggplot(consumption_path_data, aes(x = Age, y = Optimal_Wealth)) + geom_line(linewidth = 0.85) + scale_y_continuous(labels = scales::comma, breaks = scales::pretty_breaks(n = 6)) + scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + labs(x = "Age", y = "Optimal pension wealth (DKK)") + theme_bw()
optimal_wealth_plot
ggsave(filename = file.path(figure_dir, "Optimal Post-Retirement Wealth - Representative Member.png"), plot = optimal_wealth_plot, width = 8, height = 5, dpi = 300)

# SUMMARY INFORMATION FOR CHAPTER 5
consumption_summary <- data.frame(Police_Number = police_cons, Initial_Age = age0_cons, Initial_Wealth = x0_cons, Risk_Profile = profile_cons, Gamma = gamma_cons, Rho = rho_discount, Effective_Rate_Nu = nu_cons, Initial_Annuity_Factor = abar0_cons, Initial_Optimal_Consumption = x0_cons / abar0_cons, Initial_Consumption_Wealth_Ratio = 1 / abar0_cons, Terminal_Age = N_terminal)
print(consumption_summary)
write.csv(consumption_summary, file.path(results_dir, "optimal_consumption_representative_member_summary.csv"), row.names = FALSE)
# FINAL CHECKS
cat("\n==================================================\n")
cat("POST-RETIREMENT CONSUMPTION VISUALISATION COMPLETE\n")
cat("==================================================\n")
capture.output({
  cat("THEORY-PRACTICE AND WELFARE RESULTS\n")
  cat("====================================\n\n")
  cat("Cross-sectional glide-path summary:\n")
  print(as.data.frame(rq4_glidepath_summary), row.names = FALSE)
  cat("\nPre-retirement welfare summary:\n")
  print(as.data.frame(welfare_summary_by_profile), row.names = FALSE)
  cat("\nPost-retirement welfare per member:\n")
  print(post_results, row.names = FALSE)
  cat("\nPost-retirement extra-balance summary:\n")
  print(summary(post_results$extra_balance_pct))
  cat("\nRepresentative-member consumption summary:\n")
  print(consumption_summary)
  cat("\nSession information:\n")
  print(sessionInfo())}, file = file.path(results_dir, "theory_practice_key_results.txt"))

write.csv(data.frame(file = sort(list.files(results_dir)), stringsAsFactors = FALSE), file.path(results_dir, "results_manifest.csv"), row.names = FALSE)

cat("\nAll figures were written to '", figure_dir, "' and all result files were written to '", results_dir, "'.\n", sep = "")
cat("Representative age:", round(age0_cons, 2), "\n")
cat("Representative gamma:", round(gamma_cons, 4), "\n")
cat("Initial wealth:", format(round(x0_cons, 0), big.mark = ",", scientific = FALSE), "\n")
cat("Initial optimal annual consumption:", format(round(x0_cons / abar0_cons, 0), big.mark = ",", scientific = FALSE), "\n")
cat("Initial optimal consumption-to-wealth ratio:", round(100 / abar0_cons, 2), "% per year\n")
cat("Terminal optimal wealth:", round(tail(X_opt_cons, 1), 6),"\n")
cat("\nFiles created:\n")
cat("See the 'figures' directory for regenerated PNG files.\n")
cat("See the 'results' directory for CSV and key-results files.\n")
cat("==================================================\n")