# ==============================================================================
# Iron Overload and Male Hypogonadism: Random-Effects Meta-Analysis
# Transfusion-Dependent Thalassaemia (TDT) male hypogonadism prevalence
# PROSPERO: CRD420261344877
# Author: Borges JYV
# ==============================================================================
# Model: DerSimonian-Laird random effects
# Transformation: Freeman-Tukey double arcsine (PFT)
# Back-transformation: Miller (1978) with harmonic mean of sample sizes
# ==============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("meta", quietly = TRUE)) install.packages("meta")
  if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor")
  library(meta)
  library(metafor)
})

set.seed(20260318)

dat <- read.csv("data/tdt_male_hypogonadism_prevalence.csv", stringsAsFactors = FALSE)
dat$label <- paste(dat$study, dat$year)

# ---- Primary analysis: meta::metaprop (PFT, DL) ------------------------------
m <- metaprop(
  event      = events,
  n          = n_males,
  studlab    = label,
  data       = dat,
  sm         = "PFT",
  method.tau = "DL",
  method.ci  = "CP",
  random     = TRUE,
  common     = TRUE,
  title      = "TDT male hypogonadism prevalence"
)

print(summary(m))

# ---- Export pooled estimate ---------------------------------------------------
s <- summary(m)
out <- data.frame(
  model       = c("Random effects (DL)", "Fixed effects (IV)"),
  pooled_pct  = round(100 * c(s$random$TE, s$common$TE), 1),
  ci_low_pct  = round(100 * c(s$random$lower, s$common$lower), 1),
  ci_high_pct = round(100 * c(s$random$upper, s$common$upper), 1),
  I2_pct      = round(100 * m$I2, 1),
  Q           = round(m$Q, 2),
  Q_df        = m$df.Q,
  Q_p         = signif(m$pval.Q, 3),
  tau2        = round(m$tau2, 4)
)
write.csv(out, "analysis/pooled_estimate.csv", row.names = FALSE)
print(out)

# ---- Forest plot --------------------------------------------------------------
dir.create("figures", showWarnings = FALSE)
plot_forest <- function() {
  forest(m,
         xlab = "Hypogonadism prevalence (%)",
         pscale = 100,
         leftcols = c("studlab", "event", "n"),
         leftlabs = c("Study", "Events", "N"),
         rightcols = c("effect", "ci", "w.random"),
         rightlabs = c("Prevalence (%)", "95% CI", "Weight (RE)"),
         print.tau2 = TRUE, print.I2 = TRUE, print.Q = TRUE,
         col.diamond = "grey30", col.square = "grey50",
         digits = 1)
}
pdf("figures/Figure2_forest_tdt_male_hypogonadism.pdf", width = 10, height = 5); plot_forest(); dev.off()
png("figures/Figure2_forest_tdt_male_hypogonadism.png", width = 3000, height = 1500, res = 300); plot_forest(); dev.off()

# ---- Sensitivity: leave-one-out ----------------------------------------------
loo <- metainf(m, pooled = "random")
print(loo)
sink("analysis/leave_one_out.txt"); print(loo); sink()

# ---- Sensitivity: metafor cross-check ----------------------------------------
es <- escalc(measure = "PFT", xi = events, ni = n_males, data = dat)
rma_dl <- rma(yi, vi, data = es, method = "DL")
pred <- predict(rma_dl, transf = transf.ipft.hm, targs = list(ni = dat$n_males))
cat("\nmetafor cross-check (DL, PFT, harmonic-mean back-transform):\n")
print(pred)

sessionInfo()
