library(tidyverse)
library(brms)
library(cubelyr)
library(bayestestR)
library(patchwork)

d <- read.csv("data/data.csv")

# Set priors for all models
prior <- c(
  prior(normal(1, 3), class = Intercept),
  prior(normal(0, 1), class = b),
  prior(exponential(1), class = sd)
)

# Prior predictive check --------------------------------------------------
m_prior <- brm(value ~ M_Value_bdnf_s + bio_sex + (1|item) + (1|studyid),
               data = d,
               family = cumulative(),
               prior = prior,
               sample_prior = "only",
               chains = 1)

prior_pred <- ppc_sumscore(m_prior, d, probe_name = "Prior")
ggsave("figures/prior_pred.pdf", plot = prior_pred, width = 7, height = 5)

# Model fitting -----------------------------------------------------------
m_bdnf <- brm(value ~ M_Value_bdnf_s + (1 + M_Value_bdnf_s|item) + (1|studyid),
               data = d,
               family = cumulative(),
               prior = prior,
               chains = 4,
              cores = 4,
               init = "0")

m_slc6a4 <- brm(value ~ M_Value_slc6a4_s + bio_sex + (1 + M_Value_slc6a4_s|item) + (1|studyid),
               data = d,
               family = cumulative(),
               prior = prior,
               chains = 4,
               cores = 4,
               init = "0")

m_prdm8 <- brm(value ~ M_Value_prdm8_s + bio_sex + (1 + M_Value_prdm8_s + bio_sex|item) + (1|studyid),
              data = d,
              family = cumulative(),
              prior = prior,
              chains = 4,
              cores = 4,
              init = "0")

# Results/hypothesis tests ------------------------------------------------
bdnf_results <- pred_sumscore_diff(model = m_bdnf, data = d, probe_name = "bdnf")

slc6a4_results <- pred_sumscore_diff(model = m_slc6a4, data = d, probe_name = "slc6a4")

prdm8_results <- pred_sumscore_diff(model = m_prdm8, data = d, probe_name = "prdm8")

# Export contrast tables
write_csv(bind_rows(bdnf_results$preds_summary_mean, slc6a4_results$preds_summary_mean, prdm8_results$preds_summary_mean), file = "tables/model_contrasts.csv")

# save marginal effect plots
ggsave("figures/bdnf_marginal_effect.pdf", plot = bdnf_results$p_pred, width = 5, height = 5)
ggsave("figures/slc6a4_marginal_effect.pdf", plot = slc6a4_results$p_pred, width = 5, height = 5)
ggsave("figures/prdm8_marginal_effect.pdf", plot = prdm8_results$p_pred, width = 5, height = 5)

# Posterior predictive checks ---------------------------------------------
pp_bdnf <- ppc_sumscore(model = m_bdnf, d, probe_name = "bdnf", DNAm = T)
ggsave("figures/post_pred_bdnf.pdf", plot = pp_bdnf, width = 7, height = 8)

pp_slc6a4 <- ppc_sumscore(model = m_slc6a4, d, probe_name = "slc6a4", DNAm = T)
ggsave("figures/post_pred_slc6a4.pdf", plot = pp_slc6a4, width = 7, height = 8)

pp_prdm8 <- ppc_sumscore(model = m_prdm8, d, probe_name = "prdm8", DNAm = T)
ggsave("figures/post_pred_prdm8.pdf", plot = pp_slc6a4, width = 7, height = 8)
