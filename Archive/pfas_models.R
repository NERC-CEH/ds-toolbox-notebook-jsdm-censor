library(tidyverse)
library(ggplot2)
library(patchwork)
library(jsdmstan)

# read in and explore PFAS data
load("PFAS_FW_2AZZ_methods 21_24.rda")
# remotes::install_github("NERC-CEH/jsdmstan")
# View(pfas)
# limit to the 47 PFAS chemicals of interest
pfas_21 <- filter(pfas, MEAS_ANAL_METH_CODE == "21")
pfas_counts <- count(pfas_21, MEAS_DETERMINAND_CODE,DETE_SHORT_DESC)
pfas_highcounts <- filter(pfas_counts, n > 1000)
pfas_filtered <- filter(pfas_21, MEAS_DETERMINAND_CODE %in% pfas_highcounts$MEAS_DETERMINAND_CODE)
janitor::tabyl(pfas_filtered$MEAS_SIGN)
pfas_nonUnder <- filter(pfas_filtered, is.na(MEAS_SIGN))

# subset and explore data via plots
pfas_sitevisits <- pfas_filtered %>%
  select(SAMP_SMPT_USER_REFERENCE, wims_region, SAMP_ID, DATE_TIME) %>%
  distinct() %>% count(wims_region, SAMP_SMPT_USER_REFERENCE)
pfas_highsite_visits <- filter(pfas_sitevisits, n >= 20)
pfas_subset <- filter(pfas_filtered, SAMP_SMPT_USER_REFERENCE %in%
                        pfas_highsite_visits$SAMP_SMPT_USER_REFERENCE)
ggplot(pfas_subset, aes(x = DATE_TIME, y = MEAS_RESULT, colour = MEAS_DETERMINAND_CODE)) +
  geom_line() + scale_y_log10() +
  facet_wrap(~wims_region + SAMP_SMPT_USER_REFERENCE) +
  theme(legend.position = "none")


pfas_sitevisits <- pfas_filtered %>%
  filter(DATE_TIME > ymd("2023-01-01")) %>%
  select(SAMP_SMPT_USER_REFERENCE, wims_region, SAMP_ID, DATE_TIME) %>%
  distinct() %>% count(wims_region, SAMP_SMPT_USER_REFERENCE)
pfas_highsite_visits <- filter(pfas_sitevisits, n >= 20)
pfas_subset <- filter(pfas_filtered, SAMP_SMPT_USER_REFERENCE %in%
                        pfas_highsite_visits$SAMP_SMPT_USER_REFERENCE &
                        DATE_TIME > dmy("01012024") & DATE_TIME < dmy("31122024")) %>%
  group_by(SAMP_SMPT_USER_REFERENCE) %>%
  filter(sum(is.na(MEAS_SIGN))>50)
ggplot(pfas_subset, aes(x = DATE_TIME, y = MEAS_RESULT, colour = MEAS_DETERMINAND_CODE)) +
  geom_line() + scale_y_log10() +
  facet_wrap(~wims_region + SAMP_SMPT_USER_REFERENCE) +
  theme(legend.position = "none")

# pfas_subset <- filter(pfas_filtered, DATE_TIME > dmy("01012023")) %>%
#   group_by(SAMP_SMPT_USER_REFERENCE) %>%
#   filter(sum(is.na(MEAS_SIGN))>50) %>% 
#   ungroup

# get into jsdmstan format
pfas_y <- pfas_subset %>% ungroup() %>%
  select(SAMP_ID, MEAS_RESULT, MEAS_DETERMINAND_CODE,SAMP_SMPT_USER_REFERENCE, DATE_TIME, SAMP_PURPOSE_CODE, wims_region) %>%
  pivot_wider(values_from = MEAS_RESULT, names_from = MEAS_DETERMINAND_CODE)

dim(na.omit(pfas_y))
pfas_y <- na.omit(pfas_y)

# create censoring matrix (1 if left-censored, 0 if not). Treating all
# right-censored data as uncensored
Y <- pfas_y[, 6:ncol(pfas_y)]
cens_ID <- pfas_subset %>% ungroup() %>%
  select(SAMP_ID, MEAS_SIGN, MEAS_DETERMINAND_CODE,SAMP_SMPT_USER_REFERENCE, DATE_TIME, SAMP_PURPOSE_CODE, wims_region) %>%
  mutate(MEAS_SIGN = replace_na(ifelse(MEAS_SIGN == "<", 1, 0),0)) %>%
  pivot_wider(values_from = MEAS_SIGN, names_from = MEAS_DETERMINAND_CODE) %>%
  filter(SAMP_ID %in% pfas_y$SAMP_ID)
all.equal(cens_ID$SAMP_ID, pfas_y$SAMP_ID)
all.equal(cens_ID$DATE_TIME, pfas_y$DATE_TIME)
all.equal(cens_ID$SAMP_SMPT_USER_REFERENCE, pfas_y$SAMP_SMPT_USER_REFERENCE)
all.equal(colnames(cens_ID),colnames(pfas_y))
cens_ID <- cens_ID[,6:ncol(cens_ID)]

# remove chemicals that are mostly or wholly uncensored in subset for trial
# model fits
Y <- Y[, colSums(cens_ID)<(0.9*nrow(cens_ID))]
cens_ID2 <- cens_ID[, colSums(cens_ID)<(0.9*nrow(cens_ID))]
cens_ID2 <- as.matrix(cens_ID2)

# add a column that gives time as number of days since 1st January 2024
pfas_y$DATE_TIME2 <- (as.numeric(pfas_y$DATE_TIME) - 1704067201)/60/60/24
pfas_y$samp_location <- as.factor(pfas_y$SAMP_SMPT_USER_REFERENCE)
pfas_y$wims_region <- as.factor(pfas_y$wims_region)

# Limit to 2024 and not by number of measurements in site ####
pfas_subset <- filter(pfas_filtered, 
                      DATE_TIME > dmy("01012024") & DATE_TIME < dmy("31122024")) %>%
  group_by(SAMP_SMPT_USER_REFERENCE) %>%
  filter(MEAS_DETERMINAND_CODE %in% colnames(Y)) %>%
  filter(sum(is.na(MEAS_SIGN))>25) %>%
  ungroup()
ggplot(pfas_subset, aes(x = DATE_TIME, y = MEAS_RESULT, colour = MEAS_DETERMINAND_CODE)) +
  geom_line() + scale_y_log10() +
  facet_wrap(~wims_region + SAMP_SMPT_USER_REFERENCE) +
  theme(legend.position = "none")

# pfas_subset <- filter(pfas_filtered, DATE_TIME > dmy("01012023")) %>%
#   group_by(SAMP_SMPT_USER_REFERENCE) %>%
#   filter(sum(is.na(MEAS_SIGN))>50) %>% 
#   ungroup

# get into jsdmstan format
pfas_y <- pfas_subset %>% ungroup() %>%
  select(SAMP_ID, MEAS_RESULT, MEAS_DETERMINAND_CODE,SAMP_SMPT_USER_REFERENCE, DATE_TIME, SAMP_PURPOSE_CODE, wims_region) %>%
  pivot_wider(values_from = MEAS_RESULT, names_from = MEAS_DETERMINAND_CODE)

dim(na.omit(pfas_y))
pfas_y <- na.omit(pfas_y)

# create censoring matrix (1 if left-censored, 0 if not). Treating all
# right-censored data as uncensored
Y <- pfas_y[, 6:ncol(pfas_y)]
cens_ID <- pfas_subset %>% ungroup() %>%
  select(SAMP_ID, MEAS_SIGN, MEAS_DETERMINAND_CODE,SAMP_SMPT_USER_REFERENCE, DATE_TIME, SAMP_PURPOSE_CODE, wims_region) %>%
  mutate(MEAS_SIGN = replace_na(ifelse(MEAS_SIGN == "<", 1, 0),0)) %>%
  pivot_wider(values_from = MEAS_SIGN, names_from = MEAS_DETERMINAND_CODE) %>%
  filter(SAMP_ID %in% pfas_y$SAMP_ID)
all.equal(cens_ID$SAMP_ID, pfas_y$SAMP_ID)
all.equal(cens_ID$DATE_TIME, pfas_y$DATE_TIME)
all.equal(cens_ID$SAMP_SMPT_USER_REFERENCE, pfas_y$SAMP_SMPT_USER_REFERENCE)
all.equal(colnames(cens_ID),colnames(pfas_y))
cens_ID <- cens_ID[,6:ncol(cens_ID)]

# remove chemicals that are mostly or wholly uncensored in subset for trial
# model fits
Y <- Y[, colSums(cens_ID)<(0.9*nrow(cens_ID))]
cens_ID2 <- cens_ID[, colSums(cens_ID)<(0.9*nrow(cens_ID))]
cens_ID2 <- as.matrix(cens_ID2)

# add a column that gives time as number of days since 1st January 2024
pfas_y$DATE_TIME2 <- (as.numeric(pfas_y$DATE_TIME) - 1704067201)/60/60/24
pfas_y$samp_location <- as.factor(pfas_y$SAMP_SMPT_USER_REFERENCE)


pfas_mod2 <- stan_jsdm(~ s(DATE_TIME2) +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sigma = "normal(0.05,0.01)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 6000, warmup = 4000)




pfas_mod3 <- stan_jsdm(~ s(DATE_TIME2) +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sp = "normal(0,1)",
                                          sigma = "normal(0.11,0.02)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 6000, warmup = 4000, thin = 5)

pfas_mod4 <- stan_jsdm(~ s(DATE_TIME2,  bs = "cc") +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sigma = "normal(0.11,0.01)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 8000, warmup = 4000, thin = 5)
# bulk ESS warning
pfas_mod4
# Family: lognormal 
#  With parameters: sigma 
# Model type: mglmm
#   Number of species: 16
#   Number of sites: 297
#   Number of predictors: 0
# 
# Model run on 4 chains with 8000 iterations per chain (4000 warmup).
# 
# Parameters with Rhat > 1.01, or Neff/N < 0.05:
#            mean      sd       15%       85%  Rhat Bulk.ESS Tail.ESS
# sigma     0.111   0.007     0.104     0.119 1.016      292      753
# lp__  16257.337 206.186 16048.865 16474.727 1.017      265      774
saveRDS(pfas_mod4, "Archive/pfas_16chem_297obs_2024_sDateCC_sSamplocRE_sigmatightnormal_thin5.rds")

theme_set(theme_classic())
# envplot(pfas_mod3)
# ggsave("Wims region effects.png", plot = envplot(pfas_mod2), path = "Archive/",
#        width = 15, height = 11, units = "cm", dpi = 600, scale = 1.2)

filter(pfas_subset, MEAS_DETERMINAND_CODE %in% c("2960","8887","8888","8889","2959")) %>% 
  select(MEAS_DETERMINAND_CODE, DETE_SHORT_DESC) %>% distinct()
corrplot(pfas_mod4, species = c("2960","2959","8887"))
ggsave("PFAS correlations subset.png", 
       plot = corrplot(pfas_mod4, species = c("2960","2959","8887")), 
       # path = "Archive/",
       width = 15, height = 6, units = "cm", dpi = 600, scale = 1.2)

pp_check(pfas_mod4, ndraws = 100, plotfun= "ecdf_overlay") + scale_x_log10()

ggsave("PFAS pp_check.png", 
       plot = pp_check(pfas_mod4, ndraws = 100, plotfun= "ecdf_overlay") + scale_x_log10(), 
       # path = "Archive/",
       width = 10, height = 8, units = "cm", dpi = 600, scale = 1.2)

Y2 <- Y
for(i in 1:nrow(Y)){
  for(j in 1:ncol(Y)){
    Y2[i,j] <- ifelse(cens_ID2[i,j] == 1, 0.5*Y[i,j], Y[i,j])
  }
}
pfas_y2 <- pfas_y
pfas_y2$sumPFAS <- rowSums(Y2)
ggplot(pfas_y2, aes(x = wims_region, y = sumPFAS)) + 
  geom_violin() + scale_y_log10() + 
  labs(x = "WIMS region", y = "Sum of PFAS chemicals")


mod_preds <- jsdm_statsummary(pfas_mod4, post_type = "predict")
pfas_y2$predMedian <- apply(mod_preds,2, median)
pfas_y2$predQ10 <- apply(mod_preds,2, quantile, prob = 0.1)
pfas_y2$predQ90 <- apply(mod_preds,2, quantile, prob = 0.9)

ggplot(pfas_y2, aes(x = sumPFAS, y = predMedian, colour = wims_region)) +
  geom_point() +
  geom_linerange(aes(ymin = predQ10, ymax = predQ90)) +
  scale_y_log10()+
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0) +
  scale_colour_manual(values = palette.colors()[2:8]) +
  labs(x = "Naive sum of PFAS values",
       y = "Model predicted summed PFAS values")
ggsave("Sum PFAS comparison.png",# path = "Archive/",
       width = 15, height = 12, units = "cm", scale = 1.2, dpi = 600)

pfas_y2$Diff <- pfas_y2$predMedian - pfas_y2$sumPFAS
(p1 <- ggplot(filter(pfas_y2, wims_region !="wimsnepr"), aes(x = wims_region, y = Diff)) +
    geom_boxplot() +
    coord_cartesian(ylim = c(-0.005,0.005)) +
    labs(x = "Region", y = "Difference between model and data")+
    theme(panel.grid.major = element_line(colour = "grey")) +
    NULL)
ggsave("Sum PFAS comparison by Region.png", #path = "Archive/",
       width = 15, height = 12, units = "cm", scale = 1.2, dpi = 600)

Y3 <- Y
for(i in 1:nrow(Y)){
  for(j in 1:ncol(Y)){
    Y3[i,j] <- ifelse(cens_ID2[i,j] == 1, 0, Y[i,j])
  }
}
pfas_y2$sumPFAS0 <- rowSums(Y3)
pfas_y2$Diff0 <- pfas_y2$predMedian - pfas_y2$sumPFAS0
p2 <- ggplot(filter(pfas_y2, wims_region !="wimsnepr"), aes(x = wims_region, y = Diff0)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(-0.005,0.005)) +
  labs(x = "Region", y = "Difference between model and data")+
  theme(panel.grid.major = element_line(colour = "grey")) +
  NULL
p1+p2 + plot_annotation(tag_levels = "a") &
  scale_x_discrete(labels =c("an","mi","nw","so","sw","th"))
ggsave("Sum PFAS comparison by Region two panel.png", #path = "Archive/",
       width = 15, height = 8, units = "cm", scale = 1.4, dpi = 600)

wrap_plots(smoothplot(pfas_mod4))

(p1 <- smoothplot(pfas_mod4, ndraws = 100, select_smooths = list(site = 1))[[1]] +
  labs(x = "Day of Year"))
ggsave("Effect of day of year.png", width = 15, height = 12, units = "cm", dpi = 600)

(p2 <- smoothplot(pfas_mod4, ndraws = 100)[[2]] + labs(x = "Sampling location") + theme(axis.text.x = element_blank()))
ggsave("Effect of sampling location no labels.png", width = 15, height = 7, units = "cm", dpi = 600, scale = 1.5)

p1+p2 + plot_annotation(tag_levels = "a") + plot_layout(ncol = 1)
ggsave("Smooth effects PFAS model.png", path = "Archive", 
       width = 15, height = 15, units = "cm", scale = 1.25, dpi = 600)

# Older tries ####

pfas_mod2 <- stan_jsdm(~ s(DATE_TIME2) +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sigma = "cauchy(0,1)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 6000, warmup = 4000)
# maximum treedepth, R-hat, E-BFMI and ESS warnings
# maximum R-hat 1.08
# DO NOT USE

pfas_mod3 <- stan_jsdm(~ s(DATE_TIME2, k = 6) +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sigma = "normal(0.11,0.01)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 8000, warmup = 4000, thin = 5)
# E-BFMI warning (four chains) plus bulk ESS warning
pfas_mod3
# Family: lognormal 
# With parameters: sigma 
# Model type: mglmm
# Number of species: 16
# Number of sites: 297
# Number of predictors: 0
# 
# Model run on 4 chains with 8000 iterations per chain (4000 warmup).
# 
# Parameters with Rhat > 1.01, or Neff/N < 0.05:
#   mean      sd       15%       85%  Rhat Bulk.ESS Tail.ESS
# cor_species_chol[16,5]    -0.145   0.116    -0.266    -0.023 1.010      845     1768
# sigma                      0.112   0.007     0.104     0.120 1.010      352      656
# lp__                   16239.374 199.480 16027.209 16450.369 1.012      331      642
saveRDS(pfas_mod3, "pfas_16chem_297obs_2024_sDatek6_sSamplocRE_sigmatightnormal_thin5.rds")

corrplot(pfas_mod3, species = c("2960","2959","8887"))
ggsave("PFAS correlations subset.png", 
       plot = corrplot(pfas_mod3, species = c("2960","2959","8887")), 
       path = "Archive/Not_CC_spline/",
       width = 15, height = 6, units = "cm", dpi = 600, scale = 1.2)

pp_check(pfas_mod3, ndraws = 100, plotfun= "ecdf_overlay") + scale_x_log10()

ggsave("PFAS pp_check.png", 
       plot = pp_check(pfas_mod3, ndraws = 100, plotfun= "ecdf_overlay") + scale_x_log10(), 
       path = "Archive/Not_CC_spline/",
       width = 10, height = 8, units = "cm", dpi = 600, scale = 1.2)

p1 <- pp_check(pfas_mod3, plotfun = "ecdf_overlay", discrete = TRUE)
p2 <- multi_pp_check(pfas_mod3, plotfun = "ecdf_overlay", 
                     species = sample.int(16,5), discrete = TRUE)
p1+ggtitle("Sum PFAS")+pl_list + 
  plot_annotation(tag_levels = "a", tag_suffix = ")") +
  plot_layout(guides = "collect")& scale_x_log10()

ggsave("PFAS pp_check v2.png", 
       path = "Archive/Not_CC_spline/",
       width = 15, height = 8, units = "cm", dpi = 600, scale = 1.2)


Y2 <- Y
for(i in 1:nrow(Y)){
  for(j in 1:ncol(Y)){
    Y2[i,j] <- ifelse(cens_ID2[i,j] == 1, 0.5*Y[i,j], Y[i,j])
  }
}
pfas_y2 <- pfas_y
pfas_y2$sumPFAS <- rowSums(Y2)
ggplot(pfas_y2, aes(x = wims_region, y = sumPFAS)) + 
  geom_violin() + scale_y_log10() + 
  labs(x = "WIMS region", y = "Sum of PFAS chemicals")


mod_preds <- jsdm_statsummary(pfas_mod3, post_type = "predict")
pfas_y2$predMedian <- apply(mod_preds,2, median)
pfas_y2$predQ10 <- apply(mod_preds,2, quantile, prob = 0.1)
pfas_y2$predQ90 <- apply(mod_preds,2, quantile, prob = 0.9)

ggplot(pfas_y2, aes(x = sumPFAS, y = predMedian, colour = wims_region)) +
  geom_point() +
  geom_linerange(aes(ymin = predQ10, ymax = predQ90)) +
  scale_y_log10()+
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0) +
  scale_colour_manual(values = palette.colors()[2:8]) +
  labs(x = "Naive sum of PFAS values",
       y = "Model predicted summed PFAS values")
ggsave("Sum PFAS comparison.png", path = "Archive/Not_CC_spline/",
       width = 15, height = 12, units = "cm", scale = 1.2, dpi = 600)

pfas_y2$Diff <- pfas_y2$predMedian - pfas_y2$sumPFAS
(p1 <- ggplot(filter(pfas_y2, wims_region !="wimsnepr"), aes(x = wims_region, y = Diff)) +
    geom_boxplot() +
    coord_cartesian(ylim = c(-0.005,0.005)) +
    labs(x = "Region", y = "Difference between model and data")+
    theme(panel.grid.major = element_line(colour = "grey")) +
    NULL)
ggsave("Sum PFAS comparison by Region.png", path = "Archive/Not_CC_spline/",
       width = 15, height = 12, units = "cm", scale = 1.2, dpi = 600)

Y3 <- Y
for(i in 1:nrow(Y)){
  for(j in 1:ncol(Y)){
    Y3[i,j] <- ifelse(cens_ID2[i,j] == 1, 0, Y[i,j])
  }
}
pfas_y2$sumPFAS0 <- rowSums(Y3)
pfas_y2$Diff0 <- pfas_y2$predMedian - pfas_y2$sumPFAS0
p2 <- ggplot(filter(pfas_y2, wims_region !="wimsnepr"), aes(x = wims_region, y = Diff0)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(-0.005,0.005)) +
  labs(x = "Region", y = "Difference between model and data")+
  theme(panel.grid.major = element_line(colour = "grey")) +
  NULL
p1+p2 + plot_annotation(tag_levels = "a") &
  scale_x_discrete(labels =c("an","mi","nw","so","sw","th"))
ggsave("Sum PFAS comparison by Region two panel.png", path = "Archive/Not_CC_spline/",
       width = 15, height = 8, units = "cm", scale = 1.4, dpi = 600)


(p1 <- smoothplot(pfas_mod3, ndraws = 100, select_smooths = list(site = 1))[[1]] +
    labs(x = "Day of Year"))

(p2 <- smoothplot(pfas_mod4, ndraws = 100)[[2]] + labs(x = "Sampling location") + theme(axis.text.x = element_blank()))

p1+p2 + plot_annotation(tag_levels = "a") + plot_layout(ncol = 1)
ggsave("Smooth effects PFAS model.png", path = "Archive/Not_CC_spline/", 
       width = 15, height = 15, units = "cm", scale = 1.25, dpi = 600)




pfas_mod2 <- stan_jsdm(~ s(DATE_TIME2) +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)"),
                       method = "mglmm", family = "gamma",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 6000, warmup = 4000)
# 2506 divergent transitions
# Took 12.5 hours




pfas_mod2 <- stan_jsdm(~ s(DATE_TIME2) +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sigma = "cauchy(0,0.5)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 6000, warmup = 4000)
# maximum treedepth, E-BFMI, bulk and tail ESS warnings
pfas_mod2
# Family: lognormal 
# With parameters: sigma 
# Model type: mglmm
# Number of species: 16
# Number of sites: 297
# Number of predictors: 0
# 
# Model run on 4 chains with 6000 iterations per chain (4000 warmup).
# 
# Parameters with Rhat > 1.01, or Neff/N < 0.05:
#   mean      sd       15%       85%  Rhat Bulk.ESS Tail.ESS
# betas[1,2]                 -5.435   0.127    -5.568    -5.304 1.019      487      924
# betas[1,3]                 -6.273   0.127    -6.405    -6.141 1.018      480      887
# betas[1,4]                 -5.923   0.126    -6.053    -5.792 1.019      487      947
# betas[1,6]                 -5.453   0.127    -5.583    -5.323 1.018      497      803
# betas[1,7]                 -7.440   0.201    -7.649    -7.235 1.012      775     2295
# betas[1,8]                 -5.844   0.131    -5.981    -5.711 1.017      539     1053
# betas[1,9]                 -7.427   0.131    -7.562    -7.294 1.016      533     1092
# betas[1,10]                -8.756   0.180    -8.937    -8.572 1.012      815     2470
# betas[1,11]                -6.156   0.127    -6.287    -6.025 1.018      489      934
# betas[1,12]                -6.305   0.126    -6.436    -6.176 1.018      491      962
# betas[1,13]                -6.132   0.129    -6.265    -6.001 1.018      501     1013
# betas[1,14]                -8.039   0.135    -8.178    -7.899 1.016      550     1052
# betas[1,15]                -7.684   0.128    -7.816    -7.553 1.018      492      952
# betas[1,16]                -5.520   0.127    -5.653    -5.390 1.019      487      886
# sigmas_species[2]           0.382   0.022     0.358     0.405 1.013      427     1220
# sigmas_species[3]           0.304   0.021     0.282     0.325 1.017      291     1149
# sigmas_species[16]          0.358   0.020     0.337     0.379 1.017      358     1827
# z_species[3,96]             2.977   0.851     2.089     3.853 1.011      425     1619
# z_species[3,193]            4.081   0.933     3.110     5.036 1.018      278     1134
# cor_species_chol[3,2]       0.582   0.050     0.531     0.635 1.016      343      982
# cor_species_chol[3,3]       0.807   0.038     0.768     0.844 1.016      335     1069
# cor_species_chol[4,3]       0.431   0.069     0.361     0.500 1.016      351      798
# cor_species_chol[4,4]       0.885   0.036     0.851     0.919 1.015      363      779
# cor_species_chol[16,2]      0.613   0.041     0.571     0.655 1.011      474     1447
# cor_species_chol[16,4]      0.030   0.053    -0.025     0.085 1.011      976     1832
# cor_species_chol[16,16]     0.388   0.163     0.199     0.570 1.014      591     1710
# sigma                       0.111   0.011     0.100     0.123 1.046      105      194
# nfs_b[12]                   0.250   0.154     0.091     0.407 1.012      644     1357
# nfs_b[14]                  -0.545   0.160    -0.710    -0.380 1.013      737     1703
# nfs_b[17]                  -1.351   0.159    -1.514    -1.187 1.012      759     1889
# nfs_b[18]                  -0.618   0.168    -0.793    -0.445 1.011      717     1306
# nfs_b[20]                  -1.445   0.175    -1.627    -1.266 1.011      915     1759
# nfs_b[21]                  -0.550   0.175    -0.734    -0.369 1.011      614     1142
# nfs_b[22]                  -0.910   0.195    -1.109    -0.708 1.010      693     1298
# nfs_b[23]                  -0.204   0.176    -0.385    -0.020 1.010      782     1790
# nfs_b[24]                  -0.828   0.163    -0.997    -0.660 1.014      771     2028
# nfs_b[27]                  -1.366   0.163    -1.536    -1.200 1.011      688     1846
# nfs_b[28]                  -0.541   0.150    -0.694    -0.385 1.013      636     1343
# nfs_b[29]                  -0.996   0.154    -1.153    -0.838 1.014      657     1355
# nfs_b[30]                  -1.560   0.160    -1.726    -1.395 1.011      719     1523
# nfs_b[31]                  -1.688   0.164    -1.857    -1.518 1.011      759     1963
# nfs_b[32]                  -1.420   0.166    -1.591    -1.247 1.010      737     2103
# nfs_b[33]                   0.457   0.155     0.297     0.617 1.014      716     1463
# nfs_b[34]                  -0.745   0.167    -0.918    -0.571 1.012      872     2217
# nfs_b[35]                   0.725   0.151     0.568     0.884 1.012      664     1375
# nfs_b[36]                   0.355   0.161     0.187     0.521 1.012      678     1543
# nfs_b[45]                   0.456   0.153     0.296     0.613 1.012      624     1338
# nfs_b[52]                   0.669   0.147     0.517     0.821 1.016      630     1332
# nfs_b[54]                   0.307   0.149     0.153     0.461 1.012      596     1432
# nfs_b[56]                   0.661   0.153     0.501     0.822 1.013      626     1409
# nfs_b[60]                   0.364   0.157     0.201     0.527 1.012      660     1363
# nfs_b[63]                  -0.569   0.157    -0.731    -0.408 1.012      751     1611
# nfs_b[64]                  -0.426   0.158    -0.589    -0.261 1.013      657     1208
# cor_species[2,3]            0.584   0.050     0.533     0.637 1.016      340     1056
# cor_species[2,16]           0.614   0.040     0.573     0.656 1.011      464     1345
# cor_species[3,2]            0.584   0.050     0.533     0.637 1.016      340     1056
# cor_species[3,16]           0.623   0.059     0.561     0.686 1.022      316     1218
# cor_species[16,2]           0.614   0.040     0.573     0.656 1.011      464     1345
# cor_species[16,3]           0.623   0.059     0.561     0.686 1.022      316     1218
# u[193,3]                    1.130   0.166     0.957     1.301 1.013      315     1341
# lp__                    16272.442 286.698 15977.412 16571.182 1.049      103      181
# cor_species_chol[7,5]      -0.094   0.167    -0.268     0.082 1.009      372      808
# cor_species[5,5]            1.000   0.000     1.000     1.000 1.000     7398     7491
# cor_species[14,14]          1.000   0.000     1.000     1.000 1.000     8109     7615


pfas_mod3 <- stan_jsdm(~ s(DATE_TIME2, k = 6) +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sigma = "normal(0.11,0.01)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 6000, warmup = 4000)
# E-BFMI, bulk and tail ESS warnings
pfas_mod3
# Family: lognormal 
#  With parameters: sigma 
# Model type: mglmm
#   Number of species: 16
#   Number of sites: 297
#   Number of predictors: 0
# 
# Model run on 4 chains with 6000 iterations per chain (4000 warmup).
# 
# Parameters with Rhat > 1.01, or Neff/N < 0.05:
#                             mean      sd       15%       85%  Rhat Bulk.ESS Tail.ESS
# cor_species_chol[7,5]     -0.090   0.169    -0.265     0.092 1.011      441     1071
# cor_species_chol[12,8]    -0.040   0.074    -0.114     0.035 1.010     1146     2017
# cor_species_chol[13,7]     0.232   0.108     0.118     0.346 1.010      479     1389
# sigma                      0.111   0.007     0.104     0.118 1.026      202      271
# cor_species[5,7]           0.009   0.120    -0.116     0.136 1.010      448     1327
# cor_species[7,5]           0.009   0.120    -0.116     0.136 1.010      448     1327
# lp__                   16273.306 196.342 16074.357 16467.028 1.028      192      255
# cor_species[5,5]           1.000   0.000     1.000     1.000 1.000     8032     7500
# cor_species[14,14]         1.000   0.000     1.000     1.000 1.000     7557     7693

pfas_mod4 <- stan_jsdm(~ s(DATE_TIME2, k = 6, bs = "cc") +  
                         s(samp_location, bs = "re"),
                       data = pfas_y, Y = Y,
                       prior = jsdm_prior(betas = "student_t(3,-7,2)",
                                          sigma = "normal(0.11,0.01)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4,
                       control = list(adapt_delta = 0.95), 
                       iter = 8000, warmup = 4000, thin = 5)
# E-BFMI and bulk ESS warnings
pfas_mod4
# Family: lognormal 
#  With parameters: sigma 
# Model type: mglmm
#   Number of species: 16
#   Number of sites: 297
#   Number of predictors: 0
# 
# Model run on 4 chains with 8000 iterations per chain (4000 warmup).
# 
# Parameters with Rhat > 1.01, or Neff/N < 0.05:
#               mean    sd    15%    85%  Rhat Bulk.ESS Tail.ESS
# betas[1,2]  -5.423 0.126 -5.554 -5.295 1.013      515     1057
# betas[1,3]  -6.261 0.125 -6.391 -6.133 1.013      510     1304
# betas[1,4]  -5.911 0.124 -6.039 -5.784 1.013      521     1249
# betas[1,5]  -8.175 0.168 -8.349 -8.005 1.011      710     1071
# betas[1,6]  -5.440 0.125 -5.567 -5.313 1.012      533     1226
# betas[1,8]  -5.832 0.129 -5.964 -5.698 1.012      565     1236
# betas[1,9]  -7.416 0.128 -7.549 -7.287 1.013      556     1176
# betas[1,11] -6.144 0.126 -6.272 -6.017 1.015      494     1093
# betas[1,12] -6.294 0.126 -6.423 -6.166 1.014      506     1172
# betas[1,13] -6.121 0.129 -6.254 -5.991 1.013      518     1076
# betas[1,14] -8.030 0.134 -8.164 -7.895 1.014      550     1249
# betas[1,15] -7.673 0.126 -7.805 -7.545 1.015      500     1248
# betas[1,16] -5.507 0.126 -5.637 -5.380 1.013      523     1148
# nfs_b[7]     0.239 0.154  0.085  0.397 1.011      781     1444
# nfs_b[12]   -1.368 0.157 -1.530 -1.209 1.012      712     1668
# nfs_b[22]   -1.384 0.160 -1.550 -1.219 1.011      761     1352
# nfs_b[23]   -0.563 0.151 -0.720 -0.409 1.011      693     1208
# nfs_b[24]   -1.006 0.151 -1.160 -0.852 1.013      709     1514
# nfs_b[25]   -1.569 0.159 -1.733 -1.409 1.012      653     1372
# nfs_b[26]   -1.712 0.165 -1.877 -1.545 1.011      714     1615
# nfs_b[27]   -1.459 0.163 -1.621 -1.288 1.013      783     2047
# nfs_b[28]    0.436 0.155  0.274  0.599 1.011      638     1307
# nfs_b[47]    0.664 0.147  0.513  0.810 1.012      682     1517