# read in and explore PFAS data
load("~/PFAS_FW_2AZZ_methods 21_24.rda")
library(tidyverse)
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
library(ggplot2)
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

pfas_subset <- filter(pfas_filtered, DATE_TIME > dmy("01012023")) %>%
  group_by(SAMP_SMPT_USER_REFERENCE) %>%
  filter(sum(is.na(MEAS_SIGN))>50) %>% 
  ungroup

# get into jsdmstan format
pfas_y <- pfas_subset %>% ungroup() %>%
  select(SAMP_ID, MEAS_RESULT, MEAS_DETERMINAND_CODE,SAMP_SMPT_USER_REFERENCE, DATE_TIME, SAMP_PURPOSE_CODE, wims_region) %>%
  pivot_wider(values_from = MEAS_RESULT, names_from = MEAS_DETERMINAND_CODE)

dim(na.omit(pfas_y))
pfas_y <- na.omit(pfas_y)

# create censoring matrix (1 if left-censored, 0 if not). Treating all
# right-censored data as uncensored
Y <- pfas_y[, 6:52]
cens_ID <- pfas_subset %>% ungroup() %>%
  select(SAMP_ID, MEAS_SIGN, MEAS_DETERMINAND_CODE,SAMP_SMPT_USER_REFERENCE, DATE_TIME, SAMP_PURPOSE_CODE, wims_region) %>%
  mutate(MEAS_SIGN = replace_na(ifelse(MEAS_SIGN == "<", 1, 0),0)) %>%
  pivot_wider(values_from = MEAS_SIGN, names_from = MEAS_DETERMINAND_CODE) %>%
  filter(SAMP_ID %in% pfas_y$SAMP_ID)
all.equal(cens_ID$SAMP_ID, pfas_y$SAMP_ID)
all.equal(cens_ID$DATE_TIME, pfas_y$DATE_TIME)
all.equal(cens_ID$SAMP_SMPT_USER_REFERENCE, pfas_y$SAMP_SMPT_USER_REFERENCE)
all.equal(colnames(cens_ID),colnames(pfas_y))
cens_ID <- cens_ID[,6:52]

# remove chemicals that are mostly or wholly uncensored in subset for trial
# model fits
Y <- Y[, colSums(cens_ID)<423]
cens_ID2 <- cens_ID[, colSums(cens_ID)<423]
cens_ID2 <- as.matrix(cens_ID2)

# add a column that gives time as number of days since 1st January 2024
# pfas_y$DATE_TIME2 <- (as.numeric(pfas_y$DATE_TIME) - 1704067201)/60/60/24
pfas_y$samp_location <- as.factor(pfas_y$SAMP_SMPT_USER_REFERENCE)
pfas_y$wims_region <- as.factor(pfas_y$wims_region)

# model fit
pfas_mod2 <- stan_jsdm(~ wims_region +  
                        s(samp_location, bs = "re"),
                      data = pfas_y, Y = Y,
                      prior = jsdm_prior(betas = "student_t(3,0,1)",
                      #                    sp = "normal(0,1)",
                                         sigma = "normal(0.1,0.01)"),
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID2, cores = 4)
saveRDS(pfas_mod2, "Archive/pfas_mod_wimsfixef_sampranef_since2023.rds")

theme_set(theme_classic())
envplot(pfas_mod2)
ggsave("Wims region effects.png", plot = envplot(pfas_mod2), path = "Archive/",
       width = 15, height = 11, units = "cm", dpi = 600, scale = 1.2)

filter(pfas_subset, MEAS_DETERMINAND_CODE %in% c("2960","8887","8888","8889","2959")) %>% 
  select(MEAS_DETERMINAND_CODE, DETE_SHORT_DESC) %>% distinct()
corrplot(pfas_mod2, species = c("2960","2959","8887"))
ggsave("PFAS correlations subset.png", 
       plot = corrplot(pfas_mod2, species = c("2960","2959","8887")), 
       path = "Archive/",
       width = 15, height = 6, units = "cm", dpi = 600, scale = 1.2)

pp_check(pfas_mod2, ndraws = 100, plotfun= "ecdf_overlay") + scale_x_log10()

ggsave("PFAS pp_check.png", 
       plot = pp_check(pfas_mod2, ndraws = 100, plotfun= "ecdf_overlay") + scale_x_log10(), 
       path = "Archive/",
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


mod_preds <- jsdm_statsummary(pfas_mod2, post_type = "predict")
pfas_y2$predMedian <- apply(mod_preds,2, median)
pfas_y2$predQ10 <- apply(mod_preds,2, quantile, prob = 0.1)
pfas_y2$predQ90 <- apply(mod_preds,2, quantile, prob = 0.9)

ggplot(pfas_y2, aes(x = sumPFAS, y = predMedian, colour = wims_region)) +
  geom_point() +
  scale_y_log10()+
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0) +
  scale_colour_manual(values = palette.colors())

ggplot(pfas_y2, aes(x = sumPFAS, y = predMedian, colour = wims_region)) +
  geom_point() +
  geom_linerange(aes(ymin = predQ10, ymax = predQ90)) +
  scale_y_log10()+
  scale_x_log10() +
  geom_abline(slope = 1, intercept = 0) +
  scale_colour_manual(values = palette.colors()[2:8]) +
  labs(x = "Naive sum of PFAS values",
       y = "Model predicted summed PFAS values")
ggsave("Sum PFAS comparison.png", path = "Archive/",
       width = 15, height = 12, units = "cm", scale = 1.2, dpi = 600)

pfas_y2$Diff <- pfas_y2$predMedian - pfas_y2$sumPFAS
p1 <- ggplot(filter(pfas_y2, wims_region !="wimsnepr"), aes(x = wims_region, y = Diff)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(-0.005,0.005)) +
  labs(x = "WIMS region", y = "Difference between model and data")+
  theme(panel.grid.major = element_line(colour = "grey")) +
  NULL
ggsave("Sum PFAS comparison by Region.png", path = "Archive/",
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
  labs(x = "WIMS region", y = "Difference between model and data")+
  theme(panel.grid.major = element_line(colour = "grey")) +
  NULL
p1+p2 + plot_annotation(tag_levels = "a") &
  scale_x_discrete(labels =c("an","mi","nw","so","sw","th"))
ggsave("Sum PFAS comparison by Region two panel.png", path = "Archive/",
       width = 15, height = 8, units = "cm", scale = 1.4, dpi = 600)

filter(pfas_subset, MEAS_DETERMINAND_CODE %in% colnames(Y)) %>% 
  select(MEAS_DETERMINAND_CODE, DETE_SHORT_DESC, UNIT_SHORT_DESC) %>% distinct()
# pfas_mod_gamma <- stan_jsdm(~ wims_region +  
#                               s(samp_location, bs = "re"),
#                             data = pfas_y, Y = Y,
#                             prior = jsdm_prior(betas = "student_t(3,0,1)"),
#                             method = "mglmm", family = "gamma",
#                             censoring = "left", cens_ID = cens_ID2, cores = 4)


pfas_mod <- stan_jsdm(~s(DATE_TIME2) + s(DATE_TIME2, species, bs = "fs") +
                        s(DATE_TIME2, samp_location, bs = "fs"),
                      data = pfas_y, Y = Y,
                      prior = jsdm_prior(betas = "student_t(3,-6,1)",
                                         sp = "normal(0,1)",
                                         sigma = "normal(0.1,0.01)"),
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID2, cores = 4)

pfas_mod2 <- stan_jsdm(~s(DATE_TIME2, samp_location, bs = "fs"),
                      data = pfas_y, Y = Y,
                      prior = jsdm_prior(betas = "student_t(3,-6,1)",
                                         sp = "normal(0,1)",
                                         sigma = "normal(0.1,0.01)"),
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID2, cores = 4)
# model failed :(
pfas_mod
plot(pfas_mod)


pfas_mod2 <- stan_jsdm(~DATE_TIME2 + s(samp_location, bs = "re"), #+ s(DATE_TIME2, species, bs = "fs") +
                      #  s(DATE_TIME2, samp_location, bs = "fs"),
                      data = pfas_y, Y = Y,
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID2, cores = 4)
pp_check(pfas_mod2, plotfun = "ecdf_overlay")
multi_pp_check(pfas_mod2, plotfun = "ecdf_overlay")
corrplot(pfas_mod2)

pfas_mod3 <- stan_jsdm(~s(DATE_TIME2) + s(DATE_TIME2, samp_location, bs = "fs"),
                      data = pfas_y, Y = Y,
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID2, cores = 4)


pfas_mod4 <- stan_jsdm(~s(DATE_TIME2),
                       data = pfas_y, Y = Y,
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4)


pfas_mod5 <- stan_jsdm(~s(DATE_TIME2, samp_location, bs = "fs"),
                       data = pfas_y, Y = Y,
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4)

pfas_mod6 <- stan_jsdm(~s(DATE_TIME2, samp_location, bs = "fs"),
                       data = pfas_y, Y = Y, prior = jsdm_prior(sigma = "cauchy(0,1)"),
                       method = "mglmm", family = "lognormal",
                       censoring = "left", cens_ID = cens_ID2, cores = 4)

# PAHs ####
# questions for Mike: what are the SAMP_MATERIAL and SAMP_PURPOSE_CODE for?
# Loads of sites only have a subset of the six measured, what does he want us to do with that?
# Does he have a particular type of model structure in mind for this?
# load and explore PAHs data
load("~/ecw_pahs.rda")
pahs_2124 <- filter(ecw.pah2, MEAS_ANAL_METH_CODE %in% c(21,24))
mice::md.pattern(pahs_2124, rotate.names = TRUE)

# get PAHs data into format (here removing two most commonly missed PAH
# chemicals, but those could be added back in)
pas_2124_subset <- pahs_2124 %>%
  # filter(SAMP_MATERIAL == "2HZZ" & SAMP_PURPOSE_CODE == "MS") %>%
  select(SMPT_EASTING, SMPT_NORTHING, SMPT_SHORT_NAME, DATE_TIME, SAMP_ID,# MEAS_SIGN,
         # WB_NAME, WB_CAT, RBD_NAME, OPCAT_NAME,
         wims_region,
         MEAS_RESULT, MEAS_DETERMINAND_CODE) %>%
  pivot_wider(names_from = MEAS_DETERMINAND_CODE, values_from = MEAS_RESULT,
              names_prefix = "PAH") %>%
  select(-PAH8940, -PAH0772)
mice::md.pattern(pas_2124_subset, rotate.names = TRUE)

pas_2124_subset <- na.omit(pas_2124_subset)

# explore via plots - adding a GBR outline so I can see where they are
library(sf)
gbr <- st_read("N:/Data/Shapefiles/GBR_adm/GBR_adm1.shp")
eng <- filter(gbr, NAME_1 == "England")
pas_sf <- pas_2124_subset %>%
  count(wims_region, SMPT_EASTING, SMPT_NORTHING, SMPT_SHORT_NAME) %>%
  st_as_sf(coords = c("SMPT_EASTING","SMPT_NORTHING"),
                   crs = 27700)
ggplot() +
  geom_sf(data = eng) +
  geom_sf(data = pas_sf, mapping = aes(colour = n))

ggplot(pas_2124_subset, aes(x = DATE_TIME, y = SMPT_NORTHING)) +
  geom_point()

# create neighbourhood matrix (this is when I wanted to see if I could fit a mrf
# model to the wims regions. I have concluded this is almost definitely a
# terrible idea, but don't want to delete it in case I need it later)
nb_wims <- list(wimsanpr = c(2,3,7),
                wimsmipr = c(1,3,4,7),
                wimsnepr = c(1,2,4),
                wimsnwpr = c(2,3),
                wimssopr = c(6,7),
                wimsswpr = c(5,7),
                wimsthpr = c(1,2,5,6))

# get into jsdmstan format
pahs_Y <- select(pas_2124_subset, starts_with("PAH"))

# make cens_ID matrix
cens_ID_pah <- pahs_2124 %>%
  # filter(SAMP_MATERIAL == "2HZZ" & SAMP_PURPOSE_CODE == "MS") %>%
  select(SMPT_EASTING, SMPT_NORTHING, SMPT_SHORT_NAME, DATE_TIME, MEAS_SIGN, SAMP_ID,
         # WB_NAME, WB_CAT, RBD_NAME, OPCAT_NAME,
         wims_region,
         MEAS_DETERMINAND_CODE) %>%
  mutate(MEAS_SIGN = replace_na(ifelse(MEAS_SIGN == "<", 1, 0),0)) %>%
  pivot_wider(names_from = MEAS_DETERMINAND_CODE, values_from = MEAS_SIGN,
              names_prefix = "PAH") %>%
  filter(SAMP_ID %in% pas_2124_subset$SAMP_ID) %>%
  filter(SAMP_ID != 1537808 | wims_region != "wimsnepr")
all.equal(cens_ID_pah$SAMP_ID, pas_2124_subset$SAMP_ID)
all.equal(cens_ID_pah$DATE_TIME, pas_2124_subset$DATE_TIME)
all.equal(cens_ID_pah$wims_region, pas_2124_subset$wims_region)
all.equal(colnames(select(cens_ID_pah, starts_with("PAH"))),colnames(pahs_Y))

cens_ID_pah <- as.matrix(select(cens_ID_pah, starts_with("PAH")))

# model that I haven't even tried to run yet and probably never will
pahs_mod <- stan_jsdm(~s(wims_region, bs = "mrf", xt = list(nb = nb_wims)) ,
                      data = pas_2124_subset, Y = pahs_Y,
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID_pah)

pahs_mod <- stan_jsdm(~s(wims_region, bs = "mrf", xt = list(nb = nb_wims)) ,
                      data = pas_2124_subset, Y = pahs_Y,
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID_pah)


# limit to Severn

pas_2124_subset <- pahs_2124 %>%
  # filter(SAMP_MATERIAL == "2HZZ" & SAMP_PURPOSE_CODE == "MS") %>%
  select(SMPT_EASTING, SMPT_NORTHING, SMPT_SHORT_NAME, DATE_TIME, SAMP_ID,# MEAS_SIGN,
         # WB_NAME, WB_CAT, RBD_NAME, OPCAT_NAME,
         wims_region,
         MEAS_RESULT, MEAS_DETERMINAND_CODE) %>%
  pivot_wider(names_from = MEAS_DETERMINAND_CODE, values_from = MEAS_RESULT,
              names_prefix = "PAH") %>%
  select(-PAH8940, -PAH0772) %>%
  filter(grepl("SEVERN",SMPT_SHORT_NAME))
mice::md.pattern(pas_2124_subset, rotate.names = TRUE)

pas_2124_subset <- na.omit(pas_2124_subset)

# explore via plots - adding a GBR outline so I can see where they are
pas_sf <- pas_2124_subset %>%
  count(wims_region, SMPT_EASTING, SMPT_NORTHING, SMPT_SHORT_NAME) %>%
  st_as_sf(coords = c("SMPT_EASTING","SMPT_NORTHING"),
           crs = 27700)
ggplot() +
  geom_sf(data = eng) +
  geom_sf(data = pas_sf, mapping = aes(colour = n))

ggplot(pas_2124_subset, aes(x = DATE_TIME, y = SMPT_NORTHING)) +
  geom_point()

pas_2124_subset <- filter(pas_2124_subset, DATE_TIME > dmy("01012020"))

ggplot(pas_2124_subset, aes(x = DATE_TIME, y = SMPT_NORTHING)) +
  geom_point()


# get into jsdmstan format
pahs_Y <- select(pas_2124_subset, starts_with("PAH"))

# make cens_ID matrix
cens_ID_pah <- pahs_2124 %>%
  # filter(SAMP_MATERIAL == "2HZZ" & SAMP_PURPOSE_CODE == "MS") %>%
  select(SMPT_EASTING, SMPT_NORTHING, SMPT_SHORT_NAME, DATE_TIME, MEAS_SIGN, SAMP_ID,
         # WB_NAME, WB_CAT, RBD_NAME, OPCAT_NAME,
         wims_region,
         MEAS_DETERMINAND_CODE) %>%
  mutate(MEAS_SIGN = replace_na(ifelse(MEAS_SIGN == "<", 1, 0),0)) %>%
  pivot_wider(names_from = MEAS_DETERMINAND_CODE, values_from = MEAS_SIGN,
              names_prefix = "PAH")%>%
  select(-PAH8940, -PAH0772) %>%
  filter(SAMP_ID %in% pas_2124_subset$SAMP_ID)
all.equal(cens_ID_pah$SAMP_ID, pas_2124_subset$SAMP_ID)
all.equal(cens_ID_pah$DATE_TIME, pas_2124_subset$DATE_TIME)
all.equal(cens_ID_pah$wims_region, pas_2124_subset$wims_region)
all.equal(colnames(select(cens_ID_pah, starts_with("PAH"))),colnames(pahs_Y))

cens_ID_pah <- as.matrix(select(cens_ID_pah, starts_with("PAH")))
pas_2124_subset$SMPT <- as.factor(pas_2124_subset$SMPT_SHORT_NAME)
pas_2124_subset$yday <- yday(pas_2124_subset$DATE_TIME)

pahs_mod <- stan_jsdm(~s(SMPT_NORTHING, k = 4) +
                        s(yday, bs = "cc") +
                        s(SMPT, bs = "re") ,
                      data = pas_2124_subset, Y = pahs_Y,
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID_pah,
                      cores = 4)

pahs_mod2 <- stan_jsdm(~s(yday, bs = "cc") + s(SMPT, bs = "re"),
                      data = pas_2124_subset, Y = pahs_Y,
                      method = "mglmm", family = "lognormal",
                      censoring = "left", cens_ID = cens_ID_pah,
                      cores = 4)


pahs_mod3 <- stan_jsdm(~s(SMPT_NORTHING, k = 4) +
                        s(yday, bs = "cc") +
                        s(SMPT, bs = "re") ,
                      data = pas_2124_subset, Y = pahs_Y,
                      method = "mglmm", family = "lognormal",
                      prior = jsdm_prior(sigma = "normal(0.14,0.01)",
                                         sp = "normal(0,0.5)",
                                         betas = "normal(-2,1)"),
                      censoring = "left", cens_ID = cens_ID_pah,
                      cores = 4)



pahs_mod4 <- stan_jsdm(~s(SMPT_NORTHING, k = 4) +
                         s(yday, bs = "cc") +
                         s(SMPT, bs = "re") ,
                       data = pas_2124_subset, Y = pahs_Y,
                       method = "mglmm", family = "gamma",
                       # prior = jsdm_prior(sigma = "normal(0.14,0.01)",
                       #                    sp = "normal(0,0.5)",
                       #                    betas = "normal(-2,1)"),
                       censoring = "left", cens_ID = cens_ID_pah,
                       cores = 4, control = list(adapt_delta = ))


pah_data <- filter(pahs_2124, grepl("SEVERN", SMPT_SHORT_NAME) & DATE_TIME < dmy("01012018") & DATE_TIME > dmy("01012015")) %>%
  mutate(SMPT = case_match(SMPT_SHORT_NAME, "R SEVERN - BULLO PILL" ~ "BULLO PILL",
                           "R SEVERN - GATCOMBE" ~ "GATCOMBE",
                           "R SEVERN - HAYWARD ROCK" ~ "HAYWARD ROCK",
                           "R SEVERN (TIDAL) OVER BRIDGE" ~ "T OVER BRIDGE",
                           "R SEVERN (TIDAL), STONE CHUTE MINST'WRTH" ~ "T SC MINSTWRTH",
                           "R SEVERN SEDBURY CLIFFS" ~ "SEDBURY CLIFFS",
                           "RIVER SEVERN ESTUARY HOPE FIXED STATION" ~ "EST HOPE FIX ST",
                           "RIVER SEVERN ESTUARY OFF NEWPORT" ~ "EST OFF NEWPORT",
                           "RIVER SEVERN TIDAL, CHARSTON ROCKS" ~ "T CHARSTON RCKS",
                           "R SEVERN (TIDAL) 150M US CROMPTONS O/L" ~ "T US CROMPTONS",
                           "R SEVERN (TIDAL) 250M D/S LYDNEY OUTFALL" ~ "T LYDNEY OUTFALL")) %>%
  select(SMPT_NORTHING, SMPT_EASTING, SMPT, DATE_TIME, MEAS_RESULT, MEAS_DETERMINAND_CODE, MEAS_SIGN, SAMP_ID)
write.csv(pah_data, "~/Data/pah_data.csv", row.names=FALSE)
