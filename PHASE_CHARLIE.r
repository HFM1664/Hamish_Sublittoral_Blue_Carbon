###########################################################################################################################
                      ### PHASE CHARLIE; Construction of Surface area x Carbon estimates ###
###########################################################################################################################
## Author Hamish Mckee 14th to 15th of Nov 2026 
## See workflow documents for more

### Required libraries ###
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(clipr)             

### Load data ###
AFDM <- read_csv("20261117_AFDM.csv",
                 na = c("", "NA", "N/A", "n/a", "#VALUE!", "#REF!"))   

head(AFDM)
str(AFDM)


### Create a new data frame with only the relevant columns ###
AFDM_min <- AFDM %>%
  filter(!is.na(SAMPLE_ID)) %>%                 
  select(SAMPLE_ID, OTU, Surface_Area_mm2, BIOMASS, AFDM) %>%
  rename(SA = Surface_Area_mm2,                 
         BM = BIOMASS,
         LOI = AFDM) %>%
  mutate(across(c(SA, BM, LOI), as.numeric),
         org_frac = LOI / BM)                   # organic fraction (should be 0-1)

str(AFDM_min)

### Removal of any samples with organics fraction greater than 1 (likely due to measurement error) ###
AFDM_min %>% filter(org_frac > 1)
 
AFDM_min <- AFDM_min %>% filter(is.na(org_frac) | org_frac <= 1)


### Changes to make to data frame 1st update AFDM_clean on mastersheet, then update here with data - when next iteration loaded, check that the changes have been made correctly. ###
AFDM_min <- AFDM_min %>%
  mutate(BM = if_else(SAMPLE_ID == "PA13", 0.5413, BM)) ### PA13 Crucible mass was incorrect adjusted BM 13/11/2023

AFDM_min <- AFDM_min %>%
  mutate(OTU = if_else(SAMPLE_ID == "WA20", "ALG001H", OTU))
AFDM_min <- AFDM_min %>%
  mutate(OTU = if_else(SAMPLE_ID == "WA19", "ALG001H", OTU))
AFDM_min <- AFDM_min %>%
  mutate(OTU = if_else(SAMPLE_ID == "WA21", "ALG001S", OTU))
AFDM_min %>% filter(SAMPLE_ID == "WA20") %>% select(SAMPLE_ID, OTU)

AFDM_min <- AFDM_min %>% filter(SAMPLE_ID != "PA18")   # BM way to low
AFDM_min <- AFDM_min %>% filter(SAMPLE_ID != "PB58")   # LOI to high unfortunately - liquidate
AFDM_min <- AFDM_min %>% filter(SAMPLE_ID != "PA1")  ## To much uncertinity with crucible mass estimate


###########################################################################################################################
### Functions ###
###########################################################################################################################
{
### Function to produce Labelled scatter: every point named, so you can see who is odd
plot_group <- function(data, otus, x = "SA", y = "BM") {
  data %>%
    filter(OTU %in% otus) %>%
    ggplot(aes(x = .data[[x]], y = .data[[y]], colour = OTU)) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = SAMPLE_ID), size = 3, max.overlaps = 30) +
    labs(x = x, y = y) +
    theme_bw()
}

plot_LOI <- function(data, otus, x = "BM", y = "LOI") {
  data %>%
    filter(OTU %in% otus) %>%
    ggplot(aes(x = .data[[x]], y = .data[[y]], colour = OTU)) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = SAMPLE_ID), size = 3, max.overlaps = 30) +
    labs(x = x, y = y) +
    theme_bw()
}


AICc <- function(m) {
  k <- attr(logLik(m), "df"); n <- nobs(m)
  AIC(m) + 2 * k * (k + 1) / (n - k - 1)
}

### Function 4: final curve with points coloured by OTU and labelled by SAMPLE_ID (step 10)
###   chosen = your chosen model (or the ratio number for groups under 5)
###   dat    = the working data from step 5
plot_curve <- function(chosen, dat, xlab = "x", ylab = "y", title = NULL) {
  curve <- tibble(x = seq(min(dat$x), max(dat$x), length.out = 200))
  curve$y <- if (is.numeric(chosen)) chosen * curve$x else predict(chosen, newdata = curve)

  ggplot(dat, aes(x, y)) +
    geom_line(data = curve, aes(x, y), colour = "black", linewidth = 1,
              inherit.aes = FALSE) +                                    # fitted curve
    geom_point(aes(colour = OTU), size = 2) +                           # points coloured by OTU
    geom_text_repel(aes(label = SAMPLE_ID, colour = OTU), size = 3,
                    max.overlaps = 30, show.legend = FALSE) +           # sample ID labels
    labs(x = xlab, y = ylab, colour = "OTU", title = title) +
    theme_bw()
}

# 4. copy_row(): copies one phase's results to the clipboard, tab-separated (step 10)
#    chosen = chosen model (or ratio number); dat = step 5 data; aiccs = step 7 values
copy_row <- function(chosen, dat, aiccs = NULL) {
  if (is.numeric(chosen)) {                          # mean-ratio group: y = a * x
    a   <- chosen
    SEE <- sqrt(sum((dat$y - a * dat$x)^2) / (nrow(dat) - 1))
    row <- c(signif(a, 4), signif(sd(dat$y / dat$x) / sqrt(nrow(dat)), 3),
             NA, NA, NA, NA,
             signif(SEE, 3), round(100 * SEE / mean(dat$y), 1),
             NA, NA, NA,
             nrow(dat), signif(min(dat$x), 4), signif(max(dat$x), 4))
  } else {                                           # fitted model
    cf  <- coef(chosen)
    se  <- summary(chosen)$coefficients[, "Std. Error"]
    SEE <- sigma(chosen)
    d_next <- min(aiccs[aiccs != AICc(chosen)]) - AICc(chosen)
    row <- c(signif(cf[1], 4), signif(se[1], 3),     # a, a_SE
             signif(cf[2], 4), signif(se[2], 3),     # b, b_SE
             signif(cf[3], 4), signif(se[3], 3),     # c, c_SE
             signif(SEE, 3), round(100 * SEE / mean(dat$y), 1),        # SEE, SEE_pct
             round(AICc(chosen), 2), round(d_next, 2),                 # AICc, dAICc_next
             round(1 - sum(residuals(chosen)^2) /
                       sum((dat$y - mean(dat$y))^2), 3),               # R2
             nrow(dat), signif(min(dat$x), 4), signif(max(dat$x), 4))  # n, min, max
  }
  clipr::write_clip(paste(row, collapse = "\t"))
  message("Copied: ", paste(row, collapse = "  "))
}


# 3. plot_curve(): final curve, points coloured by OTU and labelled by SAMPLE_ID (step 10)
#    chosen = chosen model (or the ratio number for groups under 5); dat = step 5 data
#    band   = "SEE" (curve ± 1 SEE, default), "CI" (95% CI of the curve, lm models only) or "none"
plot_curve <- function(chosen, dat, xlab = "x", ylab = "y", title = NULL, band = "SEE") {
  curve <- tibble(x = seq(min(dat$x), max(dat$x), length.out = 200))

  if (is.numeric(chosen)) {                                   # mean ratio
    curve$y <- chosen * curve$x
    SEE     <- sqrt(sum((dat$y - chosen * dat$x)^2) / (nrow(dat) - 1))
  } else {                                                    # fitted model
    curve$y <- predict(chosen, newdata = curve)
    SEE     <- sigma(chosen)
  }

  if (band == "SEE") {                                        # curve ± 1 SEE
    curve$lo <- curve$y - SEE
    curve$hi <- curve$y + SEE
  } else if (band == "CI" && inherits(chosen, "lm")) {        # 95% CI of the curve (lm only)
    ci <- predict(chosen, newdata = curve, interval = "confidence")
    curve$lo <- ci[, "lwr"]; curve$hi <- ci[, "upr"]
  } else {
    if (band == "CI") message("CI band only works for lm models: no band drawn")
    curve$lo <- NA; curve$hi <- NA
  }

  ggplot(dat, aes(x, y)) +
    geom_ribbon(data = curve, aes(x = x, ymin = lo, ymax = hi),
                fill = "grey70", alpha = 0.4, inherit.aes = FALSE) +   # error band
    geom_line(data = curve, aes(x, y), colour = "black", linewidth = 1,
              inherit.aes = FALSE) +                                  # fitted curve
    geom_point(aes(colour = OTU), size = 2) +
    geom_text_repel(aes(label = SAMPLE_ID, colour = OTU), size = 3,
                    max.overlaps = 30, show.legend = FALSE) +
    labs(x = xlab, y = ylab, colour = "OTU", title = title,
         caption = switch(band, SEE = "Shaded: curve ± 1 SEE",
                                CI  = "Shaded: 95% confidence interval of the curve", NULL)) +
    theme_bw()
}

}

###########################################################################################################################
### OTU sets ###
###########################################################################################################################
{

alpha   <- c("POL001", "POL002", "POL007") ### Phyllodocida
beta    <- c("POL003", "POL004", "POL006", "POL008", "POL009", "POL010", "POL011") ### Sabellida
gamma   <- c("ARTH002") ### Amphipoda
delta   <- c("ARTH001") ### Pagarus bernhardus
epsilon <- c("ARTH005") ### Maja brachydactyla
zeta    <- c("ARTH003", "ARTH004", "ARTH006", "ARTH007", "ARTH008", "ARTH009", "ARTH011") ### Decapoda
eta     <- c("ARTH010", "ARTH012") ### Balanomorpha
theta   <- c("ARTH013") ### Isopoda
iota    <- c("BRY001", "BRY002", "BRY003", "BRY004", "BRY005", "BRY006", "BRY007") ### Cheilostomatida
kappa   <- c("ALG003", "ALG020") ### Chlorophyta
lambda  <- c("ASC001", "ASC003") ###Styellidae
mu      <- c("ASC004", "ASC005", "ASC006", "ASC007", "ASC008", "ASC009", "ASC010", "ASC001") ### Ascidiacea 
nu      <- c( "ANEM002") ### Actiniaria removed "ANEM001",
xi      <- c("HYD001", "HYD002", "HYD003") ### Hydrozoa
omicron <- c("ECH001") ### Marthasterias glacialis
rho     <- c("OPH001", "OPH002", "OPH003") ### Ophiuroidea
sigma   <- c("BIV001") ### Bivalvia !! RUN EACH INDIVIDUALLY 
### Run each of the following individually ### ("GAS001", "GAS002", "GAS003", "GAS005", "NUD001") Take GAS004 As face value. 
tau    <- c("GAS007", "GAS008", "GAS009") ### Nacellidae
porifera <- c("POR001", "POR002", "POR003", "POR004", "POR005", "POR006", "POR007", "POR008", "POR009")
upsilon <- c("POR004", "POR005", "POR006") ##Subertitida order
echo <- c("POR001", "POR002", "POR007") ## Encrusting globular
lima <- c("POR003") ## Elephant hide
phi     <- c("ALG001")
chi     <- c("ALG001H")
psi     <- c("ALG024")
omega   <- c("BIV002", "BIV003")
zulu    <- c("BIV004")
xray    <- c("BIV005")


PHA005 <- c()
SARG  <- c("ALG002")
RHO001 <- c("ALG006", "ALG007", "ALG008", "ALG011")
RHO002 <- c()
RHO003 <- c("ALG004", "ALG005", "ALG010", "ALG012","ALG009")
RHO004 <- c("ALG014", "ALG015", "ALG017", "ALG022", "ALG025")
RHO005 <- c("ALG023")
RHO006 <- c()
RHO007 <- c()
RHO008 <- c()

}



############################################################################################################################
           ### Working through each OTU set to produce plots and fit models for SA vs BM and BM vs LOI ###
###########################################################################################################################

###########################################################################################################################
## Alpha: Phyllodocida 
{
plot_group(AFDM_min, alpha)
# Stage 1: samples with SA and BM
AFDM_min %>% filter(OTU %in% alpha, !is.na(SA), !is.na(BM)) %>% count(OTU)

# Stage 2: samples with BM and LOI
AFDM_min %>% filter(OTU %in% alpha, !is.na(BM), !is.na(LOI)) %>% count(OTU)

ratio_dat <- AFDM_min %>% filter(OTU %in% alpha)
r <- ratio_dat$BM / ratio_dat$SA              # Stage 2: ratio_dat$LOI / ratio_dat$BM
r <- r[!is.na(r)]
mean(r)                     # the ratio:  y = ratio × x
sd(r) / sqrt(length(r))     # its standard error (for the SE column)
length(r)    
 


r <- ratio_dat$LOI / ratio_dat$BM
r <- r[!is.na(r)]
mean(r)                     # the ratio:  y = ratio × x
sd(r) / sqrt(length(r))     # its standard error (for the SE column)
length(r)      

chk   <- AFDM_min %>% filter(OTU %in% alpha, !is.na(SA), !is.na(LOI))
C_hat <- 0.5 * P2_alpha * (P1_alpha * chk$SA)
C_obs <- 0.5 * chk$LOI

C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs)          # -> C_SEE_pct (last column)


# alpha 

plot_group(AFDM_min, alpha)

## Step 1: count
AFDM_min %>% filter(OTU %in% alpha, !is.na(SA), !is.na(BM)) %>% count(OTU)    # phase 1
AFDM_min %>% filter(OTU %in% alpha, !is.na(BM), !is.na(LOI)) %>% count(OTU)   # phase 2
# n < 5 -> mean ratio for both phases

## Phase 1: SA -> BM (ratio) 
dat <- AFDM_min %>%
  filter(OTU %in% alpha, !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

P1_alpha <- mean(dat$y / dat$x)                 # store the ratio
sd(dat$y / dat$x) / sqrt(nrow(dat))             # its SE (check)

p <- plot_curve(P1_alpha, dat, "Surface area (mm²)", "Dry mass (g)", "alpha: phase 1 (ratio)")
p
ggsave("figures/alpha_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

copy_row(P1_alpha, dat)                         # type "ratio" in P1_model, paste into P1_a



## Phase 2: BM -> LOI (ratio) 
dat <- AFDM_min %>%
  filter(OTU %in% alpha, !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

P2_alpha <- mean(dat$y / dat$x)
sd(dat$y / dat$x) / sqrt(nrow(dat))

p <- plot_curve(P2_alpha, dat, "Dry mass (g)", "AFDM (g)", "alpha: phase 2 (ratio)")
p
ggsave("figures/alpha_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)

copy_row(P2_alpha, dat)                         # type "ratio" in P2_model, paste into P2_a

## Chain check -> C_SEE_pct 
chk   <- AFDM_min %>% filter(OTU %in% alpha, !is.na(SA), !is.na(LOI))
C_hat <- 0.5 * P2_alpha * (P1_alpha * chk$SA)
C_obs <- 0.5 * chk$LOI

C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs)                       # -> C_SEE_pct

}
###########################################################################################################################
### Beta: Sabellida ### 
{

AFDM_min <- AFDM_min %>% filter(SAMPLE_ID != "PB37")   # PB37 Remove SA is unworkable
plot_group(AFDM_min, beta) 
plot_group(AFDM_min, beta) + scale_x_log10() + scale_y_log10()

AFDM %>%
  filter(SAMPLE_ID == "PB37") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, Ignition_Loss, Notes)
plot_group(AFDM_min, beta) +
  facet_wrap(~ OTU, scales = "free")

pool_dat <- AFDM_min %>%
  filter(OTU %in% beta) %>%
  group_by(OTU) %>% filter(n() >= 3) %>% ungroup()

m_pooled   <- lm(log(BM) ~ log(SA),       data = pool_dat)
m_separate <- lm(log(BM) ~ log(SA) * OTU, data = pool_dat)
anova(m_pooled, m_separate)

# Stage 2: same, with BM -> LOI
m_pooled   <- lm(log(LOI) ~ log(BM),       data = pool_dat)
m_separate <- lm(log(LOI) ~ log(BM) * OTU, data = pool_dat)
anova(m_pooled, m_separate)

dat <- AFDM_min %>%
  filter(OTU %in% beta,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)     # n: decides which models you can try in step 6

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)    

# remove any model you didn't fit from this list
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)

dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))


AFDM %>%
  filter(SAMPLE_ID == "WB7") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, Ignition_Loss, Notes) ### Seems ok tbh, just very small polychaete

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)", "beta: phase 1 (power)")

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)", "GAS001: phase 1 (power)")   # ± SEE band
p
ggsave("figures/beta_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P1_beta <- chosen                 # save under its own name straight away

###Phase 2

# Stage 2: biomass on x, AFDM on y
plot_group(AFDM_min, beta, x = "BM", y = "LOI")
plot_group(AFDM_min, beta, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()

# Stage 2: biomass is x, AFDM is y
dat <- AFDM_min %>%
  filter(OTU %in% beta,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more

aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1)) ##Slight arch but will contine

dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

# Phase 2
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "beta: phase 2 (origin)", band = "SEE")
p
ggsave("figures/beta_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)

copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
copy_row(ratio, dat)              # mean-ratio group: ratio = the number from step 1

P2_beta <- chosen                 # save under its own name straight away

chk <- AFDM_min %>% filter(OTU %in% beta, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_beta, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_beta, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}

###########################################################################################################################
### Gamma: Amphipoda ###
{
plot_group(AFDM_min, gamma)
plot_group(AFDM_min, gamma) + scale_x_log10() + scale_y_log10()
AFDM %>%
  filter(SAMPLE_ID == "Pa18") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, Ignition_Loss, Notes)
dat <- AFDM_min %>%
  filter(OTU %in% gamma,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)   
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
anova(m_linear, m_quad) 
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Biomass (g)",
                "Amphipoda: phase 1 (power)", band = "SEE")
p
ggsave("figures/GAMMA_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P1_gamma <- chosen                 

##Gamma phase 2
plot_group(AFDM_min, gamma, x = "BM", y = "LOI")
plot_group(AFDM_min, gamma, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% gamma,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Biomass (g)", "AFDM (g)",
                "Amphipoda: phase 2", band = "SEE")
ggsave("figures/Amphipoda_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P2_gamma <- chosen                 # save under its own name straight away
chk <- AFDM_min %>% filter(OTU %in% gamma, !is.na(SA), !is.na(LOI))
BM_hat <- predict(P1_gamma, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_gamma, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon
1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations
}

###########################################################################################################################
### Delta: Pagarus bernhardus ##
{
plot_group(AFDM_min, delta)
plot_group(AFDM_min, delta) + scale_x_log10() + scale_y_log10()
AFDM %>%
  filter(SAMPLE_ID == "PB59") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes) ##LOI/SA looks good so keeping him in
dat <- AFDM_min %>%
  filter(OTU %in% delta,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat)     # check n dropped by 1
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1)) ## Not great but others all shit as well
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "P bernhardus: phase 1 (power)", band = "SEE")
p
ggsave("figures/Pbernhardus_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P1_delta <- chosen                 # save under its own name straight away

plot_LOI(AFDM_min, delta)+ scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% delta,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2])) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat)) 
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "P bernhardus", band = "SEE")
p
ggsave("figures/Pbern_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P2_delta <- chosen                 # save under its own name straight away

chk <- AFDM_min %>% filter(OTU %in% delta, !is.na(SA), !is.na(LOI))
BM_hat <- predict(P1_delta, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_delta, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}
###########################################################################################################################
### Epsilon: Maja brachydactyla ###
{
plot_group(AFDM_min, epsilon)
plot_group(AFDM_min, epsilon, x = "BM", y = "LOI")
## Phase 1: SA -> BM (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% epsilon, !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

P1_epsilon <- mean(dat$y / dat$x)                 # STORE the ratio (don't just print it)

p <- plot_curve(P1_epsilon, dat, "Surface area (mm²)", "Dry mass (g)", "Spider crab phase 1 (ratio)")
p
ggsave("figures/Spider crab_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

copy_row(P1_epsilon, dat)                      

## Phase 2: BM -> LOI (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% epsilon, !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

P2_epsilon <- mean(dat$y / dat$x)

p <- plot_curve(P2_epsilon, dat, "Dry mass (g)", "AFDM (g)", "spider crab: phase 2 (ratio)")
p
ggsave("figures/_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)

copy_row(P2_epsilon, dat)                         # type "ratio" in P2_model, paste into P2_a

## Chain check -> C_SEE_pct ----
chk   <- AFDM_min %>% filter(OTU %in% epsilon, !is.na(SA), !is.na(LOI))
C_hat <- 0.5 * P2_epsilon * (P1_epsilon * chk$SA)   # SA -> BM -> AFDM -> carbon
C_obs <- 0.5 * chk$LOI
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs)                       # -> C_SEE_pct (last column)

}

###########################################################################################################################
### Zeta: Decapoda ###
{
plot_group(AFDM_min, zeta) + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% zeta,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat)  
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "zeta: phase 1 (power)", band = "SEE")
p
ggsave("figures/zeta_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P1_zeta <- chosen                 # save under its own name straight away


plot_group(AFDM_min, zeta, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% zeta,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "Zeta: phase 2 (origin)", band = "SEE")
p
ggsave("figures/Zeta_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P2_zeta <- chosen  


chk <- AFDM_min %>% filter(OTU %in% zeta, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_zeta, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_zeta, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}

###########################################################################################################################
### Eta: Balanomorpha ###

AFDM_min %>% filter(OTU %in% eta, !is.na(SA), !is.na(BM)) %>% count(OTU) ## Jack shit barnacles wont count in carbon estimates


###########################################################################################################################
### Theta: Isopoda ###

AFDM_min %>% filter(OTU %in% theta, !is.na(SA), !is.na(BM)) %>% count(OTU) 
plot_group(AFDM_min, theta) ## Same as eta... no isopods.... That said very small so wont count in carbon estimates

###########################################################################################################################
### Iota: Cheilostomatida ###
{
plot_group(AFDM_min, iota) + scale_x_log10() + scale_y_log10()
plot_group(AFDM_min, iota) +
  facet_wrap(~ OTU, scales = "free")
dat <- AFDM_min %>%
  filter(OTU %in% iota,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
anova(m_linear, m_quad) 
chosen <- m_linear
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Iota: phase 1 (power)", band = "SEE")
p
ggsave("figures/Iota_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P1_iota <- chosen                 # save under its own name straight away

plot_group(AFDM_min, iota, x = "BM", y = "LOI")
plot_group(AFDM_min, iota, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()

dat <- AFDM_min %>%
  filter(OTU %in% iota,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)

p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "Iota: phase 2 (origin)", band = "SEE")
p
ggsave("figures/Iota_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P2_Iota <- chosen               # (phase 2)


chk <- AFDM_min %>% filter(OTU %in% iota, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_iota, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_Iota, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations
}
###########################################################################################################################
### Kappa: Chlorophyta ###
{
AFDM_min <- AFDM_min %>% filter(SAMPLE_ID != "WB10")     
plot_group(AFDM_min, kappa)
dat <- AFDM_min %>%
  filter(OTU %in% kappa,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)
dat <- dat %>% filter(SAMPLE_ID != "PA36")
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs      

chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
# Phase 1
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Kappa: phase 1", band = "CI")
p
ggsave("figures/Kappa_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P1_kappa <- chosen  

plot_group(AFDM_min, kappa, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% kappa,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
dat <- dat %>% filter(SAMPLE_ID != "PA36")
nrow(dat)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "kappa phase 2", band = "SEE")
p
ggsave("figures/kappa_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P2_kappa <- chosen 

chk <- AFDM_min %>% filter(OTU %in% kappa, !is.na(SA), !is.na(LOI),
                           SAMPLE_ID != "PA36")
BM_hat <- predict(P1_kappa, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_kappa, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon
1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}
###########################################################################################################################
### Lambda: Styellidae ###

AFDM_min %>% filter(OTU %in% lambda, !is.na(SA), !is.na(BM)) %>% count(OTU)
plot_group(AFDM_min, lambda) ## Only 4 samples experimenting moving ASC001 to mu DO NOT WORRY ABOUT


###########################################################################################################################
### Mu ###
{
plot_group(AFDM_min, mu)  + scale_x_log10() + scale_y_log10()
AFDM %>%
  filter(SAMPLE_ID == "PA29") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes)

plot_group(AFDM_min, mu) +
  facet_wrap(~ OTU, scales = "free")

dat <- AFDM_min %>%
  filter(OTU %in% mu,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "mu: phase 1 (power)", band = "SEE")
p
ggsave("figures/mu_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P1_mu <- chosen


plot_group(AFDM_min, mu, x = "BM", y = "LOI")
plot_group(AFDM_min, mu, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()

dat <- AFDM_min %>%
  filter(OTU %in% mu,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)     

aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
anova(m_linear, m_quad)  ## p = 4.151e-05 ***
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "mu: phase 2 (power)", band = "SEE")
p
ggsave("figures/mu_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      # fitted model (aiccs from THIS phase's step 7)
P2_mu <- chosen

chk <- AFDM_min %>% filter(OTU %in% mu, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_mu, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_mu, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)      

} 
###########################################################################################################################
### Nu ###
{
plot_group(AFDM_min, nu)

plot_group(AFDM_min, nu) + scale_x_log10() + scale_y_log10()

AFDM %>%
  filter(SAMPLE_ID == "PA28") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes)

dat <- AFDM_min %>%
  filter(OTU %in% nu,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)   
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

AFDM %>%
  filter(SAMPLE_ID == "PB33") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes)

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "nu: phase 1 (origin)", band = "SEE")
p
ggsave("figures/nu_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 

P1_nu <- chosen

plot_group(AFDM_min, nu, x = "BM", y = "LOI")
plot_group(AFDM_min, nu, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()


# Stage 2: biomass !!!!! REDO !!!!1
dat <- AFDM_min %>%
  filter(OTU %in% nu,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
anova(m_linear, m_quad)      # only if you fitted the quadratic
chosen <- m_origin

par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "nu: phase 2 (origin)", band = "SEE")
p
ggsave("figures/nu_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)   
P2_nu <- chosen

chk <- AFDM_min %>% filter(OTU %in% nu, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_nu, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_nu, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)         

}
###########################################################################################################################
### Xi ###
plot_group(AFDM_min, xi) ## Hydroids fucked so may leave
plot_group(AFDM_min, xi, x = "SA", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% xi,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
range(dat$y / dat$x)
nrow(dat)

cor.test(dat$x, dat$y, method = "spearman") ## Negative trend is actually significant.. were screwed

hyd <- AFDM_min %>%
  filter(OTU %in% xi, !is.na(SA), !is.na(LOI)) %>%
  filter(SAMPLE_ID != "WB11")              # different collection method (dropd)

r <- 0.5 * hyd$LOI / hyd$SA                # g carbon per mm², one value per sample

mean(r)                                    # carbon per mm²
sd(r) / sqrt(length(r))                    # standard error of that mean
100 * sd(r) / mean(r)                      # % variation between samples (for error on estimates)
length(r)

###########################################################################################################################
### Omicron ###
{
plot_group(AFDM_min, omicron) + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% omicron,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat)  
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Omicron: phase 1 (power)", band = "SEE")
p
ggsave("figures/Omicron_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P1_omicron <- chosen


plot_group(AFDM_min, omicron, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% omicron,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "omicron: phase 2 (origin)", band = "SEE")
p
ggsave("figures/beta_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P2_omicron <- chosen

chk <- AFDM_min %>% filter(OTU %in% omicron, !is.na(SA), !is.na(LOI))
BM_hat <- predict(P1_omicron, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_omicron, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon
1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs) 

}
###########################################################################################################################
### Rho ###
{
plot_group(AFDM_min, rho)
dat <- AFDM_min %>%
  filter(OTU %in% rho,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)   
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_linear
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1)) ## Linear Residuals are also fucked - power looks pretty good so rolling with that
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat)) ##WC24 Raised
AFDM %>%
  filter(SAMPLE_ID == "WC24") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes)
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Rho: phase 1 (power)", band = "SEE")
p
ggsave("figures/rho_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P1_rho <- chosen

dat <- AFDM_min %>%
  filter(OTU %in% rho,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "rho: phase 2 (origin)", band = "SEE")
p
ggsave("figures/rho_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      
P2_rho <- chosen

chk <- AFDM_min %>% filter(OTU %in% rho, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_rho, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_rho, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations


}

###########################################################################################################################
### Sigma ###
{
plot_group(AFDM_min, sigma) + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% sigma,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat)      

aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat)) 

AFDM %>% 
  filter(SAMPLE_ID == "PB4") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes)

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Sigma: phase 1 (power)", band = "SEE")
p
ggsave("figures/sigma_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)      
P1_sigma <- chosen

plot_group(AFDM_min, sigma, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% sigma,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2])) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "sigma: phase 2 (chosen)", band = "SEE")
p
ggsave("figures/sigma_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)     
P2_sigma <- chosen

chk <- AFDM_min %>% filter(OTU %in% sigma, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_sigma, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_sigma, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}

###########################################################################################################################
### Omega ###
{
plot_group(AFDM_min, omega)
dat <- AFDM_min %>%
  filter(OTU %in% omega,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs

anova(m_linear, m_quad) 
chosen <- m_linear
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "omega: phase 1 (power)", band = "SEE")
p
ggsave("figures/omega_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)
P1_omega <- chosen

plot_group(AFDM_min, omega, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% omega,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "omega: phase 2 (power)", band = "SEE")
p
ggsave("figures/omega_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)
P2_omega <- chosen

chk <- AFDM_min %>% filter(OTU %in% omega, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_omega, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_omega, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}
###########################################################################################################################
### Zulu ###
{
plot_group(AFDM_min, zulu)
dat <- AFDM_min %>%
  filter(OTU %in% zulu,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power ## Linear and Origin residuals bad 
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Zulu: phase 1 (power)", band = "SEE")
p
ggsave("figures/zulu_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P1_ZULU <- chosen


plot_group(AFDM_min, zulu, x = "BM", y = "LOI")
plot_group(AFDM_min, zulu, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% zulu,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "Zulu: phase 2 (origin)", band = "SEE")
p
ggsave("figures/zulu_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P2_ZULU <- chosen



chk <- AFDM_min %>% filter(OTU %in% zulu, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_ZULU, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_ZULU, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)    

}

###########################################################################################################################
### xray ###
{
plot_group(AFDM_min, xray)
dat <- AFDM_min %>%
  filter(OTU %in% xray,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)
nrow(dat)  
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b
m_quad  <- lm(y ~ x + I(x^2), data = dat) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat)) ## No samples - on the money honey
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "xray: phase 1 (power)", band = "SEE")
p
ggsave("figures/xray_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P1_xray <- chosen


plot_group(AFDM_min, xray, x = "BM", y = "LOI")+ scale_x_log10() + scale_y_log10() 
dat <- AFDM_min %>%
  filter(OTU %in% xray,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "xray: phase 2 (power)", band = "SEE")
p
ggsave("figures/xray_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P2_xray <- chosen

chk <- AFDM_min %>% filter(OTU %in% xray, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_xray, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_xray, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs) 

}

###########################################################################################################################
### Tau ### ##Removed GAS006 No examples within QUADS and clearly a different taxa size spectra
{
plot_group(AFDM_min, tau)

dat <- AFDM_min %>%
  filter(OTU %in% tau,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)

dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
AFDM %>%
  filter(SAMPLE_ID == "PC16B") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes) ## Very small sample dont see any obvious errors bar WM which also includes PC16 (the two were attached)

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "tau: phase 1 (power)", band = "SEE")
p
ggsave("figures/tau_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P1_tau <- chosen

plot_group(AFDM_min, tau, x = "BM", y = "LOI")
dat <- AFDM_min %>%
  filter(OTU %in% tau,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x
start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2])) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2) ## Curve for linear and orgin diabolical depsite lower AIC
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "tau: phase 2 (origin)", band = "SEE")
p
ggsave("figures/tau_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)   
P2_tau <- chosen


chk <- AFDM_min %>% filter(OTU %in% tau, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_tau, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_tau, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                   

}

###########################################################################################################################
### echo ###
{
plot_group(AFDM_min, echo)
dat <- AFDM_min %>%
  filter(OTU %in% echo, !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

P1_echo <- mean(dat$y / dat$x)                 # STORE the ratio (don't just print it)

p <- plot_curve(P1_echo, dat, "Surface area (mm²)", "Dry mass (g)", "echo: phase 1 (ratio)")
p
ggsave("figures/echo_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

copy_row(P1_echo, dat) 

dat <- AFDM_min %>%
  filter(OTU %in% echo, !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

P2_echo <- mean(dat$y / dat$x)

p <- plot_curve(P2_echo, dat, "Dry mass (g)", "AFDM (g)", "echo: phase 2 (ratio)")
p
ggsave("figures/echo_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)

copy_row(P2_echo, dat) 
chk   <- AFDM_min %>% filter(OTU %in% echo, !is.na(SA), !is.na(LOI))
C_hat <- 0.5 * P2_echo * (P1_echo * chk$SA)   # SA -> BM -> AFDM -> carbon
C_obs <- 0.5 * chk$LOI

C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs) 
}

###########################################################################################################################
### Upsilon ###
{
plot_group(AFDM_min, upsilon)
AFDM_min <- AFDM_min %>% filter(SAMPLE_ID != "WC18")
AFDM_min %>% filter(OTU %in% upsilon, !is.na(SA), !is.na(BM)) %>% count(OTU)
plot_group(AFDM_min, upsilon) +
  facet_wrap(~ OTU, scales = "free")

dat <- AFDM_min %>%
  filter(OTU %in% upsilon,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)


m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))

AFDM %>%
  filter(SAMPLE_ID == "WC18") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes)
dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat)) ##WC18 - Quite possible high SA so excperimenting removal from dat

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "upsilon: phase 1 (power)", band = "SEE")
p
ggsave("figures/upsilon_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P1_upsilon <- chosen

plot_group(AFDM_min, upsilon, x = "BM", y = "LOI")
plot_group(AFDM_min, upsilon, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% upsilon,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
nrow(dat)

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "upsilon: phase 2 (origin)", band = "SEE")
p
ggsave("figures/upsilon_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P2_upsilon <- chosen

ids_P1 <- dat$SAMPLE_ID        # run straight after phase 1 is finished
# ... phase 2 ...
ids_P2 <- dat$SAMPLE_ID        # run straight after phase 2 is finished

chk <- AFDM_min %>%
  filter(SAMPLE_ID %in% intersect(ids_P1, ids_P2), !is.na(SA), !is.na(LOI))

nrow(chk)      # should be 5
chk$SAMPLE_ID

chk <- AFDM_min %>% filter(OTU %in% upsilon, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_upsilon, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_upsilon, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon
1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)     

                                   # -> C_SEE_pct column in Equations

}
###########################################################################################################################
### Lima ###

{
plot_group(AFDM_min, lima)
dat <- AFDM_min %>%
  filter(OTU %in% lima,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]),
               control = nls.control(maxiter = 500))
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)        
aiccs
chosen <-  m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Lima: phase 1 (power)", band = "SEE")
p
ggsave("figures/lima_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P1_lima <- chosen

plot_group(AFDM_min, lima, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% lima,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2])) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "lima: phase 2 (power)", band = "SEE")
p
ggsave("figures/lima_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P2_lima <- chosen

chk <- AFDM_min %>% filter(OTU %in% lima, !is.na(SA), !is.na(LOI))
BM_hat <- predict(P1_lima, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_lima, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon
1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        

}
###########################################################################################################################
### Phi ###
{
plot_group(AFDM_min, phi) ### Laminaria is a complete cluster.. will require significant research to find alternatives.. 
## Phase 1: SA -> BM (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% phi, !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

P1_phi <- mean(dat$y / dat$x)                 # STORE the ratio (don't just print it)

p <- plot_curve(P1_phi, dat, "Surface area (mm²)", "Dry mass (g)", "phi: phase 1 (ratio)")
p
ggsave("figures/phi_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

copy_row(P1_phi, dat)                         # type "ratio" in P1_model, paste into P1_a

## Phase 2: BM -> LOI (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% phi, !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

P2_phi <- mean(dat$y / dat$x)

p <- plot_curve(P2_phi, dat, "Dry mass (g)", "AFDM (g)", "phi: phase 2 (ratio)")
p
ggsave("figures/phi_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)

copy_row(P2_phi, dat)                         # type "ratio" in P2_model, paste into P2_a

## Chain check -> C_SEE_pct ----
chk   <- AFDM_min %>% filter(OTU %in% phi, !is.na(SA), !is.na(LOI))
C_hat <- 0.5 * P2_phi * (P1_phi * chk$SA)   # SA -> BM -> AFDM -> carbon
C_obs <- 0.5 * chk$LOI

C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs)                       # -> C_SEE_pct (last column)
}

###########################################################################################################################
### Chi ###
{
plot_group(AFDM_min, chi)

## Phase 1: SA -> BM (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% chi, !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

P1_chi <- mean(dat$y / dat$x)                 # STORE the ratio (don't just print it)

p <- plot_curve(P1_chi, dat, "Surface area (mm²)", "Dry mass (g)", "chi: phase 1 (ratio)")
p
ggsave("figures/chi_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

copy_row(P1_chi, dat)                         # type "ratio" in P1_model, paste into P1_a

## Phase 2: BM -> LOI (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% chi, !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

P2_chi <- mean(dat$y / dat$x)

p <- plot_curve(P2_chi, dat, "Dry mass (g)", "AFDM (g)", "chi: phase 2 (ratio)")
p
ggsave("figures/chi_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)

copy_row(P2_chi, dat)                         # type "ratio" in P2_model, paste into P2_a

## Chain check -> C_SEE_pct ----
chk   <- AFDM_min %>% filter(OTU %in% chi, !is.na(SA), !is.na(LOI))
C_hat <- 0.5 * P2_chi * (P1_chi * chk$SA)   # SA -> BM -> AFDM -> carbon
C_obs <- 0.5 * chk$LOI

C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs)                       # -> C_SEE_pct (last column)
}
###########################################################################################################################
### Psi ### AKA Zostera Gonna be serious underestimate when scaled up need to id alternative scources of info - should be options online PC to carbon ideally 
{

plot_group(AFDM_min, psi)  
dat <- AFDM_min %>%
  filter(OTU %in% psi,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1)) ## Actually look good 
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "zostera: phase 1 (origin)", band = "SEE")
p
ggsave("figures/zostera_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P1_psi <- chosen

plot_group(AFDM_min, psi, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% psi,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
# Phase 2
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "zostera: phase 2 (origin)", band = "SEE")
p
ggsave("figures/zostera_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)
P2_psi <- chosen

chk <- AFDM_min %>% filter(OTU %in% psi, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_psi, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_psi, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)  

}

###########################################################################################################################
### PHA001 ###
{
PHA <- c("ALG019", "ALG016", "ALG018","ALG021") ## "ALG001","ALG001H"
plot_group(AFDM_min, PHA)
dat <- AFDM_min %>%
  filter(OTU %in% PHA,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2])) 
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2) ## Residuals awful on origin hence power chosen
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "PHA: phase 1 (power)", band = "SEE")
p
ggsave("figures/PHA_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P1_PHA <- chosen

plot_group(AFDM_min, PHA, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% PHA,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "PHA: phase 2 (power)", band = "SEE")
p
ggsave("figures/PHA_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P2_PHA <- chosen

chk <- AFDM_min %>% filter(OTU %in% PHA, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_PHA, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_PHA, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)  
}   
###########################################################################################################################
### PHA002 ### Phaeophycae corticated macrophyte
{
PHA002 <- c("ALG013", "ALG002")
plot_group(AFDM_min, PHA002)
## Phase 1: SA -> BM (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% PHA002, !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

P1_PHA002 <- mean(dat$y / dat$x)                 # STORE the ratio (don't just print it)

p <- plot_curve(P1_PHA002, dat, "Surface area (mm²)", "Dry mass (g)", "PHA002: phase 1 (ratio)")
p
ggsave("figures/PHA002a_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

copy_row(P1_PHA002, dat)                         # type "ratio" in P1_model, paste into P1_a

## Phase 2: BM -> LOI (ratio) ----
dat <- AFDM_min %>%
  filter(OTU %in% PHA002, !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

P2_PHA002 <- mean(dat$y / dat$x)

p <- plot_curve(P2_PHA002, dat, "Dry mass (g)", "AFDM (g)", "PHA002: phase 2 (ratio)")
p
ggsave("figures/PHA002_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)

copy_row(P2_PHA002, dat)                         # type "ratio" in P2_model, paste into P2_a

## Chain check -> C_SEE_pct ----
chk   <- AFDM_min %>% filter(OTU %in% PHA002, !is.na(SA), !is.na(LOI))
C_hat <- 0.5 * P2_PHA002 * (P1_PHA002 * chk$SA)   # SA -> BM -> AFDM -> carbon
C_obs <- 0.5 * chk$LOI

C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs)                       # -> C_SEE_pct (last column)
AFDM %>%
  filter(SAMPLE_ID == "ZA17") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, AFDM, Notes)
}

###########################################################################################################################
### RHO001 ### Foliose rhodophytons
{
plot_group(AFDM_min, RHO001)
dat <- AFDM_min %>%
  filter(OTU %in% RHO001,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat) 

aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
dat %>%
  mutate(std_resid = rstandard(chosen),
         cooks_D   = cooks.distance(chosen)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "Red foliose algae: phase 1 (power)", band = "SEE")
p
ggsave("figures/RHO001_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P1_RHO001 <- chosen

plot_group(AFDM_min, RHO001, x = "BM", y = "LOI")
plot_group(AFDM_min, RHO001, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()

dat <- AFDM_min %>%
  filter(OTU %in% RHO001,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat) 

aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "Red foliose: phase 2 (origin)", band = "SEE")
p
ggsave("figures/RHO001_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P2_RHO001 <- chosen

chk <- AFDM_min %>% filter(OTU %in% RHO001, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_RHO001, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_RHO001, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)    

}

###########################################################################################################################
### RHO003 ###
{
RHO003 <- c("ALG004", "ALG010","ALG009")
plot_group(AFDM_min, RHO003)
dat <- AFDM_min %>%
  filter(OTU %in% RHO003,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)   
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "RHO003: phase 1 (power)", band = "SEE")
p
ggsave("figures/RHO003_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P1_RHO003 <- chosen


plot_group(AFDM_min, RHO003, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% RHO003,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
nrow(dat)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)   
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))

p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "RHO003: phase 2 (origin)", band = "SEE")
p
ggsave("figures/RHO003_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P2_RHO003 <- chosen

chk <- AFDM_min %>% filter(OTU %in% RHO003, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_RHO003, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_RHO003, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs) 

}

###########################################################################################################################
### RHO004 ###
{
plot_group(AFDM_min, RHO004)
dat <- AFDM_min %>%
  filter(OTU %in% RHO004,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat) 
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
# Phase 1
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "RHO004: phase 1 (power)", band = "SEE")
p
ggsave("figures/RHO004_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)  
P1_RHO004 <- chosen


plot_group(AFDM_min, RHO004, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% RHO004,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_power
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "RHO004: phase 2 (power)", band = "SEE")
p
ggsave("figures/RHO004_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)
P2_RHO004 <- chosen


chk <- AFDM_min %>% filter(OTU %in% RHO004, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_RHO004, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_RHO004, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}

###########################################################################################################################
### RHO005 ### Doesnt matter none in quadrats
plot_group(AFDM_min, RHO005)


###########################################################################################################################
### GAS001 ###  
{

plot_group(AFDM_min, "GAS001")
plot_group(AFDM_min, "GAS001") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% "GAS001", 
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)     # n: decides which models you can try in step 6
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or mor

# remove any model you didn't fit from this list

aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs

anova(m_linear, m_quad)      # only if you fitted the quadratic

chosen <- m_power     # whatever you picked
plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)

dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

AFDM %>%
  filter(SAMPLE_ID == "Pa6") %>%
  select(SAMPLE_ID, OTU, Surface_Area_mm2, WET_MASS, Crucible_Mass,
         Dry_Mass, BIOMASS, Ignition_Loss, Notes)

curve   <- tibble(x = seq(min(dat$x), max(dat$x), length.out = 200))
curve$y <- predict(chosen, newdata = curve)

ggplot(dat, aes(x, y)) +
  geom_point() +
  geom_line(data = curve, colour = "red") +
  labs(x = "Surface area (mm²)", y = "Dry mass (g)") +
  theme_bw()

cf  <- coef(chosen)                                             # coefficients
se  <- summary(chosen)$coefficients[, "Std. Error"]             # their standard errors
SEE     <- sigma(chosen)                                        # standard error of the estimate (g)
SEE_pct <- 100 * SEE / mean(dat$y)                              # ... as % of the mean
d_next  <- min(aiccs[aiccs != AICc(chosen)]) - AICc(chosen)     # gap to next-best model

cat(signif(cf[1], 4), signif(se[1], 3),                         # a, a_SE
    signif(cf[2], 4), signif(se[2], 3),                         # b, b_SE
    signif(cf[3], 4), signif(se[3], 3),                         # c, c_SE
    signif(SEE, 3), round(SEE_pct, 1),                          # SEE, SEE_pct
    round(AICc(chosen), 2), round(d_next, 2),                   # AICc, dAICc_next
    round(1 - sum(residuals(chosen)^2) /
              sum((dat$y - mean(dat$y))^2), 3),                 # R2
    nrow(dat), signif(min(dat$x), 4), signif(max(dat$x), 4),    # n, min, max
    sep = "\t")

P1_GASS001 <- chosen      # Stage 1: save under its own name for step 12




### GAS001 Stage 2: LOI vs BM

plot_LOI(AFDM_min, "GAS001")               # step 2

dat <- AFDM_min %>%
  filter(OTU %in% "GAS001", 
         !is.na(SA), !is.na(BM)) %>%
  rename(x = BM, y = LOI)

# Stage 2: biomass on x, AFDM on y
plot_group(AFDM_min, "GAS001", x = "BM", y = "LOI")
plot_group(AFDM_min, "GAS001", x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()

# Stage 2: biomass is x, AFDM is y
dat <- AFDM_min %>%
  filter(OTU %in% "GAS001",
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)

m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)        # ONLY if n is 10 or more

aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs


anova(m_linear, m_quad)      # only if you fitted the quadratic
chosen <- m_power     # whatever you picked

plot(fitted(chosen), residuals(chosen)); abline(h = 0, lty = 2)

dat %>%
  mutate(std_resid = rstandard(start),
         cooks_D   = cooks.distance(start)) %>%
  filter(abs(std_resid) > 2.5 | cooks_D > 4 / nrow(dat))

curve   <- tibble(x = seq(min(dat$x), max(dat$x), length.out = 200))
curve$y <- predict(chosen, newdata = curve)

ggplot(dat, aes(x, y)) +
  geom_point() +
  geom_line(data = curve, colour = "red") +
  labs(x = "Surface area (mm²)", y = "Dry mass (g)") +   # Stage 2: "Dry mass (g)", "AFDM (g)"
  theme_bw()






# Phase 2
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)", "GAS001: phase 2 (power)")
ggsave("figures/GAS001_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)

cf  <- coef(chosen)                                             # coefficients
se  <- summary(chosen)$coefficients[, "Std. Error"]             # their standard errors
SEE     <- sigma(chosen)                                        # standard error of the estimate (g)
SEE_pct <- 100 * SEE / mean(dat$y)                              # ... as % of the mean
d_next  <- min(aiccs[aiccs != AICc(chosen)]) - AICc(chosen)     # gap to next-best model

cat(signif(cf[1], 4), signif(se[1], 3),                         # a, a_SE
    signif(cf[2], 4), signif(se[2], 3),                         # b, b_SE
    signif(cf[3], 4), signif(se[3], 3),                         # c, c_SE
    signif(SEE, 3), round(SEE_pct, 1),                          # SEE, SEE_pct
    round(AICc(chosen), 2), round(d_next, 2),                   # AICc, dAICc_next
    round(1 - sum(residuals(chosen)^2) /
              sum((dat$y - mean(dat$y))^2), 3),                 # R2
    nrow(dat), signif(min(dat$x), 4), signif(max(dat$x), 4),    # n, min, max
    sep = "\t")

P2_GAS001 <- chosen     

####
chk <- AFDM_min %>% filter(OTU %in% "GAS001", !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_GASS001, newdata = tibble(x = chk$SA))       # use your phase 1 model's name
C_hat  <- 0.5 * predict(P2_GAS001, newdata = tibble(x = BM_hat))
C_obs  <- 0.5 * chk$LOI

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)                 # R² of the chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))
100 * C_SEE / mean(C_obs)                                                # -> C_SEE_pct

}

###########################################################################################################################
### GAS002/ Orion ### If excess time incorporate with GASOO1
orion <- c("GAS002", "GAS001")
plot_group(AFDM_min, orion)



###########################################################################################################################
### icarus/ GAS003 ###
{
icarus <- c("GAS003")
plot_group(AFDM_min, "GAS003")
dat <- AFDM_min %>%
  filter(OTU %in% icarus,
         !is.na(SA), !is.na(BM)) %>%
  rename(x = SA, y = BM)

nrow(dat)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)     
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_linear
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
p <- plot_curve(chosen, dat, "Surface area (mm²)", "Dry mass (g)",
                "GAS003: phase 1 (linear)", band = "SEE")
p
ggsave("figures/GAS003_SA_vs_BM.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs) 
P1_icarus <- chosen


plot_group(AFDM_min, icarus, x = "BM", y = "LOI") + scale_x_log10() + scale_y_log10()
dat <- AFDM_min %>%
  filter(OTU %in% icarus,
         !is.na(BM), !is.na(LOI)) %>%
  rename(x = BM, y = LOI)
m_origin <- lm(y ~ 0 + x, data = dat)            # y = a·x
m_linear <- lm(y ~ x,     data = dat)            # y = a + b·x

start   <- lm(log(y) ~ log(x), data = dat)       # starting guesses for power
m_power <- nls(y ~ a * x^b, data = dat,
               start = list(a = exp(coef(start)[1]), b = coef(start)[2]))   # y = a·x^b

m_quad  <- lm(y ~ x + I(x^2), data = dat)  
aiccs <- sapply(list(origin = m_origin, linear = m_linear,
                     power = m_power, quad = m_quad), AICc)
aiccs
chosen <- m_origin
par(mfrow = c(2, 2)); plot(chosen); par(mfrow = c(1, 1))
p <- plot_curve(chosen, dat, "Dry mass (g)", "AFDM (g)",
                "GAS003: phase 2 (origin)", band = "SEE")
p
ggsave("figures/GAS003_BM_vs_LOI.png", p, width = 7, height = 5, dpi = 300)
copy_row(chosen, dat, aiccs)   
P2_icarus <- chosen


chk <- AFDM_min %>% filter(OTU %in% icarus, !is.na(SA), !is.na(LOI))

BM_hat <- predict(P1_icarus, newdata = tibble(x = chk$SA))       # stage 1
C_hat  <- 0.5 * predict(P2_icarus, newdata = tibble(x = BM_hat))  # stage 2, then x 0.5
C_obs  <- 0.5 * chk$LOI                                          # measured carbon

1 - sum((C_obs - C_hat)^2) / sum((C_obs - mean(C_obs))^2)        # R² of the full chain
C_SEE <- sqrt(sum((C_obs - C_hat)^2) / (length(C_obs) - 2))     # SEE of the full chain (g C)
100 * C_SEE / mean(C_obs)                                        # -> C_SEE_pct column in Equations

}

###########################################################################################################################
### NUD001 ###
plot_group(AFDM_min, "NUD001")##So minor dont includ in carbon

###########################################################################################################################
### GAS005 ### Low priority
plot_group(AFDM_min, "GAS004") ##So minor dont includ in carbon

    