###########################################################################################################################
                    ### PHASE FOXTROT; Quantative assessment of Carbon across le Sites ###
###########################################################################################################################

###########################################################################################################################
### PF1 Biblioteca
library(tidyverse)
library(ggplot2)
library(vegan) ## Assuming gonna be some PEAR analysis
library(ggrepel)
library(rnaturalearth) ## For plotting le sexy charts - fuck using QGIS



###########################################################################################################################
### PF2 Data insertion

setwd('/Users/hamishmckee/Documents/BAS Paper SCIDIV5000/DATA/Datasheets/AFDM')

### QUADS - SA of Quadrats - each individual SA within has corresponding OTU code
QUADS <- read.csv("20261118a_QUADS.csv", header = TRUE, na = c("", "NA", "N/A", "n/a", "#VALUE!", "#REF!"))
str(QUADS)
head(QUADS)

### OTU - list of OTUs and most importantly there corresponding relay code
OTUs <- read.csv("20261119_OTU.csv", header = TRUE, na = c("", "NA", "N/A", "n/a", "#VALUE!", "#REF!"))
head(OTUs)
str(OTUs)

OTUs %>% count(OTU) %>% filter(n > 1)
OTUs <- OTUs %>% filter(!is.na(OTU), OTU != "")
### RELAY -  Relay codes, model, coefficients to curve, SE across all stages
relay <- read.csv("20261118a_RELAY.csv", header = TRUE, na = c("", "NA", "N/A", "n/a", "#VALUE!", "#REF!"))
head(relay)


### Meta - Depth, dates and LOSTAT Data
META <- read.csv("20261119_META.csv", header = TRUE)
META <- META %>% mutate(BRG = c(340, 0, 330, 180, 90, 280, 270, 100, 130, 90))
## Forgot to add bearings will add in next iteration of meta
head(META)


###########################################################################################################################
### PF3 Execution of relay process

### QUADS insertion of Relay group

OTUs <- OTUs %>% filter(!is.na(OTU))

quadsR <- QUADS %>%
  left_join(OTUs %>% select(OTU, Relay_Group), by = "OTU", relationship = "many-to-one")

nrow(QUADS) == nrow(quadsR)      #  Is TRUE

###Randomly remove a PA transect as we have 8 not 7
set.seed(1664)                                          
pa_quads   <- unique(QUADS$Quadrat_ID[QUADS$Transect_ID == "PA"])
drop_quad  <- sample(pa_quads, 1)
drop_quad                                                 # note this in your methods
QUADS <- QUADS %>% filter(!(Transect_ID == "PA" & Quadrat_ID == drop_quad))
QUADS %>% distinct(Transect_ID, Quadrat_ID) %>% count(Transect_ID) ## Ensure 7 Quadrats across each transect


head(quadsR)
nrow(QUADS)
nrow(quadsR) ## Same number of rows so no funny buisness
quadsR %>% filter(is.na(Relay_Group)) %>% count(OTU) ## See taxa withouy a relay group


### ID Mismatches  See whch taxa from above dont have a relay group be it NA
quadsR %>% distinct(Relay_Group) %>%
  anti_join(relay, by = c("Relay_Group" = "Relay"))

###########################################################################################################################
### PF4 Biomass and carbon estimates

### Apply any relay equation:  origin/ratio, linear, power, quadratic
predict_relay <- function(model, a, b, c, x) {
  model <- tolower(model)                                   # "Power" = "power"
  case_when(
    model %in% c("origin", "ratio")   ~ a * x,
    model == "linear"                 ~ a + b * x,
    model == "power"                  ~ a * x^b,
    model %in% c("quad", "quadratic") ~ a + b * x + c * x^2,
    TRUE                              ~ NA_real_)
}

C_FACTOR <- 0.5  

quad_C <- quadsR %>%
  left_join(relay, by = c("Relay_Group" = "Relay")) %>%     # relay's group column on the right
  mutate(
    BM     = predict_relay(P1_model, P1_a, P1_b, P1_c, Area),   # phase 1: SA -> biomass
    AFDM   = predict_relay(P2_model, P2_a, P2_b, P2_c, BM),     # phase 2: biomass -> AFDM
    C      = C_FACTOR * AFDM,                                    # carbon
    BM_err = BM * P1_SEE_pct / 100,                              # ± on biomass
    C_pct  = coalesce(C_SEE_pct, sqrt(P1_SEE_pct^2 + P2_SEE_pct^2)),  # chain % error (fallback if missing)
    C_err  = C * C_pct / 100,                                    # ± on carbon
    out_of_range = P1_n > 1 & (Area < P1_minSA | Area > P1_maxSA))

head(quad_C)

### Temporary removal of Omega and GASOO3 negys will change for paper analysis or for dis if time
quad_C <- quad_C %>%
  mutate(neg_pred = BM < 0 | AFDM < 0,
         across(c(BM, AFDM, C, BM_err, C_err), ~ if_else(neg_pred %in% TRUE, NA_real_, .x)))

quad_C %>% filter(BM < 0 | C < 0) %>% nrow()   # should be 0 no neggys 




quad_summary <- quad_C %>%
  group_by(Transect_ID, Quadrat_ID) %>%                                    # your column names
  summarise(
    total_BM      = sum(BM, na.rm = TRUE),
    BM_err        = sqrt(sum(BM_err^2, na.rm = TRUE)),
    total_C       = sum(C, na.rm = TRUE),
    C_err         = sqrt(sum(C_err^2, na.rm = TRUE)),
    total_SA      = sum(Area, na.rm = TRUE),
    pct_converted = 100 * sum(Area[!is.na(C)]) / total_SA,      # % of area with a carbon estimate
    n_extrap      = sum(out_of_range, na.rm = TRUE),
    .groups = "drop")

quad_summary
nrow(quad_summary)

quad_C %>% filter(!is.na(Relay_Group), is.na(C)) %>% count(Relay_Group)   # groups that should have worked but didn't
quad_C %>% filter(BM < 0 | C < 0)                                         # impossible negatives
quad_C %>% filter(out_of_range) %>% count(Relay_Group)   

write_csv(quad_summary, "quad_summary.csv")




###########################################################################################################################
### PF5 Initial visualisation of data


quad_tot <- quad_C %>%
  group_by(Transect_ID, Quadrat_ID) %>%
  summarise(total_C = sum(C, na.rm = TRUE),
            C_err   = sqrt(sum(C_err^2, na.rm = TRUE)),
            .groups = "drop")

### Transect mean ± SE across quadrats
trans_tot <- quad_tot %>%
  group_by(Transect_ID) %>%
  summarise(mean_C = mean(total_C), se_C = sd(total_C) / sqrt(n()), n = n())
trans_tot

### Plot
p <- ggplot(quad_tot, aes(Transect_ID, total_C)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6, colour = "grey40") +
  geom_errorbar(data = trans_tot, aes(y = mean_C, ymin = mean_C - se_C, ymax = mean_C + se_C),
                width = 0.3, linewidth = 0.8) +
  geom_point(data = trans_tot, aes(y = mean_C), size = 4, shape = 18) +
  labs(x = "Transect", y = "Carbon per quadrat (g C)",
       caption = "Grey points: individual quadrats. Diamond and bars: transect mean ± SE.") +
  theme_bw()
p
ggsave("figures/carbon_per_transect.png", p, width = 7, height = 5, dpi = 300)


quad_tot %>%
  arrange(desc(total_C)) %>%
  print(n = Inf)

quad_C %>%
  filter(Quadrat_ID == "Q2_2197") %>%                  # the outlier quadrat
  select(Organism_ID, OTU, Relay_Group, Area, BM, C, out_of_range) %>%
  arrange(desc(C)) %>%
  head(10)
quad_C %>% filter(Quadrat_ID == "Q2_2197") %>% summarise(total_SA = sum(Area, na.rm = TRUE)) ### 5.1M^2 Probably an extra zero address and resolve this evening


###########################################################################################################################
### PF6 Alpha diversity assessment across le sites

com_test_ <- adonis2(community[,-c(1:2)] ~ Depth, data = meta, method = "bray") # adonis2 - PERMANOVA Depth as predictor variable  ∴ effect of depth on community composition
com_test_xdepth1 <- adonis(community[,-c(1:2)] ~ Depth, data = meta, method = "bray")


PEAR <- QUADS %>%
  filter(!is.na(OTU), OTU != "", OTU != "CHECK") %>%
  mutate(Site = case_when(substr(Transect_ID, 1, 1) == "P" ~ "Pilot",
                          substr(Transect_ID, 1, 1) == "W" ~ "Westlake",
                          substr(Transect_ID, 1, 1) == "Z" ~ "Zostera"))

head(PEAR)

comm <- PEAR %>%
  group_by(Site, Transect_ID, Quadrat_ID, OTU) %>%
  summarise(abund = sum(Area, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = OTU, values_from = abund, values_fill = 0)

comm_mat <- comm %>% select(-Site, -Transect_ID, -Quadrat_ID) %>% as.matrix()

shan <- comm %>%
  select(Site, Transect_ID, Quadrat_ID) %>%
  mutate(richness = specnumber(comm_mat),                       # number of OTUs
         shannon  = diversity(comm_mat, index = "shannon"),     # H'
         evenness = shannon / log(richness))                    

shan %>% group_by(Site) %>%
  summarise(n_quadrats = n(),
            richness = mean(richness), shannon_mean = mean(shannon),
            shannon_se = sd(shannon) / sqrt(n()), evenness = mean(evenness, na.rm = TRUE))

p <- ggplot(shan, aes(Site, shannon)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(aes(colour = Transect_ID), width = 0.15, size = 2, alpha = 0.7) +
  labs(x = "Site", y = "Shannon diversity (H′)", colour = "Transect") +
  theme_bw()
p

kruskal.test(shannon ~ Site, data = shan)                       # non-parametric, small n
pairwise.wilcox.test(shan$shannon, shan$Site, p.adjust.method = "holm")


### SWAP for GLM with transect as a random effect 



###########################################################################################################################
### PF7 Beta diversity (Bray–Curtis), NMDS and PERMANOVA (PHASE 2 of Op PEAR)

###  community matrix
keep     <- rowSums(comm_mat) > 0                       # drop any empty quadrats
meta     <- comm[keep, c("Site", "Transect_ID", "Quadrat_ID")]
comm_rel <- decostand(comm_mat[keep, ], method = "total")   # relative cover (each quadrat sums to 1)
comm_tr  <- sqrt(comm_rel)                                   # down-weight dominant OTUs

### Bray–Curtis dissimilarity
bc <- vegdist(comm_tr, method = "bray")

###  NMDS
set.seed(1664)
nmds <- metaMDS(comm_tr, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
nmds$stress                                    # < 0.2 usable, < 0.1 good
stressplot(nmds)                                    

### Plot LE NMDS
nmds_pts <- bind_cols(meta, as_tibble(scores(nmds, display = "sites")))

p <- ggplot(nmds_pts, aes(NMDS1, NMDS2, colour = Site)) +
  stat_ellipse(level = 0.95, linetype = 2) +
  geom_point(aes(shape = Transect_ID), size = 2.5) +
  scale_shape_manual(values = 0:12) +
  labs(caption = paste0("Bray–Curtis on sqrt relative cover; 2D stress = ",
                        round(nmds$stress, 3))) +
  coord_equal() +
  theme_bw()
p
ggsave("figures/nmds_bray_sites.png", p, width = 7, height = 6, dpi = 300)

### PERMANOVA
set.seed(1664)
adonis2(bc ~ Site, data = meta, permutations = 999)

### 6. Check the assumption: similar spread within each site (dispersion)
disp <- betadisper(bc, meta$Site)
anova(disp)                                              # p > 0.05 = spreads similar, PERMANOVA reliable
boxplot(disp)

### 7. Pairwise PERMANOVA (which sites differ)
pairs <- combn(unique(meta$Site), 2, simplify = FALSE)
pw <- map_dfr(pairs, function(s) {
  idx <- meta$Site %in% s
  res <- adonis2(vegdist(comm_tr[idx, ], "bray") ~ Site, data = meta[idx, ], permutations = 999)
  tibble(pair = paste(s, collapse = " vs "), R2 = res$R2[1], F = res$F[1], p = res$`Pr(>F)`[1])
})
pw %>% mutate(p_holm = p.adjust(p, "holm"))




###########################################################################################################################
### PF8 Site map ### Needs serious remedial work may infact need to use le QGIS

# install.packages(c("sf", "rnaturalearth", "rnaturalearthdata", "ggspatial"))
# install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")   # high-res coastline
library(sf)
library(rnaturalearth)
library(ggspatial)

### 1. Degrees + decimal minutes ("50°21,334") -> decimal degrees
dm_to_dd <- function(x) {
  x   <- gsub(",", ".", x)                       # decimal comma -> point
  deg <- as.numeric(sub("°.*", "", x))
  mins <- as.numeric(sub(".*°", "", x))
  deg + mins / 60
}

META <- META %>%
  mutate(lat = dm_to_dd(LAT),
         lon = -dm_to_dd(LONG))                  # west = negative

### 2. Transect direction: short line from start along the bearing (length is for display)
L <- 100                                          # metres; set to your transect length if known
META <- META %>%
  mutate(lat_end = lat + (L * cos(BRG * pi / 180)) / 111320,
         lon_end = lon + (L * sin(BRG * pi / 180)) / (111320 * cos(lat * pi / 180)))

### 3. Coastline
coast <- ne_countries(scale = "large", country = "United Kingdom", returnclass = "sf")

### 4. Map
pad <- 0.01
p <- ggplot() +
  geom_sf(data = coast, fill = "grey85", colour = "grey40") +
  geom_segment(data = META, aes(x = lon, y = lat, xend = lon_end, yend = lat_end, colour = SITE),
               arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.8) +
  geom_point(data = META, aes(lon, lat, colour = SITE), size = 2.5) +
  geom_text(data = META, aes(lon, lat, label = TRANSECT), size = 3, nudge_y = 0.0012) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", style = north_arrow_minimal) +
  coord_sf(xlim = range(META$lon) + c(-pad, pad),
           ylim = range(META$lat) + c(-pad, pad), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude", colour = "Site") +
  theme_bw()
p
ggsave("figures/site_map.png", p, width = 7, height = 7, dpi = 300)



###########################################################################################################################
### PF9 Functional group diversity


###########################################################################################################################
### PF9 Functional group diversity

### 1. Attach functional group to each organism
fg_dat <- QUADS %>%
  left_join(OTUs %>% select(OTU, FG_Code, Functional_Group), by = "OTU", relationship = "many-to-one") %>%
  filter(!is.na(FG_Code), FG_Code != "") %>%
  mutate(Site = case_when(substr(Transect_ID, 1, 1) == "P" ~ "Pilot",
                          substr(Transect_ID, 1, 1) == "W" ~ "Westlake",
                          substr(Transect_ID, 1, 1) == "Z" ~ "Zostera"))

QUADS %>% anti_join(OTUs %>% filter(!is.na(FG_Code)), by = "OTU") %>% count(OTU)   # OTUs with no FG

### 2. Quadrat x FG matrix (surface area as abundance)
fg_comm <- fg_dat %>%
  group_by(Site, Transect_ID, Quadrat_ID, FG_Code) %>%
  summarise(abund = sum(Area, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = FG_Code, values_from = abund, values_fill = 0)

fg_mat  <- fg_comm %>% select(-Site, -Transect_ID, -Quadrat_ID) %>% as.matrix()
fg_meta <- fg_comm %>% select(Site, Transect_ID, Quadrat_ID)

### 3. FG alpha diversity per quadrat
fg_alpha <- fg_meta %>%
  mutate(fg_richness = specnumber(fg_mat),
         fg_shannon  = diversity(fg_mat, index = "shannon"),
         fg_evenness = fg_shannon / log(fg_richness))

fg_alpha %>% group_by(Site) %>%
  summarise(n = n(), richness = mean(fg_richness),
            shannon = mean(fg_shannon), shannon_se = sd(fg_shannon) / sqrt(n()))

kruskal.test(fg_shannon ~ Site, data = fg_alpha)
pairwise.wilcox.test(fg_alpha$fg_shannon, fg_alpha$Site, p.adjust.method = "holm")
kruskal.test(fg_richness ~ Site, data = fg_alpha)

p1 <- ggplot(fg_alpha, aes(Site, fg_shannon)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(aes(colour = Transect_ID), width = 0.15, size = 2, alpha = 0.7) +
  labs(y = "Functional group Shannon diversity (H′)", colour = "Transect") +
  theme_bw()
p1
ggsave("figures/fg_shannon_per_site.png", p1, width = 7, height = 5, dpi = 300)

### 4. FG composition: relative cover per site (stacked bar)
fg_cover <- fg_dat %>%
  group_by(Site, Functional_Group) %>%
  summarise(area = sum(Area, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site) %>% mutate(pct = 100 * area / sum(area))

p2 <- ggplot(fg_cover, aes(Site, pct, fill = Functional_Group)) +
  geom_col() +
  labs(y = "Relative cover (%)", fill = "Functional group") +
  theme_bw()
p2
ggsave("figures/fg_cover_per_site.png", p2, width = 8, height = 5, dpi = 300)

### 5. FG composition differences: Bray–Curtis, PERMANOVA, dispersion, NMDS
keep   <- rowSums(fg_mat) > 0
fg_tr  <- sqrt(decostand(fg_mat[keep, ], "total"))
fg_bc  <- vegdist(fg_tr, "bray")
fg_m   <- fg_meta[keep, ]

set.seed(1664)
adonis2(fg_bc ~ Site, data = fg_m, permutations = 999)
anova(betadisper(fg_bc, fg_m$Site))

set.seed(1664)
fg_nmds <- metaMDS(fg_tr, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
fg_nmds$stress

fg_pts <- bind_cols(fg_m, as_tibble(scores(fg_nmds, display = "sites")))
fg_sp  <- as_tibble(scores(fg_nmds, display = "species"), rownames = "FG_Code")   # FG positions

p3 <- ggplot(fg_pts, aes(NMDS1, NMDS2)) +
  stat_ellipse(aes(colour = Site), linetype = 2) +
  geom_point(aes(colour = Site), size = 2.5) +
  geom_text(data = fg_sp, aes(label = FG_Code), size = 3, colour = "grey30") +   # which FGs drive separation
  labs(caption = paste0("FG Bray–Curtis on sqrt relative cover; stress = ", round(fg_nmds$stress, 3))) +
  coord_equal() + theme_bw()
p3
ggsave("figures/fg_nmds.png", p3, width = 7, height = 6, dpi = 300)


###########################################################################################################################
### PF10 What explains variation in carbon?

### 1. One row per quadrat: carbon + predictors
SF_CODES <- c("SHSF")      # <- FG codes that are SUSPENSION FEEDERS; edit to match your FG scheme

c_quad <- quad_C %>%
  group_by(Transect_ID, Quadrat_ID) %>%
  summarise(total_C       = sum(C, na.rm = TRUE),
            pct_converted = 100 * sum(Area[!is.na(C)]) / sum(Area, na.rm = TRUE),
            .groups = "drop")

sf_quad <- fg_dat %>%
  group_by(Transect_ID, Quadrat_ID) %>%
  summarise(prop_SF = sum(Area[FG_Code %in% SF_CODES], na.rm = TRUE) / sum(Area, na.rm = TRUE),
            total_SA = sum(Area, na.rm = TRUE), .groups = "drop")

mod_dat <- fg_alpha %>%
  left_join(shan %>% select(Quadrat_ID, richness, shannon), by = "Quadrat_ID") %>%
  left_join(c_quad,  by = c("Transect_ID", "Quadrat_ID")) %>%
  left_join(sf_quad, by = c("Transect_ID", "Quadrat_ID")) %>%
  left_join(META %>% select(TRANSECT, Avg_depth), by = c("Transect_ID" = "TRANSECT")) %>%
  mutate(log_C = log10(total_C + 0.001))          # log as in Barnes & Sands; +0.001 avoids log(0)

summary(mod_dat)                                  # check for NAs after joins

### 2. Does carbon differ between sites? (ANOVA, as in the paper)
m_site <- lm(log_C ~ Site, data = mod_dat)
anova(m_site); summary(m_site)$r.squared
TukeyHSD(aov(log_C ~ Site, data = mod_dat))

### 3. Each predictor on its own: R² and p
preds <- c("fg_richness", "fg_shannon", "richness", "shannon", "prop_SF", "Avg_depth", "total_SA")
single <- map_dfr(preds, function(v) {
  m <- lm(reformulate(v, "log_C"), data = mod_dat)
  tibble(predictor = v, slope = coef(m)[2], R2 = summary(m)$r.squared,
         p = anova(m)$`Pr(>F)`[1])
}) %>% arrange(desc(R2))
single

### 4. FGD vs site: how much does each explain, and how much is shared?
vp <- varpart(mod_dat$log_C, ~ fg_shannon, ~ Site, data = mod_dat)
vp; plot(vp, Xnames = c("FG diversity", "Site"))

### 5. Does FGD explain carbon beyond site differences?
m_full <- lm(log_C ~ Site + fg_shannon, data = mod_dat)
anova(m_site, m_full)                             # p < 0.05 -> FGD adds explanatory power
summary(m_full)

### 6. Check assumptions
par(mfrow = c(2, 2)); plot(m_full); par(mfrow = c(1, 1))

### 7. Plot: carbon vs FG diversity, by site
p <- ggplot(mod_dat, aes(fg_shannon, total_C, colour = Site)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(group = 1), method = "lm", colour = "black", se = TRUE) +
  scale_y_log10() +
  labs(x = "Functional group diversity (H′)", y = "Carbon per quadrat (g C, log scale)") +
  theme_bw()
p
ggsave("figures/carbon_vs_fgd.png", p, width = 7, height = 5, dpi = 300)