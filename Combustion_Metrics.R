#=============================================================================#
# Title: Plot Level Combustion Data
# Created by: Zach Madsen
# Created for: Paige Paulsen
# Date Created: 05/05/2026
#=============================================================================#

#==================================#
# Libraries----
#==================================#
library(tidyverse)
library(fs)

#==================================#
# Data----
#==================================#
cce.soil   <- read.csv('Field_Data/Analysis_Ready/CCE_soils_Master.csv')
cce.site   <- read_csv("Field_Data/Analysis_Ready/CCE_site_Master.csv")
cce.trees  <- read_csv("Field_Data/Analysis_Ready/CCE_combustion_Master.csv")
cce.shrubs <- read_csv("Field_Data/Analysis_Ready/CCE_shrub_Master.csv")
cce.cwd    <- read_csv("Field_Data/Analysis_Ready/CCE_cwd_Master.csv")

#=============================================================================#
# 1) Compile Soil Combustion Data----
#=============================================================================#

# Calculate Adventitious Root Height Mean
cce.soil <- cce.soil %>%
  mutate(
    ar_height.mean = case_when(
      !is.na(ar_height1) | !is.na(ar_height2) | !is.na(ar_height3) ~
        rowMeans(across(c(ar_height1, ar_height2, ar_height3)), na.rm = TRUE),
      !is.na(ar_height_avg) ~ ar_height_avg,
      TRUE ~ NA_real_
    )
  )

# Calculate burn depth and proportion of soil depth
cce.soil <- cce.soil %>%
  mutate(
    burn_depth_calc = round(ar_height.mean + 3.4, 2),
    prop_soil_depth = round(burn_depth_calc / (depth1 + burn_depth_calc), 3)
  )

# Summarize data at plot level
cce.plot.soil.sum <- cce.soil %>%
  group_by(fire_scar, site, plot) %>%
  summarise(
    depth1 = mean(depth1, na.rm = TRUE) %>% round(2),
    burn_depth_calc = mean(burn_depth_calc, na.rm = TRUE) %>% round(2),
    prop_soil_depth = mean(prop_soil_depth, na.rm = TRUE) %>% round(3),
    .groups = "drop"
  )

#=============================================================================#
# 2) Trees + Snags----
#=============================================================================#
# Plot area = 2m x 10m = 20 m2
PLOT_AREA_TREES <- 20

# Expand rows where count > 1 (groups of identical trees recorded once)
cce.trees.0 <- cce.trees %>% filter(is.na(count) | count == 0)
cce.trees.exp <- cce.trees %>%
  filter(!is.na(count) & count > 0) %>%
  { data.frame(lapply(., rep, .$count)) }
cce.trees <- bind_rows(
  cce.trees %>% filter(is.na(count)),
  cce.trees.exp,
  cce.trees.0
)


# Basal area per tree (cm2) — use 0.1 cm DBH if only BD measured
cce.trees <- cce.trees %>%
  mutate(basal.area.dbh = if_else(!is.na(dbh), pi * (dbh / 2)^2, pi * (0.1 / 2)^2))

#--------------------------------------------------------------------#
## 2a) Relative basal area by site (for unknown species allometry)----
#--------------------------------------------------------------------#
rel.ba.site <- cce.trees %>%
  filter(!tree_species %in% c("unknown", "unknown conifer", "unknown deciduous")) %>%
  group_by(tree_species, fire_scar, site) %>%
  summarise(ba = sum(basal.area.dbh) / 60, .groups = "drop") %>%
  group_by(fire_scar, site) %>%
  mutate(site.ba = sum(ba, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    rel.ba.as = if_else(tree_species == "aspen",        ba / site.ba, NA_real_),
    rel.ba.br = if_else(tree_species == "birch",        ba / site.ba, NA_real_),
    rel.ba.bs = if_else(tree_species == "black spruce", ba / site.ba, NA_real_),
    rel.ba.bp = if_else(tree_species == "poplar",       ba / site.ba, NA_real_),
    rel.ba.ws = if_else(tree_species == "white spruce", ba / site.ba, NA_real_)
  )

# Pivot to wide and fill NAs with 0
rel.ba.wide <- rel.ba.site %>%
  select(fire_scar, site, rel.ba.as, rel.ba.br, rel.ba.bs, rel.ba.bp, rel.ba.ws) %>%
  pivot_longer(cols = starts_with("rel.ba"), names_to = "sp", values_to = "val") %>%
  filter(!is.na(val)) %>%
  pivot_wider(names_from = sp, values_from = val, values_fill = 0)

# Conifer-only relative BA (for unknown conifer)
rel.ba.con <- cce.trees %>%
  filter(tree_species %in% c("black spruce", "white spruce")) %>%
  group_by(tree_species, fire_scar, site) %>%
  summarise(ba = sum(basal.area.dbh) / 60, .groups = "drop") %>%
  group_by(fire_scar, site) %>%
  mutate(con.site.ba = sum(ba, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    rel.ba.bs.c = if_else(tree_species == "black spruce", ba / con.site.ba, NA_real_),
    rel.ba.ws.c = if_else(tree_species == "white spruce", ba / con.site.ba, NA_real_)
  ) %>%
  select(fire_scar, site, rel.ba.bs.c, rel.ba.ws.c) %>%
  pivot_longer(cols = starts_with("rel.ba"), names_to = "sp", values_to = "val") %>%
  filter(!is.na(val)) %>%
  pivot_wider(names_from = sp, values_from = val, values_fill = 0)

# Join relative BA back to trees
# any_of() used because not all species may be present in every dataset
cce.trees <- cce.trees %>%
  left_join(rel.ba.wide, by = c("fire_scar", "site")) %>%
  left_join(rel.ba.con,  by = c("fire_scar", "site")) %>%
  mutate(across(any_of(c("rel.ba.as", "rel.ba.br", "rel.ba.bs", "rel.ba.bp", "rel.ba.ws",
                         "rel.ba.bs.c", "rel.ba.ws.c")), ~ replace_na(., 0)))

#--------------------------------------------------------------------#
## 2b) Biomass per tree (g/tree) — allometric equations----
#--------------------------------------------------------------------#

# -- BLACK SPRUCE (Boby et al. 2010) --
cce.trees <- cce.trees %>%
  mutate(
    boby.needles.bs = case_when(
      tree_species == "black spruce" & !is.na(dbh) & dbh >= 2.7 ~ (1.9*(dbh^2)+16.0) + (16.0*(dbh^2)+288.4),
      tree_species == "black spruce" & !is.na(dbh) & dbh < 2.7  ~ 3.9  * log(dbh),
      tree_species == "black spruce" & is.na(dbh)               ~ 18.3 * (basal_diameter^2),
      TRUE ~ NA_real_
    ),
    boby.fine.bs = case_when(
      tree_species == "black spruce" & !is.na(dbh) & dbh >= 2.7 ~ (11.3*(dbh^2)+81.2) + (0.5*(dbh^2)+1.3),
      tree_species == "black spruce" & !is.na(dbh) & dbh < 2.7  ~ 11.3 * (dbh^2),
      tree_species == "black spruce" & is.na(dbh)               ~ 7.5  * (basal_diameter^2),
      TRUE ~ NA_real_
    ),
    boby.coarse.bs = case_when(
      tree_species == "black spruce" & !is.na(dbh) & dbh >= 2.7 ~ 11.8 * (dbh^2) - 90.6,
      tree_species == "black spruce" & !is.na(dbh) & dbh < 2.7  ~ 14.4 * (dbh^2),
      tree_species == "black spruce" & is.na(dbh)               ~ 3.1  * (basal_diameter^2),
      TRUE ~ NA_real_
    ),
    # Stem biomass — Andrew's coefficients (Boby et al. 2010)
    stemwood.bark.bs = case_when(
      tree_species == "black spruce" & !is.na(dbh) & dbh >= 2.7 ~ 117.7 * (dbh^2) - 253.9,
      tree_species == "black spruce" & !is.na(dbh) & dbh < 2.7  ~ 113.5 * (dbh^2),
      tree_species == "black spruce" & is.na(dbh)               ~ 12.0  * (basal_diameter^2),
      TRUE ~ NA_real_
    )
  ) %>%
  # Total BS biomass = sum of components (per Andrew)
  rowwise() %>%
  mutate(
    boby.total.biomass.bs = if_else(
      tree_species == "black spruce",
      sum(c(boby.needles.bs, boby.fine.bs, boby.coarse.bs, stemwood.bark.bs), na.rm = TRUE),
      NA_real_
    )
  ) %>%
  ungroup()

# -- ASPEN (Alexander et al. 2012) --
cce.trees <- cce.trees %>%
  mutate(
    total.biomass.as = case_when(
      tree_species == "aspen" & !is.na(dbh) ~ 134.1  * (dbh^2.26),
      tree_species == "aspen" & is.na(dbh)  ~ 56.83  * (basal_diameter^2.49),
      TRUE ~ NA_real_
    ),
    foliage.as = case_when(
      tree_species == "aspen" & !is.na(dbh) ~ 18.98  * (dbh^1.53),
      tree_species == "aspen" & is.na(dbh)  ~ 8.73   * (basal_diameter^2),
      TRUE ~ NA_real_
    ),
    live.crown.as = case_when(
      tree_species == "aspen" & !is.na(dbh) ~ 41.74  * (dbh^1.83),
      tree_species == "aspen" & is.na(dbh)  ~ 18.16  * (basal_diameter^2.43),
      TRUE ~ NA_real_
    ),
    stemwood.bark.as = case_when(
      tree_species == "aspen" & !is.na(dbh) ~ 64.01  * (dbh^2.51),
      tree_species == "aspen" & is.na(dbh)  ~ 12.55  * (basal_diameter^2.98),
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    branches.as = if_else(tree_species == "aspen", live.crown.as - foliage.as, NA_real_),
    fine.as     = if_else(tree_species == "aspen", branches.as * 0.15, NA_real_),
    coarse.as   = if_else(tree_species == "aspen", branches.as * 0.85, NA_real_)
  )

# -- BIRCH (Alexander et al. 2012) --
cce.trees <- cce.trees %>%
  mutate(
    total.biomass.br = case_when(
      tree_species == "birch" & !is.na(dbh) ~ 164.18 * (dbh^2.29),
      tree_species == "birch" & is.na(dbh)  ~ 26.29  * (basal_diameter^2.68),
      TRUE ~ NA_real_
    ),
    foliage.br = case_when(
      tree_species == "birch" & !is.na(dbh) ~ 6.39   * (dbh^2.1),
      tree_species == "birch" & is.na(dbh)  ~ 0.63   * (basal_diameter^2.65),
      TRUE ~ NA_real_
    ),
    live.crown.br = case_when(
      tree_species == "birch" & !is.na(dbh) ~ 15.15  * (dbh^2.49),
      tree_species == "birch" & is.na(dbh)  ~ 1.38   * (basal_diameter^3.08),
      TRUE ~ NA_real_
    ),
    stemwood.bark.br = case_when(
      tree_species == "birch" & !is.na(dbh) ~ 147.96 * (dbh^2.25),
      tree_species == "birch" & is.na(dbh)  ~ 26.81  * (basal_diameter^2.62),
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    branches.br = if_else(tree_species == "birch", live.crown.br - foliage.br, NA_real_),
    fine.br     = if_else(tree_species == "birch", branches.br * 0.15, NA_real_),
    coarse.br   = if_else(tree_species == "birch", branches.br * 0.85, NA_real_)
  )

# -- POPLAR (Alexander et al. 2012) --
cce.trees <- cce.trees %>%
  mutate(
    total.biomass.bp = case_when(
      tree_species == "poplar" & !is.na(dbh) ~ 133.71 * (dbh^2.29),
      tree_species == "poplar" & is.na(dbh)  ~ 58.29  * (basal_diameter^2.44),
      TRUE ~ NA_real_
    ),
    foliage.bp = case_when(
      tree_species == "poplar" & !is.na(dbh) ~ 22.16  * (dbh^1.62),
      tree_species == "poplar" & is.na(dbh)  ~ 11.22  * (basal_diameter^1.76),
      TRUE ~ NA_real_
    ),
    live.crown.bp = case_when(
      tree_species == "poplar" & !is.na(dbh) ~ 17.24  * (dbh^2.3),
      tree_species == "poplar" & is.na(dbh)  ~ 7.36   * (basal_diameter^2.46),
      TRUE ~ NA_real_
    ),
    stemwood.bark.bp = case_when(
      tree_species == "poplar" & !is.na(dbh) ~ 98.26  * (dbh^2.32),
      tree_species == "poplar" & is.na(dbh)  ~ 43.35  * (basal_diameter^2.47),
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    branches.bp = if_else(tree_species == "poplar", live.crown.bp - foliage.bp, NA_real_),
    fine.bp     = if_else(tree_species == "poplar", branches.bp * 0.15, NA_real_),
    coarse.bp   = if_else(tree_species == "poplar", branches.bp * 0.85, NA_real_)
  )

# -- WHITE SPRUCE (Alexander et al. 2012) --
cce.trees <- cce.trees %>%
  mutate(
    total.biomass.ws = case_when(
      tree_species == "white spruce" & !is.na(dbh) ~ 96.77  * (dbh^2.4),
      tree_species == "white spruce" & is.na(dbh)  ~ 53.74  * (basal_diameter^2.45),
      TRUE ~ NA_real_
    ),
    foliage.ws = case_when(
      tree_species == "white spruce" & !is.na(dbh) ~ 25.22  * (dbh^2.04),
      tree_species == "white spruce" & is.na(dbh)  ~ 2.59   * (basal_diameter^2.58),
      TRUE ~ NA_real_
    ),
    live.crown.ws = case_when(
      tree_species == "white spruce" & !is.na(dbh) ~ 29.34  * (dbh^2.24),
      tree_species == "white spruce" & is.na(dbh)  ~ 3.01   * (basal_diameter^2.79),
      TRUE ~ NA_real_
    ),
    stemwood.bark.ws = case_when(
      tree_species == "white spruce" & !is.na(dbh) ~ 48.44  * (dbh^2.51),
      tree_species == "white spruce" & is.na(dbh)  ~ 39.13  * (basal_diameter^2.44),
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    branches.ws = if_else(tree_species == "white spruce", live.crown.ws - foliage.ws, NA_real_),
    fine.ws     = if_else(tree_species == "white spruce", branches.ws * 0.60, NA_real_),
    coarse.ws   = if_else(tree_species == "white spruce", branches.ws * 0.40, NA_real_)
  )

# -- LARCH / TAMARACK (Ker 1980b, via Ter-Mikaelian & Korzukhin 1997) --
cce.trees <- cce.trees %>%
  mutate(
    total.biomass.la = case_when(
      tree_species == "larch" & !is.na(dbh) ~ 0.0946 * (dbh^2.3572),
      TRUE ~ NA_real_
    ),
    foliage.la = case_when(
      tree_species == "larch" & !is.na(dbh) ~ 0.0061 * (dbh^1.9790),
      TRUE ~ NA_real_
    ),
    branches.la = case_when(
      tree_species == "larch" & !is.na(dbh) ~ 0.0178 * (dbh^2.1727),
      TRUE ~ NA_real_
    ),
    stemwood.bark.la = case_when(
      tree_species == "larch" & !is.na(dbh) ~ 0.0609 * (dbh^2.4472),
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    fine.la   = if_else(tree_species == "larch", branches.la * 0.15, NA_real_),
    coarse.la = if_else(tree_species == "larch", branches.la * 0.85, NA_real_)
  )

#--------------------------------------------------------------------#
## 2c) Combustion per tree----
#--------------------------------------------------------------------#
cce.trees <- cce.trees %>%
  mutate(
    # Black spruce
    foliage.bs.comb       = boby.needles.bs * (foliage_combustion / 100),
    fine.bs.comb          = boby.fine.bs    * (fine_combustion    / 100),
    coarse.bs.comb        = boby.coarse.bs  * (coarse_combustion  / 100),
    dead.branches.bs.comb = boby.total.biomass.bs * 0.088,
    cones.bs.comb         = if_else(!is.na(cone_combustion),
                                    if_else(!is.na(dbh) & dbh >= 2.7,
                                            8.6 * (dbh^2) - 32.1,
                                            11.4 * (dbh^2)) * (cone_combustion / 100),
                                    NA_real_),
    # Aspen
    foliage.as.comb       = foliage.as * (foliage_combustion / 100),
    fine.as.comb          = fine.as    * (fine_combustion    / 100),
    coarse.as.comb        = coarse.as  * (coarse_combustion  / 100),
    dead.branches.as.comb = total.biomass.as * 0.014,
    # Birch
    foliage.br.comb       = foliage.br * (foliage_combustion / 100),
    fine.br.comb          = fine.br    * (fine_combustion    / 100),
    coarse.br.comb        = coarse.br  * (coarse_combustion  / 100),
    dead.branches.br.comb = total.biomass.br * 0.014,
    # Poplar
    foliage.bp.comb       = foliage.bp * (foliage_combustion / 100),
    fine.bp.comb          = fine.bp    * (fine_combustion    / 100),
    coarse.bp.comb        = coarse.bp  * (coarse_combustion  / 100),
    dead.branches.bp.comb = total.biomass.bp * 0.014,
    # White spruce
    foliage.ws.comb       = foliage.ws * (foliage_combustion / 100),
    fine.ws.comb          = fine.ws    * (fine_combustion    / 100),
    coarse.ws.comb        = coarse.ws  * (coarse_combustion  / 100),
    dead.branches.ws.comb = total.biomass.ws * 0.088,
    # Larch (tamarack) — dead branch fraction same as deciduous (0.014)
    foliage.la.comb       = foliage.la * (foliage_combustion / 100),
    fine.la.comb          = fine.la    * (fine_combustion    / 100),
    coarse.la.comb        = coarse.la  * (coarse_combustion  / 100),
    dead.branches.la.comb = total.biomass.la * 0.014
  )

#--------------------------------------------------------------------#
## 2d) Pre-fire, combustion, post-fire totals per tree----
#--------------------------------------------------------------------#
cce.trees <- cce.trees %>%
  rowwise() %>%
  mutate(
    prefire.tree.ind = sum(boby.total.biomass.bs, total.biomass.as, total.biomass.br,
                           total.biomass.bp, total.biomass.ws, total.biomass.la, na.rm = TRUE),
    combust.tree.ind = sum(foliage.bs.comb, fine.bs.comb, coarse.bs.comb,
                           dead.branches.bs.comb, cones.bs.comb,
                           foliage.as.comb, fine.as.comb, coarse.as.comb, dead.branches.as.comb,
                           foliage.br.comb, fine.br.comb, coarse.br.comb, dead.branches.br.comb,
                           foliage.bp.comb, fine.bp.comb, coarse.bp.comb, dead.branches.bp.comb,
                           foliage.ws.comb, fine.ws.comb, coarse.ws.comb, dead.branches.ws.comb,
                           foliage.la.comb, fine.la.comb, coarse.la.comb, dead.branches.la.comb,
                           na.rm = TRUE),
    postfire.tree.ind = prefire.tree.ind - combust.tree.ind
  ) %>%
  ungroup()

#--------------------------------------------------------------------#
## 2e) Snag biomass (stem only, Alexander et al. 2012)----
#--------------------------------------------------------------------#
cce.snags <- cce.trees %>%
  filter(snag_pre_fire == "yes") %>%
  mutate(
    stemwood.bark.snag = case_when(
      tree_species == "black spruce" & !is.na(dbh) & dbh >= 2.7 ~ 117.7 * (dbh^2) - 253.9,
      tree_species == "black spruce" & !is.na(dbh) & dbh < 2.7  ~ 113.5 * (dbh^2),
      tree_species == "black spruce" & is.na(dbh)               ~ 12.0  * (basal_diameter^2),
      tree_species == "aspen"        & !is.na(dbh) ~ 64.01  * (dbh^2.51),
      tree_species == "aspen"        & is.na(dbh)  ~ 12.55  * (basal_diameter^2.98),
      tree_species == "birch"        & !is.na(dbh) ~ 147.96 * (dbh^2.25),
      tree_species == "birch"        & is.na(dbh)  ~ 26.81  * (basal_diameter^2.62),
      tree_species == "poplar"       & !is.na(dbh) ~ 98.26  * (dbh^2.32),
      tree_species == "poplar"       & is.na(dbh)  ~ 43.35  * (basal_diameter^2.47),
      tree_species == "white spruce" & !is.na(dbh) ~ 48.44  * (dbh^2.51),
      tree_species == "white spruce" & is.na(dbh)  ~ 39.13  * (basal_diameter^2.44),
      tree_species == "larch"        & !is.na(dbh) ~ 0.0609 * (dbh^2.4472),
      TRUE ~ NA_real_
    ),
    combust.snag.ind  = stemwood.bark.snag * (stem_combustion / 100),
    postfire.snag.ind = stemwood.bark.snag - combust.snag.ind
  )

#--------------------------------------------------------------------#
## 2f) Plot-level tree summaries----
#--------------------------------------------------------------------#

# Live trees only
cce.trees.live <- cce.trees %>% filter(snag_pre_fire == "no" | is.na(snag_pre_fire))

PLOTS.TREES <- cce.trees.live %>%
  group_by(fire_scar, site, plot) %>%
  summarise(
    prefire.trees    = sum(prefire.tree.ind,  na.rm = TRUE) / PLOT_AREA_TREES,
    combustion.trees = sum(combust.tree.ind,  na.rm = TRUE) / PLOT_AREA_TREES,
    postfire.trees   = sum(postfire.tree.ind, na.rm = TRUE) / PLOT_AREA_TREES,
    dens.trees       = n()                                  / PLOT_AREA_TREES,
    basal.area.trees = sum(basal.area.dbh,    na.rm = TRUE) / PLOT_AREA_TREES,
    .groups = "drop"
  ) %>%
  mutate(
    prefire.trees.carbon    = prefire.trees    * 0.5,
    combustion.trees.carbon = combustion.trees * 0.5,
    postfire.trees.carbon   = postfire.trees   * 0.5
  )

# Snags
PLOTS.SNAGS <- cce.snags %>%
  group_by(fire_scar, site, plot) %>%
  summarise(
    prefire.snags    = sum(stemwood.bark.snag, na.rm = TRUE) / PLOT_AREA_TREES,
    combustion.snags = sum(combust.snag.ind,   na.rm = TRUE) / PLOT_AREA_TREES,
    postfire.snags   = sum(postfire.snag.ind,  na.rm = TRUE) / PLOT_AREA_TREES,
    .groups = "drop"
  ) %>%
  mutate(
    prefire.snags.carbon    = prefire.snags    * 0.5,
    combustion.snags.carbon = combustion.snags * 0.5,
    postfire.snags.carbon   = postfire.snags   * 0.5
  )

#--------------------------------------------------------------------#
## 2g) Species-level density, basal area, biomass at plot level----
#--------------------------------------------------------------------#

# Species density (stems/m2)
# any_of() used in rename because not all species may be present
PLOTS.DENS <- cce.trees.live %>%
  group_by(fire_scar, site, plot, tree_species) %>%
  summarise(sp.dens = n() / PLOT_AREA_TREES, .groups = "drop") %>%
  pivot_wider(names_from = tree_species, values_from = sp.dens,
              names_glue = "{tree_species}.dens", values_fill = 0) %>%
  rename_with(~ str_replace_all(., c("black spruce" = "bs", "white spruce" = "ws",
                                     "aspen" = "as", "birch" = "br", "poplar" = "bp")),
              ends_with(".dens")) %>%
  rename(any_of(c(as.pre.dens = "as.dens", br.pre.dens = "br.dens",
                  bs.pre.dens = "bs.dens", bp.pre.dens = "bp.dens",
                  ws.pre.dens = "ws.dens")))

if (!"bp.pre.dens" %in% names(PLOTS.DENS)) PLOTS.DENS$bp.pre.dens <- 0
if (!"ws.pre.dens" %in% names(PLOTS.DENS)) PLOTS.DENS$ws.pre.dens <- 0

# Total tree density
PLOTS.DENS <- PLOTS.DENS %>%
  rowwise() %>%
  mutate(total.tree.pre.dens = sum(c_across(ends_with(".pre.dens")), na.rm = TRUE)) %>%
  ungroup()

# Species basal area (cm2/m2)
PLOTS.BA <- cce.trees.live %>%
  group_by(fire_scar, site, plot, tree_species) %>%
  summarise(sp.ba = sum(basal.area.dbh, na.rm = TRUE) / PLOT_AREA_TREES, .groups = "drop") %>%
  pivot_wider(names_from = tree_species, values_from = sp.ba,
              names_glue = "{tree_species}.ba", values_fill = 0) %>%
  rename_with(~ str_replace_all(., c("black spruce" = "bs", "white spruce" = "ws",
                                     "aspen" = "as", "birch" = "br", "poplar" = "bp")),
              ends_with(".ba")) %>%
  rename(any_of(c(as.pre.ba = "as.ba", br.pre.ba = "br.ba",
                  bs.pre.ba = "bs.ba", bp.pre.ba = "bp.ba",
                  ws.pre.ba = "ws.ba")))

# Species biomass (g/m2)
PLOTS.BIOMASS <- cce.trees.live %>%
  group_by(fire_scar, site, plot) %>%
  summarise(
    as.pre.biomass  = sum(total.biomass.as,      na.rm = TRUE) / PLOT_AREA_TREES,
    br.pre.biomass  = sum(total.biomass.br,      na.rm = TRUE) / PLOT_AREA_TREES,
    bp.pre.biomass  = sum(total.biomass.bp,      na.rm = TRUE) / PLOT_AREA_TREES,
    bs.pre.biomass  = sum(boby.total.biomass.bs, na.rm = TRUE) / PLOT_AREA_TREES,
    ws.pre.biomass  = sum(total.biomass.ws,      na.rm = TRUE) / PLOT_AREA_TREES,
    as.comb.biomass = sum(coarse.as.comb + fine.as.comb + foliage.as.comb + dead.branches.as.comb, na.rm = TRUE) / PLOT_AREA_TREES,
    br.comb.biomass = sum(coarse.br.comb + fine.br.comb + foliage.br.comb + dead.branches.br.comb, na.rm = TRUE) / PLOT_AREA_TREES,
    bp.comb.biomass = sum(coarse.bp.comb + fine.bp.comb + foliage.bp.comb + dead.branches.bp.comb, na.rm = TRUE) / PLOT_AREA_TREES,
    bs.comb.biomass = sum(coarse.bs.comb + fine.bs.comb + foliage.bs.comb + dead.branches.bs.comb, na.rm = TRUE) / PLOT_AREA_TREES,
    ws.comb.biomass = sum(coarse.ws.comb + fine.ws.comb + foliage.ws.comb + dead.branches.ws.comb, na.rm = TRUE) / PLOT_AREA_TREES,
    .groups = "drop"
  ) %>%
  mutate(
    as.post.biomass = as.pre.biomass - as.comb.biomass,
    br.post.biomass = br.pre.biomass - br.comb.biomass,
    bp.post.biomass = bp.pre.biomass - bp.comb.biomass,
    bs.post.biomass = bs.pre.biomass - bs.comb.biomass,
    ws.post.biomass = ws.pre.biomass - ws.comb.biomass,
    total.tree.pre.biomass = as.pre.biomass + br.pre.biomass + bp.pre.biomass +
      bs.pre.biomass + ws.pre.biomass
  )

#--------------------------------------------------------------------#
## 2h) ForestType----
#--------------------------------------------------------------------#
SPECIES.INDEX <- PLOTS.BIOMASS %>%
  left_join(PLOTS.DENS, by = c("fire_scar", "site", "plot")) %>%
  mutate(across(any_of(c("bs.pre.dens", "ws.pre.dens", "as.pre.dens",
                         "br.pre.dens", "bp.pre.dens")), ~ replace_na(., 0))) %>%
  mutate(
    PIMA = ((bs.pre.biomass / total.tree.pre.biomass) + (bs.pre.dens / total.tree.pre.dens)) / 2,
    PIGL = ((ws.pre.biomass / total.tree.pre.biomass) + (ws.pre.dens / total.tree.pre.dens)) / 2,
    POTR = ((as.pre.biomass / total.tree.pre.biomass) + (as.pre.dens / total.tree.pre.dens)) / 2,
    BENE = ((br.pre.biomass / total.tree.pre.biomass) + (br.pre.dens / total.tree.pre.dens)) / 2,
    POBA = ((bp.pre.biomass / total.tree.pre.biomass) + (bp.pre.dens / total.tree.pre.dens)) / 2
  ) %>%
  mutate(
    ForestType = case_when(
      PIMA >= 0.66        ~ "PIMA",
      PIGL >= 0.66        ~ "PIGL",
      BENE >= 0.66        ~ "BENE",
      POTR >= 0.66        ~ "POTR",
      POBA >= 0.66        ~ "POBA",
      PIMA + PIGL >= 0.66 ~ "MIXED.CON",
      BENE + POTR >= 0.66 ~ "MIXED.DEC",
      TRUE                ~ "MIXED"
    ),
    ForestType = replace_na(ForestType, "NON.FORESTED")
  ) %>%
  select(fire_scar, site, plot, ForestType)

#=============================================================================#
# 3) Shrubs----
#=============================================================================#
# Quadrat area = 1m x 1m, 3 quadrats per plot = 3 m2 per plot
PLOT_AREA_SHRUBS <- 3

# Filter to rows with shrubs present and valid measurements
cce.shrubs.filt <- cce.shrubs %>%
  filter(shrub_regen_present == "yes",
         !is.na(shrub_regen_count),
         shrub_regen_count > 0)

# Pivot BD measurements to long, replicate rows by count
cce.shrubs.long <- cce.shrubs.filt %>%
  pivot_longer(cols = c(s_bd1, s_bd2, s_bd3), names_to = "bd_rep", values_to = "bd") %>%
  filter(!is.na(bd)) %>%
  group_by(fire_scar, site, plot, location, shrub_regen_species) %>%
  mutate(rep_weight = shrub_regen_count / n()) %>%
  ungroup()

# Shrub relative BA by site (for unknown species allometry — Berner et al. 2015)
shrub.rel.ba <- cce.shrubs.long %>%
  filter(!shrub_regen_species %in% c("Unknown", NA)) %>%
  group_by(fire_scar, site, shrub_regen_species) %>%
  summarise(ba = sum(pi * (bd / 2)^2, na.rm = TRUE), .groups = "drop") %>%
  group_by(fire_scar, site) %>%
  mutate(site.ba = sum(ba)) %>%
  ungroup() %>%
  mutate(
    rel.ba.a = if_else(shrub_regen_species == "Alnus",  ba / site.ba, NA_real_),
    rel.ba.s = if_else(shrub_regen_species == "Salix",  ba / site.ba, NA_real_),
    rel.ba.b = if_else(shrub_regen_species == "Betula", ba / site.ba, NA_real_)
  ) %>%
  select(fire_scar, site, rel.ba.a, rel.ba.s, rel.ba.b) %>%
  pivot_longer(starts_with("rel.ba"), names_to = "sp", values_to = "val") %>%
  filter(!is.na(val)) %>%
  pivot_wider(names_from = sp, values_from = val, values_fill = 0)

cce.shrubs.long <- cce.shrubs.long %>%
  left_join(shrub.rel.ba, by = c("fire_scar", "site")) %>%
  mutate(across(any_of(c("rel.ba.a", "rel.ba.s", "rel.ba.b")), ~ replace_na(., 0)))

if (!"rel.ba.a" %in% names(cce.shrubs.long)) cce.shrubs.long$rel.ba.a <- 0
if (!"rel.ba.s" %in% names(cce.shrubs.long)) cce.shrubs.long$rel.ba.s <- 0
if (!"rel.ba.b" %in% names(cce.shrubs.long)) cce.shrubs.long$rel.ba.b <- 0

# Shrub allometry (Berner et al. 2015, g/stem)
cce.shrubs.long <- cce.shrubs.long %>%
  mutate(
    AGB.shrub = case_when(
      shrub_regen_species == "Alnus"   ~ 13.31 * (bd^3.15),
      shrub_regen_species == "Salix"   ~ 27.58 * (bd^2.36),
      shrub_regen_species == "Betula"  ~ 28.10 * (bd^2.97),
      shrub_regen_species == "Unknown" ~ rel.ba.a * 13.31 * (bd^3.15) +
        rel.ba.s * 27.58 * (bd^2.36) +
        rel.ba.b * 28.10 * (bd^2.97),
      TRUE ~ NA_real_
    ),
    stem.shrub = case_when(
      shrub_regen_species == "Alnus"   ~ 4.28  * (bd^3.68),
      shrub_regen_species == "Salix"   ~ 20.62 * (bd^2.29),
      shrub_regen_species == "Betula"  ~ 17.47 * (bd^2.36),
      shrub_regen_species == "Unknown" ~ rel.ba.a * 4.28  * (bd^3.68) +
        rel.ba.s * 20.62 * (bd^2.29) +
        rel.ba.b * 17.47 * (bd^2.36),
      TRUE ~ NA_real_
    ),
    branch.shrub = case_when(
      shrub_regen_species == "Alnus"   ~ 2.23  * (bd^3.19),
      shrub_regen_species == "Salix"   ~ 5.19  * (bd^2.37),
      shrub_regen_species == "Betula"  ~ 5.73  * (bd^4.06),
      shrub_regen_species == "Unknown" ~ rel.ba.a * 2.23  * (bd^3.19) +
        rel.ba.s * 5.19  * (bd^2.37) +
        rel.ba.b * 5.73  * (bd^4.06),
      TRUE ~ NA_real_
    ),
    ng.shrub = case_when(
      shrub_regen_species == "Alnus"   ~ 10.88 * (bd^1.55),
      shrub_regen_species == "Salix"   ~ 10.54 * (bd^1.71),
      shrub_regen_species == "Betula"  ~ 4.57  * (bd^2.45),
      shrub_regen_species == "Unknown" ~ rel.ba.a * 10.88 * (bd^1.55) +
        rel.ba.s * 10.54 * (bd^1.71) +
        rel.ba.b * 4.57  * (bd^2.45),
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    prefire.shrub.ind = stem.shrub + branch.shrub + ng.shrub
  )

# Plot-level shrub summary
PLOTS.SHRUBS <- cce.shrubs.long %>%
  group_by(fire_scar, site, plot) %>%
  summarise(
    prefire.shrubs    = sum(prefire.shrub.ind * rep_weight, na.rm = TRUE) / PLOT_AREA_SHRUBS,
    dens.shrubs       = sum(rep_weight,                     na.rm = TRUE) / PLOT_AREA_SHRUBS,
    basal.area.shrubs = sum(pi * (bd / 2)^2 * rep_weight,  na.rm = TRUE) / PLOT_AREA_SHRUBS,
    .groups = "drop"
  ) %>%
  mutate(
    prefire.shrubs.carbon = prefire.shrubs * 0.5
  )

#=============================================================================#
# 4) CWD Carbon----
#=============================================================================#
# Plot transect = 10m, diameter in cm
# Using Nalder et al. 1999 for <7cm, Ter-Mikaelian / Manies et al. 2005 for >=7cm
# Output in g C / m2 (multiply Mg/ha by 50)
cce.cwd <- cce.cwd %>%
  filter(!is.na(diameter_cm), diameter_cm >= 5) %>%
  mutate(size_class = if_else(diameter_cm < 7, "sm", "lg"))

if (!"rel.ba.bp" %in% names(rel.ba.wide)) rel.ba.wide$rel.ba.bp <- 0
if (!"rel.ba.ws" %in% names(rel.ba.wide)) rel.ba.wide$rel.ba.ws <- 0

# Join site-level relative BA for unknown CWD species
cwd.rel.ba <- rel.ba.wide %>%
  mutate(
    rel.ba.dec = rel.ba.as + rel.ba.br + rel.ba.bp,
    rel.ba.con = rel.ba.bs + rel.ba.ws
  ) %>%
  select(fire_scar, site, rel.ba.bs, rel.ba.ws, rel.ba.dec, rel.ba.con)

cce.cwd <- cce.cwd %>%
  left_join(cwd.rel.ba, by = c("fire_scar", "site"))

# Carbon per piece (g C / m2)
cce.cwd <- cce.cwd %>%
  mutate(
    postfire.cwd.piece = case_when(
      # Small (<7cm) — Nalder et al. 2000 multipliers, transect = 10m, *50 to get gC/m2
      species == "black spruce" & size_class == "sm" ~ (1 * 2.68 / 10) * 50,
      species == "white spruce" & size_class == "sm" ~ (1 * 2.30 / 10) * 50,
      species %in% c("aspen", "birch") & size_class == "sm" ~ (1 * 1.70 / 10) * 50,
      species %in% c("unknown", "unknown conifer") & size_class == "sm" ~
        ((rel.ba.dec * (1 * 1.70 / 10)) + (rel.ba.bs * (1 * 2.68 / 10)) +
           (rel.ba.ws * (1 * 2.30 / 10))) * 50,
      # Large (>=7cm) — Ter-Mikaelian et al. 2008 wood densities (g/cm3), transect length = 80 (8m x 10)
      # Decay classes: hard = class I, soft = class III, crumbly = class V (Manies et al. 2005)
      # Black spruce densities from Ter-Mikaelian Table 3
      species == "black spruce" & size_class == "lg" & decay_class == "hard"    ~ (diameter_cm^2 * 0.427) / 80,
      species == "black spruce" & size_class == "lg" & decay_class == "soft"    ~ (diameter_cm^2 * 0.289) / 80,
      species == "black spruce" & size_class == "lg" & decay_class == "crumbly" ~ (diameter_cm^2 * 0.151) / 80,
      # Birch densities from Ter-Mikaelian Table 3 (white birch)
      species == "birch" & size_class == "lg" & decay_class == "hard"    ~ (diameter_cm^2 * 0.521) / 80,
      species == "birch" & size_class == "lg" & decay_class == "soft"    ~ (diameter_cm^2 * 0.357) / 80,
      species == "birch" & size_class == "lg" & decay_class == "crumbly" ~ (diameter_cm^2 * 0.192) / 80,
      TRUE ~ NA_real_
    ),
    prefire.cwd.piece = postfire.cwd.piece * (1 / (1 - (combustion / 100))),
    comb.cwd.piece    = prefire.cwd.piece - postfire.cwd.piece
  )

# Plot-level CWD summary
PLOTS.CWD <- cce.cwd %>%
  group_by(fire_scar, site, plot) %>%
  summarise(
    prefire.cwd.carbon    = sum(prefire.cwd.piece, na.rm = TRUE),
    combustion.cwd.carbon = sum(comb.cwd.piece,    na.rm = TRUE),
    postfire.cwd.carbon   = sum(postfire.cwd.piece, na.rm = TRUE),
    .groups = "drop"
  )

#=============================================================================#
# 5) Join all plot-level data----
#=============================================================================#
PLOTS.ALL <- PLOTS.TREES %>%
  full_join(PLOTS.SNAGS,       by = c("fire_scar", "site", "plot")) %>%
  full_join(PLOTS.SHRUBS,      by = c("fire_scar", "site", "plot")) %>%
  full_join(PLOTS.CWD,         by = c("fire_scar", "site", "plot")) %>%
  full_join(PLOTS.BIOMASS,     by = c("fire_scar", "site", "plot")) %>%
  full_join(PLOTS.DENS,        by = c("fire_scar", "site", "plot")) %>%
  full_join(PLOTS.BA,          by = c("fire_scar", "site", "plot")) %>%
  full_join(SPECIES.INDEX,     by = c("fire_scar", "site", "plot")) %>%
  full_join(cce.plot.soil.sum, by = c("fire_scar", "site", "plot"))

# Fill missing snag/shrub/cwd values with 0
PLOTS.ALL <- PLOTS.ALL %>%
  mutate(across(c(prefire.snags, combustion.snags, postfire.snags,
                  prefire.snags.carbon, combustion.snags.carbon, postfire.snags.carbon,
                  prefire.shrubs, prefire.shrubs.carbon,
                  prefire.cwd.carbon, combustion.cwd.carbon, postfire.cwd.carbon),
                ~ replace_na(., 0)))

# Aboveground totals (biomass, g/m2)
PLOTS.ALL <- PLOTS.ALL %>%
  mutate(
    prefire.above    = prefire.trees  + prefire.shrubs  + prefire.snags,
    combustion.above = combustion.trees + combustion.snags,
    postfire.above   = postfire.trees + postfire.snags
  )

#=============================================================================#
# 6) Join site master and export----
#=============================================================================#
# Left join expands site master to plot level (A, B, C per site)
PLOTS.FINAL <- cce.site %>%
  distinct(fire_scar, site, .keep_all = TRUE) %>%
  left_join(PLOTS.ALL, by = c("fire_scar", "site")) %>%
  select(-any_of(c("unknown.dens", "unknown conifer.dens", "NA.dens",
                   "unknown.ba", "unknown conifer.ba",
                   "NA.ba", "bp.pre.ba")))

write_csv(PLOTS.FINAL, "Field_Data/Analysis_Ready/CCE_plot_combustion.csv")

#=============================================================================#
# End of Script
#=============================================================================#