# CCE Plot-Level Combustion

**Author:** Zach Madsen  
**Created for:** Xanthe Walker and CCE Project 

Compiles pre-fire biomass, combustion, and post-fire biomass estimates at the plot level across trees, snags, shrubs, coarse woody debris (CWD), and soil. Output is a single analysis-ready CSV joined with site metadata.

---

## Inputs

All files read from `Field_Data/Analysis_Ready/`:

- `CCE_soils_Master.csv`
- `CCE_site_Master.csv`
- `CCE_combustion_Master.csv` — trees and snags
- `CCE_shrub_Master.csv`
- `CCE_cwd_Master.csv`

## Output

`Field_Data/Analysis_Ready/CCE_plot_combustion.csv`

One row per plot. Contains pre-fire biomass, combustion, and post-fire biomass (and carbon equivalents) for each fuel type, plus species composition indices and soil burn metrics.

---

## Methods

### 1. Soil

Adventitious root height (ARH) measured on burned black spruce trees is used to reconstruct pre-fire organic soil depth (Boby et al. 2010). Where individual ARH measurements exist (up to three per point), they are averaged. Where only a pre-calculated average is recorded (`ar_height_avg`), that value is used directly. A 3.4 cm offset is then added to correct for the position of roots relative to the moss surface:

```
burn_depth = mean(ARH) + 3.4 cm
prop_soil_depth = burn_depth / (post-fire depth + burn_depth)
```

Both metrics are averaged to the plot level across all measurement points.

---

### 2. Trees

**Plot area:** 20 m² (2 × 10 m transect)

#### Allometric equations

Biomass is estimated per tree using published power equations of the form `y = a * x^b`, where x is DBH (cm) for trees ≥ 1.4 m tall or basal diameter (BD, cm) for trees < 1.4 m tall. Components (foliage, fine branches, coarse branches, stemwood/bark) are calculated separately to support combustion calculations.

**Black spruce — Boby et al. (2010), Appendix A**

For trees with DBH ≥ 2.7 cm:

| Component | Equation |
|---|---|
| New needles | `1.9 * DBH² + 16.0` |
| Old needles | `16.0 * DBH² + 288.4` |
| Fine branches | `(11.3 * DBH² + 81.2) + (0.5 * DBH² + 1.3)` |
| Coarse branches | `11.8 * DBH² − 90.6` |
| Stemwood/bark | `117.7 * DBH² − 253.9` |
| Cones | `8.6 * DBH² − 32.1` |

For trees with DBH < 2.7 cm:

| Component | Equation |
|---|---|
| Needles | `3.9 * ln(DBH)` |
| Fine branches | `11.3 * DBH²` |
| Coarse branches | `14.4 * DBH²` |
| Cones | `11.4 * DBH²` |
| Stemwood/bark | `113.5 * DBH²` |

For trees measured by basal diameter only:

| Component | Equation |
|---|---|
| Needles | `18.3 * BD²` |
| Fine branches | `7.5 * BD²` |
| Coarse branches | `3.1 * BD²` |
| Stemwood/bark | `12.0 * BD²` |

**Aspen (*Populus tremuloides*) — Alexander et al. (2012)**

| Component | DBH equation | BD equation |
|---|---|---|
| Total biomass | `134.10 * DBH^2.26` | `56.83 * BD^2.49` |
| Foliage | `18.98 * DBH^1.53` | `8.73 * BD^2.00` |
| Live crown | `41.74 * DBH^1.83` | `18.16 * BD^2.43` |
| Stemwood/bark | `64.01 * DBH^2.51` | `12.55 * BD^2.98` |

**Birch (*Betula neoalaskana*) — Alexander et al. (2012)**

| Component | DBH equation | BD equation |
|---|---|---|
| Total biomass | `164.18 * DBH^2.29` | `26.29 * BD^2.68` |
| Foliage | `6.39 * DBH^2.10` | `0.63 * BD^2.65` |
| Live crown | `15.15 * DBH^2.49` | `1.38 * BD^3.08` |
| Stemwood/bark | `147.96 * DBH^2.25` | `26.81 * BD^2.62` |

**Poplar (*Populus balsamifera*) — Alexander et al. (2012)**

| Component | DBH equation | BD equation |
|---|---|---|
| Total biomass | `133.71 * DBH^2.29` | `58.29 * BD^2.44` |
| Foliage | `22.16 * DBH^1.62` | `11.22 * BD^1.76` |
| Live crown | `17.24 * DBH^2.30` | `7.36 * BD^2.46` |
| Stemwood/bark | `98.26 * DBH^2.32` | `43.35 * BD^2.47` |

**Larch / Tamarack (*Larix laricina*) — Ker (1980b), via Ter-Mikaelian & Korzukhin (1997)**

DBH-only equations (no basal diameter equivalent available):

| Component | Equation |
|---|---|
| Total biomass | `0.0946 * DBH^2.3572` |
| Foliage | `0.0061 * DBH^1.9790` |
| Branches (total) | `0.0178 * DBH^2.1727` |
| Stemwood/bark | `0.0609 * DBH^2.4472` |

Fine/coarse branch split uses deciduous proportions (15%/85%). Dead branch fraction uses the deciduous coefficient (1.4% of total biomass). Note: no basal diameter equations are available for larch; trees recorded without DBH will return NA for biomass.

**White spruce (*Picea glauca*) — Alexander et al. (2012)**

| Component | DBH equation | BD equation |
|---|---|---|
| Total biomass | `96.77 * DBH^2.40` | `53.74 * BD^2.45` |
| Foliage | `25.22 * DBH^2.04` | `2.59 * BD^2.58` |
| Live crown | `29.34 * DBH^2.24` | `3.01 * BD^2.79` |
| Stemwood/bark | `48.44 * DBH^2.51` | `39.13 * BD^2.44` |

Branch biomass for deciduous species and white spruce is derived from live crown minus foliage. Fine and coarse fractions are split as follows: deciduous = 15%/85%; white spruce = 60%/40%.

Larch is reclassified as black spruce prior to allometry. Unknown species are apportioned by site-level relative basal area across known species.

#### Combustion

Field-recorded combustion scores (%) are applied to each biomass component. Dead branch biomass is estimated as a fixed fraction of total biomass: 8.8% for conifers, 1.4% for deciduous species (Boby et al. 2010).

#### Plot-level summary

All biomass values are divided by plot area (20 m²) to express as g/m². Carbon equivalents use a 0.5 conversion factor.

#### Forest type classification

Each plot is assigned a forest type based on a species dominance index — the mean of relative density and relative biomass. A species must exceed 66% to be classified as dominant.

| Code | Meaning |
|---|---|
| PIMA | Black spruce |
| PIGL | White spruce |
| POTR | Trembling aspen |
| BENE | Paper birch |
| POBA | Balsam poplar |
| MIXED.CON | Mixed conifer |
| MIXED.DEC | Mixed deciduous |
| MIXED | No dominant type |
| NON.FORESTED | No trees recorded |

---

### 3. Snags

Pre-fire standing dead trees are treated separately using stem-only biomass (no foliage or branch components). The same species-specific stemwood/bark equations from Alexander et al. (2012) and Boby et al. (2010) are applied. Combustion uses the field-recorded stem combustion score. Scaled to g/m² using the same 20 m² plot area.

---

### 4. Shrubs

**Plot area:** 3 m² (three 1 × 1 m quadrats)

Shrub biomass is estimated using allometric equations from Berner et al. (2015), relating basal diameter to stem, branch, and new growth biomass:

| Species | Component | Equation |
|---|---|---|
| *Alnus* | Stem | `4.28 * BD^3.68` |
| *Alnus* | Branch | `2.23 * BD^3.19` |
| *Alnus* | New growth | `10.88 * BD^1.55` |
| *Salix* | Stem | `20.62 * BD^2.29` |
| *Salix* | Branch | `5.19 * BD^2.37` |
| *Salix* | New growth | `10.54 * BD^1.71` |
| *Betula* | Stem | `17.47 * BD^2.36` |
| *Betula* | Branch | `5.73 * BD^4.06` |
| *Betula* | New growth | `4.57 * BD^2.45` |

Pre-fire shrub biomass = stem + branch + new growth. No combustion scores were collected for shrubs; combustion and post-fire estimates are not calculated for this fuel type.

---

### 5. Coarse Woody Debris (CWD)

**Small pieces (5–< 7 cm diameter):** Mass estimated using multipliers from Nalder et al. (2000) for Class III pieces in the Northwest Territories, converted to g C/m²:

| Species | Multiplier (M) | 
|---|---|
| Black spruce | 2.68 |
| White spruce | 2.30 |
| Aspen/birch | 1.70 |

Formula: `postfire carbon = (1 * M / transect_length) * 50`

**Large pieces (≥ 7 cm diameter):** Mass estimated using wood density values from Ter-Mikaelian et al. (2008), Table 3. Field decay classes (hard, soft, crumbly) correspond to Manies et al. (2005) decay classes I, III, and V respectively.

| Species | Decay class | Density (g/cm³) |
|---|---|---|
| Black spruce | hard (I) | 0.427 |
| Black spruce | soft (III) | 0.289 |
| Black spruce | crumbly (V) | 0.151 |
| Birch | hard (I) | 0.521 |
| Birch | soft (III) | 0.357 |
| Birch | crumbly (V) | 0.192 |

Formula: `postfire carbon = (diameter_cm² * density) / 80`

Pre-fire carbon is back-calculated from post-fire carbon and the field-recorded combustion score:

```
prefire = postfire / (1 − combustion/100)
combustion loss = prefire − postfire
```

---

### 6. Plot-level summaries

Each fuel type is summarized to the plot level before joining. All biomass values are divided by plot area to express as g/m². Carbon equivalents use a 0.5 conversion factor throughout.

**Trees and snags** are summarized separately. For live trees, plot-level outputs include total pre-fire biomass, combustion, and post-fire biomass (g/m²), tree density (stems/m²), and basal area (cm²/m²), plus carbon equivalents for all three biomass pools. Species-level density, basal area, and biomass are also calculated separately for each species. Snag outputs follow the same structure using stem-only biomass.

**Shrubs** are summarized to pre-fire biomass (g/m²), density (stems/m²), basal area (cm²/m²), and carbon equivalent. No combustion or post-fire values are reported.

**CWD** is summarized to pre-fire, combustion, and post-fire carbon (g C/m²) directly, as the equations output carbon rather than dry mass.

**Soil** metrics (burn depth and proportion of soil depth consumed) are averaged across all measurement points within each plot.

**Forest type** is assigned per plot based on the species dominance index described in Section 2.

### 7. Final join and export

All plot-level summaries are joined sequentially by `fire_scar`, `site`, and `plot` using `full_join`, so that plots appearing in any fuel type dataset are retained even if absent from others. Missing snag, shrub, and CWD values are zero-filled after joining, indicating those fuel types were absent rather than unsampled.

Aboveground totals are then calculated across fuel types:

```
prefire.above    = prefire.trees  + prefire.shrubs + prefire.snags
combustion.above = combustion.trees + combustion.snags
postfire.above   = postfire.trees + postfire.snags
```

Shrubs are included in the pre-fire aboveground total but excluded from combustion and post-fire totals as combustion was not estimated for this fuel type.

The site master (`CCE_site_Master.csv`) is joined last using a left join from the site master, attaching site-level metadata to each plot. Plots without a matching site entry are dropped. Unknown-species density and basal area columns are removed from the final output before export to `CCE_plot_combustion.csv`.

---

## Corrections from Previous Version

The following errors were identified and corrected during verification against source literature:

**Large CWD wood densities (Ter-Mikaelian et al. 2008):** The previous code used incorrect density values (0.4336, 0.3950 for black spruce; 0.4192, 0.3613 for aspen/birch) that could not be verified against the source paper. These have been replaced with the correct values from Table 3. Additionally, soft and crumbly decay classes were previously grouped together with a single density value; they are now assigned separate densities consistent with Manies et al. (2005) decay class definitions. Aspen was also removed from the large CWD calculations as it does not appear in the CCE field data; birch-specific densities are now used directly.

**Black spruce needle biomass (DBH < 2.7 cm):** The previous code used `3.9 * DBH²` but Boby et al. (2010) Appendix A, Table A2 specifies `3.9 * ln(DBH)` for trees in this size class. Corrected to `3.9 * log(DBH)` (`log()` in R computes the natural logarithm).

**Larch reclassification removed:** The previous code reclassified all larch (*Larix laricina*) records as black spruce prior to allometry. Larch occurs at 6 sites in the CCE dataset (31 records). Species-specific allometric equations from Ker (1980b) via Ter-Mikaelian & Korzukhin (1997) have been added and larch is now treated as its own species. Combustion and post-fire estimates for shrubs were previously calculated assuming 100% combustion, as no field combustion scores were collected. These columns have been removed. Only pre-fire shrub biomass is now reported.

**Nalder et al. citation year:** Previously cited as 1999 in code comments; corrected to 2000 (Nalder et al. 2000, *International Journal of Wildland Fire* 9:85–99).

---

## References

- Alexander, H.D., et al. 2012. Implications of increased deciduous cover on stand structure and aboveground carbon pools of Alaskan boreal forests. *Ecosphere* 3(5):45.
- Berner, L.T., et al. 2015. Biomass allometry for alder, dwarf birch, and willow in boreal forest and tundra ecosystems. *Forest Ecology and Management* 337:110–118.
- Boby, L.A., et al. 2010. Quantifying fire severity, carbon, and nitrogen emissions in Alaska's boreal forest. *Ecological Applications* 20:1633–1647.
- Ker, M.F. 1980. Tree biomass equations for ten major species in Cumberland County, Nova Scotia. Canadian Forestry Service, Maritimes Forest Research Centre, Information Report M-X-108. 26 p.
- Manies, K.L., et al. 2005. Woody debris along an upland chronosequence in boreal Manitoba and its impact on long-term carbon storage. *Canadian Journal of Forest Research* 35:472–482.
- Nalder, I.A., et al. 2000. Physical properties of dead and downed round-wood fuels in the boreal forests of western and northern Canada. *International Journal of Wildland Fire* 9:85–99.
- Ter-Mikaelian, M.T., and M.D. Korzukhin. 1997. Biomass equations for sixty-five North American tree species. *Forest Ecology and Management* 97:1–24.
- Ter-Mikaelian, M.T., et al. 2008. Amount of downed woody debris and its prediction using stand characteristics in boreal and mixedwood forests of Ontario, Canada. *Canadian Journal of Forest Research* 38:2189–2197.

---

## Dependencies

```r
library(tidyverse)
library(fs)
```
