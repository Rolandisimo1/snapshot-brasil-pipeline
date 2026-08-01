# Snapshot Brasil — Final Combined Dataset

Three files, all joined on `site_id` / `array_id`:

## final_deployments.csv (1,205 rows)
One row per camera deployment (a camera running continuously at one location
for one period). Columns: network (Snapshot / Atlantic / WI), array,
placename, deployment_id, longitude, latitude, start_date, end_date,
camera_days, array_id, site_id, cdays.

`site_id` groups deployments that represent the same physical camera location
across servicing visits (card swaps). `array_id` groups sites into the
spatial sampling unit used as the replication unit for every model in this
pipeline (>=10 cameras, clustered within 5 km, 200 m thinning applied).

## final_detections_mammals.csv (24,135 rows)
One row per independent mammal detection (60-minute independence interval
applied). Columns: network, array, placename, deployment_id, scientific_name
(as recorded), sci_mdd (harmonized to Mammal Diversity Database v2.4),
common_name, date, time, array_id, site_id.

Restricted to MAMMALS throughout the pipeline: the Atlantic Forest network's
source compilation recorded mammals only, so including birds would make every
bird species falsely "absent" on every Atlantic array -- a data artifact, not
ecology.

## final_covariates.csv (1,125 rows)
One row per site_id, with pre-extracted remote-sensing covariates (Google
Earth Engine): MapBiomas land-cover fractions at 100 m (forest_100m,
savanna_100m, pasture_100m, cropland_100m, agriculture_100m [legacy combined
pasture+cropland, superseded by the two separate columns], cropland_pct_100m
[legacy duplicate of cropland_100m]) and 1000 m (native_veg_1000m), two
global cross-checks (wc_tree_pct_100m, treecover2000_100m), canopy height
(canopy_height_m_point, canopy_height_m_100m), GHSL built-surface
(ghsl_built_100m), TerraClimate temperature and precipitation (temp_mean_C,
precip_annual_mm), MapBiomas point-level habitat class (pt_hab), WorldCover
point class (wc_class_point), and WDPA protected-area status (in_pa, 1/0).

The seven-layer FINAL MODELING SET used throughout the pipeline (after
removing redundant/collinear layers) is: forest_100m, savanna_100m,
pasture_100m, cropland_100m, native_veg_1000m, temp_mean_C, precip_annual_mm.
in_pa is tested as its own predictor, kept separate from this set.

Not every site has covariates (1,125 of 1,205 deployment-table sites) -- a
small number of sites lack a clean coordinate or fall just outside covariate
raster coverage.

## Joining the three files

```r
library(data.table)
dep <- fread("final_deployments.csv")
det <- fread("final_detections_mammals.csv")
cov <- fread("final_covariates.csv")

# effort per site
eff <- dep[, .(camera_days = sum(camera_days, na.rm=TRUE)), by = .(network, array_id, site_id)]

# join covariates onto sites
site_cov <- merge(eff, cov, by = "site_id", all.x = TRUE)

# detections already carry site_id/array_id directly
```
