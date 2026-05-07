library(lidR)
library(terra)
library(sf)

## Wolfswarte site ##

# read LAS and polygons (10 m buffer shape file)
las_ww <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250819-080319_Wolfswarte.las")
poly_train <- st_read("D:/Kiran/Harz_Results/Plot_buff10/Plot_buff10.shp") #clip from a shape file (we added 10m buffer to the original polygon)

# CHM settings
thr <- c(0, 2, 5, 10, 15)
edg <- c(0, 1.5)

# empty list to store normalized LAS
las_norm_list <- list()

# loop plots: clip -> ground classify -> normalize -> CHM
for (i in 1:2) {
  
  id <- paste0("ww_00", i)   # ww_001, ww_002
  poly_ww <- poly_train[poly_train$id == id, ] #because we have 2 plots in the same polygon, we have to clip it separately
  
  las_roi_ww      <- clip_roi(las_ww, poly_ww) #clipping two polygons within roi
  las_ww_grn      <- classify_ground(las_roi_ww, csf()) #since our las is not ground-classified, we first classify it with
  las_norm_roi_ww <- normalize_height(las_ww_grn, knnidw(k = 10, p = 2)) #normalizing height (Kriging-based method)
  
  
  # store normalized LAS for later metrics
  las_norm_list[[paste0("ww", i)]] <- las_norm_roi_ww #we use this for height percentiles and gap fraction
}  
  
##create CHM, Khosravipour et al. pitfree algorithm, pit free method avoids holes under dense canopy
chm_ww <- rasterize_canopy(las_norm_roi_ww, 0.1, pitfree(thr, edg))
writeRaster(
  chm_ww,
  paste0("D:/Kiran/Harz_Results/CHMs/CHM_WW", i, ".tif"),
  overwrite = TRUE
)


## LiDAR height percentiles (h95 and h10) at 10 cm resolution
percentiles <- list(h95 = 0.95, h10 = 0.10)

for (p_name in names(percentiles)) {
  for (ww in names(las_norm_list)) {
    
    result <- pixel_metrics(
      las_norm_list[[ww]],
      ~quantile(Z, percentiles[[p_name]]),
      res = 0.10
    )
    
    writeRaster(
      result,
      paste0(
        "D:/Kiran/Harz_Results/",
        toupper(p_name), "s/",        #create folders such as H95s and H10s
        toupper(p_name), "_", toupper(ww), ".tif"
      ),
      overwrite = TRUE
    )
  }
}

## Additional LiDAR structural metrics (Hmean, Hmedian, Hsd, IQR, Cover > 2m) at 10 cm resolution
for (ww in names(las_norm_list)) {
  
  metrics_extra <- pixel_metrics(
    las_norm_list[[ww]],
    ~list(
      Hmean   = mean(Z, na.rm = TRUE),                                # mean height 
      Hmed = median(Z, na.rm = TRUE),                                 # median height (H50) 
      Hsd     = sd(Z, na.rm = TRUE),                                  # height standard deviation 
      IQR     = quantile(Z, 0.75, na.rm = TRUE) -
        quantile(Z, 0.25, na.rm = TRUE),                     # interquartile range
      Cover_above2m = sum(Z > 2, na.rm = TRUE) / sum(Z > 0, na.rm = TRUE)   # canopy cover above 2 m 
    ),
    res = 0.10
  )
  
  out_dirs <- c(Hmean="Hmeans", Hmed="Hmeds", IQR="IQRs", Hsd="Hsds", Cover_above2m="Cover_above2m")   # c() collects key-value pairs into one object
  for (m in names(out_dirs)) writeRaster(metrics_extra[[m]],
      paste0("D:/Kiran/Harz_Results/", out_dirs[m], "/", m, "_", toupper(ww), ".tif"),
      overwrite = TRUE)

}

##calculating gap fraction, proportion of LiDAR returns below a height threshold, e.g.2 m
for (ww in names(las_norm_list)) {
  
  gap_frac <- pixel_metrics(
    las_norm_list[[ww]],
    ~mean(Z < 2),
    res = 0.1
  )
  
  writeRaster(
    gap_frac,
    paste0(
      "D:/Kiran/Harz_Results/GapFrac_below2m/GapFrac_below2m_",
      toupper(ww),
      ".tif"
    ),
    overwrite = TRUE
  )
}

rm(las_ww) #removing the las file for next steps as it consumes higher memory
gc()

##For remaining sites, comments remain the same as Wolfswarte site except plot and polygon numbers##

## Grane site ##

# read LAS
las_gos <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250829-100646_Grane_OutliersRemoved.las")

# CHM settings
thr <- c(0, 2, 5, 10, 15)
edg <- c(0, 1.5)

# store normalized LAS for later metrics
las_norm_list <- list()

# loop plots: clip -> ground classify -> normalize -> CHM
for (i in 1:2) {
  
  id <- paste0("gos_00", i)   # gos_001, gos_002
  poly_gos <- poly_train[poly_train$id == id, ] #because we have 2 plots in the same polygon, we have to clip it separately
  
  las_roi_gos      <- clip_roi(las_gos, poly_gos) 
  las_gos_grn      <- classify_ground(las_roi_gos, csf()) 
  las_norm_roi_gos <- normalize_height(las_gos_grn, knnidw(k = 10, p = 2)) 
  # save normalized LAS
  las_norm_list[[paste0("gos", i)]] <- las_norm_roi_gos 
}

##create CHM
chm_gos <- rasterize_canopy(las_norm_roi_gos, 0.1, pitfree(thr, edg))
writeRaster(
  chm_gos,
  paste0("D:/Kiran/Harz_Results/CHMs/CHM_GOS", i, ".tif"),
  overwrite = TRUE
)


##calculating LiDAR derived height percentile (h95 and h10) at 10 cm
percentiles <- list(h95 = 0.95, h10 = 0.10)

for (p_name in names(percentiles)) {
  for (gos in names(las_norm_list)) {
    
    result <- pixel_metrics(
      las_norm_list[[gos]],
      ~quantile(Z, percentiles[[p_name]]),
      res = 0.10
    )
    
    writeRaster(
      result,
      paste0(
        "D:/Kiran/Harz_Results/",
        toupper(p_name), "s/",
        toupper(p_name), "_", toupper(gos), ".tif"
      ),
      overwrite = TRUE
    )
  }
}

## Additional LiDAR structural metrics (Hmean, Hmedian, Hsd, IQR, Cover > 2m) at 10 cm resolution
for (gos in names(las_norm_list)) {
  
  metrics_extra <- pixel_metrics(
    las_norm_list[[gos]],
    ~list(
      Hmean   = mean(Z, na.rm = TRUE),                                # mean height 
      Hmed = median(Z, na.rm = TRUE),                                 # median height (H50) 
      Hsd     = sd(Z, na.rm = TRUE),                                  # height standard deviation 
      IQR     = quantile(Z, 0.75, na.rm = TRUE) -
        quantile(Z, 0.25, na.rm = TRUE),                     # interquartile range
      Cover_above2m = sum(Z > 2, na.rm = TRUE) / sum(Z > 0, na.rm = TRUE)   # canopy cover above 2 m 
    ),
    res = 0.10
  )
  
  out_dirs <- c(Hmean="Hmeans", Hmed="Hmeds", IQR="IQRs", Hsd="Hsds", Cover_above2m="Cover_above2m")   # c() collects key-value pairs into one object
  for (m in names(out_dirs)) writeRaster(metrics_extra[[m]],
                                         paste0("D:/Kiran/Harz_Results/", out_dirs[m], "/", m, "_", toupper(gos), ".tif"),
                                         overwrite = TRUE)
  
}

##calculating gap fraction (returns below 2m)
for (gos in names(las_norm_list)) {
  
  gap_frac_gos <- pixel_metrics(
    las_norm_list[[gos]],
    ~mean(Z < 2),
    res = 0.1
  )
  
  writeRaster(
    gap_frac_gos,
    paste0(
      "D:/Kiran/Harz_Results/GapFrac_below2m/GapFrac_below2m_",
      toupper(gos),
      ".tif"
    ),
    overwrite = TRUE
  )
}

rm(las_gos)
gc()

## Gelber Brink Site ##

# read las files
las_gb <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250814-065407_GelberBrink_OutliersRemoved.las")

# CHM settings
thr <- c(0, 2, 5, 10, 15)
edg <- c(0, 1.5)

# store normalized LAS for later metrics
las_norm_list <- list()

# loop plots: clip -> ground classify -> normalize -> CHM
for (i in 1:3) {
  
  id <- paste0("gb_00", i)   # gb_001, gb_002, gb_003
  poly_gb <- poly_train[poly_train$id == id, ] #because we have 3 plots in the same polygon, we have to clip it separately
  
  las_roi_gb      <- clip_roi(las_gb, poly_gb)
  las_gb_grn      <- classify_ground(las_roi_gb, csf())
  las_norm_roi_gb <- normalize_height(las_gb_grn, knnidw(k = 10, p = 2))
  
  # save normalized LAS
  las_norm_list[[paste0("gb", i)]] <- las_norm_roi_gb
}

# CHM
chm_gb <- rasterize_canopy(las_norm_roi_gb, 0.1, pitfree(thr, edg))
writeRaster(
  chm_gb,
  paste0("D:/Kiran/Harz_Results/CHMs/CHM_GB", i, ".tif"),
  overwrite = TRUE
)


## LiDAR height percentiles (h95 and h10)
percentiles <- list(h95 = 0.95, h10 = 0.10)

for (p_name in names(percentiles)) {
  for (gb in names(las_norm_list)) {
    
    result <- pixel_metrics(
      las_norm_list[[gb]],
      ~quantile(Z, percentiles[[p_name]]),
      res = 0.10
    )
    
    writeRaster(
      result,
      paste0(
        "D:/Kiran/Harz_Results/",
        toupper(p_name), "s/",
        toupper(p_name), "_", toupper(gb), ".tif"
      ),
      overwrite = TRUE
    )
  }
}

## Additional LiDAR structural metrics (Hmean, Hmedian, Hsd, IQR, Cover > 2m) at 10 cm resolution
for (gb in names(las_norm_list)) {
  
  metrics_extra <- pixel_metrics(
    las_norm_list[[gb]],
    ~list(
      Hmean   = mean(Z, na.rm = TRUE),                                # mean height 
      Hmed = median(Z, na.rm = TRUE),                                 # median height (H50) 
      Hsd     = sd(Z, na.rm = TRUE),                                  # height standard deviation 
      IQR     = quantile(Z, 0.75, na.rm = TRUE) -
        quantile(Z, 0.25, na.rm = TRUE),                     # interquartile range
      Cover_above2m = sum(Z > 2, na.rm = TRUE) / sum(Z > 0, na.rm = TRUE)   # canopy cover above 2 m 
    ),
    res = 0.10
  )
  
  out_dirs <- c(Hmean="Hmeans", Hmed="Hmeds", IQR="IQRs", Hsd="Hsds", Cover_above2m="Cover_above2m")   # c() collects key-value pairs into one object
  for (m in names(out_dirs)) writeRaster(metrics_extra[[m]],
                                         paste0("D:/Kiran/Harz_Results/", out_dirs[m], "/", m, "_", toupper(gb), ".tif"),
                                         overwrite = TRUE)
  
}

## Gap fraction (returns below 2 m)
for (gb in names(las_norm_list)) {
  
  gap_frac_gb <- pixel_metrics(
    las_norm_list[[gb]],
    ~mean(Z < 2),
    res = 0.1
  )
  
  writeRaster(
    gap_frac_gb,
    paste0(
      "D:/Kiran/Harz_Results/GapFrac_below2m/GapFrac_below2m_",
      toupper(gb),
      ".tif"
    ),
    overwrite = TRUE
  )
}

rm(las_gb)
gc()

## Christianenhaus Site ##

# read las files
las_chh <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250820-154733_Christianenhaus.las")

# CHM settings
thr <- c(0, 2, 5, 10, 15)
edg <- c(0, 1.5)

# store normalized LAS for later metrics
las_norm_list <- list()

# loop plots: clip -> ground classify -> normalize -> CHM
for (i in 1:4) {
  
  id <- paste0("chh_00", i)   # chh_001 ... chh_004
  poly_chh <- poly_train[poly_train$id == id, ] #because we have 4 plots in the same polygon, we have to clip it separately
  
  las_roi_chh      <- clip_roi(las_chh, poly_chh)
  las_chh_grn      <- classify_ground(las_roi_chh, csf())
  las_norm_roi_chh <- normalize_height(las_chh_grn, knnidw(k = 10, p = 2))
  
  # save normalized LAS
  las_norm_list[[paste0("chh", i)]] <- las_norm_roi_chh
}

# CHM
chm_chh <- rasterize_canopy(las_norm_roi_chh, 0.1, pitfree(thr, edg))
writeRaster(
  chm_chh,
  paste0("D:/Kiran/Harz_Results/CHMs/CHM_CHH", i, ".tif"),
  overwrite = TRUE
)


## LiDAR height percentiles (h95 and h10)
percentiles <- list(h95 = 0.95, h10 = 0.10)

for (p_name in names(percentiles)) {
  for (chh in names(las_norm_list)) {
    
    result <- pixel_metrics(
      las_norm_list[[chh]],
      ~quantile(Z, percentiles[[p_name]]),
      res = 0.10
    )
    
    writeRaster(
      result,
      paste0(
        "D:/Kiran/Harz_Results/",
        toupper(p_name), "s/",
        toupper(p_name), "_", toupper(chh), ".tif"
      ),
      overwrite = TRUE
    )
  }
}

## Additional LiDAR structural metrics (Hmean, Hmedian, Hsd, IQR, Cover > 2m) at 10 cm resolution
for (chh in names(las_norm_list)) {
  
  metrics_extra <- pixel_metrics(
    las_norm_list[[chh]],
    ~list(
      Hmean   = mean(Z, na.rm = TRUE),                                # mean height 
      Hmed = median(Z, na.rm = TRUE),                                 # median height (H50) 
      Hsd     = sd(Z, na.rm = TRUE),                                  # height standard deviation 
      IQR     = quantile(Z, 0.75, na.rm = TRUE) -
        quantile(Z, 0.25, na.rm = TRUE),                     # interquartile range
      Cover_above2m = sum(Z > 2, na.rm = TRUE) / sum(Z > 0, na.rm = TRUE)   # canopy cover above 2 m 
    ),
    res = 0.10
  )
  
  out_dirs <- c(Hmean="Hmeans", Hmed="Hmeds", IQR="IQRs", Hsd="Hsds", Cover_above2m="Cover_above2m")   # c() collects key-value pairs into one object
  for (m in names(out_dirs)) writeRaster(metrics_extra[[m]],
                                         paste0("D:/Kiran/Harz_Results/", out_dirs[m], "/", m, "_", toupper(chh), ".tif"),
                                         overwrite = TRUE)
  
}

## Gap fraction (returns below 2 m)
for (chh in names(las_norm_list)) {
  
  gap_frac_chh <- pixel_metrics(
    las_norm_list[[chh]],
    ~mean(Z < 2),
    res = 0.1
  )
  
  writeRaster(
    gap_frac_chh,
    paste0(
      "D:/Kiran/Harz_Results/GapFrac_below2m/GapFrac_below2m_",
      toupper(chh),
      ".tif"
    ),
    overwrite = TRUE
  )
}

rm(las_chh)
gc()

## Molkenhausschause Site ##

# read las files
las_mk <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250809-160907_Molkenhauschausee_OutliersRemoved.las")

# CHM settings
thr <- c(0, 2, 5, 10, 15)
edg <- c(0, 1.5)

# store normalized LAS for later metrics
las_norm_list <- list()

# loop plots: clip -> ground classify -> normalize -> CHM
for (i in 1:4) {
  
  id <- paste0("mk_00", i)   # mk_001 ... mk_004
  poly_mk <- poly_train[poly_train$id == id, ] #because we have 4 plots in the same polygon, we have to clip it separately
  
  las_roi_mk      <- clip_roi(las_mk, poly_mk)
  las_mk_grn      <- classify_ground(las_roi_mk, csf())
  las_norm_roi_mk <- normalize_height(las_mk_grn, knnidw(k = 10, p = 2))
  
  # save normalized LAS
  las_norm_list[[paste0("mk", i)]] <- las_norm_roi_mk
}

# CHM
chm_mk <- rasterize_canopy(las_norm_roi_mk, 0.1, pitfree(thr, edg))
writeRaster(
  chm_mk,
  paste0("D:/Kiran/Harz_Results/CHMs/CHM_MK", i, ".tif"),
  overwrite = TRUE
)


## LiDAR height percentiles (h95 and h10)
percentiles <- list(h95 = 0.95, h10 = 0.10)

for (p_name in names(percentiles)) {
  for (mk in names(las_norm_list)) {
    
    result <- pixel_metrics(
      las_norm_list[[mk]],
      ~quantile(Z, percentiles[[p_name]]),
      res = 0.10
    )
    
    writeRaster(
      result,
      paste0(
        "D:/Kiran/Harz_Results/",
        toupper(p_name), "s/",
        toupper(p_name), "_", toupper(mk), ".tif"
      ),
      overwrite = TRUE
    )
  }
}

## Additional LiDAR structural metrics (Hmean, Hmedian, Hsd, IQR, Cover > 2m) at 10 cm resolution
for (mk in names(las_norm_list)) {
  
  metrics_extra <- pixel_metrics(
    las_norm_list[[mk]],
    ~list(
      Hmean   = mean(Z, na.rm = TRUE),                                # mean height 
      Hmed = median(Z, na.rm = TRUE),                                 # median height (H50) 
      Hsd     = sd(Z, na.rm = TRUE),                                  # height standard deviation 
      IQR     = quantile(Z, 0.75, na.rm = TRUE) -
        quantile(Z, 0.25, na.rm = TRUE),                     # interquartile range
      Cover_above2m = sum(Z > 2, na.rm = TRUE) / sum(Z > 0, na.rm = TRUE)   # canopy cover above 2 m 
    ),
    res = 0.10
  )
  
  out_dirs <- c(Hmean="Hmeans", Hmed="Hmeds", IQR="IQRs", Hsd="Hsds", Cover_above2m="Cover_above2m")   # c() collects key-value pairs into one object
  for (m in names(out_dirs)) writeRaster(metrics_extra[[m]],
                                         paste0("D:/Kiran/Harz_Results/", out_dirs[m], "/", m, "_", toupper(mk), ".tif"),
                                         overwrite = TRUE)
  
}

## Gap fraction (returns below 2 m)
for (mk in names(las_norm_list)) {
  
  gap_frac_mk <- pixel_metrics(
    las_norm_list[[mk]],
    ~mean(Z < 2),
    res = 0.1
  )
  
  writeRaster(
    gap_frac_mk,
    paste0(
      "D:/Kiran/Harz_Results/GapFrac_below2m/GapFrac_below2m_",
      toupper(mk),
      ".tif"
    ),
    overwrite = TRUE
  )
}

rm(las_mk)
gc()

## Magdeburger Huette Site ##

# read las files
las_mgh1 <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250902-145550_Magdeburgerhuette_grndclass.las")
las_mgh2 <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250830-095418_Magdeburgerhuette_grandclass.las")

# CHM settings
thr <- c(0, 2, 5, 10, 15)
edg <- c(0, 1.5)

# IMPORTANT: which LAS to use for each plot (two different LAS files)
las_for_plot <- list(
  mgh_001 = las_mgh2,
  mgh_002 = las_mgh1,
  mgh_003 = las_mgh2,
  mgh_004 = las_mgh1,
  mgh_005 = las_mgh1,
  mgh_006 = las_mgh2,
  mgh_007 = las_mgh2
)

# store normalized LAS for later metrics
las_norm_list <- list()

# loop plots: clip -> normalize (already ground-classified) -> CHM
for (i in 1:7) {
  
  id <- paste0("mgh_00", i)  # mgh_001 ... mgh_007
  poly_mgh <- poly_train[poly_train$id == id, ] #because we have 7 plots in the same polygon, we have to clip it separately
  
  las_roi_mgh      <- clip_roi(las_for_plot[[id]], poly_mgh)
  ##since our two las files are already ground-classified, we do not do classify_ground step for this site.
  las_norm_roi_mgh <- normalize_height(las_roi_mgh, knnidw(k = 10, p = 2))
  
  # save normalized LAS
  las_norm_list[[paste0("mgh", i)]] <- las_norm_roi_mgh
}

# CHM
chm_mgh <- rasterize_canopy(las_norm_roi_mgh, 0.1, pitfree(thr, edg))
writeRaster(
  chm_mgh,
  paste0("D:/Kiran/Harz_Results/CHMs/CHM_MGH", i, ".tif"),
  overwrite = TRUE
)


## LiDAR height percentiles (h95 and h10)
percentiles <- list(h95 = 0.95, h10 = 0.10)

for (p_name in names(percentiles)) {
  for (mgh in names(las_norm_list)) {
    
    result <- pixel_metrics(
      las_norm_list[[mgh]],
      ~quantile(Z, percentiles[[p_name]]),
      res = 0.10
    )
    
    writeRaster(
      result,
      paste0(
        "D:/Kiran/Harz_Results/",
        toupper(p_name), "s/",
        toupper(p_name), "_", toupper(mgh), ".tif"
      ),
      overwrite = TRUE
    )
  }
}

## Additional LiDAR structural metrics (Hmean, Hmedian, Hsd, IQR, Cover > 2m) at 10 cm resolution
for (mgh in names(las_norm_list)) {
  
  metrics_extra <- pixel_metrics(
    las_norm_list[[mgh]],
    ~list(
      Hmean   = mean(Z, na.rm = TRUE),                                # mean height 
      Hmed = median(Z, na.rm = TRUE),                                 # median height (H50) 
      Hsd     = sd(Z, na.rm = TRUE),                                  # height standard deviation 
      IQR     = quantile(Z, 0.75, na.rm = TRUE) -
        quantile(Z, 0.25, na.rm = TRUE),                     # interquartile range
      Cover_above2m = sum(Z > 2, na.rm = TRUE) / sum(Z > 0, na.rm = TRUE)   # canopy cover above 2 m 
    ),
    res = 0.10
  )
  
  out_dirs <- c(Hmean="Hmeans", Hmed="Hmeds", IQR="IQRs", Hsd="Hsds", Cover_above2m="Cover_above2m")   # c() collects key-value pairs into one object
  for (m in names(out_dirs)) writeRaster(metrics_extra[[m]],
                                         paste0("D:/Kiran/Harz_Results/", out_dirs[m], "/", m, "_", toupper(mgh), ".tif"),
                                         overwrite = TRUE)
  
}

## Gap fraction (returns below 2 m)
for (mgh in names(las_norm_list)) {
  
  gap_frac_mgh <- pixel_metrics(
    las_norm_list[[mgh]],
    ~mean(Z < 2),
    res = 0.1
  )
  
  writeRaster(
    gap_frac_mgh,
    paste0(
      "D:/Kiran/Harz_Results/GapFrac_below2m/GapFrac_below2m_",
      toupper(mgh),
      ".tif"
    ),
    overwrite = TRUE
  )
}

rm(las_mgh1, las_mgh2)
gc()

## Quesenbank site ##

# read las files
las_qb <- readLAS("D:/Kiran/Harz/02_LiDAR/YS-20250806-172814_Queesenbank_OutliersRemoved.las")

# CHM settings
thr <- c(0, 2, 5, 10, 15)
edg <- c(0, 1.5)

# store normalized LAS for later metrics
las_norm_list <- list()

# loop plots: clip -> ground classify -> normalize -> CHM
for (i in 1:5) {
  
  id <- paste0("qb_00", i)   # qb_001 ... qb_005
  poly_qb <- poly_train[poly_train$id == id, ] #because we have 5 plots in the same polygon, we have to clip it separately
  
  las_roi_qb      <- clip_roi(las_qb, poly_qb)
  las_qb_grn      <- classify_ground(las_roi_qb, csf())
  las_norm_roi_qb <- normalize_height(las_qb_grn, knnidw(k = 10, p = 2))
  
  # save normalized LAS
  las_norm_list[[paste0("qb", i)]] <- las_norm_roi_qb
}

# CHM
chm_qb <- rasterize_canopy(las_norm_roi_qb, 0.1, pitfree(thr, edg))
writeRaster(
  chm_qb,
  paste0("D:/Kiran/Harz_Results/CHMs/CHM_QB", i, ".tif"),
  overwrite = TRUE
)


## LiDAR height percentiles (h95 and h10)
percentiles <- list(h95 = 0.95, h10 = 0.10)

for (p_name in names(percentiles)) {
  for (qb in names(las_norm_list)) {
    
    result <- pixel_metrics(
      las_norm_list[[qb]],
      ~quantile(Z, percentiles[[p_name]]),
      res = 0.10
    )
    
    writeRaster(
      result,
      paste0(
        "D:/Kiran/Harz_Results/",
        toupper(p_name), "s/",
        toupper(p_name), "_", toupper(qb), ".tif"
      ),
      overwrite = TRUE
    )
  }
}

## Additional LiDAR structural metrics (Hmean, Hmedian, Hsd, IQR, Cover > 2m) at 10 cm resolution
for (qb in names(las_norm_list)) {
  
  metrics_extra <- pixel_metrics(
    las_norm_list[[qb]],
    ~list(
      Hmean   = mean(Z, na.rm = TRUE),                                # mean height 
      Hmed = median(Z, na.rm = TRUE),                                 # median height (H50) 
      Hsd     = sd(Z, na.rm = TRUE),                                  # height standard deviation 
      IQR     = quantile(Z, 0.75, na.rm = TRUE) -
        quantile(Z, 0.25, na.rm = TRUE),                     # interquartile range
      Cover_above2m = sum(Z > 2, na.rm = TRUE) / sum(Z > 0, na.rm = TRUE)   # canopy cover above 2 m 
    ),
    res = 0.10
  )
  
  out_dirs <- c(Hmean="Hmeans", Hmed="Hmeds", IQR="IQRs", Hsd="Hsds", Cover_above2m="Cover_above2m")   # c() collects key-value pairs into one object
  for (m in names(out_dirs)) writeRaster(metrics_extra[[m]],
                                         paste0("D:/Kiran/Harz_Results/", out_dirs[m], "/", m, "_", toupper(qb), ".tif"),
                                         overwrite = TRUE)
  
}

## Gap fraction (returns below 2 m)
for (qb in names(las_norm_list)) {
  
  gap_frac_qb <- pixel_metrics(
    las_norm_list[[qb]],
    ~mean(Z < 2),
    res = 0.1
  )
  
  writeRaster(
    gap_frac_qb,
    paste0(
      "D:/Kiran/Harz_Results/GapFrac_below2m/GapFrac_below2m_",
      toupper(qb),
      ".tif"
    ),
    overwrite = TRUE
  )
}

rm(las_qb)
gc()