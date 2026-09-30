# ==============================================================================
# ROBUST 10-SPECIES SDM PIPELINE
# Output Directory: D:/socket2.0
# ==============================================================================

library(terra)
library(geodata)
library(maxnet)
library(pROC)
library(ggplot2)
library(sf)
library(spocc)

output_dir <- "D:/socket2.0"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n[STEP 1] Downloading and Processing Climate Rasters ...\n")
env_stack_global  <- worldclim_global(var = "bio", res = 10, path = output_dir)
env_ssp245_global <- cmip6_world(model = "CanESM5", ssp = "245", time = "2041-2060", var = "bioc", res = 10, path = output_dir)
env_ssp585_global <- cmip6_world(model = "CanESM5", ssp = "585", time = "2041-2060", var = "bioc", res = 10, path = output_dir)

# Download India spatial boundary to crop rasters (fixes memory & NA issues)
india_extent <- ext(68, 98, 6, 38)
env_stack  <- crop(env_stack_global, india_extent)
env_ssp245 <- crop(env_ssp245_global, india_extent)
env_ssp585 <- crop(env_ssp585_global, india_extent)

names(env_ssp245) <- names(env_stack)
names(env_ssp585) <- names(env_stack)

species_list <- c(
  "Withania somnifera",
  "Lantana camara",
  "Parthenium hysterophorus",
  "Aconitum heterophyllum",
  "Picrorhiza kurroa",
  "Rheum emodi",
  "Pterocarpus santalinus",
  "Garcinia indica",
  "Gymnema sylvestre",
  "Tinospora cordifolia"
)

master_summary_table <- data.frame()

cat("\n[STEP 2] Running Batch SDM Loop ...\n")

for (sp in species_list) {
  
  sp_clean_name <- gsub(" ", "_", sp)
  cat(sprintf("\n--------------------------------------------------\n"))
  cat(sprintf(" Processing Species: %s\n", sp))
  cat(sprintf("--------------------------------------------------\n"))
  
  # A. Safe GBIF Query with Retry Logic
  coords <- NULL
  tryCatch({
    cat(" -> Querying GBIF for occurrence points ...\n")
    gbif_res <- occ(query = sp, from = "gbif", limit = 1500, has_coords = TRUE)
    occ_df <- occ2df(gbif_res)
    
    if (nrow(occ_df) >= 10) {
      coords <- data.frame(
        longitude = as.numeric(occ_df$longitude),
        latitude  = as.numeric(occ_df$latitude)
      )
      coords <- na.omit(coords)
    }
  }, error = function(e) {
    cat(sprintf(" -> GBIF fetch issue for '%s': %s\n", sp, e$message))
  })
  
  if (is.null(coords) || nrow(coords) < 10) {
    cat(sprintf(" -> Skipping %s due to insufficient coordinate data.\n", sp))
    next
  }
  
  # B. Spatial Thinning
  occ_sf <- st_as_sf(coords, coords = c("longitude", "latitude"), crs = 4326)
  cell_ids <- cellFromXY(env_stack, st_coordinates(occ_sf))
  thinned_coords <- coords[!duplicated(cell_ids), ]
  thinned_coords <- na.omit(thinned_coords)
  
  write.csv(thinned_coords, file.path(output_dir, paste0(sp_clean_name, "_thinned_occ.csv")), row.names = FALSE)
  
  # C. Environmental Extraction & Background Sampling
  pres_env <- extract(env_stack, thinned_coords)[, -1]
  pres_env <- na.omit(pres_env)
  
  set.seed(42)
  bg_coords <- spatSample(env_stack, size = 10000, method = "random", xy = TRUE, na.rm = TRUE)[, c("x", "y")]
  bg_env <- extract(env_stack, bg_coords)[, -1]
  bg_env <- na.omit(bg_env)
  
  env_matrix <- rbind(pres_env, bg_env)
  p_vec <- c(rep(1, nrow(pres_env)), rep(0, nrow(bg_env)))
  
  # D. Fit Maxnet
  cat(" -> Fitting Maxnet model ...\n")
  max_mod <- maxnet(p = p_vec, data = env_matrix, f = maxnet.formula(p_vec, env_matrix))
  
  # E. Spatial Predictions (Robust Prediction Handling)
  cat(" -> Predicting spatial suitability ...\n")
  pred_curr   <- predict(env_stack,  max_mod, type = "cloglog", na.rm = TRUE)
  pred_ssp245 <- predict(env_ssp245, max_mod, type = "cloglog", na.rm = TRUE)
  pred_ssp585 <- predict(env_ssp585, max_mod, type = "cloglog", na.rm = TRUE)
  
  writeRaster(pred_curr,   file.path(output_dir, paste0(sp_clean_name, "_current_suitability.tif")), overwrite = TRUE)
  writeRaster(pred_ssp245, file.path(output_dir, paste0(sp_clean_name, "_ssp245_2050.tif")), overwrite = TRUE)
  writeRaster(pred_ssp585, file.path(output_dir, paste0(sp_clean_name, "_ssp585_2050.tif")), overwrite = TRUE)
  
  # F. Metrics Calculation
  cat(" -> Computing evaluation metrics & variable importance ...\n")
  pres_p <- predict(max_mod, env_matrix[p_vec == 1, ], type = "cloglog")
  bg_p   <- predict(max_mod, env_matrix[p_vec == 0, ], type = "cloglog")
  
  roc_o <- roc(c(rep(1, length(pres_p)), rep(0, length(bg_p))), c(pres_p, bg_p), quiet = TRUE)
  auc_v <- as.numeric(auc(roc_o))
  
  coords_r <- coords(roc_o, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  opt_thresh <- as.numeric(coords_r$threshold)
  tss_v <- coords_r$sensitivity + coords_r$specificity - 1
  
  var_names <- colnames(env_matrix)
  imp_df <- data.frame(Variable = var_names, Importance = 0)
  
  for (i in seq_along(var_names)) {
    temp_e <- env_matrix
    temp_e[, var_names[i]] <- sample(temp_e[, var_names[i]])
    perm_p <- predict(max_mod, temp_e, type = "cloglog")
    perm_r <- roc(c(rep(1, length(pres_p)), rep(0, length(bg_p))), perm_p, quiet = TRUE)
    imp_df$Importance[i] <- max(0, auc_v - as.numeric(auc(perm_r)))
  }
  
  if (sum(imp_df$Importance) > 0) {
    imp_df$Percent_Contribution <- (imp_df$Importance / sum(imp_df$Importance)) * 100
  } else {
    imp_df$Percent_Contribution <- 100 / length(var_names)
  }
  imp_df <- imp_df[order(-imp_df$Percent_Contribution), ]
  
  p_imp <- ggplot(imp_df, aes(x = reorder(Variable, Percent_Contribution), y = Percent_Contribution)) +
    geom_bar(stat = "identity", fill = "#2b5c8f", color = "black", width = 0.7) +
    coord_flip() +
    labs(
      title = paste(sp, "- Variable Importance"),
      subtitle = sprintf("AUC = %.3f | TSS = %.3f | Threshold = %.4f", auc_v, tss_v, opt_thresh),
      x = "Bioclimatic Variables", y = "Contribution (%)"
    ) +
    theme_minimal()
  
  ggsave(file.path(output_dir, paste0(sp_clean_name, "_variable_importance.png")), plot = p_imp, width = 7, height = 5, dpi = 300)
  
  # G. Area Change Calculation (km²)
  cat(" -> Quantifying land surface area change ...\n")
  bin_curr   <- pred_curr   >= opt_thresh
  bin_ssp245 <- pred_ssp245 >= opt_thresh
  bin_ssp585 <- pred_ssp585 >= opt_thresh
  
  writeRaster(bin_curr,   file.path(output_dir, paste0(sp_clean_name, "_binary_current.tif")), overwrite = TRUE)
  writeRaster(bin_ssp245, file.path(output_dir, paste0(sp_clean_name, "_binary_ssp245.tif")), overwrite = TRUE)
  writeRaster(bin_ssp585, file.path(output_dir, paste0(sp_clean_name, "_binary_ssp585.tif")), overwrite = TRUE)
  
  cell_areas  <- cellSize(pred_curr, unit = "km")
  area_curr   <- sum(values(bin_curr   * cell_areas), na.rm = TRUE)
  area_ssp245 <- sum(values(bin_ssp245 * cell_areas), na.rm = TRUE)
  area_ssp585 <- sum(values(bin_ssp585 * cell_areas), na.rm = TRUE)
  
  pct_245 <- ((area_ssp245 - area_curr) / area_curr) * 100
  pct_585 <- ((area_ssp585 - area_curr) / area_curr) * 100
  
  # H. Response Curves
  cat(" -> Generating response curves ...\n")
  top_2_vars <- head(imp_df$Variable, 2)
  for (v in top_2_vars) {
    v_seq <- seq(min(env_matrix[, v], na.rm = TRUE), max(env_matrix[, v], na.rm = TRUE), length.out = 100)
    mean_df <- as.data.frame(lapply(as.data.frame(env_matrix), mean, na.rm = TRUE))
    pred_df <- mean_df[rep(1, 100), ]
    pred_df[, v] <- v_seq
    pred_df$suitability <- predict(max_mod, pred_df, type = "cloglog")
    pred_df$var_value <- v_seq
    
    p_resp <- ggplot(pred_df, aes(x = var_value, y = suitability)) +
      geom_line(color = "#1f77b4", size = 1.2) +
      geom_hline(yintercept = opt_thresh, linetype = "dashed", color = "red", alpha = 0.7) +
      labs(title = paste(sp, "- Response Curve:", v), x = v, y = "Habitat Suitability") +
      theme_minimal()
    
    ggsave(file.path(output_dir, paste0(sp_clean_name, "_response_", v, ".png")), plot = p_resp, width = 6, height = 4, dpi = 300)
  }
  
  sp_summary <- data.frame(
    Species = sp,
    AUC = round(auc_v, 3),
    TSS = round(tss_v, 3),
    Threshold_MaxSSS = round(opt_thresh, 4),
    Top_Variable = imp_df$Variable[1],
    Current_Area_sqkm = round(area_curr, 2),
    SSP245_Area_sqkm = round(area_ssp245, 2),
    SSP245_Change_pct = round(pct_245, 2),
    SSP585_Area_sqkm = round(area_ssp585, 2),
    SSP585_Change_pct = round(pct_585, 2)
  )
  
  master_summary_table <- rbind(master_summary_table, sp_summary)
  cat(sprintf(" -> [COMPLETED SUCCESSFULLY] %s\n", sp))
}

write.csv(master_summary_table, file.path(output_dir, "ALL_10_SPECIES_SDM_SUMMARY.csv"), row.names = FALSE)
cat("\n==================================================\n")
cat(" ALL 10 SPECIES PROCESSED! Results saved in D:/socket2.0\n")
cat("==================================================\n")
