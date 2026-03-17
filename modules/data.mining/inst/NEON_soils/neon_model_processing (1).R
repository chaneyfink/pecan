
# Function to get model output of NEON sites

# This function will input
# what years of data you want (2000, 2025)
# what model output you want (like "SoilResp") *Needs to be in form that model uses
# And give you a table of that variable that for all the NEON sites with most recent batch data

#' @param start_year numeric in YYYY format: start date of the requested time window.
#' @param end_year numeric in YYYY format: end date of the requested time window.
#' @param var1 character: variable name.
#' @param var2 character: variable name.
#'
#' @return df: a df containing daily means of requested variables for all NEON sites in specified time frame
#' 
# ------------------------------------------------------

model_processing_neonsites <- function(start_year, end_year, var1, var2 = NULL){
  # First, figure out the NEON sites
  
  # get the coordinates of 6400 sites in Dongchen's output
  
  Neonplots <- read.csv("/usr4/ugrad/chaney/R/3_25/Neon_sites_terrestrial.csv") 
  
  # Add a column for site numbers in Site_Info (site_info has coordinates of each location in df)
  Site_Info <- readRDS("/usr4/ugrad/chaney/R/3_25/site.locs.rds")
  Site_Info$Site_Number <- seq(1, nrow(Site_Info)) 
  
  # Extract latitude and longitude from Neonplots and Site_Info tables
  neon_coords <- Neonplots[ ,2:3] %>%
    select(Latitude = field_latitude, Longitude = field_longitude)
  
  Site_Info_coords <- Site_Info %>%
    select(Site_Number, lat, lon)
  
  # Put Site_Info latitude/longitude in a matrix
  Site_Coordinates <- as.matrix(Site_Info[, c("lat", "lon")])
  
  # Same for Neonplots
  Neonplots_Coordinates <- as.matrix(Neonplots[, 2:3])
  
  # make list to store matched list
  matched <- data.frame(id = numeric(nrow(Neonplots)), # id in forecast
                        site_id = numeric(nrow(Neonplots))) # NEON id
  
  # Loop through each row in Neonplots and calculate each plot's closest site number
  for (i in 1:nrow(Neonplots)) {
    # Extract the coordinates for the current row in Neonplots
    plot_coords <- Neonplots_Coordinates[i, ]
    
    # Compute distances from this plot to all site coordinates (longlat means its in km)
    distances <- spDistsN1(Site_Coordinates, plot_coords, longlat = TRUE)
    
    # Find the index of the minimum distance for closest site
    closest_site_index <- which.min(distances)
    
    # Assign the corresponding Site_Number from Site_Info to the current row in Neonplots
    matched$id[i] <- Site_Info$Site_Number[closest_site_index]
    matched$site_id[i] <- Neonplots$field_site_id[i]
  }
  
  
  # Make year vector
  # Each file has all ens members for 1 year
  Years <- start_year:end_year
  
  # Process Ensemble Function
  process_ensemble <- function(year) {
    year_list <- list()
    
    # go through all the sites
    for (i in seq_len(nrow(matched))) {
      site_id <- matched$site_id[i]
      ID <- matched$id[i]
        
        file_path <- paste0(
          "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/SDA_8k_site/NA_SDA_no_Debias", "/Job_", ceiling(ID/200),"/merged_nc/",
          year, ".nc"
        )
        
        # Open file
        nc_data <- tryCatch({
          nc_open(file_path)
        }, error = function(e) {
          return(NULL)
        })
        
        # file open
        if (!is.null(nc_data)) {
          # time dimension ??
          if ("time" %in% names(nc_data$dim)) {
          
              # From dongchen's function
            
              # Extract variable for site for year
            
              time.val <- nc_data$dim$time$vals
              time.unit <- nc_data$dim$time$units
              origin <- strsplit(x = time.unit, split = "since ", fixed = TRUE)[[1]][2]
              real_time <- as.POSIXct(time.val*3600*24, origin = origin, tz = "UTC")
              # Use entire year
              time <- which(real_time >= as.Date(paste0(year, "-01-01")) & real_time <= as.Date(paste0(year, "-12-31")))
              
              # grab ensemble size.
              ensemble.size <- nc_data$dim$ensemble$len
              
              # Site ID
              ind <- which(nc_data$dim$site$vals == ID)
              
              # Get output
              res <- ncdf4::ncvar_get(nc_data, var1, start = c(ind, 1, time[1]), count = c(1, ensemble.size, length(time)))
              
              # Average over columns (ensemble members)
              out1 <- colMeans(res) 
              
              # If second variable
              if (is.null(var2) == FALSE){
                res2 <- ncdf4::ncvar_get(nc_data, var2, start = c(ind, 1, time[1]), count = c(1, ensemble.size, length(time)))
                out2 <- colMeans(res)
              }else{
                out2 <- NA
              }
              
              
              # make data frame
              df <- data.frame(Site_ID = site_id,
                              Time = substr(real_time[time], 1, 10),
                              Output1 = out1,
                              Output2 = out2)
              
              # save the results
              year_list[[length(year_list) + 1]] <- df
            } 
          
          # Close NetCDF File
          nc_close(nc_data)
        }
    }
    
    # If its a list of dfs, return list as df
    if (is.list(year_list)) {
      if (all(sapply(year_list, is.data.frame))) {
        return(rbindlist(year_list, fill = TRUE))
      }
    }
  }
  
  # Apply parallel calculation for all years at the same time
  num_cores <- detectCores() - 1
  results <- mclapply(Years, process_ensemble, mc.cores = num_cores)
  
  # Combine results from all years
  results_df <- rbindlist(results, fill = TRUE)
  
  # Make into data table
  daily_mean_data <- as.data.table(results_df)
  
  # Format Time to character for processing
  daily_mean_data[, Time := as.character(Time)]
  
  # Calculate daily mean (time resolution is 1 day)
  final_daily_data <- daily_mean_data[, .(
    Output1 = mean(Output1, na.rm = TRUE), Output2 = mean(Output2, na.rm = TRUE)), by = .(Time, Site_ID)]
  
  return(final_daily_data)
  
} 
