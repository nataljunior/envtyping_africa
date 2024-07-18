library(stringr)
library(dplyr)


filtered_data <- read.csv('full_ds.csv') 

######################################## filter data ####################################

# filter to keep only data "READY_FOR_ANALYSIS" | MACRO_REGION = ME1, ME2, ME3, ME4
filtered_data <- filtered_data[filtered_data$TRIAL_STATUS %in% c("READY_FOR_ANALYSIS"), ]
filtered_data <- filtered_data[filtered_data$MACRO_REGION %in% c("ME1","ME2", "ME3", "ME4"), ]
filtered_data <- filtered_data[filtered_data$CROP_SEASON %in% c(2018, 2019, 2020, 2021, 2022), ]
filtered_data$GRAIN_YLD_13PCT <- as.numeric(filtered_data$GRAIN_YLD_13PCT)
filtered_data$MAT_DAYS <- as.numeric(filtered_data$MAT_DAYS) 
filtered_data$R8_PL_Height <- as.numeric(filtered_data$R8_PL_Height) 
filtered_data <- subset(filtered_data, GRAIN_YLD_13PCT != 0)
filtered_data <- filtered_data[!is.na(filtered_data$GRAIN_YLD_13PCT), ]

######################################### Transform formats ##############################

df <- filtered_data

# Convert relevant columns to appropriate data types
df <- df %>%
  mutate(env_id = paste0(LOCATION_NAME, "_",
                         str_sub(CROP_SEASON, -2, -1), "_",
                         SEASON_NO)) %>%
  
  dplyr::select(env_id, TRIAL_NAME, SITE_LAT, SITE_LONG, SITE_ALT, MACRO_REGION, Country, REP_NO,
                PLANT_DATE, HARV_DATE, Crop_Variety, GRAIN_YLD_13PCT, MAT_DAYS)


df$SITE_LAT <- as.numeric(df$SITE_LAT) 
df$SITE_LONG <- as.numeric(df$SITE_LONG)
df$PLANT_DATE <- as.Date(df$PLANT_DATE, format = "%m/%d/%Y")
df$HARV_DATE <- as.Date(df$HARV_DATE, format = "%m/%d/%Y")
df$GRAIN_YLD_13PCT <- as.numeric(df$GRAIN_YLD_13PCT)

################################ Cleaning1 ##############################################
## Handling null values
colSums(is.na(df))

#site_lat = 89 null values / site_long 89 null values / Plant_date 225 / harv_date = 1418
df <- na.omit(df) # 22954 - 21357 / WITH MAT_DAYS = 22954 - 20531

##################################### envirotyping ######################################

# Select desired columns
env_t_ds <- df[, c("env_id", "SITE_LAT", "SITE_LONG", 
                   "PLANT_DATE", "HARV_DATE", "Country", "SITE_ALT", "GRAIN_YLD_13PCT")]

# Rename columns for clarity
colnames(env_t_ds) <- c("env_id", "lat", "lon", "start_date", "end_date", 
                        "Country", "alt", "yield")

#write.csv(env_t_ds, "env_t_ds2.csv")

##################################### xxxxxxxxxxxxxxxxxxxxx #############################
#rm(data)

require(EnvRtype)

data <- read.csv("dataset/env_t_ds.csv") # only 5 environments (location x year)

env.i = data$env             # environment ID, I presume that this is yours, right?
lat = data$lat               # latitude coordinates
lon = data$lon               # longitude coordinates
plant.date = data$start      # start, year-month-day
harv.date = data$end         # end, year-month-day


remotes::install_github("gcostaneto/envirotypeR",force=TRUE)


###################################################################################################
# (1) daily weather data & processed variables  ####
###################################################################################################


# add weather data
df_clim = envirotypeR::get_weather(env.id = env.i,
                                   lat = lat,lon = lon,
                                   start.day = plant.date
                                   ,end.day = harv.date,
                                   parallel = T) # parallel = T for a wide number of environments


# add altitude
df_clim <- envirotypeR::get_spatial(digital.raster = terra::rast("https://urldefense.com/v3/__https://raw.githubusercontent.com/gcostaneto/envirotypeR/main/inst/extdata/wc2.1_2.5m_elev.tif__;!!DZ3fjg!4W6ell_Vi6Y5dciPUbqQg2a3tI2uhq7CjjPjNdVj52GFmHpe1ktl7vpICJZRbxFXL8pxS4xdAe4jYHxr7E-W_EjEkA$ "),
                                    env.dataframe = df_clim,lat = 'LAT',lng = 'LON',env.id = 'env',
                                    name.feature = 'elevation', # you decide the name. you can use whatever ytou want to call
                                    merge = T ) # merge = T means that you want to combine it with the original df.clim table

write.csv(df_clim, "df_clim.csv")

# add soil

###################################################################################################
# (3) pulling soil data from SoilGrids  ####
###################################################################################################


# downloading soil data!
soil_variables =  c("bdod", # soil bulk density
                    "cec",  # conductivity eletric capacity
                    "clay",  # clay content
                    "nitrogen",  # nitrogen content
                    "ocd", # organic density
                    "ocs",  # organic content
                    "phh2o",  # ph
                    "sand", # sand content
                    "silt", # sint content
                    "soc") # soil organic carbon (total)

df.soil <- envirotypeR::get_soil(env.id = env.i,lat = lat,lon = lon,
                                 variables.names =  soil_variables )

write.csv(df.soil, "df_soil.csv")


# depth layers: example: 5_15cm from 5cm to 15cm of depth
# statistics: q_5 (quantile 5%), q_95 (quantile 95%), kind of min and max. mean = average value


# add cardinals
df_clim <- df_clim |> EnvRtype::param_temperature(Tbase1 = 8,Tbase2 = 30,Topt1 = 35,Topt2 = 45, merge = T) # remember to adjust the Tbase etc for the crop you are interested





################################################## #############################################

library(tidyr)

window_names <- c('0-20', '21-34', '35-79', '80-104', '+105')
time_windows <- c(0, 21, 35, 80, 105)

# Define variables and time windows
vars <- c("T2M", "T2M_MAX", "T2M_MIN", "T2MDEW", "ALLSKY_SFC_LW_DWN", "ALLSKY_SFC_SW_DWN",
          "ALLSKY_SFC_SW_DNI", "ALLSKY_SFC_PAR_TOT", "ALLSKY_SFC_UVA", "ALLSKY_SFC_UVB", "PRECTOT", "EVPTRNS",
          "QV2M", "RH2M", "GWETROOT", "GWETTOP", "FROST_DAYS", "WS2M",
          "P_ETP", "VPD", "N", "RTA", "n", "TH1",
          "TH2", "PAR_TEMP", "GDD", "FRUE", "T2M_RANGE")


# Ensure df.clim is a data frame
df_clim2 <- as.data.frame(x = df_clim)


# Create a new column with the window names based on the daysFromStart column
df_clim2 <- as.data.frame(df_clim2) %>%
  mutate(window_name = case_when(
    daysFromStart >= time_windows[1] & daysFromStart < time_windows[2] ~ window_names[1],
    daysFromStart >= time_windows[2] & daysFromStart < time_windows[3] ~ window_names[2],
    daysFromStart >= time_windows[3] & daysFromStart < time_windows[4] ~ window_names[3],
    daysFromStart >= time_windows[4] & daysFromStart < time_windows[5] ~ window_names[4],
    daysFromStart >= time_windows[5] ~ window_names[5],
    #daysFromStart >= time_windows[5] & daysFromStart < window_names[6] ~ window_names[5],
    #daysFromStart >= time_windows[6] ~ window_names[6],
    TRUE ~ NA_character_ # Handle cases outside the defined ranges
  ))



calculate_group_quantiles <- function(dataframe, vars, window_names) {
  # Ensure df is a data frame
  dataframe <- as.data.frame(dataframe)
  
  # Group by 'env' and 'window_name' and calculate quantiles for each variable
  group_quantiles <- dataframe %>%
    group_by(env, window_names) %>%
    summarise(across(all_of(vars), 
                     list(q10 = ~quantile(., probs = 0.1, na.rm = TRUE),
                          q50 = ~quantile(., probs = 0.5, na.rm = TRUE),
                          q90 = ~quantile(., probs = 0.9, na.rm = TRUE))), 
              .groups = "drop") %>%
    ungroup()
  
  return(group_quantiles)
}


# Call the function to calculate statistics for each variable grouped by env and window_name
group_stats <- calculate_group_quantiles(df_clim2, vars, window_names)


vars2 <- names(group_stats)[sapply(group_stats, is.numeric)]
#numeric_dataset <- group_quantile_pivot[, numeric_vars]

group_means_pivot <- group_stats %>%
  pivot_wider(names_from = "window_name",
              values_from = all_of(vars2),
              names_glue = "{.value}_{window_name}")


#########################################################################################################

# Scale the numeric columns
scaled_group_stats <- group_means_pivot %>%
  mutate(across(where(is.numeric), 
                ~scale(.) %>% as.vector()))

write.csv(scaled_group_stats, "scaled_group_stats.csv")


#######################################################################################################
### algumas localidade geram valores nulos em localidades que nao vao acima de 105 na col 105+
### vou excluir a coluna e analisar o resultado

# Remove columns with any null values
scaled_group_stats_clean <- scaled_group_stats[ , colSums(is.na(scaled_group_stats)) == 0]

# Check the dimensions of the original and cleaned dataframes
print(paste("Original dataframe dimensions:", dim(scaled_group_stats)[1], "rows,", dim(scaled_group_stats)[2], "columns"))
print(paste("Cleaned dataframe dimensions:", dim(scaled_group_stats_clean)[1], "rows,", dim(scaled_group_stats_clean)[2], "columns"))

# Check the dimensions of the original and cleaned dataframes
print(paste("Original dataframe dimensions:", dim(scaled_group_stats)[1], "rows,", dim(scaled_group_stats)[2], "columns"))
print(paste("Cleaned dataframe dimensions:", dim(scaled_group_stats_clean)[1], "rows,", dim(scaled_group_stats_clean)[2], "columns"))


# Remove rows with null values 
#cleaned_group_stats <- scaled_group_stats %>%
#  filter(!is.na(average_yield))

numerical_data <- scaled_group_stats_clean %>% 
  select(-env) %>%  # Remove the 'env' column
  scale()



################################ Adding yield avg ######################################################
env_t_ds <- data

avg_per_env <- env_t_ds %>%
  group_by(env) %>%
  summarise(average_yield = mean(yield, na.rm = TRUE),
            average_alt = mean(alt, na.rm = TRUE))

yield_merge <- group_means_pivot %>%
  left_join(avg_per_env, by = c("env" = "env_id"))

#######################################################################################################
#### trabalhar melhor essa parte, gera uma tabela estranha sem os env e com num inteiro


library(superheat)

# First, ensure all character columns are factors
cleaned_group_stats <- cleaned_group_stats %>%
  mutate(across(where(is.character), as.factor))

# Now convert all columns to numeric
numeric_group_stats <- cleaned_group_stats %>%
  mutate(across(everything(), ~ as.numeric(as.factor(.))))

# Check the structure of the new dataframe
str(numeric_group_stats)

# If you want to keep the original column names
colnames(numeric_group_stats) <- colnames(cleaned_group_stats)


w <- is.numeric(scaled_group_stats)

superheat(numeric_group_stats,
          pretty.order.rows = T,
          pretty.order.cols = T,
          bottom.label.text.angle = 90,
          bottom.label.text.size = 2,
          yt.axis.name.size = 1)


?superheat


################################ removing correlated columns ###########################

library(caret)
library(dplyr)

# Function to remove highly correlated variables
remove_correlated <- function(data, threshold = 0.99) {
  correlation_matrix <- cor(data)
  highly_correlated <- findCorrelation(correlation_matrix, cutoff = threshold, verbose = FALSE)
  return(data[, -highly_correlated])
}

# Apply the function to your dataset
uncorrelated_data <- remove_correlated(numerical_data, threshold = 0.99)

# Print the number of columns before and after
cat("Number of columns before:", ncol(numerical_data), "\n")
cat("Number of columns after:", ncol(uncorrelated_data), "\n")

# Print the names of removed columns
#removed_cols <- setdiff(names(numeric_group_stats), names(uncorrelated_data))
#cat("Removed columns:", paste(removed_cols, collapse = ", "), "\n")

# from 336 to 280 features

########################################## superheat #########################################

library(superheat)

superheat(uncorrelated_data,
          pretty.order.rows = T,
          pretty.order.cols = T,
          bottom.label.text.angle = 90,
          bottom.label.text.size = 2,
          yt.axis.name.size = 1)


################################ merging weather with grain yield #############################

merged_df <- merge(env_t_ds, scaled_group_stats, by.x = 'env_id', by.y = 'env')

##### Recursive feature selection

library(caret)
library(randomForest)

# Select features (excluding specified columns)
excluded_columns <- c("env_id", "lat", "lon", "start_date", "end_date", "Country", "alt", "yield")
features <- numeric_group_stats[, !names(numeric_group_stats) %in% excluded_columns]

# Set target variable
target <- merged_df$yield

# Ensure all features are numeric
features <- features %>% mutate_all(as.numeric)

# Remove any columns with NA or infinite values
features <- features %>% select_if(~ !any(is.na(.) | is.infinite(.)))

# Set up the control parameters for RFE
ctrl <- rfeControl(functions = rfFuncs,
                   method = "cv",
                   number = 5)  # 5-fold cross-validation

# Perform RFE
set.seed(123)  # for reproducibility
rfe_result <- rfe(x = features, 
                  y = target,
                  sizes = c(1:ncol(features)),
                  rfeControl = ctrl)

# Print the results
print(rfe_result)

# Plot the results
plot(rfe_result, type = c("g", "o"))

# Get the optimal number of features
optimal_n_features <- rfe_result$optsize

# Get the selected features
selected_features <- predictors(rfe_result)

# Print selected features
cat("Number of selected features:", length(selected_features), "\n")
cat("Selected features:", paste(selected_features, collapse = ", "), "\n")

# Create a new dataset with only the selected features and the target variable
selected_data <- cbind(uncorrelated_data[, selected_features], yield = target)


