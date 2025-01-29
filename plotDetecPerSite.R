# rm(list = ls())
# 1. Loading packages ---- 
library(dplyr)
library(sf)
library(ggplot2)

########################
# 2. Loading data ---- #
########################
# setwd("/Users/santiagogutierrezzapata/Documents/Universidad de Huelva - Doctorado/Papers/Activity patterns") # Change working directory

# Read data
all_data <- read.csv("dataActivityPatterns.csv", header = TRUE)

##############################################
# 3. get detections per sp for each site ---- #
##############################################

# Split the data
data_per_sp <- split(all_data, all_data$sp)

# Function to get detection counts per location
get_detection <- function(data) {
  
  final_data <- list()
  sp <- names(data)
  
  for (i in seq_along(data)) {
    # Access the dataframe for the current species
    data_t <- data[[i]]
    
    # Group and summarise data by site
    data_t <- data_t %>%
      group_by(site) %>%
      summarise(
        site = unique(site),             # Rename 'site' column to 'loc'
        n_det = n(),            # Count number of detections
      )
    
    # Store the summarised data in the final list
    final_data[[sp[i]]] <- data_t
  }
  
  return(final_data)
}

# Apply the function
det_per_loc <- get_detection(data_per_sp)


######################################
# 4. Get total of active days per site --- #
######################################

# Read the operation table
operation_tb <- read.csv("problemsCamtraps.csv")

loc <- as.character(unique(all_data$site))

# Crear un objeto vacío para almacenar los resultados
activity_summary <- data.frame(site = character(), 
                               total_active_days = integer(), 
                               stringsAsFactors = FALSE)

current_loc <- "11"

for (current_loc in loc) {
  
  operation_tb_filtered <- operation_tb %>%
    filter(station == current_loc) %>%
    select_if(~ any(!is.na(.) & . != "")) %>%
    mutate(across(2:ncol(.), ~ as.POSIXct(., format = "%d/%m/%Y %H:%M")))
  
  loc_activity <- seq(from = operation_tb_filtered$setup_date[1], to = operation_tb_filtered$retrieval_date[1], by = "day")
  
  activity_tb <- data.frame(
    datetime = loc_activity,
    activity = 1
  )
  
  problem_from <- operation_tb_filtered[grep("^problem.*_from", names(operation_tb_filtered))]
  problem_to <- operation_tb_filtered[grep("^problem.*_to", names(operation_tb_filtered))]
  
  # Si no hay problemas, el sitio estuvo activo todo el período
  if (length(problem_from) == 0 || length(problem_to) == 0) {
    activity_tb$activity <- 1  # Todos los minutos están activos
  } else {
    # Crear una matriz lógica para marcar inactividad
    inactive_matrix <- matrix(FALSE, nrow = length(loc_activity), ncol = length(problem_from))
    
    # Llenar la matriz de inactividad con TRUE en los periodos correspondientes
    for (i in 1:length(problem_from)) {
      inactive_matrix[, i] <- loc_activity >= problem_from[[i]] & loc_activity <= problem_to[[i]]
    }
    
    # Actualizar la columna 'activity'
    activity_tb$activity <- ifelse(rowSums(inactive_matrix) > 0, 0, 1)
  }
  # Calcular la suma total de días activos
  total_active_days <- sum(activity_tb$activity)
  
  # Almacenar el nombre del sitio y los días activos en un dataframe resumen
  activity_summary <- rbind(activity_summary, data.frame(
    site = current_loc,
    total_active_days = total_active_days
  ))
}


######################################
# 5. Divide detections per active days --- #
######################################

for (i in seq_along(det_per_loc)) {
  
  current_data <- det_per_loc[[i]]
  
  # Verifica que 'loc' exista en ambos data frames antes de hacer el merge
  if ("site" %in% names(current_data) && "site" %in% names(activity_summary)) {
    
    # Combina el data frame actual con activity_summary por la columna 'loc'
    det_per_loc[[i]] <- merge(current_data, activity_summary, by = "site", all.x = TRUE)
    
  } else {
    warning(paste("La columna 'loc' no está presente en el elemento", i))
  }
}

# Create normalized detection column
for (i in seq_along(det_per_loc)) {
  
  det_per_loc[[i]]$freq <- (det_per_loc[[i]]$n_det*1000) / det_per_loc[[i]]$total_active_days
  
}

# Create text column
for (i in seq_along(det_per_loc)) {
  
  det_per_loc[[i]]$det_days <- paste(det_per_loc[[i]]$n_det,"/", det_per_loc[[i]]$total_active_days)
  
}

# Create text column
for (i in seq_along(det_per_loc)) {
  
  det_per_loc[[i]]$text <- as.character(det_per_loc[[i]]$freq)
  
}


# Create log column 
for (i in seq_along(det_per_loc)) {
  
  det_per_loc[[i]]$log <- log10(det_per_loc[[i]]$freq*10000)
  
}

det_per_loc[[1]]

######################################
# 6. add the coordinates the the sites --- #
######################################

sites <- read.csv("allSites.csv")

for (i in seq_along(det_per_loc)) {
  
  # Asegúrate de que `data_per_sp` es una lista o un data frame accesible por índice
  current_data <- det_per_loc[[i]]
  
  # Verifica que 'loc' exista en ambos data frames antes de hacer el merge
  if ("site" %in% names(current_data) && "site" %in% names(sites)) {
    
    # Combina el data frame actual con activity_summary por la columna 'loc'
    det_per_loc[[i]] <- merge(current_data, sites, by = "site", all.x = TRUE)
    
  } else {
    warning(paste("La columna 'loc' no está presente en el elemento", i))
  }
}

##########################
# 7. Create the plot ---#
##########################

# Load the shape files for create the maps
flood_zone <- read_sf(dsn = "map/Flood zone.shp")
non_flood_zone <- read_sf(dsn = "map/Non-flood zone.shp" )

# Genetta genetta
genet <- as.character("Genetta genetta")

# Convert the data into a sf object
genet_sf <- st_as_sf(det_per_loc[[genet]],
                               coords = c("longitud", "latitud"), 
                               crs = 4326)
absence_gen <- sites[!sites$site %in% det_per_loc[[genet]]$site,]
absence_gen_sf <- st_as_sf(absence_gen,
                           coords = c("longitud", "latitud"),
                           crs = 4326)

p1 <- ggplot() + 
  geom_sf(data = flood_zone,
          alpha = 1,
          fill = "#F5F5F5") + 
  geom_sf(data = non_flood_zone,
          alpha = 1,
          fill = "white") +
  geom_sf(data = absence_gen_sf,
          color = "black",
          shape = 3) +
  geom_sf(data = genet_sf,
          aes(size = freq),
          color = "black",
          alpha = 1) + 
  geom_text(data = genet_sf, 
            aes(geometry = geometry, label = round(freq, 3)),
            stat = "sf_coordinates",
            size = 1.5, 
            color = "red") +
  # scale_size_continuous(range = c(1, 5), trans = "sqrt") +
  labs(x = "Latitude",
       y = "Longitude",
       title = genet) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "white"),
        legend.position = "none")

ggsave("newPlots/detectionsSp/detGenet.png", p1, width = 6, height = 6, dpi = 800)

# Vulpes vulpes
fox <- as.character("Vulpes vulpes")

# Convert the data into a sf object
fox_sf <- st_as_sf(det_per_loc[[fox]],
                     coords = c("longitud", "latitud"), 
                     crs = 4326)
absence_fox <- sites[!sites$site %in% det_per_loc[[fox]]$site,]
absence_fox_sf <- st_as_sf(absence_fox,
                           coords = c("longitud", "latitud"),
                           crs = 4326)

p2 <- ggplot() + 
  geom_sf(data = flood_zone,
          alpha = 1,
          fill = "#F5F5F5") + 
  geom_sf(data = non_flood_zone,
          alpha = 1,
          fill = "white") +
  geom_sf(data = absence_fox_sf,
          color = "black",
          shape = 3) +
  geom_sf(data = fox_sf,
          aes(size = freq),
          color = "black",
          alpha = 1) + 
  geom_text(data = fox_sf, 
            aes(geometry = geometry, label = round(freq, 3)),
            stat = "sf_coordinates",
            size = 1.5, 
            color = "red") +
  # scale_size_continuous(range = c(0.1, 2.9), 
  #                       breaks = seq(0.1, 2.9, by = 0.1),
  #                       name = "Log size") +
  labs(x = "Latitude",
       y = "Longitude",
       title = fox) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "white"),
        legend.position = "none")

ggsave("newPlots/detectionsSp/detFox.png", p2, width = 6, height = 6, dpi = 800)

# Meles meles
mel <- as.character("Meles meles")

# Convert the data into a sf object
mel_sf <- st_as_sf(det_per_loc[[mel]],
                   coords = c("longitud", "latitud"), 
                   crs = 4326)
absence_mel <- sites[!sites$site %in% det_per_loc[[mel]]$site,]
absence_mel_sf <- st_as_sf(absence_mel,
                           coords = c("longitud", "latitud"),
                           crs = 4326)

p3 <- ggplot() + 
  geom_sf(data = flood_zone,
          alpha = 1,
          fill = "#F5F5F5") + 
  geom_sf(data = non_flood_zone,
          alpha = 1,
          fill = "white") +
  geom_sf(data = absence_mel_sf,
          color = "black",
          shape = 3) +
  geom_sf(data = mel_sf,
          aes(size = freq),
          color = "black",
          alpha = 1) +
  geom_text(data = mel_sf, 
            aes(geometry = geometry, label = round(freq, 3)),
            stat = "sf_coordinates",
            size = 1.5, 
            color = "red") +
  labs(x = "Latitude",
       y = "Longitude",
       title = mel) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "white"),
        legend.position = "none")

ggsave("newPlots/detectionsSp/detMel.png", p3, width = 6, height = 6, dpi = 800)

# Herpestes ichneumon
her <- as.character("Herpestes ichneumon")

# Convert the data into a sf object
her_sf <- st_as_sf(det_per_loc[[her]],
                   coords = c("longitud", "latitud"), 
                   crs = 4326)
absence_her <- sites[!sites$site %in% det_per_loc[[her]]$site,]
absence_her_sf <- st_as_sf(absence_her,
                           coords = c("longitud", "latitud"),
                           crs = 4326)

p4 <- ggplot() + 
  geom_sf(data = flood_zone,
          alpha = 1,
          fill = "#F5F5F5") + 
  geom_sf(data = non_flood_zone,
          alpha = 1,
          fill = "white") +
  geom_sf(data = absence_her_sf,
          color = "black",
          shape = 3) +
  geom_sf(data = her_sf,
          aes(size = freq),
          color = "black",
          alpha = 1) + 
  geom_text(data = her_sf, 
            aes(geometry = geometry, label = round(freq, 3)),
            stat = "sf_coordinates",
            size = 1.5, 
            color = "red") +
#  scale_size_continuous(range = c(1, 5), trans = "sqrt") +
  labs(x = "Latitude",
       y = "Longitude",
       title = her) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "white"),
        legend.position = "none")

ggsave("newPlots/detectionsSp/detHer.png", p4, width = 6, height = 6, dpi = 800)


# Lynx pardinus
lyn <- as.character("Lynx pardinus")

# Convert the data into a sf object
lyn_sf <- st_as_sf(det_per_loc[[lyn]],
                   coords = c("longitud", "latitud"), 
                   crs = 4326)
absence_lyn <- sites[!sites$site %in% det_per_loc[[lyn]]$site,]
absence_lyn_sf <- st_as_sf(absence_lyn,
                           coords = c("longitud", "latitud"),
                           crs = 4326)

p5 <- ggplot() + 
  geom_sf(data = flood_zone,
          alpha = 1,
          fill = "#F5F5F5") + 
  geom_sf(data = non_flood_zone,
          alpha = 1,
          fill = "white") +
  geom_sf(data = absence_lyn_sf,
          color = "black",
          shape = 3) +
  geom_sf(data = lyn_sf,
          aes(size = freq),
          color = "black",
          alpha = 1) + 
  geom_text(data = lyn_sf, 
            aes(geometry = geometry, label = round(freq, 3)),
            stat = "sf_coordinates",
            size = 1.5, 
            color = "red") +
 # scale_size_continuous(range = c(1, 5), trans = "sqrt") +
  labs(x = "Latitude",
       y = "Longitude",
       title = lyn) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "white"),
        legend.position = "none")

ggsave("newPlots/detectionsSp/detLyn.png", p5, width = 6, height = 6, dpi = 800)


# Lepus granatensis
lep <- as.character("Lepus granatensis")

# Convert the data into a sf object
lep_sf <- st_as_sf(det_per_loc[[lep]],
                   coords = c("longitud", "latitud"), 
                   crs = 4326)
absence_lep <- sites[!sites$site %in% det_per_loc[[lep]]$site,]
absence_lep_sf <- st_as_sf(absence_lep,
                           coords = c("longitud", "latitud"),
                           crs = 4326)

p6 <- ggplot() + 
  geom_sf(data = flood_zone,
          alpha = 1,
          fill = "#F5F5F5") + 
  geom_sf(data = non_flood_zone,
          alpha = 1,
          fill = "white") +
  geom_sf(data = absence_lep_sf,
          color = "black",
          shape = 3) +
  geom_sf(data = lep_sf,
          aes(size = freq),
          color = "black",
          alpha = 1) + 
  geom_text(data = lep_sf, 
            aes(geometry = geometry, label = round(freq, 3)),
            stat = "sf_coordinates",
            size = 1.5, 
            color = "red") +
 # scale_size_continuous(range = c(1, 5), trans = "sqrt") +
  labs(x = "Latitude",
       y = "Longitude",
       title = lep) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "white"),
        legend.position = "none")

ggsave("newPlots/detectionsSp/detLep.png", p6, width = 6, height = 6, dpi = 800)


# Oryctolagus cuniculus
ory <- as.character("Oryctolagus cuniculus")

# Convert the data into a sf object
ory_sf <- st_as_sf(det_per_loc[[ory]],
                   coords = c("longitud", "latitud"), 
                   crs = 4326)
absence_ory <- sites[!sites$site %in% det_per_loc[[ory]]$site,]
absence_ory_sf <- st_as_sf(absence_ory,
                           coords = c("longitud", "latitud"),
                           crs = 4326)

p7 <- ggplot() + 
  geom_sf(data = flood_zone,
          alpha = 1,
          fill = "#F5F5F5") + 
  geom_sf(data = non_flood_zone,
          alpha = 1,
          fill = "white") +
  geom_sf(data = absence_ory_sf,
          color = "black",
          shape = 3) +
  geom_sf(data = ory_sf,
          aes(size = freq),
          color = "black",
          alpha = 1) + 
  geom_text(data = ory_sf, 
            aes(geometry = geometry, label = round(freq, 3)),
            stat = "sf_coordinates",
            size = 1.5, 
            color = "red") +
#  scale_size_continuous(range = c(1, 5), trans = "sqrt") +
  labs(x = "Latitude",
       y = "Longitude",
       title = ory) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(colour = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "white"),
        legend.position = "none")

ggsave("newPlots/detectionsSp/detOry.png", p7, width = 6, height = 6, dpi = 800)

