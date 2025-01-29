# rm(list = ls())
# 1. Loading packages ---- 
library(dplyr)
library(activity)
library(overlap)
library(data.table)

########################
# 2. Loading data ---- #
########################
# setwd("") # Change working directory

# Read data
all_data <- read.csv("dataActivityPatterns.csv", header = TRUE)

# Read the operation table
operation_tb <- read.csv("problemsCamtraps.csv")

# Define the species
vulpes <- as.character("Vulpes vulpes")
meles <-as.character("Meles meles")
genetta <- as.character("Genetta genetta")
herpestes <- as.character("Herpestes ichneumon")
lynx <- as.character("Lynx pardinus")
lepus <- as.character("Lepus granatensis")
oryctolagus <- as.character("Oryctolagus cuniculus")

##########################################################
# Split data to calculate average anchored times per specie
##########################################################
data_per_sp <- split(all_data, all_data$sp)

# Latitude and longitude of a site
lat <- 36.96693
lon <- -6.46713

# Calculate average anchored times
vulpes_solar<- solartime(data_per_sp[[vulpes]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                    format = "%Y-%m-%d %H:%M")
meles_solar <- solartime(data_per_sp[[meles]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                   format = "%Y-%m-%d %H:%M")
genetta_solar <- solartime(data_per_sp[[genetta]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                     format = "%Y-%m-%d %H:%M")
herpestes_solar <- solartime(data_per_sp[[herpestes]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                       format = "%Y-%m-%d %H:%M")
lynx_solar <- solartime(data_per_sp[[lynx]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                  format = "%Y-%m-%d %H:%M")
oryctolagus_solar <- solartime(data_per_sp[[oryctolagus]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                         format = "%Y-%m-%d %H:%M")
lepus_solar <- solartime(data_per_sp[[lepus]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                   format = "%Y-%m-%d %H:%M")

# --- Fits von Mises kernel density to time-of-day data
# Minutes of the day
t <- seq(0, 2*pi, length = 1440)

bw_vulpes <- getBandWidth(vulpes_solar$solar, kmax = 3)
pdf_vulpes <- densityFit(vulpes_solar$solar, t , bw_vulpes)

bw_meles <- getBandWidth(meles_solar$solar, kmax = 3)
pdf_meles <- densityFit(meles_solar$solar, t , bw_meles)

bw_genetta <- getBandWidth(genetta_solar$solar, kmax = 3)
pdf_genetta <- densityFit(genetta_solar$solar, t , bw_genetta)

bw_herpestes <- getBandWidth(herpestes_solar$solar, kmax =3)
pdf_herpestes <- densityFit(herpestes_solar$solar, t , bw_herpestes)

bw_lynx <- getBandWidth(lynx_solar$solar, kmax =3)
pdf_lynx <- densityFit(lynx_solar$solar, t , bw_lynx)

bw_oryctolagus <- getBandWidth(oryctolagus_solar$solar, kmax =3)
pdf_oryctolagus <- densityFit(oryctolagus_solar$solar, t , bw_oryctolagus)

bw_lepus <- getBandWidth(lepus_solar$solar, kmax =3)
pdf_lepus <- densityFit(lepus_solar$solar, t , bw_lepus)


#####################
### --- Real data ###
#####################

# Define a function to obtain real data
get_real_data <- function(spA, spB) {
  
  # Initialize an empty list to store the results
  data_list <- list()
  
  detections <- all_data %>%
    # Group by site
    group_by(site) %>%
    # Count the presence of each species at each site
    summarise(
      spA_present = any(sp == spA),
      spB_present = any(sp == spB),
      .groups = "drop"
    ) %>%
    # Filter sites where both species are present
    filter(spA_present & spB_present)
  
  # Define the sites to be analyzed
  loc <- as.vector(detections$site)
  
  # Create an activity table for each site, tracking activity at 
  # each minute during the study period. This allows us to exclude 
  # inactive periods (e.g., when the camera stopped working).
  for (current_loc in loc) {
    
    # Filter data for the current site and remove empty columns
    operation_tb_filtered <- operation_tb %>%
      filter(station == current_loc) %>%
      select_if(~ any(!is.na(.) & . != "")) %>%
      mutate(across(2:ncol(.), ~ as.POSIXct(., format = "%d/%m/%Y %H:%M")))
    
    # Generate a minute-by-minute sequence from setup date to retrieval date
    loc_activity <- seq(from = operation_tb_filtered$setup_date[1], 
                        to = operation_tb_filtered$retrieval_date[1], 
                        by = "min")
    
    # Create a table to store site activity. By default, we assume the site 
    # was active in all minutes.
    activity_tb <- data.frame(
      datetime = loc_activity,
      activity = 1
    )
    
    # Extract start dates of inactive periods
    problem_from <- operation_tb_filtered[grep("^problem.*_from", names(operation_tb_filtered))]
    # Extract end dates of inactive periods
    problem_to <- operation_tb_filtered[grep("^problem.*_to", names(operation_tb_filtered))]
    
    # If there are no inactive periods, assume the site was active throughout
    if (length(problem_from) == 0 || length(problem_to) == 0) {
      activity_tb$activity <- 1  # All minutes are active
    } else {
      # Create a logical matrix to mark inactive periods
      inactive_matrix <- matrix(FALSE, nrow = length(loc_activity), ncol = length(problem_from))
      
      # Each minute in the study period (loc_activity) is compared with 
      # the site's inactive periods.
      for (i in 1:length(problem_from)) {
        inactive_matrix[, i] <- loc_activity >= problem_from[[i]] & loc_activity <= problem_to[[i]]
      }
      
      # Update the 'activity' column. If any row has a value greater than 0 (TRUE),
      # it means the site was inactive (0). Otherwise, it was active (1).
      activity_tb$activity <- ifelse(rowSums(inactive_matrix) > 0, 0, 1)
    }
    
    # Filter real data for the current site and the specified species (spA and spB).
    # Remove seconds from the 'datetime' column to allow merging with the activity table.
    presence <- all_data %>%
      mutate(datetime = as.POSIXct(clock_date_time, 
                                   format = "%Y-%m-%d %H:%M")) %>%
      filter(site %in% current_loc,
             sp %in% c(spA, spB)) %>%
      select(datetime, site, sp)
    
    # Merge the activity table with species detections
    tabla <- merge(activity_tb, presence, by = "datetime", all.x = TRUE)
    
    # Create a list to store information for the next iteration
    listino <- vector("list", 1)
    
    for (i in 1:1) {
      tabla <- tabla %>%
        mutate(
          # Separate species detections into different columns
          sp_A = case_when(
            sp == spA ~ 1,
            activity == 1 & is.na(sp) ~ 0,
            sp == spB ~ 0,
            TRUE ~ NA_real_
          ),
          sp_B = case_when(
            sp == spB ~ 1,
            activity == 1 & is.na(sp) ~ 0,
            sp == spA ~ 0,
            TRUE ~ NA_real_
          )
        ) %>%
        # Filter to keep species detections and inactivity periods
        filter(spA == 1 | spB == 1 | is.na(spA) | is.na(spB)) %>%
        # Define interval types. If a species detection is followed by an inactive 
        # period, it is NA. Also, colaculate 'time_diff', which stores the time 
        # difference between detections of different species.
        mutate(
          interval = case_when(
            sp_A == 1 & lead(sp_B) == 1 ~ "AB",
            sp_B == 1 & lead(sp_A) == 1 ~ "BA",
            TRUE ~ NA_character_
          ),
          time_diff = as.numeric(abs(difftime(datetime, lead(datetime), units = "mins")))
        ) %>%
        # Keep only AB and BA intervals
        filter(!is.na(interval))
      
      # Store the result in the list
      listino[[i]] <- tabla
    }
    
    # Store results for each site in the main list
    data_list[[current_loc]] <- listino
  }
  return(data_list)
}


# Run function to get the real data for Rabbit vs Predators 
data_vul_ory <- get_real_data(spA = vulpes, spB = oryctolagus)
data_mel_ory <- get_real_data(spA = meles, spB = oryctolagus)
data_gen_ory <- get_real_data(spA = genetta, spB = oryctolagus)
data_her_ory <- get_real_data(spA = herpestes, spB = oryctolagus)
data_lyn_ory <- get_real_data(spA = lynx, spB = oryctolagus)

# Run function to get the real data for Hare vs Predators
data_vul_lep <- get_real_data(spA = vulpes, spB = lepus)
data_mel_lep <- get_real_data(spA = meles, spB = lepus)
data_gen_lep <- get_real_data(spA = genetta, spB = lepus)
data_her_lep <- get_real_data(spA = herpestes, spB = lepus)
data_lyn_lep <- get_real_data(spA = lynx, spB = lepus)

# In the previous function, we did not account for the time between detections when the camera was not operational.
# As a result, some locations lack time intervals. To address this, we will include the median time to encounter,
# calculated using the minimum (e.g. from detection spA to start inactive period) and maximum 
# (e.g. from detection spA to detection spB) encounter times

# Create a function to update the real data by incorporating the median time to encounter.
include_times <- function(data_list, spA, spB) {
  
  data <- list()
  
  for (i in seq_along(data_list)) {
    
    listino <- vector("list", length(data_list[[i]]))
    
    for (j in seq_along(data_list[[i]])) {
      dataframe <- data_list[[i]][[j]]
      
      # Check if the dataframe is empty. We will only modify sites that lack time intervals
      # due to a period of inactivity.
      if (nrow(dataframe) == 0) {
        # Obtain the time to encounter without considering the inactivity period.
        # We assume this represents the maximum possible time between detections.
        excluded_sites <- all_data %>%
          filter(sp %in% c(spA, spB)) %>%
          arrange(site, clock_date_time) %>%
          group_by(site) %>%
          mutate(
            interval = case_when(
              sp == spA & lead(sp) == spB ~ "AB",
              sp == spB & lead(sp) == spA ~ "BA",
              sp == spA & lead(sp) == spA ~ "AA",
              sp == spB & lead(sp) == spB ~ "BB",
              TRUE ~ NA_character_
            ),
            time_diff = abs(difftime(clock_date_time, lead(clock_date_time), units = "mins"))
          ) %>%
          filter(interval %in% c("AB", "BA"), site == names(data_list[i])) %>%
          select(site, clock_date_time, sp, interval, time_diff)
        
        # Filter the activity table for the site being processed.
        problems_tb <- operation_tb %>%
          filter(station == names(data_list[i])) %>%
          select_if(~ any(!is.na(.) & . != "")) %>%
          mutate(across(2:ncol(.), ~ as.POSIXct(., format = "%d/%m/%Y %H:%M")))
        
        # Extract the start time of the inactivity period that caused the exclusion
        # of the only available time interval for this site.
        # We assume this represents the minimum possible time between detections.
        p_from <- problems_tb[grep("^problem.*_from", names(problems_tb))]
        p_from <- as.data.frame(p_from[, p_from > excluded_sites$clock_date_time])
        p_from <- p_from[, 1]
        
        # Calculate the time difference between the first species detection and
        # the start of the inactivity period.
        difference <- abs(difftime(excluded_sites$clock_date_time, p_from, units = "mins"))
        
        # Compute the average of the minimum and maximum times between detections.
        m_difference <- mean(c(as.numeric(difference), as.numeric(excluded_sites$time_diff)), na.rm = TRUE)
        
        # Update the table that lacked time intervals.
        excluded_sites <- excluded_sites %>% mutate(time_diff = m_difference)
        dataframe <- excluded_sites
      } 
      
      # If a site already has at least 1 interval, skip the previous step.
      listino[[j]] <- dataframe
    }
    
    # Update the main list.
    data[[names(data_list[i])]] <- listino
  }
  
  return(data)
}


# Run function to update the real data for Rabbit vs Predators 
data_vul_ory <- include_times(data_vul_ory, spA = vulpes, spB = oryctolagus)
data_mel_ory <- include_times(data_mel_ory, spA = meles, spB = oryctolagus)
data_gen_ory <- include_times(data_gen_ory, spA = genetta, spB = oryctolagus)
data_her_ory <- include_times(data_her_ory, spA = herpestes, spB = oryctolagus)
data_lyn_ory <- include_times(data_lyn_ory, spA = lynx, spB = oryctolagus)

# Run function to update the real data for Hare vs Predators
data_vul_lep <- include_times(data_vul_lep, spA = vulpes, spB = lepus)
data_mel_lep <- include_times(data_mel_lep, spA = meles, spB = lepus)
data_gen_lep <- include_times(data_gen_lep, spA = genetta, spB = lepus)
data_her_lep <- include_times(data_her_lep, spA = herpestes, spB = lepus)
data_lyn_lep <- include_times(data_lyn_lep, spA = lynx, spB = lepus)


######################
### --- Simulation ###
######################
n.rep <- 100

# Define function for get the simulated data
create_simulation <- function(pdfA, pdfB, rep, spA, spB) {
  
  data_list <- list()
  
  # Realizamos el mismo proceso que con los datos reales para obtener la tabla de 
  # actividad de cada sitio para simular el mismo proceso y descartar las detecciones
  # que se vean interrumpidas por un periodo de inactividad.
  
  detections <- all_data %>%
    # Agrupar por sitio
    group_by(site) %>%
    # Contar la presencia de cada especie en cada sitio
    summarise(
      spA_present = any(sp == spA),
      spB_present = any(sp == spB),
      .groups = "drop"
    ) %>%
    # Filtrar los sitios donde ambas especies están presentes
    filter(spA_present & spB_present)
  
  # Definir el sitio que se va a analizar
  loc <- as.vector(detections$site)
  
  for (current_loc in loc) {
    
    # Filtramos por sitio y excluimos las columnas que no tienen datos
    operation_tb_filtered <- operation_tb %>%
      filter(station == current_loc) %>%
      select_if(~ any(!is.na(.) & . != "")) %>%
      mutate(across(2:ncol(.), ~ as.POSIXct(., format = "%d/%m/%Y %H:%M")))
    
    # Creamos una secuencia por minuto desde la setup date hasta la retrival date
    # de ese sitio
    loc_activity <- seq(from = operation_tb_filtered$setup_date[1], to = operation_tb_filtered$retrieval_date[1], by = "min")
    
    # Creamos la tabla donde almacenaremos la operatividad de los sitios. Asumimos
    # que el sitio estuvo activo en todos los minutos.
    activity_tb <- data.frame(
      datetime = loc_activity,
      activity = 1,
      time = format(loc_activity, "%H:%M")
    )
    
    # Extraemos las fechas de inicio de los periodos de inactividad
    problem_from <- operation_tb_filtered[grep("^problem.*_from", names(operation_tb_filtered))]
    # Extraemos las fechas de fin de los periodos de inactividad
    problem_to <- operation_tb_filtered[grep("^problem.*_to", names(operation_tb_filtered))]
    
    # Si no hay problemas, el sitio estuvo activo todo el período
    if (length(problem_from) == 0 || length(problem_to) == 0) {
      activity_tb$activity <- 1  # Todos los minutos están activos
    } else {
      # Crear una matriz lógica para marcar inactividad
      inactive_matrix <- matrix(FALSE, nrow = length(loc_activity), ncol = length(problem_from))
      
      # Llenar la matriz de inactividad con TRUE en los periodos correspondientes.
      # Comparamos cada minuto del periodo de estudio del sitio (loc_activity) con
      # los problemas del sitio.
      for (i in 1:length(problem_from)) {
        inactive_matrix[, i] <- loc_activity >= problem_from[[i]] & loc_activity <= problem_to[[i]]
      }
      
      # Actualizar la columna 'activity'. Como usamos valores lógicos en la iteración
      # anterior, si una fila es mayor a 1 (TRUE) significa inactividad (0), si no lo
      # es significa activa (1)
      activity_tb$activity <- ifelse(rowSums(inactive_matrix) > 0, 0, 1)
    }
    
    
    # Acá estoy obtenido la probabilidad de detección de cada especie en cada sitio
    # a partir del número de detecciónes dividido por el número de días en los que
    # el sitio estuvo activo. Luego escalo la esa probabilidad usando su patron de
    # actividad con el fin de obtener la probabilidad de la especie para cada minuto del
    # día.
    
    # Days of activity of the site
    days <- nrow(activity_tb[activity_tb$activity == 1, ]) / 1440
    # Detection probability of the spA
    p.spA <- as.numeric(nrow(all_data[all_data$sp == spA & all_data$site == current_loc, ]) / days)
    # Detection probability of the spB
    p.spB <- as.numeric(nrow(all_data[all_data$sp == spB & all_data$site == current_loc, ]) / days)
    # Detectionn probability of the spA scaled using its activity pattern
    p_spA_scaled <- pdfA * (p.spA / sum(pdfA))
    # Detection probability of the spB scaled using its activity pattern
    p_spB_scaled <- pdfB * (p.spB / sum(pdfB))
    
    p.species <- data.table(time = format(seq(from = as.POSIXct("00:00", format = "%H:%M"), 
                                              to = as.POSIXct("23:59", format = "%H:%M"), 
                                              by = "1 min"), "%H:%M"),
                            p.time_spA = p_spA_scaled,
                            p.time_spB = p_spB_scaled
    )
    
    
    # # En simulation_tb tengo las siguientes columnas:
    # # - datetime: secuencia de tiempo por minuto del sitio X durante el periodo de estudio.
    # # - activity: indica si la cámara estuvo activa en un minuto dado.
    # # - p.time_spA: probabilidad de detección del depredador en cada minuto del día.
    # # - p.time_spB: probabilidad de detección de la presa en cada minuto del día.
    # 
    listino <- vector("list", rep) # Lista para almacenar las simulaciones.
    
    for (i in 1:rep) {

      tabla <- simulation_tb %>%
        # Ordeno por fecha y hora
        arrange(datetime) %>%
        mutate(
          # Creo dos nuevas columnas en las que simulo la presencia de la especie A y B
          # unicamente en los minutos que la cámara estuvo activa. El n() en rbinom 
          # indica el tamaño del dataset.
          spA = ifelse(activity == 1, rbinom(n = n(), size = 1, prob = p.time_spA), NA) ,
          spB = ifelse(activity == 1, rbinom(n = n(), size = 1, prob = p.time_spB), NA),
        ) %>%
        # Retengo las filas unicamente con presencia de especies y los NA, ya que 
        # indican inactividad de la cámara
        filter(spA == 1 | spB == 1 | is.na(spA) | is.na(spB)) %>%
        mutate(
          # Acá creo los intervalos.
          interval = case_when(
            spA == 1 & lead(spB) == 1 ~ "AB", # Especie A seguida por especie B.
            spB == 1 & lead(spA) == 1 ~ "BA", # Especie B seguida por especie A.
            # Si las dos especies están presentes en el mismo minuto. 
            spA == 1 & spB == 1 ~ "ERROR", # Ambas especies detectadas en el mismo minuto
            # Descartar intervalos BB, AA o evitar calcular intervalos cuando cualquiera de
            # las especies están presentes y luego hay un periodo de inactividad. El NA_character_
            # es necesario por ...
            TRUE ~ NA_character_ # Descartar intervalos AA, BB e intervalos AB o BA que no 
            # puedo usar por inactividad de la cámara
          ),
          # Calculo la diferencia de tiempo entre detecciones consecutivas.
          time_diff = as.numeric(abs(difftime(datetime, lead(datetime), units = "mins")))
        ) %>%
        # Descarto los intervalos que no me interesan
        filter(!is.na(interval))
      
      # Me aseguré que el código estuviera usando correctamente las probabilidades de cada minuto. 
      # Lo comprobé visualizando del patrón de actividad con las detecciones simuladas. 
      # También comprobé que el número de detecciones simuladas para cada especie 
      # sea coherente con las detecciones reales.
      
      listino[[i]] <- tabla
    }
    
    data_list[[current_loc]] <- listino
  }
  
  return(data_list)  
}

# Define function to generate the simulated data
create_simulation <- function(pdfA, pdfB, rep, spA, spB) {
  
  data_list <- list()
  
  # We follow the same process as with the real data to obtain the activity table  
  # for each site, simulating the same conditions and discarding intervals  
  # that are interrupted due to inactivity periods.
  
  detections <- all_data %>%
    # Group by site
    group_by(site) %>%
    # Count the presence of each species at each site
    summarise(
      spA_present = any(sp == spA),
      spB_present = any(sp == spB),
      .groups = "drop"
    ) %>%
    # Filter sites where both species are present
    filter(spA_present & spB_present)
  
  # Define the site to be analyzed
  loc <- as.vector(detections$site)
  
  for (current_loc in loc) {
    
    # Filter by site and remove empty columns
    operation_tb_filtered <- operation_tb %>%
      filter(station == current_loc) %>%
      select_if(~ any(!is.na(.) & . != "")) %>%
      mutate(across(2:ncol(.), ~ as.POSIXct(., format = "%d/%m/%Y %H:%M")))
    
    # Create a minute-by-minute sequence from the setup date to the retrieval date
    loc_activity <- seq(from = operation_tb_filtered$setup_date[1], to = operation_tb_filtered$retrieval_date[1], by = "min")
    
    # Create the table to store site activity. We assume the site was active in all minutes.
    activity_tb <- data.frame(
      datetime = loc_activity,
      activity = 1,
      time = format(loc_activity, "%H:%M")
    )
    
    # Extract start times of inactivity periods
    problem_from <- operation_tb_filtered[grep("^problem.*_from", names(operation_tb_filtered))]
    # Extract end times of inactivity periods
    problem_to <- operation_tb_filtered[grep("^problem.*_to", names(operation_tb_filtered))]
    
    # If there are no inactivity periods, assume the site was fully active
    if (length(problem_from) == 0 || length(problem_to) == 0) {
      activity_tb$activity <- 1  # All minutes are active
    } else {
      # Create a logical matrix to mark inactivity
      inactive_matrix <- matrix(FALSE, nrow = length(loc_activity), ncol = length(problem_from))
      
      # Each minute in the study period (loc_activity) is compared with 
      # the site's inactive periods.
      for (i in 1:length(problem_from)) {
        inactive_matrix[, i] <- loc_activity >= problem_from[[i]] & loc_activity <= problem_to[[i]]
      }
      
      # Update the 'activity' column. If any row has a value greater than 0 (TRUE),
      # it means the site was inactive (0). Otherwise, it was active (1).
      activity_tb$activity <- ifelse(rowSums(inactive_matrix) > 0, 0, 1)
    }
    
    # Here, we obtain the detection probability of each species at each site  
    # based on the number of detections divided by the number of days the site was active.  
    # Then, we scale that probability using the species' activity pattern  
    # to obtain the probability for each minute of the day.
    
    # Days the site was active
    days <- nrow(activity_tb[activity_tb$activity == 1, ]) / 1440
    # Detection probability for spA
    p.spA <- as.numeric(nrow(all_data[all_data$sp == spA & all_data$site == current_loc, ]) / days)
    # Detection probability for spB
    p.spB <- as.numeric(nrow(all_data[all_data$sp == spB & all_data$site == current_loc, ]) / days)
    # Detection probability for spA scaled using its activity pattern
    p_spA_scaled <- pdfA * (p.spA / sum(pdfA))
    # Detection probability for spB scaled using its activity pattern
    p_spB_scaled <- pdfB * (p.spB / sum(pdfB))
    
    p.species <- data.table(time = format(seq(from = as.POSIXct("00:00", format = "%H:%M"), 
                                              to = as.POSIXct("23:59", format = "%H:%M"), 
                                              by = "1 min"), "%H:%M"),
                            p.time_spA = p_spA_scaled,
                            p.time_spB = p_spB_scaled
    )
    
    
    simulation_tb <- merge(activity_tb, p.species, by = "time", all.x = TRUE)
    
    # The 'simulation_tb' table contains the following columns:
    # - datetime: minute-by-minute time sequence for site X during the study period.
    # - activity: indicates whether the camera was active at a given minute.
    # - p.time_spA: predator detection probability at each minute of the day.
    # - p.time_spB: prey detection probability at each minute of the day.
    
    listino <- vector("list", rep) # List to store simulations.
    
    for (i in 1:rep) {
      
      tabla <- simulation_tb %>%
        # Sort by date and time
        arrange(datetime) %>%
        mutate(
          # Create two new columns simulating the presence of species A and B  
          # only in minutes when the camera was active. The n() in rbinom  
          # determines the dataset size.
          spA = ifelse(activity == 1, rbinom(n = n(), size = 1, prob = p.time_spA), NA) ,
          spB = ifelse(activity == 1, rbinom(n = n(), size = 1, prob = p.time_spB), NA),
        ) %>%
        # Retain only rows with species presence or NA (indicating camera inactivity)
        filter(spA == 1 | spB == 1 | is.na(spA) | is.na(spB)) %>%
        mutate(
          # Create encounter intervals.
          interval = case_when(
            spA == 1 & lead(spB) == 1 ~ "AB", # Species A followed by species B.
            spB == 1 & lead(spA) == 1 ~ "BA", # Species B followed by species A.
            # If both species are detected at the same minute.
            spA == 1 & spB == 1 ~ "ERROR", # Both species detected at the same minute
            # Discard BB and AA intervals or cases where an interval cannot  
            # be calculated due to camera inactivity. The NA_character_  
            # ensures proper handling of these cases.
            TRUE ~ NA_character_ # Discard AA, BB, and invalid AB/BA intervals due to inactivity
          ),
          # Calculate time difference between consecutive detections.
          time_diff = as.numeric(abs(difftime(datetime, lead(datetime), units = "mins")))
        ) %>%
        # Remove unwanted intervals
        filter(!is.na(interval))
      
      # Ensure that the code correctly uses probabilities for each minute.  
      # This was verified by comparing the activity pattern with simulated detections  
      # and checking that the number of simulated detections per species  
      # aligns with real detections.
      
      listino[[i]] <- tabla
    }
    
    data_list[[current_loc]] <- listino
  }
  
  return(data_list)  
}


# Run the simulation Rabbit vs Predators
sim_vul_ory <- create_simulation(spA = vulpes, spB = oryctolagus, pdfA = pdf_vulpes,
                                 pdfB = pdf_oryctolagus, rep = n.rep)
sim_mel_ory <- create_simulation(spA = meles, spB = oryctolagus, pdfA = pdf_meles,
                                 pdfB = pdf_oryctolagus, rep = n.rep)
sim_gen_ory <- create_simulation(spA = genetta, spB = oryctolagus, pdfA = pdf_genetta,
                                 pdfB = pdf_oryctolagus, rep = n.rep)
sim_her_ory <- create_simulation(spA = herpestes, spB = oryctolagus, pdfA = pdf_herpestes,
                                 pdfB = pdf_oryctolagus, rep = n.rep)
sim_lyn_ory <- create_simulation(spA = lynx, spB = oryctolagus, pdfA = pdf_lynx,
                                 pdfB = pdf_oryctolagus, rep = n.rep)

# Run the simulation Hare vs Predators
sim_vul_lep <- create_simulation(spA = vulpes, spB = lepus, pdfA = pdf_vulpes,
                                 pdfB = pdf_lepus, rep = n.rep)
sim_mel_lep <- create_simulation(spA = meles, spB = lepus, pdfA = pdf_meles,
                                 pdfB = pdf_lepus, rep = n.rep)
sim_gen_lep <- create_simulation(spA = genetta, spB = lepus, pdfA = pdf_genetta,
                                 pdfB = pdf_lepus, rep = n.rep)
sim_her_lep <- create_simulation(spA = herpestes, spB = lepus, pdfA = pdf_herpestes,
                                 pdfB = pdf_lepus, rep = n.rep)
sim_lyn_lep <- create_simulation(spA = lynx, spB = lepus, pdfA = pdf_lynx,
                                 pdfB = pdf_lepus, rep = n.rep)

###################
### Comparison ####
###################

# Create a function to generate a matrix containing both real and simulated data, 
# distinguishing between interaction intervals (AB, BA). The comparison is made 
# by calculating the average time of AB or BA intervals for each site and comparing 
# it with the mean time of AB or BA intervals for each site across all simulation repetitions.

data_matrix <- function(data_list, t_interval) {
  # Initialize the matrix
  data_matrix <- matrix(FALSE, nrow = length(data_list), ncol = length(data_list[[1]]))
  
  site_names <- names(data_list)
  
  # Fill the matrix
  for (i in seq_along(data_list)) {
    for (j in seq_along(data_list[[i]])) {
      dataframe <- data_list[[i]][[j]]
      
      dataframe <- dataframe %>%
        filter(interval %in% t_interval)
      
      # Check if the dataframe is empty
      if (nrow(dataframe) == 0) {
        data_matrix[i, j] <- NA  # If no interaction intervals occur in a simulation
      } else {
        # Calculate the mean time difference
        data_matrix[i, j] <- mean(dataframe$time_diff, na.rm = TRUE)
      }
    }
  }
  # Assign site names
  if (!is.null(site_names)) {
    rownames(data_matrix) <- site_names
  } else {
    warning("Problem: site names are missing.")
  }
  return(data_matrix)
}

######################################################
######.     Get matrix with the real data       ######
######################################################

# Get real data for intervals AB and Rabbit vs predators
mel_ory_AB <- data_matrix(data_mel_ory, "AB")
vul_ory_AB <- data_matrix(data_vul_ory, "AB")
gen_ory_AB <- data_matrix(data_gen_ory, "AB")
her_ory_AB <- data_matrix(data_her_ory, "AB")
lyn_ory_AB <- data_matrix(data_lyn_ory, "AB")

# Get real data for intervals BA and Rabbit vs predators
mel_ory_BA <- data_matrix(data_mel_ory, "BA")
vul_ory_BA <- data_matrix(data_vul_ory, "BA")
gen_ory_BA <- data_matrix(data_gen_ory, "BA")
her_ory_BA <- data_matrix(data_her_ory, "BA")
lyn_ory_BA <- data_matrix(data_lyn_ory, "BA")

# Get real data for intervals AB and Hare vs predators
mel_lep_AB <- data_matrix(data_mel_lep, "AB")
vul_lep_AB <- data_matrix(data_vul_lep, "AB")
gen_lep_AB <- data_matrix(data_gen_lep, "AB")
her_lep_AB <- data_matrix(data_her_lep, "AB")
lyn_lep_AB <- data_matrix(data_lyn_lep, "AB")

# Get real data for intervals BA and Hare vs predators
mel_lep_BA <- data_matrix(data_mel_lep, "BA")
vul_lep_BA <- data_matrix(data_vul_lep, "BA")
gen_lep_BA <- data_matrix(data_gen_lep, "BA")
her_lep_BA <- data_matrix(data_her_lep, "BA")
lyn_lep_BA <- data_matrix(data_lyn_lep, "BA")


#################################################
######.     Get simulated data matrix       #####
#################################################

# Get simulated data for intervals AB and Rabbit vs predators
s_mel_ory_AB <- data_matrix(sim_mel_ory, "AB")
s_vul_ory_AB <- data_matrix(sim_vul_ory, "AB")
s_gen_ory_AB <- data_matrix(sim_gen_ory, "AB")
s_her_ory_AB <- data_matrix(sim_her_ory, "AB")
s_lyn_ory_AB <- data_matrix(sim_lyn_ory, "AB")

# Get simulated data for intervals BA and Rabbit vs predators
s_mel_ory_BA <- data_matrix(sim_mel_ory, "BA")
s_vul_ory_BA <- data_matrix(sim_vul_ory, "BA")
s_gen_ory_BA <- data_matrix(sim_gen_ory, "BA")
s_her_ory_BA <- data_matrix(sim_her_ory, "BA")
s_lyn_ory_BA <- data_matrix(sim_lyn_ory, "BA")

# Get simulated data for intervals AB and Hare vs predators
s_mel_lep_AB <- data_matrix(sim_mel_lep, "AB")
s_vul_lep_AB <- data_matrix(sim_vul_lep, "AB")
s_gen_lep_AB <- data_matrix(sim_gen_lep, "AB")
s_her_lep_AB <- data_matrix(sim_her_lep, "AB")
s_lyn_lep_AB <- data_matrix(sim_lyn_lep, "AB")

# Get simulated data for intervals BA and Hare vs predators
s_mel_lep_BA <- data_matrix(sim_mel_lep, "BA")
s_vul_lep_BA <- data_matrix(sim_vul_lep, "BA")
s_gen_lep_BA <- data_matrix(sim_gen_lep, "BA")
s_her_lep_BA <- data_matrix(sim_her_lep, "BA")
s_lyn_lep_BA <- data_matrix(sim_lyn_lep, "BA")


###############################################

# Crear una función para comparar los datos reales con los datos simulados. Acá
# estamos teniendo en cuenta las siguientes posibilidades que se pueden dar
# (tomamos como ejemplo AB):

# Reales | Simulados
# AB     | AB        <- Lo que de el resultado
# AB     | NA        <- Simulados > Reales
# NA     | AB        <- Reales > Simulados
# NA     | NA        <- NA 

comparison <- function(real_data, simulated_data, specie, t_interval) {
  
  comparison_matrix <- matrix(FALSE, nrow = nrow(real_data), ncol = n.rep)
  
  rownames(comparison_matrix) <- rownames(real_data)
  
  for (i in 1:nrow(real_data)) {
    comparison_matrix[i, ] <- ifelse(
      is.na(real_data[i, ]) & is.na(simulated_data[i, ]), NA,
      ifelse(
        !is.na(real_data[i, ]) & is.na(simulated_data[i, ]), FALSE,
        ifelse(
          is.na(real_data[i, ]) & !is.na(simulated_data[i, ]), TRUE,
          ifelse(
            !is.na(real_data[i, ]) & !is.na(simulated_data[i, ]),
            real_data[i, ] >= simulated_data[i, ],
            NA
          )
        )
      )
    )
  }
  
# Calcular proporción de TRUE
n_total <- sum(!is.na(comparison_matrix))
# Exvluimos los NA de la comparación
n_true <- sum(comparison_matrix, na.rm = TRUE)

# Calcular intervalo de confianza y prueba de hipótesis
ic_result <- prop.test(n_true, n_total)

# Crear un dataframe inicial con los resultados
results <- data.frame(
  sp = specie,
  interval = t_interval,
  proportion = ic_result$estimate,
  lower_CI = ic_result$conf.int[1],
  upper_CI = ic_result$conf.int[2],
  p_value = ic_result$p.value
)

return(results)

}
# Create a function to compare real data with simulated data.  
# The function accounts for the following possible cases (using AB as an example):  
#  
# Is the real data greater than the simulated data?
#
# Real Data | Simulated Data  
# --------- | --------------  
# AB        | AB             <- result of the comparison  
# AB        | NA             <- Simulated > Real  
# NA        | AB             <- Real > Simulated  
# NA        | NA             <- NA (no comparison possible)  

comparison <- function(real_data, simulated_data, specie, t_interval) {
  
  comparison_matrix <- matrix(FALSE, nrow = nrow(real_data), ncol = n.rep)
  
  rownames(comparison_matrix) <- rownames(real_data)
  
  # Fill the comparison matrix based on the defined conditions
  for (i in 1:nrow(real_data)) {
    comparison_matrix[i, ] <- ifelse(
      is.na(real_data[i, ]) & is.na(simulated_data[i, ]), NA,  # Both values are NA (no comparison)
      ifelse(
        !is.na(real_data[i, ]) & is.na(simulated_data[i, ]), FALSE,  # Real > Simulated
        ifelse(
          is.na(real_data[i, ]) & !is.na(simulated_data[i, ]), TRUE,  # Simulated > Real
          ifelse(
            !is.na(real_data[i, ]) & !is.na(simulated_data[i, ]),  # Both values present
            real_data[i, ] >= simulated_data[i, ],  # Compare values
            NA
          )
        )
      )
    )
  }
  
  # Calculate the proportion of TRUE values
  n_total <- sum(!is.na(comparison_matrix))  # Exclude NA values
  n_true <- sum(comparison_matrix, na.rm = TRUE)  # Count TRUE values
  
  # Compute confidence interval and hypothesis test
  ic_result <- prop.test(n_true, n_total)
  
  # Create a dataframe to store results
  results <- data.frame(
    sp = specie,
    interval = t_interval,
    proportion = ic_result$estimate,  # Proportion of TRUE values
    lower_CI = ic_result$conf.int[1],  # Lower confidence interval
    upper_CI = ic_result$conf.int[2],  # Upper confidence interval
    p_value = ic_result$p.value  # p-value from hypothesis test
  )
  
  return(results)
}


##########################
# Results carnivores vs rabbit
##########################
r_vul_ory_AB <- comparison(vul_ory_AB, s_vul_ory_AB, vulpes, "AB")
r_vul_ory_BA <- comparison(vul_ory_BA, s_vul_ory_BA, vulpes, "BA")
r_mel_ory_AB <- comparison(mel_ory_AB, s_mel_ory_AB, meles, "AB")
r_mel_ory_BA <- comparison(mel_ory_BA, s_mel_ory_BA, meles, "BA")
r_gen_ory_AB <- comparison(gen_ory_AB, s_gen_ory_AB, genetta, "AB")
r_gen_ory_BA <- comparison(gen_ory_BA, s_gen_ory_BA, genetta, "BA")
r_her_ory_AB <- comparison(her_ory_AB, s_her_ory_AB, herpestes, "AB")
r_her_ory_BA <- comparison(her_ory_BA, s_her_ory_BA, herpestes, "BA")
r_lyn_ory_AB <- comparison(lyn_ory_AB, s_lyn_ory_AB, lynx, "AB")
r_lyn_ory_BA <- comparison(lyn_ory_BA, s_lyn_ory_BA, lynx, "BA")

rabbit_results <- rbind(
  r_lyn_ory_BA, 
  r_vul_ory_BA, 
  r_mel_ory_BA, 
  r_gen_ory_BA, 
  r_her_ory_BA,
  r_lyn_ory_AB, 
  r_vul_ory_AB, 
  r_mel_ory_AB, 
  r_gen_ory_AB, 
  r_her_ory_AB
)

##########################
# Results carnivores vs hare
##########################

r_vul_lep_AB <- comparison(vul_lep_AB, s_vul_lep_AB, vulpes, "AB")
r_vul_lep_BA <- comparison(vul_lep_BA, s_vul_lep_BA, vulpes, "BA")
r_mel_lep_AB <- comparison(mel_lep_AB, s_mel_lep_AB, meles, "AB")
r_mel_lep_BA <- comparison(mel_lep_BA, s_mel_lep_BA, meles, "BA")
r_gen_lep_AB <- comparison(gen_lep_AB, s_gen_lep_AB, genetta, "AB")
r_gen_lep_BA <- comparison(gen_lep_BA, s_gen_lep_BA, genetta, "BA")
r_her_lep_AB <- comparison(her_lep_AB, s_her_lep_AB, herpestes, "AB")
r_her_lep_BA <- comparison(her_lep_BA, s_her_lep_BA, herpestes, "BA")
r_lyn_lep_AB <- comparison(lyn_lep_AB, s_lyn_lep_AB, lynx, "AB")
r_lyn_lep_BA <- comparison(lyn_lep_BA, s_lyn_lep_BA, lynx, "BA")

hare_results <- rbind(
  r_lyn_lep_BA, 
  r_vul_lep_BA, 
  r_mel_lep_BA, 
  r_gen_lep_BA, 
  r_her_lep_BA,
  r_lyn_lep_AB, 
  r_vul_lep_AB, 
  r_mel_lep_AB, 
  r_gen_lep_AB, 
  r_her_lep_AB
)


