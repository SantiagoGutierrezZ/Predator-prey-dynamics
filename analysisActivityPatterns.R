# 1. Loading packages ---- 
library(dplyr)
library(activity)
library(overlap)

########################
# 2. Loading data ---- #
########################
setwd("/Users/santiagogutierrezzapata/Documents/Universidad de Huelva - Doctorado/Papers/Activity patterns") # Change working directory
all_data <- read.csv("dataActivityPatterns.csv", header = TRUE)

# Time as a posixct object
all_data$clock_date_time <- as.POSIXct(paste(all_data$clock_date_time), "%Y-%m-%d %H:%M", tz = "GMT")

#################################
# 3. Average anchored times ---- #
#################################
# Split data to calculate average achored time per specie
data_per_sp <- split(all_data, all_data$sp)

# Latitude and longitude of a site
lat <- 36.96693
lon <- -6.46713

# Calculate averege achored times
vulpes <- solartime(data_per_sp[["Vulpes vulpes"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                    format = "%Y-%m-%d %H:%M")
meles <- solartime(data_per_sp[["Meles meles"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                  format = "%Y-%m-%d %H:%M")
genetta <- solartime(data_per_sp[["Genetta genetta"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                   format = "%Y-%m-%d %H:%M")
herpestes <- solartime(data_per_sp[["Herpestes ichneumon"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                       format = "%Y-%m-%d %H:%M")
lynx <- solartime(data_per_sp[["Lynx pardinus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                       format = "%Y-%m-%d %H:%M")
oryctolagus <- solartime(data_per_sp[["Oryctolagus cuniculus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                       format = "%Y-%m-%d %H:%M")
lepus <- solartime(data_per_sp[["Lepus granatensis"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                       format = "%Y-%m-%d %H:%M")

# Fit circular kernel to the average-anchored times
vulpes_s <- fitact(vulpes$solar, sample = "data")
meles_s <- fitact(meles$solar, sample = "data")
genetta_s <- fitact(genetta$solar, sample = "data")
herpestes_s <- fitact(herpestes$solar, sample = "data")
lynx_s <- fitact(lynx$solar, sample = "data")
oryctolagus_s <- fitact(oryctolagus$solar, sample = "data")
lepus_s <- fitact(lepus$solar, sample = "data")

#######################################################
# Get the detection probability of the species per min
#######################################################





# Graficar la distribución resultante
plot(t, pdf_l, type = "l", col = "blue", lwd = 2,
     xlab = "Tiempo (minutos)", ylab = "Probabilidad",
     main = "Distribución multimodal ajustada de probabilidad p")
abline(h=0.00004, lty=2)
grid()

#################################################
# 4. Set up boundaries for sunrise and sunset  ##
#    boxes and ticks marks in x-axis           ##
#################################################
# Set up boundaries for sunrise and sunset boxes during summer and autumn seasons

# Calculate the latest sunrise time in radians for summer and autumn
max_sunrise <- max(all_data$sunrise[all_data$season == "summer" | all_data$season == "autumn"])
# Calculate the earliest sunrise time in radians for summer and autumn
min_sunrise <- min(all_data$sunrise[all_data$season == "summer" | all_data$season == "autumn"])
# Calculate the latest sunset time in radians for summer and autumn
max_sunset <- max(all_data$sunset[all_data$season == "summer" | all_data$season == "autumn"])
# Calculate the earliest sunset time in radians for summer and autumn
min_sunset <- min(all_data$sunset[all_data$season == "summer" | all_data$season == "autumn"])

# Function to draw daytime polygons (representing sunrise and sunset ranges)
draw_day_polygons <- function() {
  # Draw polygon from 0 to max sunrise time, covering y-range from 0 to 10
  polygon(x = c(0, max_sunrise, max_sunrise, 0), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  # Draw polygon from 0 to min sunrise time
  polygon(x = c(0, min_sunrise, min_sunrise, 0), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  # Draw polygon from min sunset to 2*pi (full day), representing the end of the day
  polygon(x = c(min_sunset, (2 * pi), (2 * pi), min_sunset), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  # Draw polygon from max sunset to 2*pi, marking the latest sunset period
  polygon(x = c(max_sunset, (2 * pi), (2 * pi), max_sunset), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  box()
}

# Function to draw nighttime polygons
draw_night_polygons <- function() {
  # Draw polygon for nighttime starting at 0 up to max sunrise time
  polygon(x = c(0, max_sunrise, max_sunrise, 0), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  # Draw polygon for nighttime from 0 to min sunrise time
  polygon(x = c(0, min_sunrise, min_sunrise, 0), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  # Draw polygon starting at min sunset shifted by -2*pi, to 0, to represent early night period
  polygon(x = c(min_sunset - 2 * pi, 0, 0, min_sunset - 2 * pi), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  # Draw polygon starting at max sunset shifted by -2*pi to 0
  polygon(x = c(max_sunset - 2 * pi, 0, 0, max_sunset - 2 * pi), y = c(0, 0, 10, 10), 
          col = rgb(0.3, 0.3, 0.3, 0.3), border = NA)
  box()
}

# Define tick marks in radians for a daytime x-axis
hours_in_radians_day <- c(0, 3, 6, 9, 12, 15, 18, 21, 24) * (pi / 12)
draw_day_xaxis <- function() {
  axis(1, at = hours_in_radians_day, 
       labels = c("0", "3:00", "6:00", "9:00", "12:00", "15:00", "18:00", "21:00", "24:00"), 
       tick = TRUE, las = 1)
}

# Define tick marks in radians for a nighttime x-axis centered on midnight
hours_in_radians_night <- c(-12, -9, -6, -3, 0, 3, 6, 9, 12) * (pi / 12)
draw_night_xaxis <- function() {
  axis(1, at = hours_in_radians_night,
       labels = c("12:00", "15:00", "18:00", "21:00", "24:00", "3:00", "6:00", "9:00", "12:00"),
       tick = TRUE, las = 1)
}

#############################################################################
# 5. Define functions for to draw the diel activity pattern for each specie #
#############################################################################
col_meles <- "#00BFB2"
draw_plot_meles <- function(centre = "",
                            add = FALSE) {
  plot(meles_s, 
       data = "none",  
       xunit="radians", yunit="density",  
       ylim=c(0,0.6), 
       centre = centre, 
       ylab = "", xlab = "", 
       tline = list(col = col_meles, lwd = 3),
       dline = list(col = "transparent"), 
       cline = list(col = "transparent"),
       xaxis=list(xaxt="n"),
       add = add)
}
col_vulpes <- "#FF8C00"
draw_plot_vulpes <- function(centre = "",
                             add = FALSE) {
  plot(vulpes_s,
       data = "none",
       xunit="radians", yunit="density", 
       ylim=c(0,0.6), 
       centre = centre, 
       ylab="", xlab="", 
       tline = list(col = col_vulpes, lwd = 3), 
       dline = list(col = "transparent"),
       cline = list(col = "transparent"),
       xaxis=list(xaxt="n"),
       add = add)
}
col_genetta <- "#163070"
draw_plot_genetta <- function(centre = "",
                              add = FALSE) {
  plot(genetta_s, 
       data = "none",
       xunit="radians", yunit="density", 
       ylim=c(0,0.6), 
       centre = centre, 
       ylab="", xlab="", 
       tline = list(col = col_genetta, lwd = 3),
       dline = list(col = "transparent"),
       cline = list(col = "transparent"),
       xaxis=list(xaxt="n"),
       add = add)
}
col_herpestes <- "#DB0D44"
draw_plot_herpestes <- function(centre = "",
                                add = FALSE) {
  plot(herpestes_s, 
       data = "none",
       xunit="radians", yunit="density", 
       ylim=c(0,0.6), 
       centre = centre, 
       ylab="", xlab="", 
       tline = list(col = col_herpestes, lwd = 3),
       dline = list(col = "transparent"),
       cline = list(col = "transparent"),
       xaxis=list(xaxt="n"),
       add = add)
}
col_lynx <- "#7423D9"
draw_plot_lynx <- function(centre = "",
                           add = FALSE) {
  plot(lynx_s, 
       data = "none",
       xunit="radians", yunit="density", 
       ylim=c(0,0.6), 
       centre = centre, 
       ylab="", xlab="", 
       tline = list(col = col_lynx, lwd = 3),
       dline = list(col = "transparent"),
       cline = list(col = "transparent"),
       xaxis=list(xaxt="n"),
       add = add)
}

col_lepus <- "#4BEB74"
draw_plot_lepus <- function(centre = "",
                            add = FALSE) {
  plot(lepus_s, 
       data = "none",
       xunit="radians", yunit="density", 
       ylim=c(0,0.8), 
       centre = centre, 
       ylab="", xlab="", 
       tline = list(col = col_lepus, lwd = 3),
       dline = list(col = "transparent"),
       cline = list(col = "transparent"),
       xaxis=list(xaxt="n"),
       add = add)
}
col_oryctolagus <- "#FF33A8"
draw_plot_oryctolagus <- function(centre = "",
                                  add = FALSE) {
  plot(oryctolagus_s, 
       data = "none",
       xunit="radians", yunit="density", 
       ylim=c(0,0.8), 
       centre = centre, 
       ylab="", xlab="", 
       tline = list(col = col_oryctolagus, lwd = 3),
       dline = list(col = "transparent"),
       cline = list(col = "transparent"),
       xaxis=list(xaxt="n"),
       add = add)
}

####################################################
###                                              ###
### 5. Tempora overlap between predator/predator ###
###                                              ###
####################################################
# Size of legend
size_legend <- 2

# Vulpes vulpes
png("newPlots/species/allSp2.png", width = 3000, height = 1200)
par(family="serif", cex.lab = 4, cex.axis = 4)
layout(matrix(1:8, ncol = 4, nrow = 2))
# png("newPlots/species/fox.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_vulpes(centre = "night", add = FALSE)
draw_night_xaxis()
lines(rep(cmean(data_per_sp[["Vulpes vulpes"]]$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(data_per_sp[["Vulpes vulpes"]]$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Red fox"), col = c(col_vulpes), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

# lines(rep(mean(data_per_sp[["Vulpes vulpes"]]$sunset [data_per_sp[["Vulpes vulpes"]]$season=="summer" | data_per_sp[["Vulpes vulpes"]]$season=="autumn" ]),2) - 2*pi, c(0,2500), col="red")
# lines(rep(mean(data_per_sp[["Vulpes vulpes"]]$sunrise [data_per_sp[["Vulpes vulpes"]]$season=="summer" | data_per_sp[["Vulpes vulpes"]]$season=="autumn" ]),2), c(0,2500), col="red")
# dev.off()

# Meles meles
# png("newPlots/species/badger.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_meles(centre = "night", add = FALSE)
draw_night_xaxis()
lines(rep(cmean(data_per_sp[["Meles meles"]]$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(data_per_sp[["Meles meles"]]$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Badger"), col = c(col_meles), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
# dev.off()

# Genetta genetta
# png("newPlots/species/genet.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_genetta(centre = "night", add = FALSE)
draw_night_xaxis()
lines(rep(cmean(data_per_sp[["Genetta genetta"]]$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(data_per_sp[["Genetta genetta"]]$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Genet"), col = c(col_genetta), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
# dev.off()

# Herpestes ichneumon
# png("newPlots/species/mongoose.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_herpestes(centre = "day", add = FALSE)
draw_day_xaxis()
lines(rep(cmean(data_per_sp[["Herpestes ichneumon"]]$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(data_per_sp[["Herpestes ichneumon"]]$sunset),2), c(0,2500), col="red", lwd = 2)
legend("top", c("Mongoose"), col = c(col_herpestes), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
# dev.off()

# Lynx pardinus
# png("newPlots/species/lynx.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_lynx(centre = "night", add = FALSE)
draw_night_xaxis()
lines(rep(cmean(data_per_sp[["Lynx pardinus"]]$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(data_per_sp[["Lynx pardinus"]]$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Lynx"), col = c(col_lynx), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
# dev.off()

# Oryctolagus cuniculus
# png("newPlots/species/rabbit.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_oryctolagus(centre = "night", add = FALSE)
draw_night_xaxis()
lines(rep(cmean(data_per_sp[["Oryctolagus cuniculus"]]$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(data_per_sp[["Oryctolagus cuniculus"]]$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Rabbit"), col = c(col_oryctolagus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
# dev.off()

# Lepus granatensis
# png("newPlots/species/hare.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_lepus(centre = "night", add = FALSE)
draw_night_xaxis()
lines(rep(cmean(data_per_sp[["Lepus granatensis"]]$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(data_per_sp[["Lepus granatensis"]]$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Hare"), col = c(col_lepus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
layout(1)
dev.off()

#####################################################
###                                               ###
### 6. Temporal overlap between predator/predator ###
###                                               ###
#####################################################

size_legend <- 4
size_delta <- 4

png("newPlots/pairs/allPairs2.png", width = 2000, height = 2800)
par(family="serif", cex.lab = 4, cex.axis = 4)
layout(matrix(1:10, ncol = 2, nrow = 5))

# Vulpes vulpes vs Meles meles
# png("newPlots/pairs/fox-badger.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_vulpes(centre = "night", add = FALSE)
draw_plot_meles(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Red fox", "Badger"), col = c(col_vulpes, col_meles), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_vul_mel <- overlapEst(vulpes$solar, meles$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_vul_mel <- bootstrap(vulpes$solar, meles$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_vul_mel)
conf <- bootCI(ove_vul_mel, boost_vul_mel)
delta <- paste0("∆4: ", round(ove_vul_mel, 3)," ", "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Vulpes vulpes vs Genetta genetta
# png("newPlots/pairs/fox-genet.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_vulpes(centre = "night", add = FALSE)
draw_plot_genetta(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Red fox", "Genet"), col = c(col_vulpes, col_genetta), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_vul_gen <- overlapEst(vulpes$solar, genetta$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_vul_gen <- bootstrap(vulpes$solar, genetta$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_vul_gen)
conf <- bootCI(ove_vul_gen, boost_vul_gen)
delta <- paste0("∆4: ", round(ove_vul_gen, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Vulpes vulpes vs Herpestes ichneumon
# png("newPlots/pairs/fox-mongoose.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_vulpes(centre = "day", add = FALSE)
draw_plot_herpestes(centre = "day", add = TRUE)
draw_day_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2), c(0,2500), col="red", lwd = 2)
legend("topleft", c("Red fox", "Mongoose"), col = c(col_vulpes, col_herpestes), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_vul_her <- overlapEst(vulpes$solar, herpestes$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_vul_her <- bootstrap(vulpes$solar, herpestes$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_vul_her)
conf <- bootCI(ove_vul_her, boost_vul_her)
delta <- paste0("∆4: ", round(ove_vul_her, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 3.14, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Vulpes vulpes vs Lynx pardinus
# png("newPlots/pairs/fox-lynx.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_vulpes(centre = "night", add = FALSE)
draw_plot_lynx(centre = "day", add = FALSE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2 * pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Red fox", "Lynx"), col = c(col_vulpes, col_lynx), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_vul_lyn <- overlapEst(vulpes$solar, lynx$solar, type = "Dhat1") # Overlapping coefficient Dht4
# boost_vul_lyn <- bootstrap(vulpes$solar, lynx$solar, 999, smooth = TRUE, type = "Dhat1")
mean_over <- mean(boost_vul_lyn)
conf <- bootCI(ove_vul_lyn, boost_vul_lyn)
delta <- paste0("∆1: ", round(ove_vul_lyn, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Meles meles vs Genetta genetta
# png("newPlots/pairs/badger-genet.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_meles(centre = "night", add = FALSE)
draw_plot_genetta(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2 * pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Badger", "Genet"), col = c(col_meles, col_genetta), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_mel_gen <- overlapEst(meles$solar, genetta$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_mel_gen <- bootstrap(meles$solar, genetta$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_mel_gen)
conf <- bootCI(ove_mel_gen, boost_mel_gen)
delta <- paste0("∆4: ", round(ove_mel_gen, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Meles meles vs Lynx pardinus
# png("newPlots/pairs/badger-lynx.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_meles(centre = "night", add = FALSE)
draw_plot_lynx(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2 * pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Badger", "Lynx"), col = c(col_meles, col_lynx), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_mel_lyn <- overlapEst(meles$solar, lynx$solar, type = "Dhat1") # Overlapping coefficient Dht4
# boost_mel_lyn <- bootstrap(meles$solar, lynx$solar, 999, smooth = TRUE, type = "Dhat1")
mean_over <- mean(boost_mel_lyn)
conf <- bootCI(ove_mel_lyn, boost_mel_lyn)
delta <- paste0("∆1: ", round(ove_mel_lyn, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Meles meles vs Herpestes ichneumon
# png("newPlots/pairs/badger-mongoose.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_meles(centre = "day", add = FALSE)
draw_plot_herpestes(centre = "day", add = TRUE)
draw_day_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2), c(0,2500), col="red", lwd = 2)
legend("topleft", c("Badger", "Mongoose"), col = c(col_meles, col_herpestes), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_mel_her <- overlapEst(meles$solar, herpestes$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_mel_her <- bootstrap(meles$solar, herpestes$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_mel_her)
conf <- bootCI(ove_mel_her, boost_mel_her)
delta <- paste0("∆4: ", round(ove_mel_her, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 3.14, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Genetta genetta vs Lynx pardinus
# png("newPlots/pairs/genet-lynx.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_genetta(centre = "night", add = FALSE)
draw_plot_lynx(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2 * pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Genet", "Lynx"), col = c(col_genetta, col_lynx), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_gen_lyn <- overlapEst(genetta$solar, lynx$solar, type = "Dhat1") # Overlapping coefficient Dht4
# boost_gen_lyn <- bootstrap(genetta$solar, lynx$solar, 999, smooth = TRUE, type = "Dhat1")
mean_over <- mean(boost_gen_lyn)
conf <- bootCI(ove_gen_lyn, boost_gen_lyn)
delta <- paste0("∆1: ", round(ove_gen_lyn, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Genetta genetta vs Herpestes ichneumon
# png("newPlots/pairs/genet-mongoose.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_genetta(centre = "day", add = FALSE)
draw_plot_herpestes(centre = "day", add = TRUE)
draw_day_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2), c(0,2500), col="red", lwd = 2)
legend("topleft", c("Genet", "Mongoose"), col = c(col_genetta, col_herpestes), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_gen_her <- overlapEst(genetta$solar, herpestes$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_gen_her <- bootstrap(genetta$solar, herpestes$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_gen_her)
conf <- bootCI(ove_gen_her, boost_gen_her)
delta <- paste0("∆4: ", round(ove_gen_her, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 3.14, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Lynx pardinus vs Herpestes ichneumon
# png("newPlots/pairs/lynx-mongoose.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_lynx(centre = "day", add = FALSE)
draw_plot_herpestes(centre = "day", add = TRUE)
draw_day_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2), c(0,2500), col="red", lwd = 2)
legend("topleft", c("Lynx", "Mongoose"), col = c(col_lynx, col_herpestes), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_lyn_her <- overlapEst(lynx$solar, herpestes$solar, type = "Dhat1") # Overlapping coefficient Dht4
# boost_lyn_her <- bootstrap(lynx$solar, herpestes$solar, 999, smooth = TRUE, type = "Dhat1")
mean_over <- mean(boost_lyn_her)
conf <- bootCI(ove_lyn_her, boost_lyn_her)
delta <- paste0("∆1: ", round(ove_lyn_her, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 3.14, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
layout(1)
dev.off()

#####################################################
###                                               ###
### 7. Temporal overlap between predator/prey     ###
###                                               ###
#####################################################

png("newPlots/predator-prey/predator-prey.png", width = 2000, height = 2800)
par(family="serif", cex.lab = 4, cex.axis = 4)
layout(matrix(1:10, ncol = 2, nrow = 5))

# Vulpes vulpes vs Oryctolagus cuniculus
# png("newPlots/predator-prey/fox-rabbit.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_vulpes(centre = "night", add = FALSE)
draw_plot_oryctolagus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Red fox", "Rabbit"), col = c(col_vulpes, col_oryctolagus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_ory_vul <- overlapEst(vulpes$solar, oryctolagus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_ory_vul <- bootstrap(vulpes$solar, oryctolagus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_ory_vul)
conf <- bootCI(ove_ory_vul, boost_ory_vul)
delta <- paste0("∆4: ", round(ove_ory_vul, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Meles meles vs Oryctolagus cuniculus
# png("newPlots/predator-prey/badger-rabbit.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_meles(centre = "night", add = FALSE)
draw_plot_oryctolagus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Badger", "Rabbit"), col = c(col_meles, col_oryctolagus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_ory_mel <- overlapEst(meles$solar, oryctolagus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_ory_mel <- bootstrap(meles$solar, oryctolagus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_ory_mel)
conf <- bootCI(ove_ory_mel, boost_ory_mel)
delta <- paste0("∆4: ", round(ove_ory_mel, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Genetta genetta vs Oryctolagus cuniculus
# png("newPlots/predator-prey/genet-rabbit.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_genetta(centre = "night", add = FALSE)
draw_plot_oryctolagus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Genet", "Rabbit"), col = c(col_genetta, col_oryctolagus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)


ove_ory_gen <- overlapEst(genetta$solar, oryctolagus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_ory_gen <- bootstrap(genetta$solar, oryctolagus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_ory_gen)
conf <- bootCI(ove_ory_gen, boost_ory_gen)
delta <- paste0("∆4: ", round(ove_ory_gen, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Herpestes ichneumon vs Oryctolagus cuniculus
# png("newPlots/predator-prey/mongoose-rabbit.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_herpestes(centre = "day", add = FALSE)
draw_plot_oryctolagus(centre = "day", add = TRUE)
draw_day_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2), c(0,2500), col="red", lwd = 2)
legend("topleft", c("Mongoose", "Rabbit"), col = c(col_herpestes, col_oryctolagus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_ory_her <- overlapEst(herpestes$solar, oryctolagus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_ory_her <- bootstrap(herpestes$solar, oryctolagus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_ory_her)
conf <- bootCI(ove_ory_her, boost_ory_her)
delta <- paste0("∆4: ", round(ove_ory_her, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 3.14, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Lynx pardinus vs Oryctolagus cuniculus
# png("newPlots/predator-prey/lynx-rabbit.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_lynx(centre = "night", add = FALSE)
draw_plot_oryctolagus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Lynx", "Rabbit"), col = c(col_lynx, col_oryctolagus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_ory_lyn <- overlapEst(lynx$solar, oryctolagus$solar, type = "Dhat1") # Overlapping coefficient Dht4
# boost_ory_lyn <- bootstrap(lynx$solar, oryctolagus$solar, 999, smooth = TRUE, type = "Dhat1")
mean_over <- mean(boost_ory_lyn)
conf <- bootCI(ove_ory_lyn, boost_ory_lyn)
delta <- paste0("∆1: ", round(ove_ory_lyn, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Vulpes vulpes vs Lepus granatensis
# png("newPlots/predator-prey/fox-hare.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_vulpes(centre = "night", add = FALSE)
draw_plot_lepus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Red fox", "Hare"), col = c(col_vulpes, col_lepus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_lep_vul <- overlapEst(vulpes$solar, lepus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_lep_vul <- bootstrap(vulpes$solar, lepus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_lep_vul)
conf <- bootCI(ove_lep_vul, boost_lep_vul)
delta <- paste0("∆4: ", round(ove_lep_vul, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Meles meles vs Lepus granatensis
# png("newPlots/predator-prey/badger-hare.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_meles(centre = "night", add = FALSE)
draw_plot_lepus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Badger", "Hare"), col = c(col_meles, col_lepus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_lep_mel <- overlapEst(meles$solar, lepus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_lep_mel <- bootstrap(meles$solar, lepus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_lep_mel)
conf <- bootCI(ove_lep_mel, boost_lep_mel)
delta <- paste0("∆4: ", round(ove_lep_mel, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Genetta genetta vs Lepus granatensis
# png("newPlots/predator-prey/genet-hare.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_genetta(centre = "night", add = FALSE)
draw_plot_lepus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Genet", "Hare"), col = c(col_genetta, col_lepus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_lep_gen <- overlapEst(genetta$solar, lepus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_lep_gen <- bootstrap(genetta$solar, lepus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_lep_gen)
conf <- bootCI(ove_lep_gen, boost_lep_gen)
delta <- paste0("∆4: ", round(ove_lep_gen, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Herpestes ichneumon vs Lepus granatensis
# png("newPlots/predator-prey/mongoose-hare.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_herpestes(centre = "day", add = FALSE)
draw_plot_lepus(centre = "day", add = TRUE)
draw_day_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2), c(0,2500), col="red", lwd = 2)
legend("topleft", c("Mongoose", "Hare"), col = c(col_herpestes, col_lepus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_lep_her <- overlapEst(herpestes$solar, lepus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_lep_her <- bootstrap(herpestes$solar, lepus$solar, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_lep_her)
conf <- bootCI(ove_lep_her, boost_lep_her)
delta <- paste0("∆4: ", round(ove_lep_her, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 3.14, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
# dev.off()

# Lynx pardinus vs Lepus granatensis
# png("newPlots/predator-prey/lynx-hare.png", width = 800, height = 600)
# par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_lynx(centre = "night", add = FALSE)
draw_plot_lepus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Lynx", "Hare"), col = c(col_lynx, col_lepus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)

ove_lep_lyn <- overlapEst(lynx$solar, lepus$solar, type = "Dhat1") # Overlapping coefficient Dht4
# boost_lep_lyn <- bootstrap(lynx$solar, lepus$solar, 999, smooth = TRUE, type = "Dhat1")
mean_over <- mean(boost_lep_lyn)
conf <- bootCI(ove_lep_lyn, boost_lep_lyn)
delta <- paste0("∆1: ", round(ove_lep_lyn, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)
layout(1)
dev.off()

#####################################################
###                                               ###
### 8. Temporal overlap between prey/prey         ###
###                                               ###
#####################################################
size_legend <- 2
size_delta <- 2

# Oryctolagus cuniculus vs Lepus granatensis
png("newPlots/rabbit-hare.png", width = 800, height = 600)
par(family="serif", cex.lab = 2.2, cex.axis = 2.2)
draw_plot_oryctolagus(centre = "night", add = FALSE)
draw_plot_lepus(centre = "night", add = TRUE)
draw_night_xaxis()
lines(rep(cmean(all_data$sunrise),2), c(0,2500), col="red", lwd = 2)
lines(rep(cmean(all_data$sunset),2) - 2*pi, c(0,2500), col="red", lwd = 2)
legend("topleft", c("Rabbit", "Hare"), col = c(col_oryctolagus, col_lepus), 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
# OVERLAP COEFFICIENTS AND BOOSTRAP CONFIDENC INTERVALS
ove_lep_ory <- overlapEst(oryctolagus$solar, lepus$solar, type = "Dhat4") # Overlapping coefficient Dht4
# boost_lep_ory <- bootstrap(oryctolagus$solar, lepus$solar, 999, smooth = TRUE, type = "Dhat4") 
mean_over <- mean(boost_lep_ory)
conf <- bootCI(ove_lep_ory, boost_lep_ory)
delta <- paste0("∆4: ", round(ove_lep_ory, 3), "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.65, labels = delta, cex = size_delta, font = 2, pos = 3)
dev.off()

##########################################################################
##                                                                      ##
##   SPATIOTEMPORAL ANALYSIS        ##
##                                                                      ##
##########################################################################

# We are going to obtain the data for the spatiotemporal analysis. For doing that, we need to select two
# species, A & B and find type of interval between consecutive detections in a specific location 
# (AB, BA, AA, BB) and the difftime between these intervals.

# List of predators
predators <- c("Meles meles", "Genetta genetta", "Herpestes ichneumon", "Lynx pardinus", "Vulpes vulpes")
rabbit <- "Oryctolagus cuniculus"
hare <- "Lepus granatensis"

# Define the function to process spatiotemporal data for a predator-rabbit pair
rabbit_spatiotemporal_data <- function(predator, rabbit, data) {
  data %>%
    # Filter the dataset to include only records of Predator and rabbit
    filter(sp %in% c(predator, rabbit)) %>%
    # Arrange the dataset by site and time of detection in ascending order
    arrange(site, clock_date_time) %>%
    # Group the data by site to process interactions within each location independently
    group_by(site) %>%
    # Create new columns:
    mutate(
      # Define interaction intervals based on consecutive species detections
      interval = case_when(
        sp == predator & lead(sp) == rabbit ~ "AB",
        sp == rabbit & lead(sp) == predator ~ "BA",
        sp == predator & lead(sp) == predator ~ "AA",
        sp == rabbit & lead(sp) == rabbit ~ "BB",
        TRUE ~ NA_character_
      ),
      # Calculate the time difference (in minutes) between consecutive detections
      time_diff = abs(difftime(clock_date_time, lead(clock_date_time), units = "mins"))
    ) %>%
    # Filter the dataset to keep only interactions where a fox and a rabbit are involved
    filter(interval %in% c("AB", "BA")) %>%
    # Select only the relevant columns for the final dataset
    select(site, clock_date_time, sp, interval, time_diff)
}

# Apply the function to each predator and store the results in a list
rabbit_predators_data <- lapply(predators, function(predator) {
  rabbit_spatiotemporal_data(predator, rabbit, all_data)
})

# Assign names to the list for clarity
names(rabbit_predators_data) <- predators

# Define the function to process spatiotemporal data for a predator-hare pair
hare_spatiotemporal_data <- function(predator, hare, data) {
  data %>%
    # Filter the dataset to include only records of Predator and hare
    filter(sp %in% c(predator, hare)) %>%
    # Arrange the dataset by site and time of detection in ascending order
    arrange(site, clock_date_time) %>%
    # Group the data by site to process interactions within each location independently
    group_by(site) %>%
    # Create new columns:
    mutate(
      # Define interaction intervals based on consecutive species detections
      interval = case_when(
        sp == predator & lead(sp) == hare ~ "AB",
        sp == hare & lead(sp) == predator ~ "BA",
        sp == predator & lead(sp) == predator ~ "AA",
        sp == hare & lead(sp) == hare ~ "BB",
        TRUE ~ NA_character_
      ),
      # Calculate the time difference (in minutes) between consecutive detections
      time_diff = abs(difftime(clock_date_time, lead(clock_date_time), units = "mins"))
    ) %>%
    # Filter the dataset to keep only interactions where a fox and a hare are involved
    filter(interval %in% c("AB", "BA")) %>%
    # Select only the relevant columns for the final dataset
    select(site, clock_date_time, sp, interval, time_diff)
}

# Apply the function to each predator and store the results in a list
hare_predators_data <- lapply(predators, function(predator) {
  hare_spatiotemporal_data(predator, hare, all_data)
})

# Assign names to the list for clarity
names(hare_predators_data) <- predators


lynx_mongoose2 <- all_data %>%
  # Filter the dataset to include only records of Predator and hare
  filter(sp %in% c("Lynx pardinus", "Herpestes ichneumon")) %>%
  # Arrange the dataset by site and time of detection in ascending order
  arrange(site, clock_date_time) %>%
  # Group the data by site to process interactions within each location independently
  group_by(site) %>%
  # Create new columns:
  mutate(
    # Define interaction intervals based on consecutive species detections
    interval = case_when(
      sp == "Lynx pardinus" & lead(sp) == "Herpestes ichneumon" ~ "AB",
      sp == "Herpestes ichneumon" & lead(sp) == "Lynx pardinus" ~ "BA",
      sp == "Lynx pardinus" & lead(sp) == "Lynx pardinus" ~ "AA",
      sp == "Herpestes ichneumon" & lead(sp) == "Herpestes ichneumon" ~ "BB",
      TRUE ~ NA_character_
    ),
    # Calculate the time difference (in minutes) between consecutive detections
    time_diff = abs(difftime(clock_date_time, lead(clock_date_time), units = "mins"))
  ) %>%
  # Filter the dataset to keep only interactions where a fox and a hare are involved
  # filter(interval %in% c("AB", "BA")) %>%
  # Select only the relevant columns for the final dataset
  select(site, clock_date_time, sp, interval, time_diff)

lynx_mongoose2 <- lynx_mongoose2 %>%
  mutate(time_diff = abs(difftime(clock_date_time, lead(clock_date_time), units = "mins")))

predators_spatiotemporal_data <- list(lynx_badger = lynx_badger,
                                      lynx_fox = lynx_fox)


predators_spatiotemporal_data[["lynx_badger"]]

save(all_data,
     hare_predators_data,
     rabbit_predators_data,
     lynx_badger, lynx_fox, lynx_genet, lynx_mongoose,
     badger_fox, badger_genet, badger_lynx, badger_mongoose,
     genet_badger, genet_fox, genet_lynx, genet_mongoose,
     fox_badger, fox_genet, fox_lynx, fox_mongoose,
     mongoose_badger, mongoose_fox, mongoose_lynx, mongoose_genet, file = "spatiotemporalData.RData")