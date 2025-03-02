# 1. Loading packages ---- 
library(dplyr)
library(activity)
library(overlap)

########################
# 2. Loading data ---- #
########################
setwd("/Users/santiagogutierrezzapata/Documents/Universidad de Huelva - Doctorado/Papers/Activity patterns") # Change working directory
all_data <- read.csv("dataActivityPatterns.csv", header = TRUE)
load("seasonAnalysisActivity.RData")
# Time as a posixct object
all_data$clock_date_time <- as.POSIXct(paste(all_data$clock_date_time), "%Y-%m-%d %H:%M", tz = "GMT")

#################################
# 3. get data per each season ---- #
#################################

# Spring
spring_data <- all_data[all_data$season == "spring",]
spring_per_sp <- split(spring_data, spring_data$sp)

# Summer
summer_data <- all_data[all_data$season == "summer",]
summer_per_sp <- split(summer_data, summer_data$sp)

# autumn
autumn_data <- all_data[all_data$season == "autumn",]
autumn_per_sp <- split(autumn_data, autumn_data$sp)

# Winter
winter_data <- all_data[all_data$season == "winter",]
winter_per_sp <- split(winter_data, winter_data$sp)


#####################
##  Vulpes vulpes  ##
#####################

# Latitude and longitude of a site
lat <- 36.96693
lon <- -6.46713

# Calculate averege achored times
vulpes_spring <- solartime(spring_per_sp[["Vulpes vulpes"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                    format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
vulpes_summer <- solartime(summer_per_sp[["Vulpes vulpes"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                           format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
vulpes_autumn <- solartime(autumn_per_sp[["Vulpes vulpes"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                           format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
vulpes_winter <- solartime(winter_per_sp[["Vulpes vulpes"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                           format = "%Y-%m-%d %H:%M")

spring_vul <- fitact(vulpes_spring$clock, sample = "data")
summer_vul <- fitact(vulpes_summer$clock, sample = "data")
autumn_vul <- fitact(vulpes_autumn$clock, sample = "data")
winter_vul <- fitact(vulpes_winter$clock, sample = "data")


###############################
##  Meles meles  ##
#####################

# Calculate averege achored times
meles_spring <- solartime(spring_per_sp[["Meles meles"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
meles_summer <- solartime(summer_per_sp[["Meles meles"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
meles_autumn <- solartime(autumn_per_sp[["Meles meles"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
meles_winter <- solartime(winter_per_sp[["Meles meles"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")

spring_mel <- fitact(meles_spring$clock, sample = "data")
summer_mel <- fitact(meles_summer$clock, sample = "data")
autumn_mel <- fitact(meles_autumn$clock, sample = "data")
winter_mel <- fitact(meles_winter$clock, sample = "data")

###############################
##  Genetta genetta  ##
#####################

# Calculate averege achored times
genetta_spring <- solartime(spring_per_sp[["Genetta genetta"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
genetta_summer <- solartime(summer_per_sp[["Genetta genetta"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
genetta_autumn <- solartime(autumn_per_sp[["Genetta genetta"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
genetta_winter <- solartime(winter_per_sp[["Genetta genetta"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")

spring_gen <- fitact(genetta_spring$clock, sample = "data")
summer_gen <- fitact(genetta_summer$clock, sample = "data")
autumn_gen <- fitact(genetta_autumn$clock, sample = "data")
winter_gen <- fitact(genetta_winter$clock, sample = "data")


###############################
##  Herpestes ichneumon  ##
#####################

# Calculate averege achored times
herpestes_spring <- solartime(spring_per_sp[["Herpestes ichneumon"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
herpestes_summer <- solartime(summer_per_sp[["Herpestes ichneumon"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
herpestes_autumn <- solartime(autumn_per_sp[["Herpestes ichneumon"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
herpestes_winter <- solartime(winter_per_sp[["Herpestes ichneumon"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                            format = "%Y-%m-%d %H:%M")

spring_her <- fitact(herpestes_spring$clock, sample = "data")
summer_her <- fitact(herpestes_summer$clock, sample = "data")
autumn_her <- fitact(herpestes_autumn$clock, sample = "data")
winter_her <- fitact(herpestes_winter$clock, sample = "data")


###############################
##  Lynx pardinus  ##
#####################

# Calculate averege achored times
lynx_spring <- solartime(spring_per_sp[["Lynx pardinus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                              format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
lynx_summer <- solartime(summer_per_sp[["Lynx pardinus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                              format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
lynx_autumn <- solartime(autumn_per_sp[["Lynx pardinus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                              format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
lynx_winter <- solartime(winter_per_sp[["Lynx pardinus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                              format = "%Y-%m-%d %H:%M")

spring_lyn <- fitact(lynx_spring$clock, sample = "data")
summer_lyn <- fitact(lynx_summer$clock, sample = "data")
autumn_lyn <- fitact(lynx_autumn$clock, sample = "data")
winter_lyn <- fitact(lynx_winter$clock, sample = "data")

###############################
##  Lepus granatensis ##
#####################

# Calculate averege achored times
lepus_spring <- solartime(spring_per_sp[["Lepus granatensis"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                         format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
lepus_summer <- solartime(summer_per_sp[["Lepus granatensis"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                         format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
lepus_autumn <- solartime(autumn_per_sp[["Lepus granatensis"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                         format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
lepus_winter <- solartime(winter_per_sp[["Lepus granatensis"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                         format = "%Y-%m-%d %H:%M")

spring_lep <- fitact(lepus_spring$clock, sample = "data")
summer_lep <- fitact(lepus_summer$clock, sample = "data")
autumn_lep <- fitact(lepus_autumn$clock, sample = "data")
winter_lep <- fitact(lepus_winter$clock, sample = "data")

###############################
##  Oryctolagus cuniculus ##
#####################

# Calculate averege achored times
oryctolagus_spring <- solartime(spring_per_sp[["Oryctolagus cuniculus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
oryctolagus_summer <- solartime(summer_per_sp[["Oryctolagus cuniculus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
oryctolagus_autumn <- solartime(autumn_per_sp[["Oryctolagus cuniculus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")
# Calculate averege achored times
oryctolagus_winter <- solartime(winter_per_sp[["Oryctolagus cuniculus"]]$clock_date_time, lat = lat, lon = lon, tz = 1,
                          format = "%Y-%m-%d %H:%M")

spring_ory <- fitact(oryctolagus_spring$clock, sample = "data")
summer_ory <- fitact(oryctolagus_summer$clock, sample = "data")
autumn_ory <- fitact(oryctolagus_autumn$clock, sample = "data")
winter_ory <- fitact(oryctolagus_winter$clock, sample = "data")



###############################
##  Oryctolagus cuniculus ##
#####################

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

col_spring <- "green"
col_summer <- "orange"
col_autumn <- "red"
col_winter <- "blue"

###############################
##  Vulpes vulpes ##
#####################

png("newPlots/seasons/vulpes.png", width = 800, height = 600)
par(family="serif", cex.lab = 1, cex.axis = 1)
plot(spring_vul, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_spring, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_vul, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_vul, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_vul, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", c("Red fox", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()


###############################
##  Meles meles ##
#####################

png("newPlots/seasons/meles.png", width = 800, height = 600)
par(family="serif", cex.lab = 1, cex.axis = 1)
plot(spring_mel, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_spring, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_mel, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_mel, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_mel, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", c("Badger", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()

###############################
##  Genetta genetta ##
#####################

png("newPlots/seasons/genet.png", width = 800, height = 600)
par(family="serif", cex.lab = 1, cex.axis = 1)
plot(spring_gen, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_spring, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_gen, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_gen, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_gen, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", c("Genet", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()


###############################
##  Herpestes ichneumon ##
#####################

png("newPlots/seasons/mongoose.png", width = 800, height = 600)
par(family="serif", cex.lab = 1, cex.axis = 1)
plot(spring_her, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "day", ylab="", xlab="", tline = list(col = col_spring, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_her, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "day", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_her, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "day", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_her, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "day", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_day_xaxis()
legend("topleft", c("Mongoose", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()


###############################
##  Lynx pardinus ##
#####################

png("newPlots/seasons/lynx.png", width = 800, height = 600)
par(family="serif", cex.lab = 1, cex.axis = 1)
plot(spring_lyn, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_spring, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_lyn, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_lyn, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_lyn, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", c("Lynx", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()

###############################
##  Oryctolagus cuniculus ##
#####################

png("newPlots/seasons/rabbit.png", width = 800, height = 600)
par(family="serif", cex.lab = 1, cex.axis = 1)
plot(spring_ory, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_spring, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_ory, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_ory, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_ory, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", c("Rabbit", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()

###############################
##  Lepus granatensis ##
#####################

png("newPlots/seasons/hare.png", width = 800, height = 600)
par(family="serif", cex.lab = 1, cex.axis = 1)
plot(spring_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_spring, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", c("Hare", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()


###############################
##  Vulpes vulpes /  ##
#####################
col_vulpes <- "#FF8C00"
col_meles <- "#00BFB2"
col_genetta <- "#163070"
col_herpestes <- "#DB0D44"
col_lynx <- "#7423D9"
col_lepus <- "#4BEB74"
col_oryctolagus <- "#FF33A8"

size_legend <- 1
size_delta <- 1

png("newPlots/seasons/fox-preys.png", width = 2000, height = 2800)
par(family="serif", cex.lab = 4, cex.axis = 4)
layout(matrix(1:8, ncol = 2, nrow = 4))

plot(spring_vul, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_vulpes, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(spring_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_lepus, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", "spring", col = "white", 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
vul_lep_spr <- overlapEst(vulpes_spring$clock, lepus_spring$clock, type = "Dhat4") # Overlapping coefficient Dht4
# boost_vul_lep <- bootstrap(vulpes_spring$clock, lepus_spring$clock, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_vul_lep)
conf <- bootCI(vul_lep_spr, boost_vul_lep)
delta <- paste0("∆4: ", round(vul_lep_spr, 3)," ", "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)



plot(summer_vul, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_vulpes, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = FALSE)
plot(summer_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_lepus, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", "summer", col = "white", 
       box.col = "transparent", bg = "transparent", cex = size_legend, lty = 1)
vul_lep_sum <- overlapEst(vulpes_summer$clock, lepus_summer$clock, type = "Dhat4") # Overlapping coefficient Dht4
# boost_vul_lep_sum <- bootstrap(vulpes_summer$clock, lepus_summer$clock, 999, smooth = TRUE, type = "Dhat4")
mean_over <- mean(boost_vul_lep_sum)
conf <- bootCI(vul_lep_spr, boost_vul_lep_sum)
delta <- paste0("∆4: ", round(vul_lep_sum, 3)," ", "(", round(conf[1,1], 2), " ", round(conf[1,2], 2), ")")
text(x = 0, y = 0.5, labels = delta, cex = size_delta, font = 2, pos = 3)


plot(summer_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_summer, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(autumn_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_autumn, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
plot(winter_lep, data = "none", xunit="radians", yunit="density", ylim=c(0,0.8), 
     centre = "night", ylab="", xlab="", tline = list(col = col_winter, lwd = 3),
     dline = list(col = "transparent"),  cline = list(col = "transparent"), xaxis=list(xaxt="n"),
     add = TRUE)
draw_night_xaxis()
legend("topleft", c("Hare", "spring", "summer", "autumn", "winter"), col = c("white", col_spring, col_summer, col_autumn, col_winter),
       box.col = "transparent", bg = "transparent", cex = 1, lty = 1)

dev.off()
