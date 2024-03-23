
################################################################
# LIBRARIES
################################################################
library(stringr)
library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)
library(lme4)
library(sjPlot)
library(car)
library(ggplot2)
library(ggpubr)

# Folder path containing logger data
logger.path <- "~/Desktop/GCS/GCS_git/Loggers Data/"
setwd("~/Desktop/GCS/GCS_git/Loggers Data/")

################################################################
# FUNCTIONS
################################################################

read.loggers <- function(fn) {			# fn = filename
  print(paste("Reading ", fn))			# print which file is being read
  
  ###############################################################
  # Header
  ###############################################################
  # Use header from first line in original file to create new header names
  suppressMessages(header.part1 <- read_excel(fn, sheet=1, n_max = 2))
  
  # Provide ports (i.e., ports begin at 2nd column); results in "port1" etc.
  port <- paste0("port", str_extract(names(header.part1)[-1], "\\d"))
  
  # measurement (choosing last word of the measurement type, so "% VWC" will be VWC)
  measurement <- str_extract(header.part1[2, -1], '\\w+$')
  
  # Combine port with measurement, and use this as headers for the file:
  final.header <- paste(port, measurement, sep = ".")
  
  ###############################################################
  # Import
  ###############################################################
  # Create helper data set to establish how many columns are going to be read in
  suppressMessages(first.row <- read_excel(fn, sheet=1, col_names = F, skip = 11, n_max = 15, na = "#N/A!"))
  
  # Read in file; specify column types, with first column being "date" and everything else (i.e., length(first.row)-1) is numeric
  suppressMessages(data <- read_excel(fn, sheet=1, col_names = F, skip = 12, na = "#N/A!", col_types = c("date", rep("numeric", ncol(first.row)-1))))
  
  ###############################################################
  # Processing
  ###############################################################
  # Add new header names ("final.header") to data file
  names(data) <- c("date.time", final.header)
  
  # Convert %VWC to m3/m3 if present
  columns.to.convert <- any(grepl("%", header.part1[2,]))
  if(columns.to.convert) { # if %VWC is detected, divide those columns by 100
    cols <- grep("%", header.part1[2,])
    data[,cols] <- data[,cols]/100
  }
  
  # Remove any columns that have no data (NA) or only contain 0's
  new_data <- data %>% select_if(~ !(all(is.na(.) | . == 0)))
  
  # Sometimes previous year's environmental data is still present, so use only current year's data
  current.year <- sub(".*(\\d{4}).*", "\\1", basename(fn))

  # Select data from current year onwards (say, starting January 1):
  if(current.year != "2023") {
    new_data <- filter(new_data, date.time>= paste0(current.year, "-01-01") & date.time <= paste0(current.year, "-12-31"))
  } 
  # 2023 loggers were still in place from 2022, so 2022 needs to be included
  if(current.year == "2023") {
    new_data <- filter(new_data, date.time>= paste0(as.numeric(current.year)-1, "-01-01"))
  }
  
  # Melt data to long format
  data.melted <- 
    new_data %>%
    pivot_longer(!date.time, names_to = "sensor", values_to = "value", values_drop_na = TRUE)
 
  # Add site id to the long data
  data.melted$logger.name <- basename(fn)
  
  return(data = data.melted)
}

###############################################################
# DATA IMPORT
###############################################################
# List all the file names under the data folder "DataForPaper" (see path above)
file.names <- list.files(logger.path, pattern="*.xls*", recursive = T, full.names=T, include.dirs = T)

# Test individual files
#pos <- which(basename(file.names)=="BE191-2019.xls")
#pos <- which(basename(file.names)=="BE1-2021.xlsx")
#fn <- file.names[pos]

# Read in first file, the use a loop to read in the rest
df <- read.loggers(file.names[1])

for (i in 2:length(file.names)) {			# For each file name, read the data
  single.file <- read.loggers(file.names[i])	# Use FUN 2 to read the data within
  df <- rbind(df, single.file)
  # Keep track of what has been imported:
  if(i %% 10 == 0) {
    print(i)
  }
}

test <- df # copy to test, so I don't have to rerun the loop again (just in case)
###############################################################
# DATA PROCESSING
###############################################################
# Create columns for port number and sensor type
df$port <- gsub("\\..*","",df$sensor)
df$measurement <- gsub(".*\\.","",df$sensor)

# Remove battery percentage, pressure and voltage
df <- droplevels(filter(df, !(measurement %in% c("Percent", "Pressure", "Voltage"))))

# Rename VWC and Content to be just VWC
df$measurement <- factor(df$measurement)
levels(df$measurement) <- c("VWC", "EC", "Temperature", "Temperature", "VWC")

# Remove bad data (e.g., negative values, extreme high values, only zeroes or NA's)
df <- 
  df %>%
  filter(!(measurement == "EC")) %>%
  filter(!(measurement == "VWC" & (value <= 0 | value >= 0.9))) %>%
  filter(!(measurement == "Temperature" & (value <= -20 | value >= 70))) %>%
  filter(port != "port8") %>% # exclude port 8 (logger temperature)
  filter(!(is.na(value) | value == 0)) # Remove any zeros and NAs


###############################################################
# MERGE LOGGER DATA WITH FIELD INFO AND SOIL DEPTH
###############################################################
sites <-read_excel("~/Desktop/GCS/master.xlsx", sheet="Sites")
depth <-read_excel("~/Desktop/GCS/master.xlsx", sheet="sensors")
management <- read_excel("~/Desktop/GCS/master.xlsx", sheet="fields")

# Pre-processing:
# Re-name, then reorder factor levels
management$tillage <- factor(management$tillage)
management$irrigated <- factor(management$irrigated)
management$fertilizer <- factor(management$fertilizer)
management$crop.rotation <- factor(management$crop.rotation)
management$cover.crop <- factor(management$cover.crop)
management$residue <- factor(management$residue)

unique(management$tillage)

levels(management$tillage) <- c("minimal till","no-till", "tillage")
levels(management$irrigated) <- c("Dryland", "Irrigated")
levels(management$fertilizer) <- c("No fertilizer", "Fertilized")
levels(management$crop.rotation) <- c("Continuous cotton", "Crop rotation")
levels(management$cover.crop) <- c("Fallow", "Winter cover crop")
levels(management$residue) <- c("Medium-High residue", "None-Low residue")

sites$perc.clay.discr <- character(length = nrow(sites))
sites$perc.clay.discr[sites$perc.clay<=12] = "<12% clay"
sites$perc.clay.discr[(sites$perc.clay>12) & (sites$perc.clay<=25) ] = "12-25% clay"
sites$perc.clay.discr[sites$perc.clay>25] = ">25%clay"
sites$perc.clay.discr <- factor(sites$perc.clay.discr, levels = c("<12% clay", "12-25% clay", ">25%clay"))

# Select only important columns of sites and fields
sites <- 
  sites %>% 
  select(sample, logger.name, farm, field, year, perc.clay, perc.clay.discr)

management <- 
  management %>%
  select(farm:tillage, irrigated, fertilizer, crop.rotation, residue, cover.crop)

# Remove environmental data from "...-2023.xlsx" files
df.2017.2022 <- 
  df %>% 
  filter(!grepl("-2023", logger.name))

# select only rows that have a sensor name
sites <- filter(sites, !(is.na(logger.name)))

# MERGING:
master <- left_join(df.2017.2022, depth)
master <- left_join(master, sites)
master <- left_join(master, management)


###############################################################
# DAILY, MONTHLY AND ANNUAL TEMPERATURE & VWC AVERAGES
###############################################################

# Calculate average daily fluctuation for the year of each sensor
# (Note: port 8 is logger temperature, not soil temperature)

# Add year, month, day, then group by such that daily mean, min, max and range are calculated
daily.avg <- 
  master %>%
  filter(port != "port8") %>%
  mutate(year = year(date.time), month = month(date.time), day = day(date.time)) %>%
  group_by(logger.name, field, port, depth, measurement, year, month, day, perc.clay.discr, tillage, irrigated, fertilizer, crop.rotation, residue, cover.crop) %>%
  summarise(mean.value = mean(value, na.rm = T), min.value = min(value, na.rm = T), max.value = max(value, na.rm = T), range = max.value - min.value)

monthly.avg <- 
  daily.avg %>%
  group_by(logger.name, field, port, depth, measurement, year, month, perc.clay.discr, tillage, irrigated, fertilizer, crop.rotation, residue, cover.crop) %>%
  summarise(mean.value = mean(mean.value), min.value = mean(min.value), max.value = mean(max.value), range = mean(range), nDays = length(unique(day))) %>%
  arrange(year)

annual.avg <- 
  daily.avg %>%
  group_by(logger.name, field, port, depth, measurement, year, perc.clay.discr, tillage, irrigated, fertilizer, crop.rotation, residue, cover.crop) %>%
  summarise(mean.value = mean(mean.value), min.value = mean(min.value), max.value = mean(max.value), range = mean(range)) %>%
  arrange(year)

#write.csv(annual.avg, paste0(logger.path,"checking_ports.csv"), row.names = F)

###############################################################
# PREP FOR ANALYSIS
###############################################################
# Have temperature and VWC next to each other (so VWC can be used as a predictor)
# For July 26: normalize the monthly data for each year (to get rid of year-effects)
daily.avg.scaled <- 
  daily.avg %>%
  filter(month %in% 8:9) %>%
  group_by(depth, measurement, year) %>%
  mutate(mean.value = scale(mean.value), 
         min.value = scale(min.value), 
         max.value = scale(max.value), 
         range = scale(range))
  
# Test if all values are very close to zero
daily.avg.scaled %>%
  group_by(depth, measurement, year) %>%
  summarize(mean.sc = mean(mean.value), min.sc = mean(min.value), max.sc = mean(max.value), range.sc = mean(range), range.sd = sd(range)) %>%
  View()


monthly.avg.scaled <- 
  daily.avg.scaled %>%
  group_by(logger.name, field, port, depth, measurement, year, month, perc.clay.discr, tillage, irrigated, fertilizer, crop.rotation, residue, cover.crop) %>%
  summarise(mean.value = mean(mean.value), min.value = mean(min.value), max.value = mean(max.value), range = mean(range), nDays = length(unique(day)))

monthly.scaled <- 
  monthly.avg.scaled %>% 
  pivot_wider(names_from = measurement,
              values_from = c(mean.value, min.value, max.value, range, nDays)) %>%
  drop_na()

# Shorten the column names (e.g. min.value_Temperature becomes min.Temperature)
names(monthly.scaled) <- gsub("value_", "", names(monthly.scaled))

# Add unique field id
monthly.scaled$field.year <- with(monthly.scaled, paste(field, year, sep = "_"))

###############################################################
# BEST PREDICTORS
###############################################################

# Mixed-effects model for July - September 
  # 1) mean daily T (scaled)
  growing.season.meanT.15cm <- lmer(mean.Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                   data = subset(monthly.scaled,  (depth == 15)))
  growing.season.meanT.0cm <- lmer(mean.Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                   data = subset(monthly.scaled,  (depth == 0)))
  # 2) max daily T (scaled)
  growing.season.maxT.15cm <- lmer(max.Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                   data = subset(monthly.scaled,  (depth == 15)))
  growing.season.maxT.0cm <- lmer(max.Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                   data = subset(monthly.scaled,  (depth == 0)))
  # 3) minimum daily T (scaled)
  growing.season.minT.15cm <- lmer(min.Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                   data = subset(monthly.scaled,  (depth == 15)))
  growing.season.minT.0cm <- lmer(min.Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                  data = subset(monthly.scaled,  (depth == 0)))
  # 4) DTR (daily temperature range)
  growing.season.DTR.15cm <- lmer(range_Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                  data = subset(monthly.scaled,  (depth == 15)))
  growing.season.DTR.0cm <- lmer(range_Temperature ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                 data = subset(monthly.scaled,  (depth == 0)))

# Table S2
tab_model(
  growing.season.meanT.15cm, 
  growing.season.meanT.0cm, 
  growing.season.maxT.15cm, 
  growing.season.maxT.0cm,
  growing.season.minT.15cm, 
  growing.season.minT.0cm,
  growing.season.DTR.15cm,
  growing.season.DTR.0cm,
  auto.label = F)

Anova(growing.season.meanT.15cm)
summary(growing.season.meanT.15cm)
Anova(growing.season.meanT.0cm)
summary(growing.season.meanT.0cm)
Anova(growing.season.maxT.15cm)
Anova(growing.season.maxT.0cm)
Anova(growing.season.minT.15cm)
Anova(growing.season.minT.0cm)
Anova(growing.season.DTR.15cm)
summary(growing.season.DTR.15cm)
Anova(growing.season.DTR.0cm)
summary(growing.season.DTR.0cm)

# Mixed-effects model for July - September 
# 1) mean soil moisture
growing.season.vwc.15cm <- lmer(mean.VWC ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop  + perc.clay.discr + (1|field.year), 
                                  data = subset(monthly.scaled,  (depth == 15)))

growing.season.vwc.0cm <- lmer(mean.VWC ~ tillage + irrigated + fertilizer + crop.rotation + residue + cover.crop + perc.clay.discr + (1|field.year), 
                                 data = subset(monthly.scaled,  (depth == 0)))

tab_model(growing.season.vwc.15cm, growing.season.vwc.0cm)
Anova(growing.season.vwc.15cm)
summary(growing.season.vwc.15cm)
Anova(growing.season.vwc.0cm)
summary(growing.season.vwc.0cm)

# plots
# custom theme for plots
theme <- theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.background = element_rect(fill = NA, linewidth = 1),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        panel.border = element_rect(linewidth  = 1))

# VWC by tillage 

monthly.avg.VWC <- monthly.avg %>% 
  filter(measurement == "VWC") %>% 
  filter(month %in% c(6, 7,8, 9))

VWC.tillage <- monthly.avg.VWC %>% 
  filter(tillage != "NA") %>% 
  group_by(tillage, depth, month) %>% 
  summarise(mean.VWC = mean(mean.value))

VWC.tillage$tillage <-factor( VWC.tillage$tillage,
                              levels =  c("no-till", "minimal till", "tillage"))

(a <- ggplot(VWC.tillage, aes(month, mean.VWC)) +
  geom_point(aes(shape = tillage), size = 3) +
  geom_line(aes(lty = tillage), linewidth = 1) +
  facet_wrap(~ depth, ncol = 1,
             labeller = labeller(depth = c("0" = "0 cm", "15" = "15 cm"))) +
  scale_shape_manual(labels =  c("No-till", "Minimal-till","Till"),
                     values = c(15, 16, 17)) +
  scale_linetype_manual(labels =  c("No-till", "Minimal-till","Till"),
                        values = c(1,3,5)) +
  labs(x = "Months", y = expression("VWC ("*m^3/m^3*")")) +
  scale_x_continuous(labels = c("Jun", "July", "Aug", "Sep")) +
  scale_y_continuous(limits = c(0.05, 0.25)) +
  theme +
  theme(legend.position = "top"))

# VWC by crop rotation

VWC.cr <- monthly.avg.VWC %>% 
  filter(crop.rotation != "NA") %>% 
  group_by(crop.rotation, depth, month) %>% 
  summarise(mean.VWC = mean(mean.value))


(b <- ggplot(VWC.cr, aes(month, mean.VWC, color = crop.rotation)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_wrap(~ depth, ncol = 1,
             labeller = labeller(depth = c("0" = "0 cm", "15" = "15 cm"))) +
  scale_color_manual(labels = c("Continuous Cotton", "Crop Rotation"),
                     values = c("#606060", "#FF9933")) +
  labs(x = "Months", y = expression("VWC ("*m^3/m^3*")")) +
  scale_x_continuous(labels = c("Jun", "July", "Aug", "Sep")) +
  scale_y_continuous(limits = c(0.05, 0.25)) +
  theme +
  theme(legend.position = "top"))

# VWC by irrigation

VWC.irrigation <- monthly.avg.VWC %>% 
  filter(irrigated != "NA") %>% 
  group_by(irrigated, depth, month) %>% 
  summarise(mean.VWC = mean(mean.value))

(c <- ggplot(VWC.irrigation, aes(month, mean.VWC, color = irrigated)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_line(linewidth = 1, alpha = 0.8) +
  facet_wrap(~ depth, ncol =1,
             labeller = labeller(depth = c("0" = "0 cm", "15" = "15 cm"))) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                     values = c("#FF0000", "#0000FF")) +
  labs(x = "Months", y = expression("VWC ("*m^3/m^3*")")) +
  scale_x_continuous(labels = c("Jun", "July", "Aug", "Sep")) +
  scale_y_continuous(limits = c(0.05, 0.25)) +
  theme +
  theme(legend.position = "top"))


(VWC_combined <- ggarrange(a + rremove("xlab"),
          c + rremove("xlab") + rremove("ylab") + rremove("y.text"), 
          b + rremove("xlab") + rremove("ylab") + rremove("y.text"),
          labels = "auto",
          ncol = 3,
          align = "hv"))



# DTR vs moisture graph


(Dtr <- ggplot(monthly.scaled, aes(x = mean.VWC, y = range_Temperature, color = irrigated)) +
  geom_point(size = 3, alpha = 1) +
  geom_smooth(method = "lm") +
  facet_wrap(~depth,
             labeller = labeller(depth = c("0" = "0 cm", "15" = "15 cm"))) +
  scale_color_manual(labels = c("Dryland", "Irrigated"),
                     values = c("#FF0000", "#0000FF")) +
  labs(x = expression("VWC ("*m^3/m^3*")"), 
       y = expression("DTR ("*~degree*C*")"))+
  theme)


# Saving the graphs

setwd("~/Desktop/GCS/graphs/")

Fig_3 <- ggsave("Fig3.pdf", VWC_combined, dpi = 300, width = 15, height = 6 ,
                units = "in", device = "pdf")
Fig_S1 <- ggsave("S1.pdf", Dtr, dpi = 300, width = 12, height = 6 ,
                units = "in", device = "pdf")





