#This script processes GC data and graphs the trial data
#To run this script:
# 1. Add this file to the source folder with scripts for each day
# 2. The metadata is in the same colder
# 3. input contains all the fles and matches to the file names
# 4. output will conatin the combined data that is worked from

# 1. Setup ----

# 1.1 Libraries----
install.packages("devtools")
install.packages("rstudioapi")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("patchwork")
install.packages("openxlsx")
library(devtools)
library(rstudioapi)
library(ggplot2)
library(dplyr)
library(patchwork)
library(openxlsx)

# 1.2 Set wd and folders----
#set working directory to same folder as script
setwd(dirname(getActiveDocumentContext()$path))
output <- file.path("/Users/abisese/Desktop/Capstone/GC CODE/output")
input <- file.path("/Users/abisese/Desktop/Capstone/GC CODE/input")
gc_meta <- read.xlsx("/Users/abisese/Desktop/Capstone/GC CODE/GC_metadata.xlsx")

#list csvs in Folder
csv_files <- list.files(path = input, pattern = "\\.csv$", full.names = TRUE)
#get csv basenames
csv_names <- sub("\\.csv$", "", basename(csv_files))

# 1.3 Check Meta/Files (CHECKPOINT)----
for (i in 1:nrow(gc_meta)) {
  fname <- gc_meta$file.name[i]
  if (fname %in% csv_names) {
    message(paste("Match found for:", fname))} 
  else {
    message(paste("No match found for:", fname))}
}

# 1.4 Combining Files and Populating Expt Details----
comb.data <- bind_rows(lapply(1:nrow(gc_meta), function(i) {
  fname <- gc_meta$file.name[i]
  if (!(fname %in% csv_names)) return(NULL)
  temp <- read.csv(csv_files[csv_names == fname])
  #check for common columns
  needed <- c("Peak", "Time", "Area", "Height")
  missing <- setdiff(needed, names(temp))
  # add missing columns as NA (avoid error)
  if (length(missing) > 0) temp[missing] <- NA
  #add needed columns form gc_meta
  temp[, needed] %>% mutate(
      type      = gc_meta$type[i],
      voc.ul    = gc_meta$voc.ul[i],
      method    = gc_meta$method[i],
      mof.ex.s  = gc_meta$mof.ex.s[i],
      mof.mg    = gc_meta$mof.mg[i],
      species   = gc_meta$mof.type[i],
      activation = gc_meta$mof.active[i],
      mof.trial = ifelse(grepl("MOF[A-Za-z]+(\\d|BLANK)", fname), sub(".*MOF([A-Za-z]+)(\\d|BLANK).*", "\\1", fname), NA),
      file.name = fname
    )
}))

for(i in 1:nrow(comb.data[comb.data$type == "blank", ])) {
  new.row <- comb.data[comb.data$type == "blank", ][i, ]
  new.row$type <- "mof"
  comb.data <- bind_rows(comb.data, new.row)
}

#add 0.5 uL to the calibration data
comb.data$voc.ul <- comb.data$voc.ul + 0.55
#add moles of etOH column
comb.data$mol.voc <- (0.0544012/(15*10^-3))*comb.data$voc.ul
#comb.data <- comb.data[comb.data$mof.ex.s <= 1500,]
comb.data <- comb.data[comb.data$mof.mg <= 14,]

#write.csv(comb.data, file = file.path(output, "combined_gc_data.csv"), row.names = FALSE)

#comb.data <- read.csv("/Users/abisese/Desktop/Capstone/GC CODE/output/combined_gc_data.csv")

# 1.5 Individual Peak Files ----
EtOH <- comb.data[comb.data$Peak == 1,]
#write.csv(EtOH, file = file.path(output, "combined_etoh_peak.csv"), row.names = FALSE)
Hex <- comb.data[comb.data$Peak == 2,]
#write.csv(Hex, file = file.path(output, "combined_hexane_peak.csv"), row.names = FALSE)
Benz <- comb.data[comb.data$Peak == 3,]
#write.csv(Benz, file = file.path(output, "combined_benzene_peak.csv"), row.names = FALSE)

# 2. Gas Calibration Curves by Peak ----
# 2.1 EtOH Peak----
gas.EtOH.std <- EtOH[EtOH$type == "std" & EtOH$method == "gas", ]

#filter df and calculate SD and SE
cal.EtOH.gas <- gas.EtOH.std %>%
  dplyr::group_by(mol.voc) %>%
  dplyr::summarise(
    mean_area = mean(Area, na.rm = TRUE),
    n         = dplyr::n(),
    area_sd   = if (n > 1) sd(Area, na.rm = TRUE) else 0,
    area_se   = area_sd / sqrt(n),
    .groups   = "drop"
  )


#get regression line, R^2, and eqn
model <- lm(mean_area ~ mol.voc, data = cal.EtOH.gas)
eq <- paste0("y = ", round(coef(model)[2], 4), "x + ", round(coef(model)[1], 4), "\nR² = ", round(summary(model)$r.squared, 4))

#plot data
ggplot(cal.EtOH.gas, aes(x = mol.voc, y = mean_area)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_area - area_se,
                    ymax = mean_area + area_se),
                    width = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text",
           x = max(cal.EtOH.gas$mol.voc)*0.6,
           y = max(cal.EtOH.gas$mean_area)*0.9,
           label = eq,
           hjust = 0)+
  scale_x_continuous(
    breaks = seq(min(cal.EtOH.gas$mol.voc), max(cal.EtOH.gas$mol.voc), by = 4),
    minor_breaks = seq(min(cal.EtOH.gas$mol.voc), max(cal.EtOH.gas$mol.voc), by = 2)) +
  scale_y_continuous(
    breaks = seq(0, 1000, by = 150),
    minor_breaks = seq(0, max(cal.EtOH.gas$mean_area)+150, by = 50)) +
  labs(
    title = "Calibration Curve for EtOH Peak",
    x = "Moles Injected (mol)",
    y = "Peak Area") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line()
  )

ggplot(cal.EtOH.gas, aes(x = mol.voc, y = mean_area)) +
  geom_point(size = 3, color = "#5A102B") +
  geom_errorbar(aes(ymin = mean_area - area_se,
                    ymax = mean_area + area_se),
                width = 0.3) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1, color = "#CC5500") +
  annotate("text",
           x = max(cal.EtOH.gas$mol.voc) * 0.6,
           y = max(cal.EtOH.gas$mean_area) * 0.9,
           label = eq,
           hjust = 0) +
  scale_x_continuous(
    breaks = seq(min(cal.EtOH.gas$mol.voc),
                 max(cal.EtOH.gas$mol.voc), by = 4),
    minor_breaks = seq(min(cal.EtOH.gas$mol.voc),
                       max(cal.EtOH.gas$mol.voc), by = 2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(
    breaks = seq(0, 1000, by = 150),
    minor_breaks = seq(0,
                       max(cal.EtOH.gas$mean_area) + 150, by = 50),
    expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Calibration Curve for EtOH Peak",
       x = "Headspace Moles (mol)",
       y = "Peak Area") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )

# 2.2 Hexane Peak ----
gas.hex.std <- Hex[Hex$type == "std" & Hex$method == "gas", ]

#filter df and calculate SD and SE
cal.hex.gas <- gas.hex.std %>%
  dplyr::group_by(mol.voc) %>%
  dplyr::summarise(
    mean_area = mean(Area, na.rm = TRUE),
    n         = dplyr::n(),
    area_sd   = if (n > 1) sd(Area, na.rm = TRUE) else 0,
    area_se   = area_sd / sqrt(n),
    .groups   = "drop"
  )


#get regression line, R^2, and eqn
model <- lm(mean_area ~ mol.voc, data = cal.hex.gas)
eq <- paste0("y = ", round(coef(model)[2], 4), "x + ", round(coef(model)[1], 4), "\nR² = ", round(summary(model)$r.squared, 4))

#plot data
ggplot(cal.hex.gas, aes(x = mol.voc, y = mean_area)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_area - area_se,
                    ymax = mean_area + area_se),
                width = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text",
           x = max(cal.hex.gas$mol.voc)*0.6,
           y = max(cal.hex.gas$mean_area)*0.9,
           label = eq,
           hjust = 0)+
  scale_x_continuous(
    breaks = seq(min(cal.hex.gas$mol.voc), max(cal.hex.gas$mol.voc), by = 4),
    minor_breaks = seq(min(cal.hex.gas$mol.voc), max(cal.hex.gas$mol.voc), by = 2)) +
  scale_y_continuous(
    breaks = seq(0, 3000, by = 300),
    minor_breaks = seq(0, max(cal.hex.gas$mean_area), by = 50)) +
  labs(
    title = "Calibration Curve for Hexane Peak",
    x = "Moles Injected (mol)",
    y = "Peak Area") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()
  )

ggplot(cal.hex.gas, aes(x = mol.voc, y = mean_area)) +
  geom_point(size = 3, color = "#5A102B") +
  geom_errorbar(aes(ymin = mean_area - area_se, ymax = mean_area + area_se), width = 0.3) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1, color = "#CC5500") +
  annotate("text",
           x = max(cal.hex.gas$mol.voc) * 0.6,
           y = max(cal.hex.gas$mean_area) * 0.9,
           label = eq,
           hjust = 0) +
  scale_x_continuous(
    breaks = seq(min(cal.hex.gas$mol.voc),
                 max(cal.hex.gas$mol.voc), by = 4),
    minor_breaks = seq(min(cal.hex.gas$mol.voc),
                       max(cal.hex.gas$mol.voc), by = 2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(
    breaks = seq(0, 3500, by = 250),
    minor_breaks = seq(0,
                       max(cal.hex.gas$mean_area) + 150, by = 50),
    expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Calibration Curve for Hexane Peak",
       x = "Headspace Moles (mol)",
       y = "Peak Area") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )

# 2.2 Benzene Peak ----
gas.benz.std <- Benz[Benz$type == "std" & Benz$method == "gas", ]

#filter df and calculate SD and SE
cal.benz.gas <- gas.benz.std %>%
  dplyr::group_by(mol.voc) %>%
  dplyr::summarise(
    mean_area = mean(Area, na.rm = TRUE),
    n         = dplyr::n(),
    area_sd   = if (n > 1) sd(Area, na.rm = TRUE) else 0,
    area_se   = area_sd / sqrt(n),
    .groups   = "drop"
  )


#get regression line, R^2, and eqn
model <- lm(mean_area ~ mol.voc, data = cal.benz.gas)
eq <- paste0("y = ", round(coef(model)[2], 4), "x + ", round(coef(model)[1], 4), "\nR² = ", round(summary(model)$r.squared, 4))

#plot data
ggplot(cal.benz.gas, aes(x = mol.voc, y = mean_area)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_area - area_se,
                    ymax = mean_area + area_se),
                width = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text",
           x = max(cal.benz.gas$mol.voc)*0.6,
           y = max(cal.benz.gas$mean_area)*0.9,
           label = eq,
           hjust = 0)+
  scale_x_continuous(
    breaks = seq(min(cal.benz.gas$mol.voc), max(cal.benz.gas$mol.voc), by = 4),
    minor_breaks = seq(min(cal.benz.gas$mol.voc), max(cal.benz.gas$mol.voc), by = 2)) +
  scale_y_continuous(
    breaks = seq(0, 4000, by = 500),
    minor_breaks = seq(0, max(cal.benz.gas$mean_area), by = 50)) +
  labs(
    title = "Calibration Curve for Benzene Peak",
    x = "Moles Injected (nmol)",
    y = "Peak Area") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()
  )

ggplot(cal.benz.gas, aes(x = mol.voc, y = mean_area)) +
  geom_point(size = 3, color = "#5A102B") +
  geom_errorbar(aes(ymin = mean_area - area_se, ymax = mean_area + area_se), width = 0.3) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 1, color = "#CC5500") +
  annotate("text",
           x = max(cal.benz.gas$mol.voc) * 0.6,
           y = max(cal.benz.gas$mean_area) * 0.9,
           label = eq,
           hjust = 0) +
  scale_x_continuous(
    breaks = seq(min(cal.benz.gas$mol.voc),
                 max(cal.benz.gas$mol.voc), by = 4),
    minor_breaks = seq(min(cal.benz.gas$mol.voc),
                       max(cal.benz.gas$mol.voc), by = 2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(
    breaks = seq(0, 4000, by = 500),
    minor_breaks = seq(0,
                       max(cal.hex.gas$mean_area) + 150, by = 50),
    expand = expansion(mult = c(0.05, 0.05))) +
  labs(title = "Calibration Curve for Benzene Peak",
       x = "Headspace Moles (mol)",
       y = "Peak Area") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 1)
  )

# 4. MOF Capturing Over Time ----

# 4.1.1. EtOH Blanking ----
#get EtOH peak
# remove NA trials
EtOH <- EtOH[EtOH$method == "gas",]
EtOH.mof <- EtOH[!is.na(EtOH$mof.trial), ]
#EtOH.mof <- EtOH.mof[EtOH.mof$species == "amin",]


#Use "blank" and subtract form Area column
EtOH.mof.b <- EtOH.mof %>% group_by(mof.trial) %>% mutate(blank_val = dplyr::first(Area[type == "blank"], default = 0),
    Area = dplyr::if_else(type == "mof" & !is.na(Area), Area - blank_val, Area)) %>%
  select(-blank_val) %>%
  ungroup()

#remove blanks
EtOH.mof.b <- EtOH.mof.b[EtOH.mof.b$type == "mof",]
EtOH.mof.b <- EtOH.mof.b[EtOH.mof.b$species != "blank",]
EtOH.mof.b <- EtOH.mof.b[EtOH.mof.b$voc.ul == 2.55,]

# 4.1.2 EtOH Quantification of Volume in Headspace ----
EtOH.mof.b$headspace.mol <- (EtOH.mof.b$Area - (10.0441))/(45.6562)
EtOH.mof.b$abs.mol <- (((EtOH.mof.b$mol.voc)) - EtOH.mof.b$headspace.mol)
EtOH.mof.b$norm.mol <- EtOH.mof.b$abs.mol/EtOH.mof.b$mof.mg
EtOH.mof.b$abs.per <- EtOH.mof.b$abs.mol/(EtOH.mof.b$mol.voc)*100

# 4.1.2. Turn false blanks to zero values ----
for(i in 1:nrow(EtOH.mof.b)) {
  if(grepl("BLANK", EtOH.mof.b$file.name[i])) {
    EtOH.mof.b$abs.mol[i] <- 0
    EtOH.mof.b$norm.mol[i] <- 0}}

EtOH.mof.b <- EtOH.mof.b %>%filter(abs.mol >= 0)
# 4.1.3. EtOH Graphing ----

ggplot(EtOH.mof.b %>% filter(activation == "n"),
       aes(x = mof.ex.s, y = norm.mol, color = species, fill = species)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE,
              alpha = 0.10, linewidth = 1.3) +
  annotate("text", x = 1200, y = 0.55, label = "n = 23",
           color = "#5A102B", fontface = "bold", size = 5) +
  annotate("text", x = 1200, y = 0.91, label = "n = 41",
           color = "#CC5500", fontface = "bold", size = 5) +
  coord_cartesian(
    xlim = c(0,
             max(EtOH.mof.b$mof.ex.s[EtOH.mof.b$species == "amin" &
                                       EtOH.mof.b$activation == "n"])),
    ylim = c(0, 1.25),
    expand = FALSE) +
  scale_x_continuous(
    breaks = seq(0,
                 max(EtOH.mof.b$mof.ex.s[EtOH.mof.b$species == "amin" & EtOH.mof.b$activation == "n"]),
                 by = 100),
    expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1.25),
    breaks = seq(0, 1.25, by = 0.25),
    expand = c(0, 0)) +
  scale_color_manual(
    name = "Linker:",
    values = c("amin" = "#5A102B", "unmod" = "#CC5500"),
    labels = c("amin" = "Amine", "unmod" = "Native")) +
  scale_fill_manual(
    name = "Linker:",
    values = c("amin" = "#861F41", "unmod" = "#E5751F"),
    labels = c("amin" = "Amine", "unmod" = "Native")) +
  labs(
    title = "Ethanol Absorbance",
    x = "Time (s)",
    y = "EtOH Absorbed (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

ggplot(EtOH.mof.b %>% filter(activation == "AO"), aes(x = mof.ex.s, y = norm.mol, color = species, fill = species)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10, linewidth = 1.3) +
  annotate("text", x = 1200, y = 0.7, label = "n = 22", color = "#5A102B", fontface = "bold", size = 5) +
  annotate("text", x = 1200, y = 1.04, label = "n = 24", color = "#CC5500", fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0, max(EtOH.mof.b$mof.ex.s[EtOH.mof.b$species == "amin" & EtOH.mof.b$activation == "n"])),
                  ylim = c(0, 1.25),
                  expand = FALSE) +
  scale_x_continuous(breaks = seq(0,
                                  max(EtOH.mof.b$mof.ex.s[EtOH.mof.b$species == "amin" & EtOH.mof.b$activation == "n"]),
                                  by = 100),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.25),
                     breaks = seq(0, 1.25, by = 0.25),
                     expand = c(0, 0)) +
  scale_color_manual(
    name = "Linker",
    values = c("amin" = "#5A102B","unmod" = "#CC5500"),
    labels = c("amin" = "Amine","unmod" = "Native")) +
  scale_fill_manual(
    name = "Linker",
    values = c("amin" = "#861F41","unmod" = "#E5751F"),
    labels = c("amin" = "Amine","unmod" = "Native")) +
  labs(
    title = "Activated Ethanol Absorbance",
    x = "Time (s)",
    y = "EtOH Absorbed (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,face = "bold",size = 15),
    axis.title = element_text(face = "bold",size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold",size = 12,hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black",fill = "white",linewidth = 0.6),
    legend.box.background = element_rect(color = "black",fill = NA,linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

# 4.2.1. Hexane Blanking ----
#get Hexane peak
# remove NA trials
Hex <- Hex[Hex$method == "gas",]
Hex.mof <- Hex[!is.na(Hex$mof.trial), ]
#Hex.mof <- Hex.mof[Hex.mof$species == "amin",]

#Use "blank" and subtract form Area column
Hex.mof.b <- Hex.mof %>% group_by(mof.trial) %>% mutate(blank_val = dplyr::first(Area[type == "blank"], default = 0),
  Area = dplyr::if_else(type == "mof" & !is.na(Area), Area - blank_val, Area)) %>%
  select(-blank_val) %>%
  ungroup()

#remove blanks
Hex.mof.b <- Hex.mof.b[Hex.mof.b$type == "mof",]
Hex.mof.b <- Hex.mof.b[Hex.mof.b$species != "blank",]
Hex.mof.b <- Hex.mof.b[Hex.mof.b$voc.ul == 2.55,]

# 4.2.2 Hexane Quantification of Volume in Headspace ----
Hex.mof.b$headspace.mol <- (Hex.mof.b$Area - (15.3605))/(134.9891)
Hex.mof.b$abs.mol <- (((Hex.mof.b$mol.voc)) - Hex.mof.b$headspace.mol)
Hex.mof.b$norm.mol <- Hex.mof.b$abs.mol/Hex.mof.b$mof.mg
Hex.mof.b$abs.per <- Hex.mof.b$abs.mol/(Hex.mof.b$mol.voc)*100

# 4.1.3. Turn false blanks to zero values ----
for(i in 1:nrow(Hex.mof.b)) {
  if(grepl("BLANK", Hex.mof.b$file.name[i])) {
    Hex.mof.b$abs.mol[i] <- 0
    Hex.mof.b$norm.mol[i] <- 0}}

Hex.mof.b <- Hex.mof.b %>%filter(abs.mol >= 0)
# 4.2.4. Hexane Graphing ----

ggplot(Hex.mof.b %>% filter(activation == "AO"), aes(x = mof.ex.s, y = norm.mol, color = species, fill = species)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10, linewidth = 1.3) +
  annotate("text", x = 1200, y = 0.2, label = "n = 22", color = "#5A102B", fontface = "bold", size = 5) +
  annotate("text", x = 1200, y = 0.89, label = "n = 24", color = "#CC5500", fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0, max(Hex.mof.b$mof.ex.s[Hex.mof.b$species == "amin" & Hex.mof.b$activation == "n"])),
                  ylim = c(0, 1.25),
                  expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,
                                  max(Hex.mof.b$mof.ex.s[Hex.mof.b$species == "amin" & Hex.mof.b$activation == "n"]),
                                  by = 100),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1.25),
                     breaks = seq(0, 1.25, by = 0.25),
                     expand = c(0, 0)) +
  scale_color_manual(
    name = "Linker",
    values = c("amin" = "#5A102B", "unmod" = "#CC5500"),
    labels = c("amin" = "Amine", "unmod" = "Native")) +
  scale_fill_manual(
    name = "Linker",
    values = c("amin" = "#861F41", "unmod" = "#E5751F"),
    labels = c("amin" = "Amine", "unmod" = "Native")) +
  labs(
    title = "Activated Hexane Absorbance",
    x = "Time (s)",
    y = "Hexane Absorbed (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,face = "bold",size = 15),
    axis.title = element_text(face = "bold",size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold",size = 12,hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black",fill = "white",linewidth = 0.6),
    legend.box.background = element_rect(color = "black",fill = NA,linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

ggplot(Hex.mof.b %>% filter(activation == "n"), aes(x = mof.ex.s, y = norm.mol, color = species, fill = species)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10, linewidth = 1.3, fullrange = FALSE) +
  annotate("text", x = 1200, y = 0.10, label = "n = 20", color = "#5A102B", fontface = "bold", size = 5) +
  annotate("text", x = 1200, y = 0.79, label = "n = 41", color = "#CC5500", fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0, max(Hex.mof.b$mof.ex.s[Hex.mof.b$species == "amin" & Hex.mof.b$activation == "n"])),
                  ylim = c(0, 1.25),
                  expand = FALSE) +
  scale_x_continuous(breaks = seq(0,
                                  max(Hex.mof.b$mof.ex.s[Hex.mof.b$species == "amin" & Hex.mof.b$activation == "n"]),
                                  by = 100),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1.25),
                     breaks = seq(0, 1.25, by = 0.25),
                     expand = c(0, 0)) +
  scale_color_manual(
    name = "Linker",
    values = c("amin" = "#5A102B", "unmod" = "#CC5500"),
    labels = c("amin" = "Amine", "unmod" = "Native")) +
  scale_fill_manual(
    name = "Linker",
    values = c("amin" = "#861F41", "unmod" = "#E5751F"),
    labels = c("amin" = "Amine", "unmod" = "Native")) +
  labs(
    title = "Hexane Absorbance",
    x = "Time (s)",
    y = "Hexane Absorbed (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,face = "bold",size = 15),
    axis.title = element_text(face = "bold",size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold",size = 12,hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black",fill = "white",linewidth = 0.6),
    legend.box.background = element_rect(color = "black",fill = NA,linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

# 4.3.1. Benzene Blanking ----
#get Benz peak
# remove NA trials
Benz <- Benz[Benz$method == "gas",]
Benz.mof <- Benz[!is.na(Benz$mof.trial), ]

#Use "blank" and subtract form Area column
Benz.mof.b <- Benz.mof %>% group_by(mof.trial) %>% mutate(blank_val = dplyr::first(Area[type == "blank"], default = 0),
  Area = dplyr::if_else(type == "mof" & !is.na(Area), Area - blank_val, Area)) %>%
  select(-blank_val) %>%
  ungroup()

#remove blanks
Benz.mof.b <- Benz.mof.b[Benz.mof.b$type == "mof",]
Benz.mof.b <- Benz.mof.b[Benz.mof.b$voc.ul == 2.55,]
Benz.mof.b <- Benz.mof.b[Benz.mof.b$species != "blank",]

# 4.3.2. Benzene Quantification of Volume in Headspace ----
Benz.mof.b$headspace.mol <- (Benz.mof.b$Area - (18.3582))/(186.3952)
Benz.mof.b$abs.mol <- (((Benz.mof.b$mol.voc)) - Benz.mof.b$headspace.mol)
Benz.mof.b$norm.mol <- Benz.mof.b$abs.mol/Benz.mof.b$mof.mg
Benz.mof.b$abs.per <- Benz.mof.b$abs.mol/(Benz.mof.b$mol.voc)*100

# 4.1.3. Turn false blanks to zero values ----
for(i in 1:nrow(Benz.mof.b)) {
  if(grepl("BLANK", Benz.mof.b$file.name[i])) {
    Benz.mof.b$abs.mol[i] <- 0
    Benz.mof.b$norm.mol[i] <- 0}}

Benz.mof.b <- Benz.mof.b %>%filter(abs.mol >= 0)
# 4.3.4. Benzene Graphing ----

ggplot(Benz.mof.b %>% filter(activation == "AO"), aes(x = mof.ex.s, y = norm.mol, color = species, fill = species)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10, linewidth = 1.3) +
  annotate("text", x = 1200, y = 0.51, label = "n = 22", color = "#5A102B", fontface = "bold", size = 5) +
  annotate("text", x = 1200, y = 0.96, label = "n = 24", color = "#CC5500", fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0, max(Benz.mof.b$mof.ex.s[Benz.mof.b$species == "amin" & Benz.mof.b$activation == "AO"])),
                  ylim = c(0, 1.25),
                  expand = FALSE) +
  scale_x_continuous(breaks = seq(0,
                                  max(Benz.mof.b$mof.ex.s[Benz.mof.b$species == "amin" & Benz.mof.b$activation == "AO"]),
                                  by = 100),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.5),
                     breaks = seq(0, 1.25, by = 0.25),
                     expand = c(0, 0)) +
  scale_color_manual(
    name = "Linker",
    values = c("amin" = "#5A102B","unmod" = "#CC5500"),
    labels = c("amin" = "Amine","unmod" = "Native")) +
  scale_fill_manual(
    name = "Linker",
    values = c("amin" = "#861F41","unmod" = "#E5751F"),
    labels = c("amin" = "Amine","unmod" = "Native")) +
  labs(
    title = "Activated Benzene Absorbance",
    x = "Time (s)",
    y = "Benzene Absorbed (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,face = "bold",size = 15),
    axis.title = element_text(face = "bold",size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold",size = 12,hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black",fill = "white",linewidth = 0.6),
    legend.box.background = element_rect(color = "black",fill = NA,linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

ggplot(Benz.mof.b %>% filter(activation == "n"), aes(x = mof.ex.s, y = norm.mol, color = species, fill = species)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10, linewidth = 1.3) +
  annotate("text", x = 1200, y = 0.22, label = "n = 23", color = "#5A102B", fontface = "bold", size = 5) +
  annotate("text", x = 1200, y = 0.86, label = "n = 41", color = "#CC5500", fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0, max(Benz.mof.b$mof.ex.s[Benz.mof.b$species == "amin" & Benz.mof.b$activation == "n"])),
                  ylim = c(0, 1.25),
                  expand = FALSE) +
  scale_x_continuous(breaks = seq(0,
                                  max(Benz.mof.b$mof.ex.s[Benz.mof.b$species == "amin" & Benz.mof.b$activation == "n"]),
                                  by = 100),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1.25),
                     breaks = seq(0, 1.25, by = 0.25),
                     expand = c(0, 0)) +
  scale_color_manual(
    name = "Linker",
    values = c("amin" = "#5A102B","unmod" = "#CC5500"),
    labels = c("amin" = "Amine","unmod" = "Native")) +
  scale_fill_manual(
    name = "Linker",
    values = c("amin" = "#861F41","unmod" = "#E5751F"),
    labels = c("amin" = "Amine","unmod" = "Native")) +
  labs(
    title = "Benzene Absorbance",
    x = "Time (s)",
    y = "Benzene Absorbed (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5,face = "bold",size = 15),
    axis.title = element_text(face = "bold",size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold",size = 12,hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black",fill = "white",linewidth = 0.6),
    legend.box.background = element_rect(color = "black",fill = NA,linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black",fill = NA,linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

# 5. Activation Changes ----
# 5.1.1 EtOH Peak ----

ggplot(EtOH.mof.b %>% filter(species == "unmod"), aes(x = mof.ex.s, y = norm.mol, colour = activation, fill = activation)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10, linewidth = 1.3) +
  # annotate("text", x = 1050, y = 0.18, label = "Not Activated", color = "#5A102B", fontface = "bold", size = 5) +
  # annotate("text", x = 1450, y = 0.88, label = "Oven Activated", color = "#CC5500", fontface = "bold", size = 5) +
  # annotate("text", x = 1250, y = 0.55, label = "Vacuum Oven Activated", color = "#1F5AA6", fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0, 1200),
                  ylim = c(0, 1.25),
                  expand = FALSE) +
  scale_x_continuous(breaks = seq(0,1300,by = 100),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1.25),
                     breaks = seq(0, 1.25, by = 0.25),
                     expand = c(0, 0)) +
  scale_colour_manual(
    name = "Activation",
    values = c("n" = "#5A102B", "AO" = "#006400", "AV" = "#1F5AA6"),
    labels = c("n" = "Not Activated", "AO" = "Oven Activated", "AV" = "Vacuum Oven Activated")) +
  scale_fill_manual(
    name = "Activation",
    values = c("n" = "#861F41", "AO" = "#007D00", "AV" = "#4F83CC"),
    labels = c("n" = "Not Activated", "AO" = "Oven Activated", "AV" = "Vacuum Oven Activated")) +
  labs(
    title = "Native Aborbance of EtOH",
    x = "Time (s)",
    y = "Absorbance (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

ggplot(EtOH.mof.b %>% filter(species == "amin"), aes(x = mof.ex.s, y = norm.mol, colour = activation, fill = activation)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10, linewidth = 1.3) +
  # annotate("text", x = 1050, y = 0.18, label = "Not Activated", color = "#5A102B", fontface = "bold", size = 5) +
  # annotate("text", x = 1450, y = 0.88, label = "Oven Activated", color = "#CC5500", fontface = "bold", size = 5) +
  # annotate("text", x = 1250, y = 0.55, label = "Vacuum Oven Activated", color = "#1F5AA6", fontface = "bold", size = 5) +
  coord_cartesian(xlim = c(0, 1300),
                  ylim = c(0, 1.25),
                  expand = FALSE) +
  scale_x_continuous(breaks = seq(0,1300,by = 100),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1, 1.25),
                     breaks = seq(0, 1.25, by = 0.25),
                     expand = c(0, 0)) +
  scale_colour_manual(
    name = "Activation",
    values = c("n" = "#5A102B", "AO" = "#006400", "AV" = "#1F5AA6"),
    labels = c("n" = "Not Activated", "AO" = "Oven Activated", "AV" = "Vacuum Oven Activated")) +
  scale_fill_manual(
    name = "Activation",
    values = c("n" = "#861F41", "AO" = "#007D00", "AV" = "#4F83CC"),
    labels = c("n" = "Not Activated", "AO" = "Oven Activated", "AV" = "Vacuum Oven Activated")) +
  labs(
    title = "Amine Absorbance of EtOH",
    x = "Time (s)",
    y = "Absorbance (mol/mg)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.6),
    legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

# 5.2.1 Hexane Peak ----

ggplot(Hex.mof.b %>% filter(species == "unmod"), aes(x = mof.ex.s, y = norm.nmol, color = activation)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10,linewidth = 1.3) +
  labs(
    title = "Absorbance Over Time of Hexane",
    x = "Time (s)",
    y = "Absorbed Volume (uL/mg)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()
  )


# 5.3.1 Benzene Peak ----
Benz.mof.b <- Benz.mof.b[Benz.mof.b$species == "amin",]

ggplot(Benz.mof.b %>% filter(species == "unmod"), aes(x = mof.ex.s, y = norm.nmol, colour = activation)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.10,linewidth = 1.3) +
  labs(
    title = "Absorbance Over Time of Benzene",
    x = "Time (s)",
    y = "Absorbed Volume (uL/mg)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line()
  )

# 6.1. Overlaid Signals for Presentation ----
ggplot(signals, aes(x = cal.t, y = cal.int)) +
  geom_line(color = "darkblue") +
  scale_x_continuous(limits = c(1.7, 4.3),labels = scales::label_number(accuracy = 0.1),expand = c(0, 0)) +
  labs(
    title = NULL,
    x = "Time (min)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )

# 7.1. Percentages for Cool Graph----

EtOH.per.aa <- EtOH.mof.b %>% filter(EtOH.mof.b$mof.ex.s >= 500 & EtOH.mof.b$species == "amin" & EtOH.mof.b$activation == "AO")
av.mol.abs <- mean(EtOH.per.aa$abs.mol)

EtOH.per.an <- EtOH.mof.b %>% filter(EtOH.mof.b$mof.ex.s >= 500 & EtOH.mof.b$species == "amin" & EtOH.mof.b$activation == "n")
av.mol.abs <- mean(EtOH.per.an$abs.mol)

EtOH.per.na <- EtOH.mof.b %>% filter(EtOH.mof.b$mof.ex.s >= 500 & EtOH.mof.b$species == "unmod" & EtOH.mof.b$activation == "AO")
av.mol.abs <- mean(EtOH.per.na$abs.mol)

EtOH.per.nn <- EtOH.mof.b %>% filter(EtOH.mof.b$mof.ex.s >= 500 & EtOH.mof.b$species == "unmod" & EtOH.mof.b$activation == "n")
av.mol.abs <- mean(EtOH.per.nn$abs.mol)

Hex.per.aa <- Hex.mof.b %>% filter(Hex.mof.b$mof.ex.s >= 500 & Hex.mof.b$species == "amin" & Hex.mof.b$activation == "AO")
av.mol.abs <- mean(Hex.per.aa$abs.mol)

Hex.per.an <- Hex.mof.b %>% filter(Hex.mof.b$mof.ex.s >= 500 & Hex.mof.b$species == "amin" & Hex.mof.b$activation == "n")
av.mol.abs <- mean(Hex.per.an$abs.mol)

Hex.per.na <- Hex.mof.b %>% filter(Hex.mof.b$mof.ex.s >= 500 & Hex.mof.b$species == "unmod" & Hex.mof.b$activation == "AO")
av.mol.abs <- mean(Hex.per.na$abs.mol)

Hex.per.nn <- Hex.mof.b %>% filter(Hex.mof.b$mof.ex.s >= 500 & Hex.mof.b$species == "unmod" & Hex.mof.b$activation == "n")
av.mol.abs <- mean(Hex.per.nn$abs.mol)

Benz.per.aa <- Benz.mof.b %>% filter(Benz.mof.b$mof.ex.s >= 500 & Benz.mof.b$species == "amin" & Benz.mof.b$activation == "AO")
av.mol.abs <- mean(Benz.per.aa$abs.mol)

Benz.per.an <- Benz.mof.b %>% filter(Benz.mof.b$mof.ex.s >= 500 & Benz.mof.b$species == "amin" & Benz.mof.b$activation == "n")
av.mol.abs <- mean(Benz.per.an$abs.mol)

Benz.per.na <- Benz.mof.b %>% filter(Benz.mof.b$mof.ex.s >= 500 & Benz.mof.b$species == "unmod" & Benz.mof.b$activation == "AO")
av.mol.abs <- mean(Benz.per.na$abs.mol)

Benz.per.nn <- Benz.mof.b %>% filter(Benz.mof.b$mof.ex.s >= 500 & Benz.mof.b$species == "unmod" & Benz.mof.b$activation == "n")
av.mol.abs <- mean(Benz.per.nn$abs.mol)

  
