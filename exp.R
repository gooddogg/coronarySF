####################################################################################################################################
#################################################EXPERIMENT 1#######################################################################
###############################################DATA ACQUISITION#####################################################################
####################################################################################################################################
#Preparations ####
#Remove riff raff from R
rm(list = ls())
#Set working directory
setwd("C:/Users/nickl/OneDrive/Skrivbord/Uni/MASTER THESIS/R")
#Load relevant libraries
#install.packages("plyr")
#install.packages("tidyverse")
#install.packages("here")
#install.packages("lubridate")
#install.packages("gridExtra")
#install.packages("corrplot")
#install.packages("DescTools")
#install.packages("kimisc")
library(plyr)
library(tidyverse)
detach("package:here", unload = TRUE)
library(here)
library(lubridate)
library(grid)
library(gridExtra)
library(corrplot)
library(DescTools)
library(kimisc)

if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}
# ####

#Load fish data ####
fish_dat <- read_csv2("fish_dat.csv") #basic data
fish_timeline <- read_csv2("fish_timeline.csv") #timeline of experimental protocols
fish_meta <- read_csv2("fish_meta.csv") #metadata
#filter(fish_meta, variable == "fish_dat$fishID") #replace the string within citation marks with the variable you want to check out
# ####

#Load smhi data on atmospheric pressure and convert the data into suitable formats; load constants; add day of week info####
#smhi data on atmospheric pressure
smhi_atm <- read_csv2("smhi_atm.csv")
smhi_atm$date_num <- as.numeric(smhi_atm$date) #Turn the dates into a numeric value using lubridate package
smhi_atm$atm_pr_kpa <- smhi_atm$atm_pr_hpa / 10 #Convert the atmospheric pressure into kPa
smhi_atm$time_uct_num <- as.numeric(smhi_atm$time_uct) #transform time in smhi_atm into a continuous variable as seconds into the day

#Load constants and turn them into objects for calculations of oxygen solubility later
cons <- read_csv2(here("mo2r2/constants.csv"))
cons_1 <- cons %>% filter(constant == "cons_1") %>% pull(value)
cons_2 <- cons %>% filter(constant == "cons_2") %>% pull(value)
cons_3 <- cons %>% filter(constant == "cons_3") %>% pull(value)
cons_4 <- cons %>% filter(constant == "cons_4") %>% pull(value)
cons_5 <- cons %>% filter(constant == "cons_5") %>% pull(value)
cons_6 <- cons %>% filter(constant == "cons_6") %>% pull(value)
cons_7 <- cons %>% filter(constant == "cons_7") %>% pull(value)
cons_8 <- cons %>% filter(constant == "cons_8") %>% pull(value)
cons_9 <- cons %>% filter(constant == "cons_9") %>% pull(value)
cons_10 <- cons %>% filter(constant == "cons_10") %>% pull(value)
cons_11 <- cons %>% filter(constant == "cons_11") %>% pull(value)
cons_12 <- cons %>% filter(constant == "cons_12") %>% pull(value)
cons_13 <- cons %>% filter(constant == "cons_13") %>% pull(value)
cons_14 <- cons %>% filter(constant == "cons_14") %>% pull(value)
cons_15 <- cons %>% filter(constant == "cons_15") %>% pull(value)
cons_16 <- cons %>% filter(constant == "cons_16") %>% pull(value)
cons_17 <- cons %>% filter(constant == "cons_17") %>% pull(value)
cons_18 <- cons %>% filter(constant == "cons_18") %>% pull(value)
cons_19 <- cons %>% filter(constant == "cons_19") %>% pull(value)
cons_20 <- cons %>% filter(constant == "cons_20") %>% pull(value)
cons_21 <- cons %>% filter(constant == "cons_21") %>% pull(value)
cons_22 <- cons %>% filter(constant == "cons_22") %>% pull(value)
cons_23 <- cons %>% filter(constant == "cons_23") %>% pull(value)
cons_24 <- cons %>% filter(constant == "cons_24") %>% pull(value)
cons_25 <- cons %>% filter(constant == "cons_25") %>% pull(value)
cons_26 <- cons %>% filter(constant == "cons_26") %>% pull(value)
cons_27 <- cons %>% filter(constant == "cons_27") %>% pull(value)
cons_28 <- cons %>% filter(constant == "cons_28") %>% pull(value)
cons_29 <- cons %>% filter(constant == "cons_29") %>% pull(value)
cons_30 <- cons %>% filter(constant == "cons_30") %>% pull(value)
cons_31 <- cons %>% filter(constant == "cons_31") %>% pull(value)
cons_32 <- cons %>% filter(constant == "cons_32") %>% pull(value)
cons_kpa <- 101.3
cons_mmhg <- 760

#Add what day and week from the start the experiments started
fish_dat$date_start_num <- as.numeric(fish_dat$date_start) #turn dates into numeric values
fish_dat$day_start <- fish_dat$date_start_num - min(pull(fish_dat, date_start_num), na.rm = T) + 1 #calculate which day from the first day each experiment started
fish_dat$wk <- floor(1 + (fish_dat$day_start / 7)) #calculate which week the experiment was performed
# ####

#Calculate morphological, hematological and general variables####
fish_dat <- mutate(fish_dat, bodymass = bodymass / 1000) #body mass in kg
fish_dat <- mutate(fish_dat, cf = 1000 * bodymass / (length^3) * 100) #condition factor
fish_dat <- mutate(fish_dat, relwetvent = (wetventmass / (1000* bodymass)) * 100) #relative wet ventricular mass
fish_dat <- mutate(fish_dat, dryventmass = compact + spongy) #dry ventricular mass
fish_dat <- mutate(fish_dat, relcomp = (compact / dryventmass) * 100) #relative compact mass
fish_dat <- mutate(fish_dat, relspleen = (spleenmass / (1000* bodymass)) * 100) #relative spleen mass
fish_dat <- fish_dat %>% rowwise() %>% mutate(hb_av = mean(c(hb1, hb2), na.rm = T)) #mean [Hb]
fish_dat <- mutate(fish_dat, hbcorr = (0.815 * hb_av) - 2.198) #corrected [Hb]
drabkins <- function(x){
  Whb <- 64458
  Fd <- 2.51 / 0.01
  Ce <- 44
  d <- 1
  
  x <- mutate(x, hb_drab = (abs * Whb * Fd) / (Ce * d * 1000))
  
  return(x)
} #Function to calculate [Hb] by drabkins method
fish_dat <- drabkins(fish_dat)
fish_dat <- fish_dat %>% rowwise() %>% mutate(hct_av = mean(c(hct1, hct2), na.rm = T)) #mean Hct
fish_dat <- fish_dat %>% rowwise() %>% mutate(osm_av = mean(c(osm1, osm2), na.rm = T)) #mean osmolality
fish_dat <- mutate(fish_dat, reldryvent = (dryventmass / (1000* bodymass)) * 100) #relative dry ventricular mass
fish_dat <- mutate(fish_dat, respfish_ratio = 10 / bodymass) #respirometer volume:fish mass ratio
fish_dat <- mutate(fish_dat, respvol_net = 10 - bodymass) #net respirometer volume
fish_dat <- fish_dat %>% rowwise() %>% mutate(cK_av = mean(c(cK1, cK2 / c2Diss), na.rm = T)) #average [K] of the two samples
fish_dat <- fish_dat %>% rowwise() %>% mutate(cNa_av = mean(c(cNa1, cNa2 / c2Diss), na.rm = T)) #average [Na] of the two samples
fish_dat <- fish_dat %>% rowwise() %>% mutate(cCl_av = mean(c(cCl1, cCl2 / c2Diss), na.rm = T)) #average [Cl] of the two samples
fish_dat <- fish_dat %>% rowwise() %>% mutate(cCa_av = mean(c(cCa1, cCa2 / c2Diss), na.rm = T)) #average [Ca] of the two samples
fish_dat <- fish_dat %>% rowwise() %>% mutate(cpH_av = mean(c(cpH1, cpH2 / c2Diss), na.rm = T)) #not so sure about this one, but supposed to be average pH of the two samples
# ####

#Load transonic flow probe calibration data and calculate calibration equations####
#Calculate the required data
trans_cal <- read_csv2("trans_cal.csv")
trans_cal <- trans_cal %>% rowwise() %>% mutate(base_av = mean(c(base_pre, base_post), na.rm = T))
trans_cal <- mutate(trans_cal, trans_corr = trans_flow - base_av)
trans_cal <- mutate(trans_cal, grav_minute = grav_flow / (duration_s / 60))
trans_cal <- mutate(trans_cal, grav_minute_rough = grav_flow / (60 / 60))

#Get the slopes and intercept of a linear regression of gravimetric flow ~ transonic flow (corrected) for each probe
lm_t732 <- summary(lm(filter(trans_cal, transonic == "t732")$grav_minute ~ filter(trans_cal, transonic == "t732")$trans_corr))
lm_t792 <- summary(lm(filter(trans_cal, transonic == "t792")$grav_minute ~ filter(trans_cal, transonic == "t792")$trans_corr))
lm_t794 <- summary(lm(filter(trans_cal, transonic == "t794")$grav_minute ~ filter(trans_cal, transonic == "t794")$trans_corr))
lm_t795 <- summary(lm(filter(trans_cal, transonic == "t795")$grav_minute ~ filter(trans_cal, transonic == "t795")$trans_corr))
lm_t771 <- summary(lm(filter(trans_cal, transonic == "t771")$grav_minute ~ filter(trans_cal, transonic == "t771")$trans_corr))
lm_t766 <- summary(lm(filter(trans_cal, transonic == "t766")$grav_minute ~ filter(trans_cal, transonic == "t766")$trans_corr))

#Gather the slopes and intercept in a new object
trans_eq <- setNames(data.frame(matrix(ncol = 3, nrow = length(unique(trans_cal$transonic)))), c("transonic", "slope", "intercept"))
trans_eq$transonic <- unique(trans_cal$transonic)
trans_eq$slope <- c(lm_t732$coefficients[2,1], lm_t792$coefficients[2,1], lm_t794$coefficients[2,1], lm_t795$coefficients[2,1], lm_t771$coefficients[2,1], lm_t766$coefficients[2,1])
trans_eq$intercept <- c(lm_t732$coefficients[1,1], lm_t792$coefficients[1,1], lm_t794$coefficients[1,1], lm_t795$coefficients[1,1], lm_t771$coefficients[1,1], lm_t766$coefficients[1,1])
# ####

#Add functions for exploring raw data and computing MO2, CO, HR, SV, EO2 for each fish and adding them into fish_dat####
sem <- function(x, na.rm = TRUE){
  sd(x, na.rm = na.rm) / sqrt(length(x))
} #function for calculating standard error of the mean

xplr <- function(x){
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  treat <- fish_dat[deparse(substitute(x)), "treatment"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  x <- filter(x, resp != "bkg_pre" & resp != "bkg_post")
  x$slope_airsat_inv <- -1 * x$slope_airsat #inverse the slope to make it more intuitive to look at
  
  #Have a look at the data to get a clear view of which values are out of line
  co <- arrange(x, co)$co
  hr <- arrange(x, hr)$hr
  slope_airsat_inv <- arrange(x, slope_airsat_inv)$slope_airsat_inv
  lotohi <- data.frame(co, hr, slope_airsat_inv)
  
  #Display max and min values (without background respiration)
  mami <- select(x, co, hr, slope_airsat_inv) %>% summarise_all(.funs = list(Min = ~ min(., na.rm = T), Max = ~ max(., na.rm = T), SEM = ~ sem(., na.rm = T)))
  mami$treat <- treat
  
  #Plot histograms of the raw data (without background respiration) to see if anything looks fishy
  histo <- ggplot(gather(select(x, co, hr, slope_airsat_inv), key = "var", value = "value"), aes(value)) +
    geom_histogram(bins = 50, fill = "black") +
    facet_wrap(~var, scales = 'free') +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 1), plot.title = element_text(hjust = 0.5))
  
  display_lotohi <- readline(prompt = "Display lotohi (Y/N)? ")
  
  if(tolower(display_lotohi) == "y"){
    View(lotohi) #displays the raw variables each ordered from low to high; but only if you input "y" or "Y" at the prompt
  }
  
  return(list(mami, histo))
} #to explore raw data in "fX" to see if there's anything fishy with the data

find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
} #to create colored tables in xplm()

mfits <- function(x){
  fishID <- deparse(substitute(x))
  x_sub <- filter(select(mutate(x, time_h_ca = round_any(time_min / 60, 0.05), fishID = fishID, time_post = time_min - min(filter(x, resp == "post_chase")$time_min)),
                         fishID, resp, time_h_ca, time_min, time_post, co_kg, eo2, do2_c, hr, sv), 
                  time_h_ca > 0)
  
  return(x_sub)
} #subset the data - for usage in xplm, tlco, tleo, tlmo, tlhr & tlsv

xplm <- function(x){
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  fishID <- fish_dat[deparse(substitute(x)), 0]
  treat <- fish_dat[deparse(substitute(x)), "treatment"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  chase_h <- min(filter(mutate(x, time_h = time_min / 60), resp == "post_chase")$time_h)
  
  co_plot <- ggplot(mutate(x, time_h = time_min / 60), aes(x = time_h, y = co_kg)) +
    geom_point(size = 0.25) +
    geom_point(data = filter(mutate(x, time_h = time_min / 60), co_kg == max(pull(x, co_kg), na.rm = T)), aes(x = time_h, y = co_kg), color = "red", size = 3) +
    geom_point(data = select(filter(mutate(x, time_h = time_min / 60), do2_c == max(pull(filter(x, resp == "post_chase"), do2_c), na.rm = T)), time_h, co_kg), color = "yellow", size = 2) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), time_post == 0 | 
                                      time_post == Closest(mfits(x)$time_post, 10) | 
                                      time_post == Closest(mfits(x)$time_post, 20) | 
                                      time_post == Closest(mfits(x)$time_post, 30) | 
                                      time_post == Closest(mfits(x)$time_post, 60) | 
                                      time_post == Closest(mfits(x)$time_post, 120) |
                                      time_post == Closest(mfits(x)$time_post, 300) |
                                      time_post == Closest(mfits(x)$time_post, 600) |
                                      time_post == Closest(mfits(x)$time_post, 900) ), co_kg, time_h), color = "magenta1", size = 1) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), co_kg, time_h), color = "skyblue1", size = 0.25) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200 & do2_c < quantile(pull(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), do2_c), probs = 0.1, na.rm = T)), co_kg, time_h), color = "blue", size = 1) +
    geom_vline(xintercept = min(filter(mutate(x, time_h = time_min / 60), resp == "post_chase")$time_h))
  
  eo2_plot <- ggplot(mutate(x, time_h = time_min / 60), aes(x = time_h, y = eo2)) +
    geom_point(size = 0.25) +
    geom_point(data = filter(mutate(x, time_h = time_min / 60), eo2 == max(pull(x, eo2), na.rm = T)), aes(x = time_h, y = eo2), color = "red", size = 3) +
    geom_point(data = select(filter(mutate(x, time_h = time_min / 60), do2_c == max(pull(filter(x, resp == "post_chase"), do2_c), na.rm = T)), time_h, eo2), color = "yellow", size = 2) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), time_post == 0 | 
                                      time_post == Closest(mfits(x)$time_post, 10) | 
                                      time_post == Closest(mfits(x)$time_post, 20) | 
                                      time_post == Closest(mfits(x)$time_post, 30) | 
                                      time_post == Closest(mfits(x)$time_post, 60) | 
                                      time_post == Closest(mfits(x)$time_post, 120) |
                                      time_post == Closest(mfits(x)$time_post, 300) |
                                      time_post == Closest(mfits(x)$time_post, 600) |
                                      time_post == Closest(mfits(x)$time_post, 900) ), eo2, time_h), color = "magenta1", size = 1) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), eo2, time_h), color = "skyblue1", size = 0.25) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200 & do2_c < quantile(pull(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), do2_c), probs = 0.1, na.rm = T)), eo2, time_h), color = "blue", size = 1) +
    geom_vline(xintercept = min(filter(mutate(x, time_h = time_min / 60), resp == "post_chase")$time_h))
  
  mo_plot <- ggplot(mutate(x, time_h = time_min / 60), aes(x = time_h, y = do2_c)) +
    geom_point(size = 0.25) +
    geom_point(data = filter(mutate(x, time_h = time_min / 60), do2_c == max(pull(x, do2_c), na.rm = T)), aes(x = time_h, y = do2_c), color = "red", size = 3) +
    geom_point(data = select(filter(mutate(x, time_h = time_min / 60), do2_c == max(pull(filter(x, resp == "post_chase"), do2_c), na.rm = T)), time_h, do2_c), color = "yellow", size = 2) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), time_post == 0 | 
                                      time_post == Closest(mfits(x)$time_post, 10) | 
                                      time_post == Closest(mfits(x)$time_post, 20) | 
                                      time_post == Closest(mfits(x)$time_post, 30) | 
                                      time_post == Closest(mfits(x)$time_post, 60) | 
                                      time_post == Closest(mfits(x)$time_post, 120) |
                                      time_post == Closest(mfits(x)$time_post, 300) |
                                      time_post == Closest(mfits(x)$time_post, 600) |
                                      time_post == Closest(mfits(x)$time_post, 900) ), do2_c, time_h), color = "magenta1", size = 1) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), do2_c, time_h), color = "skyblue1", size = 0.25) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200 & do2_c < quantile(pull(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), do2_c), probs = 0.1, na.rm = T)), do2_c, time_h), color = "blue", size = 1) +
    geom_vline(xintercept = min(filter(mutate(x, time_h = time_min / 60), resp == "post_chase")$time_h))
  
  hr_plot <- ggplot(mutate(x, time_h = time_min / 60), aes(x = time_h, y = hr)) +
    geom_point(size = 0.25) +
    geom_point(data = filter(mutate(x, time_h = time_min / 60), hr == max(pull(x, hr), na.rm = T)), aes(x = time_h, y = hr), color = "red", size = 3) +
    geom_point(data = select(filter(mutate(x, time_h = time_min / 60), do2_c == max(pull(filter(x, resp == "post_chase"), do2_c), na.rm = T)), time_h, hr), color = "yellow", size = 2) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), time_post == 0 | 
                                      time_post == Closest(mfits(x)$time_post, 10) | 
                                      time_post == Closest(mfits(x)$time_post, 20) | 
                                      time_post == Closest(mfits(x)$time_post, 30) | 
                                      time_post == Closest(mfits(x)$time_post, 60) | 
                                      time_post == Closest(mfits(x)$time_post, 120) |
                                      time_post == Closest(mfits(x)$time_post, 300) |
                                      time_post == Closest(mfits(x)$time_post, 600) |
                                      time_post == Closest(mfits(x)$time_post, 900) ), hr, time_h), color = "magenta1", size = 1) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), hr, time_h), color = "skyblue1", size = 0.25) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200 & do2_c < quantile(pull(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), do2_c), probs = 0.1, na.rm = T)), hr, time_h), color = "blue", size = 1) +
    geom_vline(xintercept = min(filter(mutate(x, time_h = time_min / 60), resp == "post_chase")$time_h))
  
  sv_plot <- ggplot(mutate(x, time_h = time_min / 60), aes(x = time_h, y = sv)) +
    geom_point(size = 0.25) +
    geom_point(data = filter(mutate(x, time_h = time_min / 60), sv == max(pull(x, sv), na.rm = T)), aes(x = time_h, y = sv), color = "red", size = 3) +
    geom_point(data = select(filter(mutate(x, time_h = time_min / 60), do2_c == max(pull(filter(x, resp == "post_chase"), do2_c), na.rm = T)), time_h, sv), color = "yellow", size = 2) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), time_post == 0 | 
                                      time_post == Closest(mfits(x)$time_post, 10) | 
                                      time_post == Closest(mfits(x)$time_post, 20) | 
                                      time_post == Closest(mfits(x)$time_post, 30) | 
                                      time_post == Closest(mfits(x)$time_post, 60) | 
                                      time_post == Closest(mfits(x)$time_post, 120) |
                                      time_post == Closest(mfits(x)$time_post, 300) |
                                      time_post == Closest(mfits(x)$time_post, 600) |
                                      time_post == Closest(mfits(x)$time_post, 900) ), sv, time_h), color = "magenta1", size = 1) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), sv, time_h), color = "skyblue1", size = 0.25) +
    geom_point(data = select(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200 & do2_c < quantile(pull(filter(mutate(mfits(x), time_h = time_min / 60), resp == "pre_chase" & time_post > -1200), do2_c), probs = 0.1, na.rm = T)), sv, time_h), color = "blue", size = 1) +
    geom_vline(xintercept = min(filter(mutate(x, time_h = time_min / 60), resp == "post_chase")$time_h))
  
  variable <- c("fishID", "treat", "max of all", "max at max mo2 post chase", "timeline", "smr", "20h pre chase")
  value_color <- c(deparse(substitute(x)), treat, "red", "yellow", "purple", "darkblue", "skyblue")
  checkit <- tableGrob(data.frame(variable, value_color))
  checkit$grobs[find_cell(checkit, 4, 3, "core-bg")][[1]][["gp"]] <- gpar(fill = "red")
  checkit$grobs[find_cell(checkit, 5, 3, "core-bg")][[1]][["gp"]] <- gpar(fill = "yellow")
  checkit$grobs[find_cell(checkit, 6, 3, "core-bg")][[1]][["gp"]] <- gpar(fill = "magenta1")
  checkit$grobs[find_cell(checkit, 7, 3, "core-bg")][[1]][["gp"]] <- gpar(fill = "blue")
  checkit$grobs[find_cell(checkit, 8, 3, "core-bg")][[1]][["gp"]] <- gpar(fill = "skyblue")
  
  plots <- grid.arrange(co_plot, eo2_plot, mo_plot, checkit, hr_plot, sv_plot, layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6)))
  
  return(plots)
} #to explore timeline data (MO2, CO, HR, SV) from start of experiment to finish

rsq <- function(x, y){
  cor(x, y) ^ 2
} # this function calculates r2 from two vectors, and is used in the function below: get_r2().

fick_it <- function(x){
  #Add a column of salinity from the value in fish_dat
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  trans_eq <- column_to_rownames(trans_eq, "transonic")
  sal <- fish_dat[deparse(substitute(x)), "salinity"]
  bodymass <- fish_dat[deparse(substitute(x)), "bodymass"]
  respvol_net <- fish_dat[deparse(substitute(x)), "respvol_net"]
  trans_slope <- trans_eq[fish_dat[deparse(substitute(x)), "transonic"], "slope"]
  trans_intercept <- trans_eq[fish_dat[deparse(substitute(x)), "transonic"], "intercept"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  trans_eq <- rownames_to_column(trans_eq, "transonic")
  x$sal <- sal
  x$bodymass <- bodymass
  x$respvol_net <- respvol_net
  x$trans_slope <- trans_slope
  x$trans_intercept <- trans_intercept
  
  #Calculate time since start in minutes
  x <- mutate(x, time_num = as.numeric(time) / 60 ) #minutes into the day
  x <- mutate(x, time_num_s = time_num * 60) #seconds into the day (for atmospheric pressure)
  x$date_num <- as.numeric(x$date)
  x <- mutate(x, day = date_num - min(pull(x, date_num), na.rm = T))
  x <- mutate(x, day_min = day * 1440) #one day = 1440 minutes
  x <- mutate(x, time_cont = day_min + time_num)
  x <- mutate(x, time_min = time_cont - min(pull(x, time_cont), na.rm = T))
  x <- mutate(x, datetime = ymd_hms(paste(x$date, x$time))) 
  x <- mutate(x, interval = interval(start = datetime, end = as_datetime(as.numeric(datetime) + dur)))
  
  #Calculate mean atmospheric pressure in kpa from smhi data
  x$atm <- mean(pull(filter(smhi_atm, date_num == min(x$date_num, na.rm = T) & time_uct_num > min(filter(x, date_num == min(x$date_num, na.rm = T))$time_num_s, na.rm = T) | date_num == min(filter(x, resp == "post_chase")$date_num, na.rm = T) | date_num == max(pull(x, date_num), na.rm = T) & time_uct_num < max(filter(x, date_num == max(x$date_num, na.rm = T))$time_num_s, na.rm = T)), atm_pr_kpa))
  
  #Calculate O2-solubility
  x <- mutate(x, po2_mmhg = (atm-exp(cons_10*(1-(cons_11/(temp+cons_25)))+(cons_12*(1-exp(cons_13*(1-(temp+cons_25)/cons_11))))-(cons_14*(1-exp(cons_15*(1-(cons_11/(temp+cons_25))))))+(cons_16*log((cons_11/(temp+cons_25)))))*atm)*cons_26*cons_mmhg/cons_kpa)
  x <- mutate(x, po2_kpa = abs(po2_mmhg) / cons_mmhg * cons_kpa)
  x <- mutate(x, sol_1 = exp((cons_17+(cons_18/(temp+cons_25))+(cons_19*log(temp+cons_25))+(cons_20*(temp+cons_25)))-((sal-cons_28)/cons_27)*((cons_21+(cons_22/(cons_25+temp))+(cons_23*log(cons_25+temp))+(cons_24*(cons_25+temp)))))/cons_29*cons_30/cons_mmhg*cons_31)
  x <- mutate(x, sol_2 = abs(sol_1) * cons_mmhg / cons_kpa)
  x <- mutate(x, o2_sol = sol_2 * cons_32) #Solubility of O2 at salinity K and temp T; unit mgO2/L/kPa
  
  #Calculate background respiration
  x_x <- c(min(pull(x, time_min), na.rm = T), max(pull(x, time_min), na.rm = T))
  x_y <- c(x %>% filter(resp == "bkg_pre") %>% pull(slope_airsat) %>% mean(), x %>% filter(resp == "bkg_post") %>% pull(slope_airsat) %>% mean())
  x_lm <- summary(lm(x_y ~ x_x)) #linear regression to find bkg respiration slope and intercept
  x <- mutate(x, bkg_slope_airsat = time_min * x_lm$coefficients[2,1] + x_lm$coefficients[1,1]) #slope air saturation (#x_lm$coefficients[1,1] = intercept; x_lm$coefficients[2,1] = slope)
  x <- mutate(x, bkg_slope_kpa = bkg_slope_airsat / 100 * po2_kpa)
  x <- mutate(x, bkg_slope_o2 = bkg_slope_kpa*60*o2_sol)
  x <- mutate(x, bkg_percent_of_mo2 = bkg_slope_airsat / slope_airsat * 100)
  
  #Calculate MO2
  x <- mutate(x, dpo2 = slope_airsat / 100 * po2_kpa) #kpa/s
  x <- mutate(x, do2_a = dpo2*60*o2_sol) #mg/l/min
  x <- mutate(x, do2_b = ((do2_a * respvol_net) - (bkg_slope_o2 * 10)) / (-1*bodymass)) #mg/min/kg
  x <- mutate(x, do2_c = do2_b * 60) #mg/kg/h
  x <- mutate(x, do2_d = do2_c / 1000) #mg/g/h
  
  #Calculate CO & SV per kg
  x <- mutate(x, co_corr = co * trans_slope + trans_intercept) #corrected CO
  x <- mutate(x, co_kg = co_corr / bodymass) #co per kg; ml/min/kg
  x <- mutate(x, sv = co_kg / hr) #sv (per kg)
  
  #Calculate Eo2
  x <- mutate(x, eo2 = do2_c / (co_kg *60)) #calculate eo2, requires calculated mo2 (do2_c) and co_kg
  
  #Remove variables that are not useful anymore
  x <- select(x, -bodymass, -respvol_net, -trans_slope, -trans_intercept)
  
  return(x)
} #for computing MO2, CO, HR, SV and EO2 in "fX"

get_r2 <- function(x){
  r2_all <- read.csv2("r2_all.csv")
  string <- paste(deparse(substitute(x)), "r2", sep = "")
  
  if(string %in% colnames(r2_all)){
    r2 <- drop_na(select(r2_all, all_of(string)))
    
  } else{
    fxrx <- read.delim2(here(paste("mo2r2/", string, ".txt", sep = "")))
    
    fxrx <- mutate(fxrx, date = ymd(date))
    fxrx <- mutate(fxrx, datetime = ymd_hms(paste(fxrx$date, seconds.to.hms(fxrx$second))))
    fxrx$seconds <- seq_along(fxrx$second)
    
    r2 <- vector("double", nrow(x))
    
    for (i in 1:nrow(x)){
      r2[[i]] <- rsq(select(filter(fxrx, datetime %within% x[i, "interval"] == TRUE), seconds),
                     select(filter(fxrx, datetime %within% x[i, "interval"] == TRUE), airsat))
    }
    
  }
  
  return(r2)
} # this function does the trick of calculating r2 provided it has access to raw data files (fxr2.txt). Unless they have already been saved in r2_all.csv

get_ficked <- function(x){
  #Add quantiles, max and scope of CO, HR, SV, MO2 and EO2 to fish_dat
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  
  fish_dat[deparse(substitute(x)), "smr_10"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.1, na.rm = T)), do2_c))
  fish_dat[deparse(substitute(x)), "smr_15"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.15, na.rm = T)), do2_c))
  fish_dat[deparse(substitute(x)), "smr_20"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.2, na.rm = T)), do2_c))
  fish_dat[deparse(substitute(x)), "mmr"] <- ifelse(nrow(drop_na(select(x, do2_c))) > 0, max(pull(filter(x, resp == "post_chase"), do2_c), na.rm = TRUE), max(pull(filter(x, resp == "post_chase"), do2_c), na.rm = FALSE))
  fish_dat[deparse(substitute(x)), "as"] <- fish_dat[deparse(substitute(x)), "mmr"] - fish_dat[deparse(substitute(x)), "smr_10"]
  
  old_smr <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.1, na.rm = T)), do2_c))
  mutated <- filter(mutate(x, time_h = time_min / 60, do2_c = do2_c - old_smr), resp == "post_chase")
  fish_dat[deparse(substitute(x)), "epoc"] <- AUC(mutated$time_h, 
                                                  mutated$do2_c, 
                                                  from = min(filter(mutated, resp == "post_chase")$time_h, na.rm = TRUE), 
                                                  to = filter(mutated, resp == "post_chase" & time_h == Closest(filter(mutated, resp == "post_chase")$time_h, 15))$time_h, 
                                                  method = "trapezoid")
  
  fish_dat[deparse(substitute(x)), "co_10"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.1, na.rm = T)), co_kg))
  fish_dat[deparse(substitute(x)), "co_15"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.15, na.rm = T)), co_kg))
  fish_dat[deparse(substitute(x)), "co_20"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.2, na.rm = T)), co_kg))
  fish_dat[deparse(substitute(x)), "co_max"] <- ifelse(nrow(drop_na(select(x, co_kg))) > 0, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = TRUE)$co_kg, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = FALSE)$co_kg)
  fish_dat[deparse(substitute(x)), "co_scope"] <- fish_dat[deparse(substitute(x)), "co_max"] - fish_dat[deparse(substitute(x)), "co_10"]
  
  fish_dat[deparse(substitute(x)), "hr_10"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.1, na.rm = T)), hr))
  fish_dat[deparse(substitute(x)), "hr_15"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.15, na.rm = T)), hr))
  fish_dat[deparse(substitute(x)), "hr_20"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.2, na.rm = T)), hr))
  fish_dat[deparse(substitute(x)), "hr_max"] <- ifelse(nrow(drop_na(select(x, co_kg))) > 0, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = TRUE)$hr, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = FALSE)$hr)
  fish_dat[deparse(substitute(x)), "hr_scope"] <- fish_dat[deparse(substitute(x)), "hr_max"] - fish_dat[deparse(substitute(x)), "hr_10"]
  
  fish_dat[deparse(substitute(x)), "sv_10"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.1, na.rm = T)), sv))
  fish_dat[deparse(substitute(x)), "sv_15"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.15, na.rm = T)), sv))
  fish_dat[deparse(substitute(x)), "sv_20"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.2, na.rm = T)), sv))
  fish_dat[deparse(substitute(x)), "sv_max"] <- ifelse(nrow(drop_na(select(x, co_kg))) > 0, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = TRUE)$sv, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = FALSE)$sv)
  fish_dat[deparse(substitute(x)), "sv_scope"] <- fish_dat[deparse(substitute(x)), "sv_max"] - fish_dat[deparse(substitute(x)), "sv_10"]
  
  fish_dat[deparse(substitute(x)), "eo2_10"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.1, na.rm = T)), eo2))
  fish_dat[deparse(substitute(x)), "eo2_15"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.15, na.rm = T)), eo2))
  fish_dat[deparse(substitute(x)), "eo2_20"] <- mean(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(x, resp == "pre_chase" & time_min > (min(filter(x, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.2, na.rm = T)), eo2))
  fish_dat[deparse(substitute(x)), "eo2_max"] <- ifelse(nrow(drop_na(select(x, co_kg))) > 0, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = TRUE)$eo2, filter(x, do2_c == max(pull(filter(x, resp == "post_chase"), do2_c)), na.rm = FALSE)$eo2)
  fish_dat[deparse(substitute(x)), "eo2_scope"] <- fish_dat[deparse(substitute(x)), "eo2_max"] - fish_dat[deparse(substitute(x)), "eo2_10"]
  
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  
  return(fish_dat)
} #to gather the calculated values from "fX" into fish_dat.

ctime <- function(fX){
  fishID <- deparse(substitute(fX))
  chase_h <- min(filter(mutate(fX, time_h = time_min / 60), resp == "post_chase")$time_h)
  max_co_h <- filter(mutate(fX, time_h = time_min / 60), co_kg == max(pull(fX, co_kg), na.rm = T))$time_h
  max_eo2_h <- filter(mutate(fX, time_h = time_min / 60), eo2 == max(pull(fX, eo2), na.rm = T))$time_h
  max_mo_h <- filter(mutate(fX, time_h = time_min / 60), do2_c == max(pull(fX, do2_c), na.rm = T))$time_h
  max_hr_h <- filter(mutate(fX, time_h = time_min / 60), hr == max(pull(fX, hr), na.rm = T))$time_h
  max_sv_h <- filter(mutate(fX, time_h = time_min / 60), sv == max(pull(fX, sv), na.rm = T))$time_h
  
  time_dat <- rbind(time_dat, cbind(fishID, chase_h, max_co_h, max_eo2_h, max_mo_h, max_hr_h, max_sv_h)) #if I want this
  
  time_dat <- mutate(time_dat, chase_h = as.numeric(chase_h), max_co_h = as.numeric(max_co_h), max_eo2_h = as.numeric(max_eo2_h),
                     max_mo_h = as.numeric(max_mo_h), max_hr_h = as.numeric(max_hr_h), max_sv_h = as.numeric(max_sv_h))
  
  return(time_dat)
} #to calculate the at what times (after starting the experimental run) max values occured

time_post <- function(fX){
  chase_h <- min(filter(mutate(fX, time_h = time_min / 60), resp == "post_chase")$time_h)
  max_h <- filter(mutate(fX, time_h = time_min / 60), time_min == max(pull(filter(fX, resp == "post_chase"), time_min), na.rm = T))$time_h
  time_after <- max_h - chase_h
  
  return(time_after)
} #returns hour of last measurement post-chase in "fX"

tlco <- function(x){
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  treat <- fish_dat[deparse(substitute(x)), "treatment"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  
  fishID <- deparse(substitute(x))
  
  x <- mfits(x)
  
  post_int <- filter(x, time_post == 0 | 
                       time_post == Closest(x$time_post, 10) | 
                       time_post == Closest(x$time_post, 20) | 
                       time_post == Closest(x$time_post, 30) | 
                       time_post == Closest(x$time_post, 60) | 
                       time_post == Closest(x$time_post, 120) |
                       time_post == Closest(x$time_post, 300) |
                       time_post == Closest(x$time_post, 600) |
                       time_post == Closest(x$time_post, 900))
  
  rco <- select(filter(x, resp == "pre_chase" & time_post > -180), co_kg, do2_c)
  sco <- select(filter(rco, do2_c < quantile(pull(rco, do2_c), probs = 0.2, na.rm = T)), co_kg, do2_c)
  rco <- mean(rco$co_kg)
  sco <- mean(sco$co_kg)
  chase <- select(mutate(post_int, chase = co_kg), chase)[1,1]
  post_10 <- select(mutate(post_int, post_10 = co_kg), post_10)[2,1]
  post_20 <- select(mutate(post_int, post_20 = co_kg), post_20)[3,1]
  post_30 <- select(mutate(post_int, post_30 = co_kg), post_30)[4,1]
  post_60 <- select(mutate(post_int, post_60 = co_kg), post_60)[5,1]
  post_120 <- select(mutate(post_int, post_120 = co_kg), post_120)[6,1]
  post_300 <- select(mutate(post_int, post_300 = co_kg), post_300)[7,1]
  post_600 <- select(mutate(post_int, post_600 = co_kg), post_600)[8,1]
  post_900 <- select(mutate(post_int, post_900 = co_kg), post_900)[9,1]
  co_all <- mutate(cbind(fishID, treat, rco, sco, chase, post_10, post_20, post_30, post_60, post_120, post_300, post_600, post_900), fishID = as.character(fishID), treat = as.character(treat))
  
  co_time <- rbind(co_time, co_all)
  
  return(co_time)
} #get timeline data for co in fish X

tleo <- function(x){
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  treat <- fish_dat[deparse(substitute(x)), "treatment"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  
  fishID <- deparse(substitute(x))
  
  x <- mfits(x)
  
  post_int <- filter(x, time_post == 0 | 
                       time_post == Closest(x$time_post, 10) | 
                       time_post == Closest(x$time_post, 20) | 
                       time_post == Closest(x$time_post, 30) | 
                       time_post == Closest(x$time_post, 60) | 
                       time_post == Closest(x$time_post, 120) |
                       time_post == Closest(x$time_post, 300) |
                       time_post == Closest(x$time_post, 600) |
                       time_post == Closest(x$time_post, 900))
  
  reo <- select(filter(x, resp == "pre_chase" & time_post > -180), eo2, do2_c)
  seo <- select(filter(reo, do2_c < quantile(pull(reo, do2_c), probs = 0.2, na.rm = T)), eo2, do2_c)
  reo <- mean(reo$eo2)
  seo <- mean(seo$eo2)
  chase <- select(mutate(post_int, chase = eo2), chase)[1,1]
  post_10 <- select(mutate(post_int, post_10 = eo2), post_10)[2,1]
  post_20 <- select(mutate(post_int, post_20 = eo2), post_20)[3,1]
  post_30 <- select(mutate(post_int, post_30 = eo2), post_30)[4,1]
  post_60 <- select(mutate(post_int, post_60 = eo2), post_60)[5,1]
  post_120 <- select(mutate(post_int, post_120 = eo2), post_120)[6,1]
  post_300 <- select(mutate(post_int, post_300 = eo2), post_300)[7,1]
  post_600 <- select(mutate(post_int, post_600 = eo2), post_600)[8,1]
  post_900 <- select(mutate(post_int, post_900 = eo2), post_900)[9,1]
  eo_all <- mutate(cbind(fishID, treat, reo, seo, chase, post_10, post_20, post_30, post_60, post_120, post_300, post_600, post_900), fishID = as.character(fishID), treat = as.character(treat))
  
  eo2_time <- rbind(eo2_time, eo_all)
  
  return(eo2_time)
} #get timeline data for eo2 in fish X

tlmo <- function(x){
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  treat <- fish_dat[deparse(substitute(x)), "treatment"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  
  fishID <- deparse(substitute(x))
  
  x <- mfits(x)
  
  post_int <- filter(x, time_post == 0 | 
                       time_post == Closest(x$time_post, 10) | 
                       time_post == Closest(x$time_post, 20) | 
                       time_post == Closest(x$time_post, 30) | 
                       time_post == Closest(x$time_post, 60) | 
                       time_post == Closest(x$time_post, 120) |
                       time_post == Closest(x$time_post, 300) |
                       time_post == Closest(x$time_post, 600) |
                       time_post == Closest(x$time_post, 900))
  
  rmr <- select(filter(x, resp == "pre_chase" & time_post > -180), do2_c)
  smr <- select(filter(rmr, do2_c < quantile(pull(rmr, do2_c), probs = 0.2, na.rm = T)), do2_c)
  rmr <- mean(rmr$do2_c)
  smr <- mean(smr$do2_c)
  chase <- select(mutate(post_int, chase = do2_c), chase)[1,1]
  post_10 <- select(mutate(post_int, post_10 = do2_c), post_10)[2,1]
  post_20 <- select(mutate(post_int, post_20 = do2_c), post_20)[3,1]
  post_30 <- select(mutate(post_int, post_30 = do2_c), post_30)[4,1]
  post_60 <- select(mutate(post_int, post_60 = do2_c), post_60)[5,1]
  post_120 <- select(mutate(post_int, post_120 = do2_c), post_120)[6,1]
  post_300 <- select(mutate(post_int, post_300 = do2_c), post_300)[7,1]
  post_600 <- select(mutate(post_int, post_600 = do2_c), post_600)[8,1]
  post_900 <- select(mutate(post_int, post_900 = do2_c), post_900)[9,1]
  mo2_all <- mutate(cbind(fishID, treat, rmr, smr, chase, post_10, post_20, post_30, post_60, post_120, post_300, post_600, post_900), fishID = as.character(fishID), treat = as.character(treat))
  
  mo2_time <- rbind(mo2_time, mo2_all)
  
  return(mo2_time)
} #get timeline data for mo2 in fish X

tlhr <- function(x){
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  treat <- fish_dat[deparse(substitute(x)), "treatment"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  
  fishID <- deparse(substitute(x))
  
  x <- mfits(x)
  
  post_int <- filter(x, time_post == 0 | 
                       time_post == Closest(x$time_post, 10) | 
                       time_post == Closest(x$time_post, 20) | 
                       time_post == Closest(x$time_post, 30) | 
                       time_post == Closest(x$time_post, 60) | 
                       time_post == Closest(x$time_post, 120) |
                       time_post == Closest(x$time_post, 300) |
                       time_post == Closest(x$time_post, 600) |
                       time_post == Closest(x$time_post, 900))
  
  rhr <- select(filter(x, resp == "pre_chase" & time_post > -180), hr, do2_c)
  shr <- select(filter(rhr, do2_c < quantile(pull(rhr, do2_c), probs = 0.2, na.rm = T)), hr, do2_c)
  rhr <- mean(rhr$hr)
  shr <- mean(shr$hr)
  chase <- select(mutate(post_int, chase = hr), chase)[1,1]
  post_10 <- select(mutate(post_int, post_10 = hr), post_10)[2,1]
  post_20 <- select(mutate(post_int, post_20 = hr), post_20)[3,1]
  post_30 <- select(mutate(post_int, post_30 = hr), post_30)[4,1]
  post_60 <- select(mutate(post_int, post_60 = hr), post_60)[5,1]
  post_120 <- select(mutate(post_int, post_120 = hr), post_120)[6,1]
  post_300 <- select(mutate(post_int, post_300 = hr), post_300)[7,1]
  post_600 <- select(mutate(post_int, post_600 = hr), post_600)[8,1]
  post_900 <- select(mutate(post_int, post_900 = hr), post_900)[9,1]
  hr_all <- mutate(cbind(fishID, treat, rhr, shr, chase, post_10, post_20, post_30, post_60, post_120, post_300, post_600, post_900), fishID = as.character(fishID), treat = as.character(treat))
  
  hr_time <- rbind(hr_time, hr_all)
  
  return(hr_time)
} #get timeline data for hr in fish X

tlsv <- function(x){
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  treat <- fish_dat[deparse(substitute(x)), "treatment"]
  fish_dat <- rownames_to_column(fish_dat, "fishID")
  
  fishID <- deparse(substitute(x))
  
  x <- mfits(x)
  
  post_int <- filter(x, time_post == 0 | 
                       time_post == Closest(x$time_post, 10) | 
                       time_post == Closest(x$time_post, 20) | 
                       time_post == Closest(x$time_post, 30) | 
                       time_post == Closest(x$time_post, 60) | 
                       time_post == Closest(x$time_post, 120) |
                       time_post == Closest(x$time_post, 300) |
                       time_post == Closest(x$time_post, 600) |
                       time_post == Closest(x$time_post, 900))
  
  rsv <- select(filter(x, resp == "pre_chase" & time_post > -180), sv, do2_c)
  ssv <- select(filter(rsv, do2_c < quantile(pull(rsv, do2_c), probs = 0.2, na.rm = T)), sv, do2_c)
  rsv <- mean(rsv$sv)
  ssv <- mean(ssv$sv)
  chase <- select(mutate(post_int, chase = sv), chase)[1,1]
  post_10 <- select(mutate(post_int, post_10 = sv), post_10)[2,1]
  post_20 <- select(mutate(post_int, post_20 = sv), post_20)[3,1]
  post_30 <- select(mutate(post_int, post_30 = sv), post_30)[4,1]
  post_60 <- select(mutate(post_int, post_60 = sv), post_60)[5,1]
  post_120 <- select(mutate(post_int, post_120 = sv), post_120)[6,1]
  post_300 <- select(mutate(post_int, post_300 = sv), post_300)[7,1]
  post_600 <- select(mutate(post_int, post_600 = sv), post_600)[8,1]
  post_900 <- select(mutate(post_int, post_900 = sv), post_900)[9,1]
  sv_all <- mutate(cbind(fishID, treat, rsv, ssv, chase, post_10, post_20, post_30, post_60, post_120, post_300, post_600, post_900), fishID = as.character(fishID), treat = as.character(treat))
  
  sv_time <- rbind(sv_time, sv_all)
  
  return(sv_time)
} #get timeline data for sv in fish X

diff_time <- find_diff <- function(x){
  fishID <- deparse(substitute(x))
  
  time_diff <- select(mutate(select(filter(mfits(x), 
                                           time_post == Closest(mfits(x)$time_post, 10) | 
                                             time_post == Closest(mfits(x)$time_post, 20) | 
                                             time_post == Closest(mfits(x)$time_post, 30) | 
                                             time_post == Closest(mfits(x)$time_post, 60) | 
                                             time_post == Closest(mfits(x)$time_post, 120) |
                                             time_post == Closest(mfits(x)$time_post, 300) |
                                             time_post == Closest(mfits(x)$time_post, 600) |
                                             time_post == Closest(mfits(x)$time_post, 900)),
                                    fishID, time_post),
                             timepoint = c(10, 20, 30, 60, 120, 300, 600, 900), timediff = time_post - timepoint), timediff)
  
  diff_10 <- select(mutate(time_diff, diff_10 = timediff), diff_10)[1, 1]
  diff_20 <- select(mutate(time_diff, diff_20 = timediff), diff_20)[2, 1]
  diff_30 <- select(mutate(time_diff, diff_30 = timediff), diff_30)[3, 1]
  diff_60 <- select(mutate(time_diff, diff_60 = timediff), diff_60)[4, 1]
  diff_120 <- select(mutate(time_diff, diff_120 = timediff), diff_120)[5, 1]
  diff_300 <- select(mutate(time_diff, diff_300 = timediff), diff_300)[6, 1]
  diff_600 <- select(mutate(time_diff, diff_600 = timediff), diff_600)[7, 1]
  diff_900 <- select(mutate(time_diff, diff_900 = timediff), diff_900)[8, 1]
  diff_all <- mutate(cbind(fishID, diff_10, diff_20, diff_30, diff_60, diff_120, diff_300, diff_600, diff_900), fishID = as.character(fishID))
  
  diff_time <- rbind(diff_time, diff_all)
  
  return(diff_time)
} #get the time deviance of the measurements from the timepoint we are looking for

rbind.na <- function (..., deparse.level = 1){
  na <- nargs() - (!missing(deparse.level))
  deparse.level <- as.integer(deparse.level)
  stopifnot(0 <= deparse.level, deparse.level <= 2)
  argl <- list(...)
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }    
  if (na == 0) 
    return(NULL)
  if (na == 1) {
    if (isS4(..1)) 
      return(rbind2(..1))
    else return(matrix(..., nrow = 1)) ##.Internal(rbind(deparse.level, ...)))
  }
  if (deparse.level) {
    symarg <- as.list(sys.call()[-1L])[1L:na]
    Nms <- function(i) {
      if (is.null(r <- names(symarg[i])) || r == "") {
        if (is.symbol(r <- symarg[[i]]) || deparse.level == 
            2) 
          deparse(r)
      }
      else r
    }
  }
  
  ## deactivated, otherwise no fill in with two arguments
  if (na == 0) {
    r <- argl[[2]]
    fix.na <- FALSE
  }
  else {
    nrs <- unname(lapply(argl, ncol))
    iV <- sapply(nrs, is.null)
    fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
    ## deactivated, otherwise data will be recycled
    #if (fix.na) {
    #    nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
    #    argl[[na]] <- rbind(rep(argl[[na]], length.out = nr), 
    #        deparse.level = 0)
    #}
    if (deparse.level) {
      if (fix.na) 
        fix.na <- !is.null(Nna <- Nms(na))
      if (!is.null(nmi <- names(argl))) 
        iV <- iV & (nmi == "")
      ii <- if (fix.na) 
        2:(na - 1)
      else 2:na
      if (any(iV[ii])) {
        for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i))) 
          names(argl)[i] <- nmi
      }
    }
    
    ## filling with NA's to maximum occuring ncols
    nCol <- as.numeric(sapply(argl, function(x) if (is.null(ncol(x))) length(x)
                              else ncol(x)))
    maxCol <- max(nCol, na.rm = TRUE)  
    argl <- lapply(argl, function(x)  if (is.null(ncol(x))) c(x, rep(NA, maxCol - length(x)))
                   else cbind(x, matrix(, nrow(x), maxCol - ncol(x))))  
    
    ## create a common name vector from the
    ## column names of all 'argl' items
    namesVEC <- rep(NA, maxCol)  
    for (i in 1:length(argl)) {
      CN <- colnames(argl[[i]])          
      m <- !(CN %in% namesVEC)
      namesVEC[m] <- CN[m]          
    }  
    
    ## make all column names from common 'namesVEC'
    for (j in 1:length(argl)) {    
      if (!is.null(ncol(argl[[j]]))) colnames(argl[[j]]) <- namesVEC
    }
    
    r <- do.call(rbind, c(argl[-1L], list(deparse.level = deparse.level)))        
  }
  
  d2 <- dim(r)
  
  ## make all column names from common 'namesVEC'
  colnames(r) <- colnames(argl[[1]])
  
  r <- rbind2(argl[[1]], r)
  
  if (deparse.level == 0) 
    return(r)
  ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
  ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
  if (ism1 && ism2) 
    return(r)
  Nrow <- function(x) {
    d <- dim(x)
    if (length(d) == 2L) 
      d[1L]
    else as.integer(length(x) > 0L)
  }
  nn1 <- !is.null(N1 <- if ((l1 <- Nrow(..1)) && !ism1) Nms(1))
  nn2 <- !is.null(N2 <- if (na == 2 && Nrow(..2) && !ism2) Nms(2))
  if (nn1 || nn2 || fix.na) {
    if (is.null(rownames(r))) 
      rownames(r) <- rep.int("", nrow(r))
    setN <- function(i, nams) rownames(r)[i] <<- if (is.null(nams)) 
      ""
    else nams
    if (nn1) 
      setN(1, N1)
    if (nn2) 
      setN(1 + l1, N2)
    if (fix.na) 
      setN(nrow(r), Nna)
  }
  r
}

cbind.na <- function (..., deparse.level = 1){
  na <- nargs() - (!missing(deparse.level))    
  deparse.level <- as.integer(deparse.level)
  stopifnot(0 <= deparse.level, deparse.level <= 2)
  argl <- list(...)   
  while (na > 0 && is.null(argl[[na]])) {
    argl <- argl[-na]
    na <- na - 1
  }
  if (na == 0) 
    return(NULL)
  if (na == 1) {         
    if (isS4(..1)) 
      return(cbind2(..1))
    else return(matrix(...))  ##.Internal(cbind(deparse.level, ...)))
  }
  if (deparse.level) {       
    symarg <- as.list(sys.call()[-1L])[1L:na]
    Nms <- function(i) {
      if (is.null(r <- names(symarg[i])) || r == "") {
        if (is.symbol(r <- symarg[[i]]) || deparse.level == 
            2) 
          deparse(r)
      }
      else r
    }
  }   
  ## deactivated, otherwise no fill in with two arguments
  if (na == 0) {
    r <- argl[[2]]
    fix.na <- FALSE
  }
  else {
    nrs <- unname(lapply(argl, nrow))
    iV <- sapply(nrs, is.null)
    fix.na <- identical(nrs[(na - 1):na], list(NULL, NULL))
    ## deactivated, otherwise data will be recycled
    #if (fix.na) {
    #    nr <- max(if (all(iV)) sapply(argl, length) else unlist(nrs[!iV]))
    #    argl[[na]] <- cbind(rep(argl[[na]], length.out = nr), 
    #        deparse.level = 0)
    #}       
    if (deparse.level) {
      if (fix.na) 
        fix.na <- !is.null(Nna <- Nms(na))
      if (!is.null(nmi <- names(argl))) 
        iV <- iV & (nmi == "")
      ii <- if (fix.na) 
        2:(na - 1)
      else 2:na
      if (any(iV[ii])) {
        for (i in ii[iV[ii]]) if (!is.null(nmi <- Nms(i))) 
          names(argl)[i] <- nmi
      }
    }
    
    ## filling with NA's to maximum occuring nrows
    nRow <- as.numeric(sapply(argl, function(x) NROW(x)))
    maxRow <- max(nRow, na.rm = TRUE)  
    argl <- lapply(argl, function(x)  if (is.null(nrow(x))) c(x, rep(NA, maxRow - length(x)))
                   else rbind.na(x, matrix(, maxRow - nrow(x), ncol(x))))
    r <- do.call(cbind, c(argl[-1L], list(deparse.level = deparse.level)))
  }
  d2 <- dim(r)
  r <- cbind2(argl[[1]], r)
  if (deparse.level == 0) 
    return(r)
  ism1 <- !is.null(d1 <- dim(..1)) && length(d1) == 2L
  ism2 <- !is.null(d2) && length(d2) == 2L && !fix.na
  if (ism1 && ism2) 
    return(r)
  Ncol <- function(x) {
    d <- dim(x)
    if (length(d) == 2L) 
      d[2L]
    else as.integer(length(x) > 0L)
  }
  nn1 <- !is.null(N1 <- if ((l1 <- Ncol(..1)) && !ism1) Nms(1))
  nn2 <- !is.null(N2 <- if (na == 2 && Ncol(..2) && !ism2) Nms(2))
  if (nn1 || nn2 || fix.na) {
    if (is.null(colnames(r))) 
      colnames(r) <- rep.int("", ncol(r))
    setN <- function(i, nams) colnames(r)[i] <<- if (is.null(nams)) 
      ""
    else nams
    if (nn1) 
      setN(1, N1)
    if (nn2) 
      setN(1 + l1, N2)
    if (fix.na) 
      setN(ncol(r), Nna)
  }
  r
}

discuss_epoc <- function(x){
  current_fish <- x
  fish_name <- deparse(substitute(x))
  
  fish_dat <- column_to_rownames(fish_dat, "fishID")
  epoc_data <- fish_dat[deparse(substitute(x)), "epoc"]
  
  old_smr <- 1.1*mean(pull(filter(current_fish, resp == "pre_chase" & time_min > (min(filter(current_fish, resp == "post_chase")$time_min) - 1200) & do2_c < quantile(pull(filter(current_fish, resp == "pre_chase" & time_min > (min(filter(current_fish, resp == "post_chase")$time_min) - 1200) ), do2_c), probs = 0.1, na.rm = T)), do2_c))
  
  mutated <- filter(mutate(current_fish, time_h = time_min / 60, do2_c = do2_c - old_smr), resp == "post_chase")
  mutated <- mutate(mutated, time_h = time_h - min(mutated$time_h))
  
  #model it with a generalized additive model
  thegam <- mgcv::gam(do2_c ~ s(time_h, bs = "cs"), data = mutated)
  y_axis_gam <- predict(thegam)
  
  #model it with a fifth order polynomial linear regression
  thelm <- lm(do2_c ~ poly(time_h, 5), data = mutated)
  y_axis_lm <- predict(thelm)
  
  #model it as exponential decay
  f <- function(x,a,b) {a * exp(b * x)}
  fm_lm <- lm(do2_c ~ time_h, mutated)
  st <- list(a = exp(coef(fm_lm)[1]), b = coef(fm_lm)[2])
  thedec <- nls(do2_c ~ f(time_h, a, b), mutated, start = st)
  y_axis_dec <- predict(thedec)
  
  #Calculate EPOC from the smoothed line
  epoc_gam <- AUC(mutated$time_h, 
                  y_axis_gam, 
                  from = min(mutated$time_h, na.rm = TRUE),
                  to = filter(mutated, time_h == Closest(mutated$time_h, 15))$time_h,
                  method = "trapezoid")
  
  epoc_lm <- AUC(mutated$time_h,
                 y_axis_lm,
                 from = min(mutated$time_h, na.rm = TRUE),
                 to = filter(mutated, time_h == Closest(mutated$time_h, 15))$time_h,
                 method = "trapezoid")
  
  epoc_dec <- AUC(mutated$time_h, y_axis_dec,
                  from = min(mutated$time_h, na.rm = TRUE),
                  to = filter(mutated, time_h == Closest(mutated$time_h, 15))$time_h,
                  method = "trapezoid")
  
  #Plot the data and overlay model fits and a legend that shows the calculated epoc
  plot_og <- plot(mutated$time_h, mutated$do2_c, ylim = c(-50, 260), ylab = "mo2 - (1.1 * smr)", xlab = "hours post chase", 
                  main = deparse(substitute(x)))
  
  output <- as_tibble(data.frame(epoc_data, epoc_gam))
  
  plot_og
  lines(mutated$time_h, mutated$do2_c, col = "black", lwd = 1, lty = 2)
  lines(mutated$time_h, predict(thegam), col = "red", lwd = 1)
  lines(mutated$time_h, predict(thelm), col = "blue", lwd = 1)
  lines(mutated$time_h, predict(thedec), col = "purple", lwd = 1)
  abline(h = 0, lwd = 2, lty = 2)
  legend(8, 260, legend = c(paste("data epoc: ",  as.character(round(epoc_data, digits = 2)), ""), paste("gam epoc: ", as.character(round(epoc_gam, digits = 2)), ""), paste("5poly epoc: ", as.character(round(epoc_lm, digits = 2)), ""), paste("exp epoc: ", as.character(round(epoc_dec, digits = 2)), "")),
         col = c("black", "red", "blue", "purple"), lty=c(2,1,1,1), cex=1)
  
} #Function to model data from post-chase to compare model fits and epoc-estimates
# ####

#Calculate and add MO2, EPOC, CO, HR, SV and EO2 for each fish into fish_dat and explore them####
f2 <- read_csv2(here("mo2r2/f2.csv"))
f2 <- fick_it(f2)
#xplr(f2)
#xplm(f2)
f2$r2 <- get_r2(f2)
fish_dat <- get_ficked(f2)

f3 <- read_csv2(here("mo2r2/f3.csv"))
f3 <- fick_it(f3)
#xplr(f3)
#xplm(f3)
f3$r2 <- get_r2(f3)
fish_dat <- get_ficked(f3)

f4 <- read_csv2(here("mo2r2/f4.csv"))
f4 <- fick_it(f4)
#xplr(f4)
#xplm(f4)
f4$r2 <- get_r2(f4)
fish_dat <- get_ficked(f4)

f5 <- read_csv2(here("mo2r2/f5.csv"))
f5 <- fick_it(f5)
#xplr(f5)
#xplm(f5)
f5$r2 <- get_r2(f5)
fish_dat <- get_ficked(f5)

f6 <- read_csv2(here("mo2r2/f6.csv"))
f6 <- fick_it(f6)
#xplr(f6)
#xplm(f6)
f6$r2 <- get_r2(f6)
fish_dat <- get_ficked(f6)

f7 <- read_csv2(here("mo2r2/f7.csv"))
f7 <- fick_it(f7)
#xplr(f7)
#xplm(f7)
f7$r2 <- get_r2(f7)
fish_dat <- get_ficked(f7)

f8 <- read_csv2(here("mo2r2/f8.csv"))
f8 <- fick_it(f8)
#xplr(f8)
#xplm(f8)
f8$r2 <- get_r2(f8)
fish_dat <- get_ficked(f8)

f9 <- read_csv2(here("mo2r2/f9.csv"))
f9 <- fick_it(f9)
#xplr(f9)
#xplm(f9)
f9$r2 <- get_r2(f9)
fish_dat <- get_ficked(f9)

f10 <- read_csv2(here("mo2r2/f10.csv"))
f10 <- fick_it(f10)
#xplr(f10)
#xplm(f10)
f10$r2 <- get_r2(f10)
fish_dat <- get_ficked(f10)

f11 <- read_csv2(here("mo2r2/f11.csv"))
f11 <- fick_it(f11)
#xplr(f11)
#xplm(f11)
f11$r2 <- get_r2(f11)
fish_dat <- get_ficked(f11)

f12 <- read_csv2(here("mo2r2/f12.csv"))
f12 <- fick_it(f12)
#xplr(f12)
#xplm(f12)
f12$r2 <- get_r2(f12)
fish_dat <- get_ficked(f12)

f13 <- read_csv2(here("mo2r2/f13.csv"))
f13 <- fick_it(f13)
#xplr(f13)
#xplm(f13)
f13$r2 <- get_r2(f13)
fish_dat <- get_ficked(f13)

f14 <- read_csv2(here("mo2r2/f14.csv"))
f14 <- fick_it(f14)
#xplr(f14)
#xplm(f14)
f14$r2 <- get_r2(f14)
fish_dat <- get_ficked(f14)

f15 <- read_csv2(here("mo2r2/f15.csv"))
f15 <- fick_it(f15)
#xplr(f15)
#xplm(f15)
f15$r2 <- get_r2(f15)
fish_dat <- get_ficked(f15)

f16 <- read_csv2(here("mo2r2/f16.csv"))
f16 <- fick_it(f16)
#xplr(f16)
#xplm(f16)
f16$r2 <- get_r2(f16)
fish_dat <- get_ficked(f16)

f17 <- read_csv2(here("mo2r2/f17.csv"))
f17 <- fick_it(f17)
#xplr(f17)
#xplm(f17)
f17$r2 <- get_r2(f17)
fish_dat <- get_ficked(f17)

f18 <- read_csv2(here("mo2r2/f18.csv"))
f18 <- fick_it(f18)
#xplr(f18)
#xplm(f18)
f18$r2 <- get_r2(f18)
fish_dat <- get_ficked(f18)

f19 <- read_csv2(here("mo2r2/f19.csv"))
f19 <- fick_it(f19)
#xplr(f19)
#xplm(f19)
f19$r2 <- get_r2(f19)
fish_dat <- get_ficked(f19)

f20 <- read_csv2(here("mo2r2/f20.csv"))
f20 <- fick_it(f20)
#xplr(f20)
#xplm(f20)
f20$r2 <- get_r2(f20)
fish_dat <- get_ficked(f20)

f21 <- read_csv2(here("mo2r2/f21.csv"))
f21 <- fick_it(f21)
#xplr(f21)
#xplm(f21)
f21$r2 <- get_r2(f21)
fish_dat <- get_ficked(f21)

f24 <- read_csv2(here("mo2r2/f24.csv"))
f24 <- fick_it(f24)
#xplr(f24)
#xplm(f24)
f24$r2 <- get_r2(f24)
fish_dat <- get_ficked(f24)

f25 <- read_csv2(here("mo2r2/f25.csv"))
f25 <- fick_it(f25)
#xplr(f25)
#xplm(f25)
f25$r2 <- get_r2(f25)
fish_dat <- get_ficked(f25)

f26 <- read_csv2(here("mo2r2/f26.csv"))
f26 <- fick_it(f26)
#xplr(f26)
#xplm(f26)
f26$r2 <- get_r2(f26)
fish_dat <- get_ficked(f26)

f27 <- read_csv2(here("mo2r2/f27.csv"))
f27 <- fick_it(f27)
#xplr(f27)
#xplm(f27)
f27$r2 <- get_r2(f27)
fish_dat <- get_ficked(f27)

f28 <- read_csv2(here("mo2r2/f28.csv"))
f28 <- fick_it(f28)
#xplr(f28)
#xplm(f28)
f28$r2 <- get_r2(f28)
fish_dat <- get_ficked(f28)

f29 <- read_csv2(here("mo2r2/f29.csv"))
f29 <- fick_it(f29)
#xplr(f29)
#xplm(f29)
f29$r2 <- get_r2(f29)
fish_dat <- get_ficked(f29)

f30 <- read_csv2(here("mo2r2/f30.csv"))
f30 <- fick_it(f30)
#xplr(f30)
#xplm(f30)
f30$r2 <- get_r2(f30)
fish_dat <- get_ficked(f30)

f31 <- read_csv2(here("mo2r2/f31.csv"))
f31 <- fick_it(f31)
#xplr(f31)
#xplm(f31)
f31$r2 <- get_r2(f31)
fish_dat <- get_ficked(f31)

f33 <- read_csv2(here("mo2r2/f33.csv"))
f33 <- fick_it(f33)
#xplr(f33)
#xplm(f33)
f33$r2 <- get_r2(f33)
fish_dat <- get_ficked(f33)

f34 <- read_csv2(here("mo2r2/f34.csv"))
f34 <- fick_it(f34)
#xplr(f34)
#xplm(f34)
f34$r2 <- get_r2(f34)
fish_dat <- get_ficked(f34)

f35 <- read_csv2(here("mo2r2/f35.csv"))
f35 <- fick_it(f35)
#xplr(f35)
#xplm(f35)
f35$r2 <- get_r2(f35)
fish_dat <- get_ficked(f35)

f36 <- read_csv2(here("mo2r2/f36.csv"))
f36 <- fick_it(f36)
#xplr(f36)
#xplm(f36)
f36$r2 <- get_r2(f36)
fish_dat <- get_ficked(f36)

f38 <- read_csv2(here("mo2r2/f38.csv"))
f38 <- fick_it(f38)
#xplr(f38)
#xplm(f38)
f38$r2 <- get_r2(f38)
fish_dat <- get_ficked(f38)

# ####

#Compute timeline data####
time_dat <- data.frame(fishID = character(), chase_h = double(), max_co_h = double(), max_eo2_h = double(), max_mo_h = double(), max_hr_h = double(), max_sv_h = double(), stringsAsFactors = FALSE) #data frame for collecting values from time_dat <- ctime()
co_time <- data.frame(fishID = character(), treat = character(), rco = double(), sco = double(), chase = double(), post_10 = double(), post_20 = double(), post_30 = double(), post_60 = double(), post_120 = double(), stringsAsFactors = FALSE) #data frame for collecting values from co_time <- tlco()
eo2_time <- data.frame(fishID = character(), treat = character(), rco = double(), sco = double(), chase = double(), post_10 = double(), post_20 = double(), post_30 = double(), post_60 = double(), post_120 = double(), stringsAsFactors = FALSE) #data frame for collecting values from co_time <- tleo()
mo2_time <- data.frame(fishID = character(), treat = character(), rco = double(), sco = double(), chase = double(), post_10 = double(), post_20 = double(), post_30 = double(), post_60 = double(), post_120 = double(), stringsAsFactors = FALSE) #data frame for collecting values from co_time <- tlmo()
hr_time <- data.frame(fishID = character(), treat = character(), rhr = double(), shr = double(), chase = double(), post_10 = double(), post_20 = double(), post_30 = double(), post_60 = double(), post_120 = double(), stringsAsFactors = FALSE) #data frame for collecting values from co_time <- tlhr()
sv_time <- data.frame(fishID = character(), treat = character(), rco = double(), sco = double(), chase = double(), post_10 = double(), post_20 = double(), post_30 = double(), post_60 = double(), post_120 = double(), stringsAsFactors = FALSE) #data frame for collecting values from co_time <- tlsv()
diff_time <- data.frame(fishID = character(), diff_10 = double(), diff_20 = double(), diff_30 = double(), diff_60 = double(), diff_120 = double(), stringsAsFactors = FALSE) #data frame for collecting values from co_time <- tlsv()

time_dat <- ctime(f2)
co_time <- tlco(f2)
eo2_time <- tleo(f2)
mo2_time <- tlmo(f2)
hr_time <- tlhr(f2)
sv_time <- tlsv(f2)
diff_time <- find_diff(f2)

time_dat <- ctime(f3)
co_time <- tlco(f3)
eo2_time <- tleo(f3)
mo2_time <- tlmo(f3)
hr_time <- tlhr(f3)
sv_time <- tlsv(f3)
diff_time <- find_diff(f3)

time_dat <- ctime(f4)
co_time <- tlco(f4)
eo2_time <- tleo(f4)
mo2_time <- tlmo(f4)
hr_time <- tlhr(f4)
sv_time <- tlsv(f4)
diff_time <- find_diff(f4)

time_dat <- ctime(f5)
co_time <- tlco(f5)
eo2_time <- tleo(f5)
mo2_time <- tlmo(f5)
hr_time <- tlhr(f5)
sv_time <- tlsv(f5)
diff_time <- find_diff(f5)

time_dat <- ctime(f6)
co_time <- tlco(f6)
eo2_time <- tleo(f6)
mo2_time <- tlmo(f6)
hr_time <- tlhr(f6)
sv_time <- tlsv(f6)
diff_time <- find_diff(f6)

time_dat <- ctime(f7)
co_time <- tlco(f7)
eo2_time <- tleo(f7)
mo2_time <- tlmo(f7)
hr_time <- tlhr(f7)
sv_time <- tlsv(f7)
diff_time <- find_diff(f7)

time_dat <- ctime(f8)
co_time <- tlco(f8)
eo2_time <- tleo(f8)
mo2_time <- tlmo(f8)
hr_time <- tlhr(f8)
sv_time <- tlsv(f8)
diff_time <- find_diff(f8)

time_dat <- ctime(f9)
co_time <- tlco(f9)
eo2_time <- tleo(f9)
mo2_time <- tlmo(f9)
hr_time <- tlhr(f9)
sv_time <- tlsv(f9)
diff_time <- find_diff(f9)

time_dat <- ctime(f10)
co_time <- tlco(f10)
eo2_time <- tleo(f10)
mo2_time <- tlmo(f10)
hr_time <- tlhr(f10)
sv_time <- tlsv(f10)
diff_time <- find_diff(f10)

time_dat <- ctime(f11)
co_time <- tlco(f11)
eo2_time <- tleo(f11)
mo2_time <- tlmo(f11)
hr_time <- tlhr(f11)
sv_time <- tlsv(f11)
diff_time <- find_diff(f11)

time_dat <- ctime(f12)
co_time <- tlco(f12)
eo2_time <- tleo(f12)
mo2_time <- tlmo(f12)
hr_time <- tlhr(f12)
sv_time <- tlsv(f12)
diff_time <- find_diff(f12)

time_dat <- ctime(f13)
co_time <- tlco(f13)
eo2_time <- tleo(f13)
mo2_time <- tlmo(f13)
hr_time <- tlhr(f13)
sv_time <- tlsv(f13)
diff_time <- find_diff(f13)

time_dat <- ctime(f14)
co_time <- tlco(f14)
eo2_time <- tleo(f14)
mo2_time <- tlmo(f14)
hr_time <- tlhr(f14)
sv_time <- tlsv(f14)
diff_time <- find_diff(f14)

time_dat <- ctime(f15)
co_time <- tlco(f15)
eo2_time <- tleo(f15)
mo2_time <- tlmo(f15)
hr_time <- tlhr(f15)
sv_time <- tlsv(f15)
diff_time <- find_diff(f15)

time_dat <- ctime(f16)
co_time <- tlco(f16)
eo2_time <- tleo(f16)
mo2_time <- tlmo(f16)
hr_time <- tlhr(f16)
sv_time <- tlsv(f16)
diff_time <- find_diff(f16)

time_dat <- ctime(f17)
co_time <- tlco(f17)
eo2_time <- tleo(f17)
mo2_time <- tlmo(f17)
hr_time <- tlhr(f17)
sv_time <- tlsv(f17)
diff_time <- find_diff(f17)

time_dat <- ctime(f18)
co_time <- tlco(f18)
eo2_time <- tleo(f18)
mo2_time <- tlmo(f18)
hr_time <- tlhr(f18)
sv_time <- tlsv(f18)
diff_time <- find_diff(f18)

time_dat <- ctime(f19)
co_time <- tlco(f19)
eo2_time <- tleo(f19)
mo2_time <- tlmo(f19)
hr_time <- tlhr(f19)
sv_time <- tlsv(f19)
diff_time <- find_diff(f19)

time_dat <- ctime(f20)
co_time <- tlco(f20)
eo2_time <- tleo(f20)
mo2_time <- tlmo(f20)
hr_time <- tlhr(f20)
sv_time <- tlsv(f20)
diff_time <- find_diff(f20)

time_dat <- ctime(f21)
co_time <- tlco(f21)
eo2_time <- tleo(f21)
mo2_time <- tlmo(f21)
hr_time <- tlhr(f21)
sv_time <- tlsv(f21)
diff_time <- find_diff(f21)

time_dat <- ctime(f24)
co_time <- tlco(f24)
eo2_time <- tleo(f24)
mo2_time <- tlmo(f24)
hr_time <- tlhr(f24)
sv_time <- tlsv(f24)
diff_time <- find_diff(f24)

time_dat <- ctime(f25)
co_time <- tlco(f25)
eo2_time <- tleo(f25)
mo2_time <- tlmo(f25)
hr_time <- tlhr(f25)
sv_time <- tlsv(f25)
diff_time <- find_diff(f25)

time_dat <- ctime(f26)
co_time <- tlco(f26)
eo2_time <- tleo(f26)
mo2_time <- tlmo(f26)
hr_time <- tlhr(f26)
sv_time <- tlsv(f26)
diff_time <- find_diff(f26)

time_dat <- ctime(f27)
co_time <- tlco(f27)
eo2_time <- tleo(f27)
mo2_time <- tlmo(f27)
hr_time <- tlhr(f27)
sv_time <- tlsv(f27)
diff_time <- find_diff(f27)

time_dat <- ctime(f28)
co_time <- tlco(f28)
eo2_time <- tleo(f28)
mo2_time <- tlmo(f28)
hr_time <- tlhr(f28)
sv_time <- tlsv(f28)
diff_time <- find_diff(f28)

time_dat <- ctime(f29)
co_time <- tlco(f29)
eo2_time <- tleo(f29)
mo2_time <- tlmo(f29)
hr_time <- tlhr(f29)
sv_time <- tlsv(f29)
diff_time <- find_diff(f29)

time_dat <- ctime(f30)
co_time <- tlco(f30)
eo2_time <- tleo(f30)
mo2_time <- tlmo(f30)
hr_time <- tlhr(f30)
sv_time <- tlsv(f30)
diff_time <- find_diff(f30)

time_dat <- ctime(f31)
co_time <- tlco(f31)
eo2_time <- tleo(f31)
mo2_time <- tlmo(f31)
hr_time <- tlhr(f31)
sv_time <- tlsv(f31)
diff_time <- find_diff(f31)

time_dat <- ctime(f33)
co_time <- tlco(f33)
eo2_time <- tleo(f33)
mo2_time <- tlmo(f33)
hr_time <- tlhr(f33)
sv_time <- tlsv(f33)
diff_time <- find_diff(f33)

time_dat <- ctime(f34)
co_time <- tlco(f34)
eo2_time <- tleo(f34)
mo2_time <- tlmo(f34)
hr_time <- tlhr(f34)
sv_time <- tlsv(f34)
diff_time <- find_diff(f34)

time_dat <- ctime(f35)
co_time <- tlco(f35)
eo2_time <- tleo(f35)
mo2_time <- tlmo(f35)
hr_time <- tlhr(f35)
sv_time <- tlsv(f35)
diff_time <- find_diff(f35)

time_dat <- ctime(f36)
co_time <- tlco(f36)
eo2_time <- tleo(f36)
mo2_time <- tlmo(f36)
hr_time <- tlhr(f36)
sv_time <- tlsv(f36)
diff_time <- find_diff(f36)

time_dat <- ctime(f38)
co_time <- tlco(f38)
eo2_time <- tleo(f38)
mo2_time <- tlmo(f38)
hr_time <- tlhr(f38)
sv_time <- tlsv(f38)
diff_time <- find_diff(f38)

# ####

#Look at the data after sorting it according to treatment and then calculate means for each variable in each treatment####
#Sort data and exclude fish that are not relevant
fish_sub <- arrange(fish_dat, treatment) #sort according to treatment

#Discard fish
fish_sub <- filter(fish_sub, fishID != "f22" & fishID != "f23" & fishID != "f30" & fishID != "f3" & fishID != "f12")

#Split treatment column into salinity and surgical treatment
fish_sub <- rbind(mutate(filter(fish_sub, treatment == "fwc" | treatment == "fwl"), sal = "fw"),
                  mutate(filter(fish_sub, treatment == "swc" | treatment == "swl"), sal = "sw"))
fish_sub <- rbind(mutate(filter(fish_sub, treatment == "fwc"), surg = "control"),
                  mutate(filter(fish_sub, treatment == "fwl"), surg = "ligated"),
                  mutate(filter(fish_sub, treatment == "swc"), surg = "control"),
                  mutate(filter(fish_sub, treatment == "swl"), surg = "ligated" ))

#Means + SD + SEM
fish_mean <- fish_sub %>% 
  group_by(treatment) %>% 
  summarise_at(vars(-fishID, -salinity, -transonic, -date_start, -setup, -respirometer, -hb1, -hb2, -hct1, -hct2, -osm1, -osm2, -cK1, -cNa1, -cCl1, -cCa1, -cpH1, -cK2, 
                    -cNa2, -cCl2, -cCa2, -cpH2, -c2Diss, -date_start_num, -day_start, -wk, - respfish_ratio, -respvol_net, -hb_av), 
               list(~ mean(., na.rm = TRUE), ~sem(., na.rm = TRUE))) #means of each treatment
fish_mean$n <- count(fish_sub, vars = treatment) %>% pull(n) #add number of animals in each group (=n)

#Create a csv with the means
write.csv2(fish_mean, here("output/fish_mean.csv"), row.names = FALSE)
# ####

#Set up some objects referencing explanatory and response variables ####
vent_vars <- c("wetventmass", "compact", "spongy", "relwetvent", "dryventmass", "relcomp", "reldryvent")
spleen_vars <- c("spleenmass", "relspleen")
respir_vars <- c("setup", "respirometer")
time_vars <- c("date_start", "date_start_num", "day_start", "wk")
resting_vars <- c("smr_10", "co_10", "hr_10", "sv_10", "eo2_10")
max_vars <- c("mmr", "co_max", "hr_max", "sv_max", "eo2_max")
scope_vars <- c("as", "co_scope", "hr_scope", "sv_scope", "eo2_scope")
hema_vars <- c("ph", "osm_av", "hbcorr", "hb_drab", "hct_av")
ion_vars <- c("cK_av", "cNa_av", "cCl_av", "cCa_av", "cpH_av")
morph_vars <- c("bodymass", "length", "cf")
exper_vars <- c("fishID", "treatment", "salinity")

resp_vars <- c(vent_vars, spleen_vars, resting_vars, max_vars, scope_vars, hema_vars, ion_vars)
expl_vars <- c(respir_vars, time_vars, morph_vars, exper_vars)
all_vars <- c(resp_vars, expl_vars)
# ####

# Collect r2 in r2_all for quicker computation of R2 later on ####
r2_all <- read.csv2("r2_all.csv")

f33r2 <- f33$r2
f34r2 <- f34$r2
f35r2 <- f35$r2
f36r2 <- f36$r2
f38r2 <- f38$r2

r2_all <- cbind.na(r2_all, f33r2, f34r2, f35r2, f36r2, f38r2)
view(r2_all)

write.csv2(r2_all, here("r2_all.csv"), row.names = FALSE)
# ####


####################################################################################################################################
#################################################EXPERIMENT 1#######################################################################
######################################DATA EXPLORATION AND STATISTICS###############################################################
####################################################################################################################################
#Inspect and visualize data
#See how many of the r2-values are left when filtering out the ones <0.90####
fish_list <- list(f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f24, f25, 
                  f26, f27, f28, f29, f30, f31, f33, f34, f35, f36, f38)

#Create data frame with the fishID, respirometer, nr of observations and some r2-diagnostics
fishID <- filter(fish_dat, respirometer != "NA" & fishID != "f22" & fishID != "f23" & fishID != "f37")$fishID
resp <- filter(fish_dat, respirometer != "NA" & fishID != "f22" & fishID != "f23" & fishID != "f37")$respirometer
observations <- sapply(fish_list, function(x) nrow(filter(x, r2 != "NA")))
nr_ok_r2 <- sapply(fish_list, function(x) nrow(filter(x, r2 >= 0.90)))
nr_low_r2 <- sapply(fish_list, function(x) nrow(filter(x, r2 < 0.90)))
percent_ok_r2 <- sapply(fish_list, function(x) nrow(filter(x, r2 >= 0.90)) / nrow(x) * 100)
percent_low_r2 <- sapply(fish_list, function(x) (nrow(x) - nrow(filter(x, r2 >= 0.90))) / nrow(x) * 100)
r2_diag <- data.frame(fishID, resp, observations, nr_ok_r2, nr_low_r2, percent_ok_r2, percent_low_r2)

r2_diag
write.csv2(r2_diag, file = here("output/r2_diag.csv"), row.names = FALSE)
# ####

#Data diagnostics ####
fish_diag <- select(fish_sub, fishID, treatment, sal, surg, transonic, setup, respirometer, smr_10, mmr, co_10, co_max)
fish_diag

#METABOLIC OXYGEN CONSUMPTION
fpd <- function(data, x, y, color){
  x <- enquo(x)
  y <- enquo(y)
  color <- enquo(color)
  
  ggp <- ggplot(data, aes(x = !! x, y = !! y, fill = !! color)) +
    geom_boxplot() +
    scale_fill_manual(values = c("chocolate", "darkolivegreen", "darkorange", "darkmagenta"))
  
  return(ggp)
} #plotting function

#smr_10 boxplot colored by respirometer
fpd(fish_diag, treatment, smr_10, respirometer)

#mmr boxplot colored by respirometer
fpd(fish_diag, treatment, mmr, respirometer)

#smr_10 boxplot colored by transonic
fpd(fish_diag, treatment, smr_10, transonic)

#mmr boxplot colored by transonic
fpd(fish_diag, treatment, mmr, transonic)

#Count number of fishes within each treatment that were run in each respective respirometer
count(filter(fish_diag, treatment == "fwc"), vars = respirometer)
count(filter(fish_diag, treatment == "fwl"), vars = respirometer)
count(filter(fish_diag, treatment == "swc"), vars = respirometer)
count(filter(fish_diag, treatment == "swl"), vars = respirometer)

#Order them from high to low
arrange(fish_diag, -smr_10)
arrange(fish_diag, -mmr)

#Check the data of each respirometer
filter(fish_diag, respirometer == "r1")
filter(fish_diag, respirometer == "r2")
filter(fish_diag, respirometer == "r3")
filter(fish_diag, respirometer == "r4")

#CONCLUSION MO2
#No consistent pattern in over/under-estimations of MO2 in certain respirometers
#f45 might be an outlier as it has consistently higher values on both MO2 and CO
#f48 might be an outlier as well but it has a "normal" smr


#CARDIAC OUTPUT
#co_10 boxplot colored by flow probe
fpd(fish_diag, treatment, co_10, transonic)

#co_max boxplot colored by flow probe
fpd(fish_diag, treatment, co_max, transonic)

#co_10 boxplot colored by respirometer
fpd(fish_diag, treatment, co_10, respirometer)

#co_max boxplot colored by respirometer
fpd(fish_diag, treatment, co_max, respirometer)

#Count number of fishes within each treatment that were run with each respective transonic flow probe
count(filter(fish_diag, treatment == "fwc"), vars = transonic)
count(filter(fish_diag, treatment == "fwl"), vars = transonic)
count(filter(fish_diag, treatment == "swc"), vars = transonic)
count(filter(fish_diag, treatment == "swl"), vars = transonic)

#Order them from high to low
arrange(fish_diag, -co_10)
arrange(fish_diag, -co_max)

#Check the probes out individually
filter(fish_diag, transonic == "t732")
filter(fish_diag, transonic == "t792")
filter(fish_diag, transonic == "t794")
filter(fish_diag, transonic == "t795")

#Check the calibration equations for the flow probes
trans_eq

#CONCLUSION CARDIAC OUTPUT
#Flow probe 794 consistently gives higher CO-values than other flow probes (especially for co_10), 
#and we used this probe for 50% of fwc, whereas it was used for 22% fwl, 22% swc and 33% swl
#The CO-values are also highest in fwc in respirometer 1+2, and corresponds to when we used flow probe 794
#So, I think this flow probe could be the cause of measurement errors
# ####

#Flow probe diagnostics ####
fish_diag <- select(fish_sub, fishID, date_start, treatment, sal, surg, transonic, setup, respirometer, smr_10, mmr, co_10, co_max, hr_10, hr_max, eo2_10, eo2_max)
fish_diag

#We used this probe before the calibrations and it gave high values back then as well
arrange(fish_diag, date_start)

#Plot the variables colored by flow probe
#Pretty obviously flow probe t794 is up to no good
fpd <- function(data, x, y){
  x <- enquo(x)
  y <- enquo(y)
  ggp <- ggplot(data, aes(x = !! x, y = !! y, fill = transonic)) +
    geom_boxplot() +
    scale_fill_manual(values = c("chocolate", "darkolivegreen", "darkorange", "darkmagenta"))
  
  return(ggp)
} #plotting function

#smr_10 #mmr #co_10 #co_max #hr_10 #hr_max #eo2_10 #eo2_max
grid.arrange(fpd(fish_diag, treatment, smr_10), fpd(fish_diag, treatment, mmr),
             fpd(fish_diag, treatment, co_10), fpd(fish_diag, treatment, co_max),
             fpd(fish_diag, treatment, hr_10), fpd(fish_diag, treatment, hr_max),
             fpd(fish_diag, treatment, eo2_10), fpd(fish_diag, treatment, eo2_max))
# ####

#Check ALL the means out graphically ####
fish_mean_only <- select(fish_mean, treatment, 
                         wetventmass_mean, relwetvent_mean, compact_mean, spongy_mean, dryventmass_mean, reldryvent_mean, relcomp_mean,
                         spleenmass_mean, relspleen_mean,
                         ph_mean, hbcorr_mean, hb_drab_mean, hct_av_mean, osm_av_mean, 
                         cK_av_mean, cNa_av_mean, cCl_av_mean, cCa_av_mean, cpH_av_mean,
                         smr_10_mean, mmr_mean, as_mean, epoc_mean,
                         co_10_mean, co_max_mean, co_scope_mean,
                         hr_10_mean, hr_max_mean, hr_scope_mean,
                         sv_10_mean, sv_max_mean, sv_scope_mean,
                         eo2_10_mean, eo2_max_mean, eo2_scope_mean)
fish_mean_plot <- gather(fish_mean_only, key = "var", value = "value", -treatment)

fish_se_only <- select(fish_mean,  treatment, 
                       wetventmass_sem, relwetvent_sem, compact_sem, spongy_sem, dryventmass_sem, reldryvent_sem, relcomp_sem,
                       spleenmass_sem, relspleen_sem,
                       ph_sem, hbcorr_sem, hb_drab_sem, hct_av_sem, osm_av_sem, 
                       cK_av_sem, cNa_av_sem, cCl_av_sem, cCa_av_sem, cpH_av_sem,
                       smr_10_sem, mmr_sem, as_sem, epoc_sem,
                       co_10_sem, co_max_sem, co_scope_sem,
                       hr_10_sem, hr_max_sem, hr_scope_sem,
                       sv_10_sem, sv_max_sem, sv_scope_sem,
                       eo2_10_sem, eo2_max_sem, eo2_scope_sem)
fish_se_plot <- gather(fish_se_only, key = "var", value = "se", -treatment)
fish_mean_plot$se <- fish_se_plot$se

#PLOTS
all_vars
ggplot(fish_mean_plot, aes(x = treatment, y = value, fill = treatment)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = value - (1.96 * se), ymax = value + (1.96 * se)), width = 0.2) +
  scale_fill_manual(values = c("darkblue", "blue", "darkred", "red")) +
  facet_wrap(~var, scales = 'free') +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

vent_vars
ggplot(data = cbind(gather(fish_mean_only, key = "var_mean", value = "value", wetventmass_mean, compact_mean, spongy_mean, relwetvent_mean, dryventmass_mean, relcomp_mean, reldryvent_mean),
                    select(gather(fish_se_only, key = "var_se", value = "se", wetventmass_sem, compact_sem, spongy_sem, relwetvent_sem, dryventmass_sem, relcomp_sem, reldryvent_sem), -treatment)), 
       aes(x = treatment, y = value, fill = treatment)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = value - (1.96 * se), ymax = value + (1.96 * se)), width = 0.2) +
  scale_fill_manual(values = c("darkblue", "blue", "darkred", "red")) +
  facet_wrap(~var_mean, scales = 'free') +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

spleen_vars
ggplot(data = cbind(gather(fish_mean_only, key = "var_mean", value = "value", spleenmass_mean, relspleen_mean),
                    select(gather(fish_se_only, key = "var_se", value = "se", spleenmass_sem, relspleen_sem), -treatment)), 
       aes(x = treatment, y = value, fill = treatment)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = value - (1.96 * se), ymax = value + (1.96 * se)), width = 0.2) +
  scale_fill_manual(values = c("darkblue", "blue", "darkred", "red")) +
  facet_wrap(~var_mean, scales = 'free') +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

c(resting_vars, max_vars, scope_vars)
#See "Stacking columns for restin - scope - max" below

hema_vars
ggplot(data = cbind(gather(fish_mean_only, key = "var_mean", value = "value", ph_mean, osm_av_mean, hbcorr_mean, hb_drab_mean, hct_av_mean),
                    select(gather(fish_se_only, key = "var_se", value = "se", ph_sem, osm_av_sem, hbcorr_sem, hb_drab_sem, hct_av_sem), -treatment)), 
       aes(x = treatment, y = value, fill = treatment)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = value - (1.96 * se), ymax = value + (1.96 * se)), width = 0.2) +
  scale_fill_manual(values = c("darkblue", "blue", "darkred", "red")) +
  facet_wrap(~var_mean, scales = 'free') +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ion_vars
ggplot(data = cbind(gather(fish_mean_only, key = "var_mean", value = "value", cK_av_mean, cNa_av_mean, cCl_av_mean, cCa_av_mean, cpH_av_mean),
                    select(gather(fish_se_only, key = "var_se", value = "se", cK_av_sem, cNa_av_sem, cCl_av_sem, cCa_av_sem, cpH_av_sem), -treatment)), 
       aes(x = treatment, y = value, fill = treatment)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = value - (1.96 * se), ymax = value + (1.96 * se)), width = 0.2) +
  scale_fill_manual(values = c("darkblue", "blue", "darkred", "red")) +
  facet_wrap(~var_mean, scales = 'free') +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

#epoc & scopes
ggplot(data = cbind(gather(fish_mean_only, key = "var_mean", value = "value", epoc_mean, as_mean, co_scope_mean, hr_scope_mean, sv_scope_mean, eo2_scope_mean),
                    select(gather(fish_se_only, key = "var_se", value = "se", epoc_sem, as_sem, co_scope_sem, hr_scope_sem, sv_scope_sem, eo2_scope_sem), -treatment)),
       aes(x = treatment, y = value, fill = treatment)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = value - (1.96 * se), ymax = value + (1.96 * se)), width = 0.2) +
  scale_fill_manual(values = c("darkblue", "blue", "darkred", "red")) +
  facet_wrap(~var_mean, scales = 'free') +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

# ####

#Stacking columns for resting - scope - max ####
stacks <- select(fish_mean, treatment, 
                 smr_10_mean, smr_10_sem, as_mean, mmr_sem, 
                 co_10_mean, co_10_sem, co_scope_mean, co_max_sem, 
                 hr_10_mean, hr_10_sem, hr_scope_mean, hr_max_sem, 
                 sv_10_mean, sv_10_sem, sv_scope_mean, sv_max_sem, 
                 eo2_10_mean, eo2_10_sem, eo2_scope_mean, eo2_max_sem)

#Add resting values and scope values for stacking columns (resting + scope = max). Also add resting SEM and maximum SEM (not scope SEM)
stacks_mo2 <- cbind(select(gather(stacks, key = mo2, value = mo2v, smr_10_mean, as_mean), treatment, mo2, mo2v), 
                    select(gather(stacks, key = mo2_se, value = mo2_sem, smr_10_sem, mmr_sem), mo2_sem))
stacks_co <- cbind(select(gather(stacks, key = co, value = cov, co_10_mean, co_scope_mean), co, cov), 
                   select(gather(stacks, key = co_se, value = co_sem, co_10_sem, co_max_sem), co_sem))
stacks_hr <- cbind(select(gather(stacks, key = hr, value = hrv, hr_10_mean, hr_scope_mean), hr, hrv), 
                   select(gather(stacks, key = hr_se, value = hr_sem, hr_10_sem, hr_max_sem), hr_sem))
stacks_sv <- cbind(select(gather(stacks, key = sv, value = svv, sv_10_mean, sv_scope_mean), sv, svv),
                   select(gather(stacks, key = sv_se, value = sv_sem, sv_10_sem, sv_max_sem), sv_sem))
stacks_eo2 <- cbind(select(gather(stacks, key = eo2, value = eo2v, eo2_10_mean, eo2_scope_mean), eo2, eo2v),
                    select(gather(stacks, key = eo2_se, value = eo2_sem, eo2_10_sem, eo2_max_sem), eo2_sem))
stacks <- cbind(stacks_mo2, rbind(select(filter(select(mutate(stacks_mo2, mo2adj = mo2v - filter(stacks_mo2, mo2 == "smr_10_mean")$mo2v), mo2adj, mo2), mo2 == "smr_10_mean"), -mo2),
                                  select(filter(select(mutate(stacks_mo2, mo2adj = mo2v - mo2v + filter(stacks_mo2, mo2 == "smr_10_mean")$mo2v), mo2adj, mo2), mo2 == "as_mean"), -mo2)),
                stacks_co, rbind(select(filter(select(mutate(stacks_co, coadj = cov - filter(stacks_co, co == "co_10_mean")$cov), coadj, co), co == "co_10_mean"), -co),
                                 select(filter(select(mutate(stacks_co, coadj = cov - cov + filter(stacks_co, co == "co_10_mean")$cov), coadj, co), co == "co_scope_mean"), -co)),
                stacks_hr, rbind(select(filter(select(mutate(stacks_hr, hradj = hrv - filter(stacks_hr, hr == "hr_10_mean")$hrv), hradj, hr), hr == "hr_10_mean"), -hr),
                                 select(filter(select(mutate(stacks_hr, hradj = hrv - hrv + filter(stacks_hr, hr == "hr_10_mean")$hrv), hradj, hr), hr == "hr_scope_mean"), -hr)),
                stacks_sv, rbind(select(filter(select(mutate(stacks_sv, svadj = svv - filter(stacks_sv, sv == "sv_10_mean")$svv), svadj, sv), sv == "sv_10_mean"), -sv),
                                 select(filter(select(mutate(stacks_sv, svadj = svv - svv + filter(stacks_sv, sv == "sv_10_mean")$svv), svadj, sv), sv == "sv_scope_mean"), -sv)),
                stacks_eo2, rbind(select(filter(select(mutate(stacks_eo2, eo2adj = eo2v - filter(stacks_eo2, eo2 == "eo2_10_mean")$eo2v), eo2adj, eo2), eo2 == "eo2_10_mean"), -eo2),
                                  select(filter(select(mutate(stacks_eo2, eo2adj = eo2v - eo2v + filter(stacks_eo2, eo2 == "eo2_10_mean")$eo2v), eo2adj, eo2), eo2 == "eo2_scope_mean"), -eo2)))

mo2p <- ggplot(stacks, aes(x = treatment, y = mo2v, fill = factor(mo2, levels = c("as_mean", "smr_10_mean")))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymin = mo2v + mo2adj - (1.96 * mo2_sem), ymax = mo2v + mo2adj + (1.96 * mo2_sem), width = 0.2)) +
  labs(y = "SMR + MMR ", fill = "mo2") +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

cop <- ggplot(stacks[rev(order(stacks$co)),], aes(x = treatment, y = cov, fill = factor(co, levels = c("co_scope_mean", "co_10_mean")))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymin = cov + coadj - (1.96 * co_sem), ymax = cov + coadj + (1.96 * co_sem), width = 0.2)) +
  labs(y = "CO Rest + Max", fill = "co") +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

hrp <- ggplot(stacks, aes(x =treatment, y = hrv, fill = factor(hr, levels = c("hr_scope_mean", "hr_10_mean")))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymin = hrv + hradj - (1.96 * hr_sem), ymax = hrv + hradj + (1.96 * hr_sem), width = 0.2)) +
  labs(y = "HR Rest + Max", fill = "hr") +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

svp <- ggplot(stacks, aes(x = treatment, y = svv, fill = factor(sv, levels = c("sv_scope_mean", "sv_10_mean")))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymin = svv + svadj - (1.96 * sv_sem), ymax = svv +svadj + (1.96 * sv_sem), width = 0.2)) +
  labs(y = "SV Rest + Max", fill = "sv") +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

eo2p <- ggplot(stacks, aes(x = treatment, y = eo2v, fill = factor(eo2, levels = c("eo2_scope_mean", "eo2_10_mean")))) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymin = eo2v + eo2adj - (1.96 * eo2_sem), ymax = eo2v +eo2adj + (1.96 * eo2_sem), width = 0.2)) +
  labs(y = "EO2 Rest + Max", fill = "eo2") +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

#Print plots side by side
grid.arrange(mo2p, cop, hrp, svp, eo2p)
# ####

#Check correlations ####
pairs(select(fish_dat, date_start_num, day_start, wk))
pairs(select(fish_dat, bodymass, length, cf))
pairs(select(fish_dat, resting_vars))
corrplot(cor(select(fish_dat, resting_vars), use = "complete.obs"), method = "number")
pairs(select(time_dat, 2:7))
corrplot(cor(select(time_dat, 2:7), use = "complete.obs"), method = "number")
pairs(select(f1, co_kg, eo2, do2_c, hr, sv))
# ####

#Plot histograms of the variables to check for normality of response variables ####
hist_dat <- gather(fish_dat, key = "var", value = "value", all_of(c(vent_vars, spleen_vars, resting_vars, max_vars, scope_vars, hema_vars, ion_vars)), bodymass, length, cf)
hist_dat$wk <- as.character(hist_dat$wk)

ggplot(hist_dat, aes(x = value, fill = setup)) +
  geom_histogram(bins = 50) +
  facet_wrap(~var, scales = "free")

ggplot(hist_dat, aes(x = value, fill = respirometer)) +
  geom_histogram(bins = 50) +
  facet_wrap(~var, scales = "free")

ggplot(hist_dat, aes(x = value, fill = wk)) +
  geom_histogram(bins = 50) +
  facet_wrap(~var, scales = "free")
# ####

#Plot ####
#mo2
ggplot(gather(fish_sub, key = "var", value = "val", smr_10, mmr, as), aes(x = sal, y = val, color = surg, fill = respirometer)) +
  geom_boxplot(size = 1) +
  facet_wrap(~var, scales = "free") +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_fill_manual(values = c("chocolate", "darkolivegreen", "darkorange", "darkmagenta"))

#co
ggplot(gather(fish_sub, key = "var", value = "val", co_10, co_max, co_scope), aes(x = sal, y = val, color = surg, fill = respirometer)) +
  geom_boxplot(size = 1) +
  facet_wrap(~var, scales = "free") +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_fill_manual(values = c("chocolate", "darkolivegreen", "darkorange", "darkmagenta"))

#hr
ggplot(gather(fish_sub, key = "var", value = "val", hr_10, hr_max, hr_scope), aes(x = sal, y = val, color = surg, fill = respirometer)) +
  geom_boxplot(size = 1) +
  facet_wrap(~var, scales = "free") +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_fill_manual(values = c("chocolate", "darkolivegreen", "darkorange", "darkmagenta"))

#sv
ggplot(gather(fish_sub, key = "var", value = "val", sv_10, sv_max, sv_scope), aes(x = sal, y = val, color = surg, fill = respirometer)) +
  geom_boxplot(size = 1) +
  facet_wrap(~var, scales = "free") +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_fill_manual(values = c("chocolate", "darkolivegreen", "darkorange", "darkmagenta"))

#eo2
ggplot(gather(fish_sub, key = "var", value = "val", eo2_10, eo2_max, eo2_scope), aes(x = sal, y = val, color = surg, fill = respirometer)) +
  geom_boxplot(size = 1) +
  facet_wrap(~var, scales = "free") +
  scale_color_manual(values = c("darkred", "darkblue")) +
  scale_fill_manual(values = c("chocolate", "darkolivegreen", "darkorange", "darkmagenta"))

#hematological variables
ggplot(gather(filter(fish_sub, hct_av >= 20), key = "var", value = "val", hema_vars), aes(x = sal, y = val, fill = surg)) +
  geom_boxplot() +
  facet_wrap(~var, scales = "free")

#morphological variables
ggplot(gather(fish_sub, key = "var", value = "val", morph_vars), aes(x = sal, y = val, fill = surg)) +
  geom_boxplot() +
  facet_wrap(~var, scales = "free")
# ####

#Plot timeline data means for each treatment ####
#hours after chase
order_of_fish <- arrange(fish_sub, date_start, fishID)$fishID
h_after <- data.frame(fishID = order_of_fish, hours_after = rbind(time_post(f2), time_post(f3), time_post(f4), time_post(f5), time_post(f6), time_post(f7), time_post(f8), 
                                                                  time_post(f9), time_post(f10), time_post(f11), time_post(f12), time_post(f13), time_post(f14), time_post(f15), 
                                                                  time_post(f16), time_post(f17), time_post(f18), time_post(f19), time_post(f20), time_post(f21), time_post(f24), time_post(f25),
                                                                  time_post(f26), time_post(f27), time_post(f28), time_post(f29), time_post(f30), time_post(f31)))
mean(h_after$hours_after)
min(filter(h_after, hours_after > min(h_after$hours_after))$hours_after) #=1189.5 minutes
h_after

#Conclusion: We have around 15 hours of recording post-chase

tlco_mean <- select(co_time, -rco) %>%
  group_by(treat) %>% 
  summarise_at(vars(-fishID,), 
               list(~ mean(., na.rm = TRUE), ~sem(., na.rm = TRUE))) #means of each treatment

tleo_mean <- select(eo2_time, -reo) %>%
  group_by(treat) %>% 
  summarise_at(vars(-fishID,), 
               list(~ mean(., na.rm = TRUE), ~sem(., na.rm = TRUE))) #means of each treatment

tlmo_mean <- select(mo2_time, -rmr) %>%
  group_by(treat) %>% 
  summarise_at(vars(-fishID,), 
               list(~ mean(., na.rm = TRUE), ~sem(., na.rm = TRUE))) #means of each treatment

tlhr_mean <- select(hr_time, -rhr) %>%
  group_by(treat) %>% 
  summarise_at(vars(-fishID,), 
               list(~ mean(., na.rm = TRUE), ~sem(., na.rm = TRUE))) #means of each treatment

tlsv_mean <- select(sv_time, -rsv) %>%
  group_by(treat) %>% 
  summarise_at(vars(-fishID,), 
               list(~ mean(., na.rm = TRUE), ~sem(., na.rm = TRUE))) #means of each treatment

tl_plot <- cbind(select(gather(tlco_mean, key = "time", value = "mean_co", 2:11), treat, time, mean_co),
                 select(gather(tlco_mean, key = "time", value = "sem_co", 12:21), sem_co),
                 select(gather(tleo_mean, key = "time", value = "mean_eo2", 2:11), mean_eo2),
                 select(gather(tleo_mean, key = "time", value = "sem_eo2", 12:21), sem_eo2),
                 select(gather(tlmo_mean, key = "time", value = "mean_mo2", 2:11), mean_mo2),
                 select(gather(tlmo_mean, key = "time", value = "sem_mo2", 12:21), sem_mo2),
                 select(gather(tlhr_mean, key = "time", value = "mean_hr", 2:11), mean_hr),
                 select(gather(tlhr_mean, key = "time", value = "sem_hr", 12:21), sem_hr),
                 select(gather(tlsv_mean, key = "time", value = "mean_sv", 2:11), mean_sv),
                 select(gather(tlsv_mean, key = "time", value = "sem_sv", 12:21), sem_sv))

tl_plot <- rbind(mutate(filter(tl_plot, time == "sco_mean"), time_min = -30),
                 mutate(filter(tl_plot, time == "chase_mean"), time_min = 0),
                 mutate(filter(tl_plot, time == "post_10_mean"), time_min = 10),
                 mutate(filter(tl_plot, time == "post_20_mean"), time_min = 20),
                 mutate(filter(tl_plot, time == "post_30_mean"), time_min = 30),
                 mutate(filter(tl_plot, time == "post_60_mean"), time_min = 60),
                 mutate(filter(tl_plot, time == "post_120_mean"), time_min = 120),
                 mutate(filter(tl_plot, time == "post_300_mean"), time_min = 300),
                 mutate(filter(tl_plot, time == "post_600_mean"), time_min = 600),
                 mutate(filter(tl_plot, time == "post_900_mean"), time_min = 900))

tl_fwc <- filter(tl_plot, treat == "fwc")
tl_fwl <- filter(tl_plot, treat == "fwl")
tl_swc <- filter(tl_plot, treat == "swc")
tl_swl <- filter(tl_plot, treat == "swl")

#CO FW
cofw <- ggplot(filter(tl_plot, treat == "fwc" | treat == "fwl"), aes(x = time_min, y = mean_co, color = treat)) +
  geom_line(data = tl_fwc, aes(x = time_min, y = mean_co)) + geom_point(data = tl_fwc, aes(x = time_min, y = mean_co)) +
  geom_line(data = tl_fwl, aes(x = time_min, y = mean_co)) + geom_point(data = tl_fwl, aes(x = time_min, y = mean_co)) +
  geom_errorbar(aes(ymin = mean_co - sem_co, ymax = mean_co + sem_co))

#CO SW
cosw <- ggplot(filter(tl_plot, treat == "swc" | treat == "swl"), aes(x = time_min, y = mean_co, color = treat)) +
  geom_line(data = tl_swc, aes(x = time_min, y = mean_co)) + geom_point(data = tl_swc, aes(x = time_min, y = mean_co)) +
  geom_line(data = tl_swl, aes(x = time_min, y = mean_co)) + geom_point(data = tl_swl, aes(x = time_min, y = mean_co)) +
  geom_errorbar(aes(ymin = mean_co - sem_co, ymax = mean_co + sem_co))

#EO2 FW
eofw <- ggplot(filter(tl_plot, treat == "fwc" | treat == "fwl"), aes(x = time_min, y = mean_eo2, color = treat)) +
  geom_line(data = tl_fwc, aes(x = time_min, y = mean_eo2)) + geom_point(data = tl_fwc, aes(x = time_min, y = mean_eo2)) +
  geom_line(data = tl_fwl, aes(x = time_min, y = mean_eo2)) + geom_point(data = tl_fwl, aes(x = time_min, y = mean_eo2)) +
  geom_errorbar(aes(ymin = mean_eo2 - sem_eo2, ymax = mean_eo2 + sem_eo2))

#EO2 sW
eosw <- ggplot(filter(tl_plot, treat == "swc" | treat == "swl"), aes(x = time_min, y = mean_eo2, color = treat)) +
  geom_line(data = tl_swc, aes(x = time_min, y = mean_eo2)) + geom_point(data = tl_swc, aes(x = time_min, y = mean_eo2)) +
  geom_line(data = tl_swl, aes(x = time_min, y = mean_eo2)) + geom_point(data = tl_swl, aes(x = time_min, y = mean_eo2)) +
  geom_errorbar(aes(ymin = mean_eo2 - sem_eo2, ymax = mean_eo2 + sem_eo2))

#MO2 FW
mofw <- ggplot(filter(tl_plot, treat == "fwc" | treat == "fwl"), aes(x = time_min, y = mean_mo2, color = treat)) +
  geom_line(data = tl_fwc, aes(x = time_min, y = mean_mo2)) + geom_point(data = tl_fwc, aes(x = time_min, y = mean_mo2)) +
  geom_line(data = tl_fwl, aes(x = time_min, y = mean_mo2)) + geom_point(data = tl_fwl, aes(x = time_min, y = mean_mo2)) +
  geom_errorbar(aes(ymin = mean_mo2 - sem_mo2, ymax = mean_mo2 + sem_mo2))

#MO2 SW
mosw <- ggplot(filter(tl_plot, treat == "swc" | treat == "swl"), aes(x = time_min, y = mean_mo2, color = treat)) +
  geom_line(data = tl_swc, aes(x = time_min, y = mean_mo2)) + geom_point(data = tl_swc, aes(x = time_min, y = mean_mo2)) +
  geom_line(data = tl_swl, aes(x = time_min, y = mean_mo2)) + geom_point(data = tl_swl, aes(x = time_min, y = mean_mo2)) +
  geom_errorbar(aes(ymin = mean_mo2 - sem_mo2, ymax = mean_mo2 + sem_mo2))

#HR FW
hrfw <- ggplot(filter(tl_plot, treat == "fwc" | treat == "fwl"), aes(x = time_min, y = mean_hr, color = treat)) +
  geom_line(data = tl_fwc, aes(x = time_min, y = mean_hr)) + geom_point(data = tl_fwc, aes(x = time_min, y = mean_hr)) +
  geom_line(data = tl_fwl, aes(x = time_min, y = mean_hr)) + geom_point(data = tl_fwl, aes(x = time_min, y = mean_hr)) +
  geom_errorbar(aes(ymin = mean_hr - sem_hr, ymax = mean_hr + sem_hr))

#HR SW
hrsw <- ggplot(filter(tl_plot, treat == "swc" | treat == "swl"), aes(x = time_min, y = mean_hr, color = treat)) +
  geom_line(data = tl_swc, aes(x = time_min, y = mean_hr)) + geom_point(data = tl_swc, aes(x = time_min, y = mean_hr)) +
  geom_line(data = tl_swl, aes(x = time_min, y = mean_hr)) + geom_point(data = tl_swl, aes(x = time_min, y = mean_hr)) +
  geom_errorbar(aes(ymin = mean_hr - sem_hr, ymax = mean_hr + sem_hr))

#SV FW
svfw <- ggplot(filter(tl_plot, treat == "fwc" | treat == "fwl"), aes(x = time_min, y = mean_sv, color = treat)) +
  geom_line(data = tl_fwc, aes(x = time_min, y = mean_sv)) + geom_point(data = tl_fwc, aes(x = time_min, y = mean_sv)) +
  geom_line(data = tl_fwl, aes(x = time_min, y = mean_sv)) + geom_point(data = tl_fwl, aes(x = time_min, y = mean_sv)) +
  geom_errorbar(aes(ymin = mean_sv - sem_sv, ymax = mean_sv + sem_sv))

#SV SW
svsw <- ggplot(filter(tl_plot, treat == "swc" | treat == "swl"), aes(x = time_min, y = mean_sv, color = treat)) +
  geom_line(data = tl_swc, aes(x = time_min, y = mean_sv)) + geom_point(data = tl_swc, aes(x = time_min, y = mean_sv)) +
  geom_line(data = tl_swl, aes(x = time_min, y = mean_sv)) + geom_point(data = tl_swl, aes(x = time_min, y = mean_sv)) +
  geom_errorbar(aes(ymin = mean_sv - sem_sv, ymax = mean_sv + sem_sv))

grid.arrange(cofw, cosw, eofw, eosw, mofw, mosw, hrfw, hrsw, svfw, svsw, layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10)))

pdf(file = "figures/tlplot.pdf", width = 8, height = 8)
grid.arrange(cofw, cosw)
grid.arrange(eofw, eosw)
grid.arrange(mofw, mosw)
grid.arrange(hrfw, hrsw)
grid.arrange(svfw, svsw)
dev.off()
# ####