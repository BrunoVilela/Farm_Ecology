# In this script we use a sensitivity analysis to identify the buffer 
# that maximizes the relationship between farming and domesticated species richness
# It requires the previous SDMs.R results and the mydata file from the analysis folder.

# Packages
library(rgeos)
library(raster)
library(maptools)
library(betareg)
library(fields)
library(mgcv)
library(gnm)
library(sampler)
library(letsR)
library(dismo)
library(lme4)
library(parallel)
library(tidyverse)

# Load
load("Data/eval.RData")
load("Data/model.RData")
load("Data/pred_hist.RData")
mydata <- read.scv("Data/mydata.csv")
origins <- rgdal::readOGR('Data/Origins_updated.shp')
species.list <- read.csv('Data/species_list.csv')
sp.names <- names(pred.hist)

# Farming data binary
mydata$farming01 <- (mydata$farming - min(mydata$farming)) / 
  (max(mydata$farming) - min(mydata$farming))

# Change origins names
origins <- origins[order(origins@data$CONTINENT), ]
origins@data$CONTINENT <- c("Central/Southern Andes",
                            "Chinese loess plateau",
                            "Eastern North America",
                            "Ethiopian plateau",
                            "Fertile Crescent",
                            "Ganges and eastern Indian plains",
                            "Japanese islands",
                            "Lower-Middle Yangtze",
                            "Mesoamerica, lowlands and highlands",
                            "Northern Lowland South America",
                            "New Guinea",
                            "Northwestern Lowland South America",
                            "South India",
                            "Savannahs of Western India",
                            "Lingnan (tropical south China)",
                            "Southwestern Amazonia",
                            "Sudanic Savannah",
                            "West African Savannah/Sahe",
                            "Western Yunnan/Eastern Tibet",
                            "West African tropical forest")
# Societies coordinates 
coords <- cbind(mydata$Longitude, mydata$Latitude)


# Get the thresholds values for each species 
thresholds <- lapply(eval, dismo::threshold)
pallete <- colorRampPalette(c("#edf8fb", "#b2e2e2", "#66c2a4",
                              "#2ca25f", "#006d2c"))
vals <- NULL
buf <- seq(10, 40000, 10)
n <- length(buf)
thresholds2 <- lapply(thresholds, function(x, y){c(as.numeric(x), y)}, y = vals)
thresholds2 <- lapply(thresholds2, function(x){ifelse(x < 0, 0, x)})

# thresh function
# Function for the threshold
threshold2 <- function(x, val = .9, thre = NULL) {
  valores <- values(x)
  if (is.null(thre)) {
    thre <- quantile(x = valores, probs = val, na.rm = TRUE) # I AM HERE
  }
  valores[valores > thre] <- 1
  valores[valores < thre] <- 0
  values(x) <- valores
  return(x)
}

# Define the original projection of the origins object
proj4string(origins) <- CRS("+proj=longlat +datum=WGS84")


# Create buffers before
origins.extract <- extract(pred.hist[[1]], origins, cellnumbers = TRUE)
xy <- xyFromCell(pred.hist[[1]], 1:ncell(pred.hist[[1]]))
rownames(xy) <- 1:ncell(pred.hist[[1]])



func_dist <- function(origins.extract, xy) {
  fields::rdist.earth(xy, 
                      xy[origins.extract[, 1], ],
                      miles = FALSE)
}
cl <- makeCluster(6)
dist.geo <- parLapply(cl, origins.extract, func_dist, xy)
stopCluster(cl)


# Save the results because it takes a while to run (few minutes)
save(dist.geo, file = "Data/dist.geo.RData")
# load("Data/dist.geo.RData")

cl <- makeCluster(6)
dist.geo.min <- parLapply(cl, dist.geo, function(x){apply(x, 1, min)})
stopCluster(cl)
rm(dist.geo)
save(dist.geo.min, file = "Data/dist_geo_min.RData")
# load("Data/dist_geo_min.RData")


# Empty lists for the different models and maps
rich.buf.I <- model.results.lm <- model.results.lmer <- model.data.list <- list()

# Loop
count <- 0 # counter

# Keep track in a log file

for (j in 1:n) {
  print(n - j)
  masked.pred <- list()
  for (i in 1:length(pred.hist)) {
    subregions <- as.character(species.list[sp.names[i] == species.list[, 1], 3])
    which.sub <- which(origins@data$CONTINENT %in% subregions)
    if (length(which.sub) == 1) {
      origins.sp <- dist.geo.min[[which.sub]]
    } else {
      origins.sp <- apply(do.call(rbind, dist.geo.min[which.sub]), 2, min)
    }
    origins.sp.n <- c(as.numeric(names(origins.sp)[origins.sp > buf[j]]))
    masked.pred[[i]] <- pred.hist[[i]]
    values(masked.pred[[i]])[origins.sp.n] <- NA
  }
  for (k in 1:length(thresholds2[[1]])) {
    count <- count + 1
    # Threshold - buffer
    r.threshold <- list()
    for (i in 1:length(masked.pred)) {
      r.threshold[[i]] <- threshold2(masked.pred[[i]], 
                                     val = thresholds2[[i]][k]) # change
    }
    # Make pam
    pam <- stack(r.threshold)
    xy <- xyFromCell(pam, cell = 1:ncell(pam))
    pam.ma <- cbind(xy, values(pam))
    rich.pam <- rowSums(pam.ma[, -(1:2)], na.rm = TRUE)
    rich <- raster(pam, 1)
    values(rich) <- rich.pam 
    colnames(pam.ma) <- c("Longitude(x)", "Latitude(y)", sp.names)
    pam.ma[is.na(pam.ma)] <- 0
    PAM.sp <- list("Presence_and_Absence_Matrix" = pam.ma, "Richness_Raster" =  rich,
                   "Species_name" = sp.names)
    class(PAM.sp) <- "PresenceAbsence"
    save(PAM.sp, file = paste0("Data/PAM_sp_", count, ".RData"))
    rich.buf.I[[count]] <- mask(rich, pred.hist[[1]])
    suitability <- raster::extract(rich.buf.I[[count]], coords)
    rem <- is.na(suitability)
    model.data <- data.frame(farming = mydata$farming01[!rem], 
                             suitability = suitability[!rem],
                             family = as.factor(mydata$Language_family[!rem]))
    model.data.list[[count]] <- model.data
    model.results.lmer[[count]] <- lmer(farming ~ suitability + (1 | family), data = model.data)
    model.results.lm[[count]] <- lm(farming ~ suitability, data = model.data)
    if (count %% 100 == 0)
      save(model.data.list, model.results.lm, model.results.lmer, rich.buf.I, file = "Data/model_buffer_results.RData")
  }
}

r2s <- sapply(model.results.lm, function(x){summary(x)[["r.squared"]]})
names(r2s) <- paste(rep(buf, each = length(thresholds2[[1]])), 
                    rep(c(names(thresholds[[1]]), vals), length(buf)),
                    sep = "_")[1:length(r2s)]

nome <- gsub(".*_", "", names(which.max(r2s)))
r2s_no <- r2s[pos_nos]
choosen <- which.max(r2s_no)
r2s_no[choosen]


ts <- sapply(model.results.lm, function(x){summary(x)$coefficients[6]})
betas <- sapply(model.results.lm, function(x){summary(x)$coefficients[2]})

buf_sizes <- as.numeric(gsub("_.*", "", names(r2s_no)))
data_plot <- tibble(buf = buf_sizes,
                    `R^2` = r2s_no, 
                    `t-values` = ts[pos_nos],
                    Coefficient = betas[pos_nos])

# Plot analysis
g <- data_plot %>%
  pivot_longer(-buf) %>%
  mutate(name = factor(name)) %>%
  ggplot(aes(x = buf, y = value)) +
  geom_line() +
  xlab("Buffer size (km)") +
  ylab("") +
  theme_classic() +
  facet_grid(name ~., scales = "free_y", switch = "y", labeller = label_parsed) +
  theme(text = element_text(size = 16),
        strip.placement = "outside",
        strip.background = element_blank())
g
ggsave("FigureS2.tiff", g, width = 10, height = 7)


