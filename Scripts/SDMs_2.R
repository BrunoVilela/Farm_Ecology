# Package
library(rgeos)
library(maptools)
library(maps)
library(dismo)

# Load
load("Data/pred_hist.RData")
load("Data/eval.RData")
load("Data/model.RData")
load("Data/Buffer_size.RData")

# Data
sp.names <- list.files('occurrences/')
sp.names <- gsub(".csv", "", sp.names)
species.list <- read.csv('Data/species_list.csv')
origins <- readShapePoly('Data/Origins/origins4.shp')
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

map()
plot(origins, add = TRUE, col = "red")
polygonsLabel(origins, origins@data$CONTINENT,
              method = "centroid", cex=.8)


masked.pred <- list()
buffer = 50
# Buffers
for (i in 1:length(pred.hist)) {
  print(i)
  subregions <- as.character(species.list[sp.names[i] == species.list[, 1], 3])
  origins.sp <- origins[origins@data$CONTINENT %in% subregions, ]
  buffers <- gBuffer(origins.sp, width = buffer) # Buffer size that maximize our analysis
  masked.pred[[i]] <- mask(pred.hist[[i]], buffers)
}

# No sut
nosut <- function(x){
  v <- values(x)
  v[!is.na(v)] <- 1
  values(x) <- v
  return(x)
}
masked.pred.nosut <- lapply(masked.pred, nosut)
rich.nosut <- stack(masked.pred.nosut)
rich.nosut <- sum(rich.nosut, na.rm = TRUE)
save(rich.nosut, file = "data/rich.nosut.RData")


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


# Threshold - buffer
thresholds <- lapply(eval, threshold)
r.threshold <- list()
for (i in 1:length(masked.pred)) {
  r.threshold[[i]] <-threshold2(masked.pred[[i]], 
                                thre = thresholds[[i]]$no_omission)
}

rich <- do.call(raster::stack, r.threshold)
rich <- sum(rich, na.rm = TRUE)

pallete <- colorRampPalette(c("#edf8fb", "#b2e2e2", "#66c2a4",
                              "#2ca25f", "#006d2c"))
rich <- mask(rich, pred.hist[[1]])
plot(rich, col = c("gray80", pallete(100)), axes = FALSE, box = FALSE)



# Threshold - no-buffer
r.threshold <- list()
for (i in 1:length(pred.hist)) {
  r.threshold[[i]] <-threshold2(pred.hist[[i]], 
                                thre = thresholds[[i]]$no_omission)
}
rich.nobuf <- do.call(raster::stack, r.threshold)
rich.nobuf <- sum(rich.nobuf, na.rm = TRUE)
pallete <- colorRampPalette(c("#edf8fb", "#b2e2e2", "#66c2a4",
                              "#2ca25f", "#006d2c"))
rich.nobuf <- mask(rich.nobuf, pred.hist[[1]])
plot(rich.nobuf, col = c("gray80", pallete(100)), axes = FALSE, box = FALSE)

data("wrld_simpl")
rich_plot <- mask(rich, wrld_simpl)
pdf(file = "Map_potential.pdf", width = 25, height = 20)
plot(rich_plot, col = c("gray80", pallete(100)), 
     axes = FALSE, box = FALSE, 
     asp = c(1, 1), ylim = c(-60, 90))
dev.off()
save(rich, file = "data/rich.Rdata")
save(rich.nobuf, file = "data/rich.nobuf.Rdata")

# Making the data into a PresenceAbsence object
library(letsR)
pam <- stack(r.threshold)
xy <- xyFromCell(pam, cell = 1:ncell(pam))
pam.ma <- cbind(xy, values(pam))
rich.pam <- rowSums(values(pam))
r.pam <- raster(pam, 1)
values(r.pam) <- rich.pam 
colnames(pam.ma) <- c("Longitude(x)", "Latitude(y)", species.names)
PAM.sp <- list("Presence_and_Absence_Matrix" = pam.ma, "Richness_Raster" = r.pam,
               "Species_name" = species.names)
class(PAM.sp) <- "PresenceAbsence"
plot(PAM.sp)
summary(PAM.sp)
save(PAM.sp, file = "PAM_sp.Rdata")
