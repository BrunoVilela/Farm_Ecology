# In this script we apply the ENMval to run the models for every 
# first domesticated species

# Packages
library(dismo)
library(maptools)
library(vegan)
library(ENMeval)

# Projection data (1900-1949)
env.path <- list.files("Data/bio # CCSM_Historical(1900-1949)", 
                       pattern = ".bil", full.names = TRUE)
env.path <- env.path[seq(1, 38, 2)]
proj.data <- stack(env.path)
names(proj.data) <- gsub("bio_._CCSM_Historical.1900.1949._", "", 
                         names(proj.data))
data("wrld_simpl")
proj.data <- mask(proj.data, wrld_simpl)


# Soil
path <- "/Users/bvilela/Box Sync/Bruno/WUSTL projects/culture diversity (now domestication)/Dataset"
soils <- stack(list.files(path,
                          pattern = "sq", full.names = T))
plot(soils)
soils <- aggregate(soils, fact = res(proj.data)[1]/res(soils)[1], na.rm = TRUE, mean)
soils <- resample(soils, proj.data)
elev <- raster(paste0(path, "/GloElev_5min.asc"))
elev <- aggregate(elev, fact = res(proj.data)[1]/res(elev)[1], na.rm = TRUE, mean)
elev <- resample(elev, proj.data)
plot(elev)
proj.data <- stack(proj.data, soils, elev)

# Projection data (1950-1999)
env.path <- list.files("Data/bio # CCSM_Modern(1950-1999)", 
                       pattern = ".bil", full.names = TRUE)
env.path <- env.path[seq(1, 38, 2)]
proj.data2 <- stack(env.path)
names(proj.data2) <- gsub("bio_._CCSM_Modern.1950.1999._", "", 
                          names(proj.data2))
proj.data2 <- stack(proj.data2, soils, elev)
proj.data2 <- mask(proj.data2, wrld_simpl)


# Subset
vars <- paste0("bio", 12:19)
proj.data <- stack(proj.data, layers = which(!names(proj.data) %in% vars))
proj.data2 <- stack(proj.data2, layers = which(!names(proj.data2) %in% vars))


# List of files
files <- list.files('occurrences/', full.names = TRUE)
# Empty lists
n <- length(files)

pred.modern <- pred.hist <- pred.holo <- pred.cont <- list()
model <- test <- eval <- AUCs_enmval <- list()
x <- 0
numbers.list <- c()
xy_bg <- xyFromCell(raster(proj.data2, 1), which(!is.na(values(raster(proj.data2, 1)))))
for (i in 1:n) {
  x <- x + 1
  print(n - i)
  occ.data <- read.csv(files[i])
  numbers.list <- c(numbers.list, i)
  pres <- na.omit(occ.data)
  pres <- pres[!is.na(extract(raster(proj.data2, 1), pres)), ]
  test.few <- nrow(pres) > 10
  bg <- xy_bg[sample(1:nrow(xy_bg), 10000, replace = TRUE), ]
  if (test.few) {
    test[[x]] <- ENMevaluate(pres, proj.data2, bg.coords = bg, method = "block", RMvalues = 1:10,
                             fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                             bin.output = F)
  } else {
    test[[x]] <- ENMevaluate(pres, proj.data2, bg.coords = bg, method = "jackknife", RMvalues = 1:10,
                             fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                             bin.output = F)
  }
  qual <- which.max(test[[x]]@results$w.AIC)
  
  AUCs_enmval[[x]] <- test[[x]]@results$Mean.AUC[[qual]]
  
  max <- test[[x]]@models[[qual]] 
  model[[x]] <- max
  pred.hist[[x]] <- predict(max, proj.data)
  test_env <- unlist(extract(proj.data2, pres))
  eval[[x]] <- evaluate(max@presence, max@absence, max)
}

# Run evaluations
evaluations <- sapply(eval, function(x){x@auc})

# Histogram AUCs
aucs <-unlist(AUCs_enmval)
tiff("Figures/histogram_auc.tiff", 30, 20, res = 300, units = "cm")
par(mar = c(5, 5, 2, 2))
hist(aucs, 5, xlab = "AUC", main = "", col = "gray80", cex.lab = 2, xlim = c(0, 1),
     ylim = c(0, length(aucs)))
abline(v = 0.7, lty = 2)
dev.off()

# Name the lists
spp <- gsub("occurrences//", "", files)
spp <- gsub(".csv", "", spp)
names(aucs) <- names(evaluations) <- names(pred.hist) <- names(model) <- names(eval) <- spp

save(eval, file = "Data/eval.RData")
save(model, file = "Data/model.RData")
save(pred.hist, file = "Data/pred_hist.RData")
save(evaluations, file = "Data/evaluations.RData")
save(files, file = "Data/species_name_sdm.RData")
save(aucs, file = "Data/auc_values.RData")
# Next buffer decision