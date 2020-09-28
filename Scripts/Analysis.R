# Analysis

# Packages
library(vegan)
library(maptools)
library(car)
library(spdep)
library(letsR)
library(lme4)
library(lmerTest)
library(sjstats)
library(letsR)

# Source this script
source("Scripts/Fun_Rnd_Categories.R")

# Data 
load("Data/Data4Bruno.2.Rdata")
load("Plots_data.RData")

rm(list = ls()[-which(ls() %in% c('MammalDiv_raster', 'VascPlant_CoKrig_KreftJetz2007_raster',
                                  'mydata', 'rich_plot', 'rich.nobuf'))])


# To remove australia and california run this code
# Shapefiles australia and california

# data(wrld_simpl)
# aus <- wrld_simpl[wrld_simpl@data$NAME == "Australia", ]
# ca <- shapefile("Data/ca-state-boundary/CA_State_TIGER2016.shp")
# 
# ca <- spTransform(ca, projection(aus))
# 
# plot(rich_plot)
# plot(aus, add = T)
# plot(ca, add = T)
# 
# coords <- data.frame(Longitude = mydata$Longitude, Latitude = mydata$Latitude)
# coordinates(coords) <- ~ Longitude + Latitude
# proj4string(coords) <- proj4string(aus)
# plot(aus)
# points(coords)
# aus_pnts <- apply(over(coords, aus), 1, function(x){!all(is.na(x))})
# cas_pnts <-  apply(over(coords, ca), 1, function(x){!all(is.na(x))})
# 
# plot(wrld_simpl)
# points(coords, pch = 20, col = "gray")
# points(coords[aus_pnts, ], pch = 20, col = "blue")
# points(coords[cas_pnts, ], pch = 20, col = "blue")
# 
# mydata <- mydata[!aus_pnts & !cas_pnts, ]


# Extract diversity data for societies
myrasters<-c('MammalDiv_raster', 'VascPlant_CoKrig_KreftJetz2007_raster', 'rich_plot', 'rich.nobuf')
prj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 +units=km")

for (i in 1:dim(mydata)[1]){
  mycoords <- as.data.frame(mydata[i,c("Longitude","Latitude")])
  mysp <- SpatialPoints(mycoords, proj4string = prj )
  
  for (j in 1:length(myrasters)) {
    mydata[i,gsub('_raster','',
                  x=myrasters[j])] <- eval(parse(text=c("mean( unlist(extract(", myrasters[j], ", mysp)), na.rm = T)")))
  }
  
  print(c(i, "out of", dim(mydata)[1], "samples"))
  
}

###############################
# explore best neighborhood structure


resu <- matrix(ncol = 2, nrow = 100)
for (j in 1:1) {
  print(j)
  kn <- knearneigh(cbind(mydata$Longitude, mydata$Latitude),
                   k = j, longlat = TRUE)
  nbs <- knn2nb(kn, sym = FALSE) # 10 neighbors
  nbs.sym <- knn2nb(kn, sym = TRUE) # 10 neighbors
  
  n <- length(nbs)
  neighProd.sym <- neighProd <- numeric(n)
  
  for(i in 1:n) {
    # Farming around
    neighProd[i] <- mean(mydata$farming[nbs[[i]]])
    neighProd.sym[i] <- mean(mydata$farming[nbs.sym[[i]]])
  }
  mydata$neighProd <- neighProd
  mydata$neighProd.sym <- neighProd.sym
  
  mymodlin_producers_spa.1 <- lm(farming ~ neighProd,
                                 data = mydata)
  mymodlin_producers_spa.2 <- lm(farming ~ neighProd.sym,
                                 data = mydata)
  
  resu[j, 1] <- summary(mymodlin_producers_spa.1)$r.square
  resu[j, 2] <- summary(mymodlin_producers_spa.2)$r.square
  
}
par(mfrow = c(1, 1))
plot(y = resu[, 1], x = 1:100, col = "red", type = "l", ylim = c(.5, max(resu)))
lines(y = resu[, 2], x = 1:100, col = "blue")
which(resu[, 1] == max(resu[, 1]))

par(mfrow = c(1, 2))
dist.coords <- lets.distmat( xy = mydata[,c('Longitude', 'Latitude')])
(correl <- lets.correl(mymodlin_producers_spa.1$residuals, dist.coords, 12, equidistant = T, plot = TRUE))
(correl <- lets.correl(mymodlin_producers_spa.2$residuals, dist.coords, 12, equidistant = T, plot = TRUE))

# best neighborhood structure includes 7 neighbors
kn <- knearneigh(cbind(mydata$Longitude, mydata$Latitude),
                 k = 7, longlat = TRUE)
nbs <- knn2nb(kn, sym = TRUE) 
n <- length(nbs)
for(i in 1:n) {
  # Farming around
  mydata$neighborhoodScore[i] <- mean(mydata$farming[nbs[[i]]])
}


plot(mydata$neighborhoodScore, mydata$neighProd)

plot(mydata$neighborhoodScore)
plot(mydata$neighProd)
###############################


mydata$rich_plot <- scale(mydata$rich_plot)
mydata$rich.nobuf <- scale(mydata$rich.nobuf)
mydata$MammalDiv <- scale(mydata$MammalDiv)
mydata$VascPlant_CoKrig_KreftJetz2007 <- scale(mydata$VascPlant_CoKrig_KreftJetz2007)
mydata$neighborhoodScore <- scale(mydata$neighborhoodScore)

##### VARIANCE COMPONENTS ANALYSES WITH DISPERSAL BUFFER #####
# note we use language as fixed in order to estimate variance components in a strandard fashion
mydata <- na.omit(mydata)
mymod.CurrEO_HistEO_NEIGH_LAN <- lm(farming01 ~ rich_plot + MammalDiv + 
                                      VascPlant_CoKrig_KreftJetz2007 + 
                                      neighborhoodScore +
                                      Language_family, data = mydata)
drop1(mymod.CurrEO_HistEO_NEIGH_LAN)

dist.coords <- as.matrix(lets.distmat(mydata[, c("Longitude", "Latitude")]))
lets.correl(residuals(mymod.CurrEO_HistEO_NEIGH_LAN), dist.coords, 12, equidistant = T, plot = TRUE)

R2.CurrEO_HistEO_NEIGH_LAN <- summary(lm(farming01 ~ rich_plot + 
                                           MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                           neighborhoodScore +
                                           Language_family, 
                                         data = mydata))$r.squared

R2.CurrEO_HistEO_LAN <- summary(lm(farming01 ~ rich_plot + 
                                     MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                     Language_family, 
                                   data = mydata))$r.squared

R2.HistEO_NEIGH_LAN <- summary(lm(farming01 ~ rich_plot + 
                                    neighborhoodScore +
                                    Language_family, 
                                  data = mydata))$r.squared

R2.CurrEO_NEIGH_LAN <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                    neighborhoodScore +
                                    Language_family, 
                                  data = mydata))$r.squared

R2.CurrEO_LAN <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                              Language_family, 
                            data = mydata))$r.squared

R2.HistEO_LAN <- summary(lm(farming01 ~ rich_plot + 
                              Language_family, 
                            data = mydata))$r.squared

R2.NEIGH_LAN <- summary(lm(farming01 ~ neighborhoodScore +
                             Language_family, 
                           data = mydata))$r.squared

R2.LAN <- summary(lm(farming01 ~ Language_family, 
                     data = mydata))$r.squared

R2.CurrEO_HistEO_NEIGH <- summary(lm(farming01 ~ rich_plot + 
                                       MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                       neighborhoodScore, 
                                     data = mydata))$r.squared

R2.CurrEO_HistEO <- summary(lm(farming01 ~ rich_plot + 
                                 MammalDiv + VascPlant_CoKrig_KreftJetz2007, 
                               data = mydata))$r.squared

R2.CurrEO_NEIGH <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                neighborhoodScore, 
                              data = mydata))$r.squared

R2.HistEO_NEIGH <- summary(lm(farming01 ~ rich_plot + 
                                neighborhoodScore, 
                              data = mydata))$r.squared

R2.CurrEO <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007, 
                        data = mydata))$r.squared

R2.HistEO <- summary(lm(farming01 ~ rich_plot , 
                        data = mydata))$r.squared

R2.NEIGH <- summary(lm(farming01 ~ neighborhoodScore, 
                       data = mydata))$r.squared


a <- R2.CurrEO_HistEO_NEIGH_LAN - R2.CurrEO_HistEO_NEIGH # language family
b <- R2.CurrEO_HistEO_NEIGH_LAN - R2.CurrEO_HistEO_LAN # neighbors
c <- R2.CurrEO_HistEO_NEIGH_LAN - R2.HistEO_NEIGH_LAN # Current eco opportunity
d <- R2.CurrEO_HistEO_NEIGH_LAN - R2.CurrEO_NEIGH_LAN # historic Eco Opp

e_f_g_h_i_j_k_l_m_n_o <- R2.CurrEO_HistEO_NEIGH_LAN - a - b - c - d
e_f_g_i_j_k_l_m_n_o <- R2.CurrEO_NEIGH - b - c
e_g_h_i_j_k_l_m_n_o <- R2.HistEO_LAN - a - d
f_g_h_i_j_k_l_m_n_o <- R2.CurrEO_HistEO - c - d
e_f_h_i_j_k_l_m_n_o <- R2.HistEO_NEIGH - b - d
e_f_g_h_j_k_l_m_n_o <- R2.CurrEO_LAN - a - c
e_f_g_h_i_k_l_m_n_o <- R2.CurrEO_NEIGH_LAN - c - b - a

e <- e_f_g_h_i_j_k_l_m_n_o - f_g_h_i_j_k_l_m_n_o
f <- e_f_g_h_i_j_k_l_m_n_o - e_g_h_i_j_k_l_m_n_o
g <- e_f_g_h_i_j_k_l_m_n_o - e_f_h_i_j_k_l_m_n_o
h <- e_f_g_h_i_j_k_l_m_n_o - e_f_g_i_j_k_l_m_n_o
i <- e_f_g_h_i_j_k_l_m_n_o - e_f_g_h_j_k_l_m_n_o
j <- e_f_g_h_i_j_k_l_m_n_o - e_f_g_h_i_k_l_m_n_o

k <- R2.CurrEO_NEIGH_LAN - R2.CurrEO - b - e - a - h
l <- R2.CurrEO_HistEO_LAN - R2.HistEO - f - b - e - a - g
m <- R2.CurrEO_HistEO_LAN - R2.LAN - f - c - j - d - i
n <- R2.CurrEO_HistEO_NEIGH - R2.NEIGH - g - c- j - d - h
o <- e_f_g_h_i_j_k_l_m_n_o - e - f - g - h - i - j - k - l - m - n

p <- 1 - R2.CurrEO_HistEO_NEIGH_LAN 

resu_par <- c(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
show <- round(resu_par, 4) * 100
show[show < 0] <- 0
show <- paste(show, " %", sep = "")


par(mfrow = c(1, 1))
showvarparts(
  4,
  labels = show,
  Xnames = c("Vert", "Hor", "CEO", "HEO"),
  bg = c("darkorange3", "red", "forestgreen",'lightblue')
)
dev.off()


dist.coords <- as.matrix(lets.distmat(mydata[, c("Longitude", "Latitude")]))
(correl <- lets.correl(mymod.CurrEO_HistEO_NEIGH_LAN$residuals, dist.coords, 12, plot = FALSE))
(correl2 <- lets.correl(mydata$farming, dist.coords, 12, plot = FALSE))

tiff(file="Figures/Correlogramresidual.tiff", width = 20,height = 20, units = "cm",res = 600)
.plotcorrel(correl[, 5], correl[, 1], correl[, 2], correl[, 3], 20, pch = 19)
.plotcorrel(correl2[, 5], correl2[, 1], correl2[, 2], correl2[, 3], 20, add = TRUE, 
            col = "gray40", pch = 17)
legend(x = 8000, y = 0.8, c("Observed", "Residuals"), col = c("gray40", "black"),
       pch = c(17, 19), lty = c(3, 3), cex = 1)
dev.off()

##############################################################
##### ANALYSES WITHOUT DISPERSAL BUFFER #####
mymod.CurrEO_HistEO_NEIGH_LAN <- lm(farming01 ~ rich.nobuf + MammalDiv + VascPlant_CoKrig_KreftJetz2007 + neighborhoodScore +
                                      Language_family, data = mydata)
dist.coords <- as.matrix(lets.distmat(mydata[, c("Longitude", "Latitude")]))
lets.correl(residuals(mymod.CurrEO_HistEO_NEIGH_LAN), dist.coords, 12, equidistant = T, plot = TRUE)

R2.CurrEO_HistEO_NEIGH_LAN <- summary(lm(farming01 ~ rich.nobuf + 
                                           MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                           neighborhoodScore +
                                           Language_family, 
                                         data = mydata))$r.squared

R2.CurrEO_HistEO_LAN <- summary(lm(farming01 ~ rich.nobuf + 
                                     MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                     Language_family, 
                                   data = mydata))$r.squared

R2.HistEO_NEIGH_LAN <- summary(lm(farming01 ~ rich.nobuf + 
                                    neighborhoodScore +
                                    Language_family, 
                                  data = mydata))$r.squared

R2.CurrEO_NEIGH_LAN <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                    neighborhoodScore +
                                    Language_family, 
                                  data = mydata))$r.squared

R2.CurrEO_LAN <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                              Language_family, 
                            data = mydata))$r.squared

R2.HistEO_LAN <- summary(lm(farming01 ~ rich.nobuf + 
                              Language_family, 
                            data = mydata))$r.squared

R2.NEIGH_LAN <- summary(lm(farming01 ~ neighborhoodScore +
                             Language_family, 
                           data = mydata))$r.squared

R2.LAN <- summary(lm(farming01 ~ Language_family, 
                     data = mydata))$r.squared

R2.CurrEO_HistEO_NEIGH <- summary(lm(farming01 ~ rich.nobuf + 
                                       MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                       neighborhoodScore, 
                                     data = mydata))$r.squared

R2.CurrEO_HistEO <- summary(lm(farming01 ~ rich.nobuf + 
                                 MammalDiv + VascPlant_CoKrig_KreftJetz2007, 
                               data = mydata))$r.squared

R2.CurrEO_NEIGH <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007 + 
                                neighborhoodScore, 
                              data = mydata))$r.squared

R2.HistEO_NEIGH <- summary(lm(farming01 ~ rich.nobuf + 
                                neighborhoodScore, 
                              data = mydata))$r.squared

R2.CurrEO <- summary(lm(farming01 ~ MammalDiv + VascPlant_CoKrig_KreftJetz2007, 
                        data = mydata))$r.squared

R2.HistEO <- summary(lm(farming01 ~ rich.nobuf , 
                        data = mydata))$r.squared

R2.NEIGH <- summary(lm(farming01 ~ neighborhoodScore, 
                       data = mydata))$r.squared


a <- R2.CurrEO_HistEO_NEIGH_LAN - R2.CurrEO_HistEO_NEIGH # language family
b <- R2.CurrEO_HistEO_NEIGH_LAN - R2.CurrEO_HistEO_LAN # neighbors
c <- R2.CurrEO_HistEO_NEIGH_LAN - R2.HistEO_NEIGH_LAN # Current eco opportunity
d <- R2.CurrEO_HistEO_NEIGH_LAN - R2.CurrEO_NEIGH_LAN # historic Eco Opp

e_f_g_h_i_j_k_l_m_n_o <- R2.CurrEO_HistEO_NEIGH_LAN - a - b - c - d
e_f_g_i_j_k_l_m_n_o <- R2.CurrEO_NEIGH - b - c
e_g_h_i_j_k_l_m_n_o <- R2.HistEO_LAN - a - d
f_g_h_i_j_k_l_m_n_o <- R2.CurrEO_HistEO - c - d
e_f_h_i_j_k_l_m_n_o <- R2.HistEO_NEIGH - b - d
e_f_g_h_j_k_l_m_n_o <- R2.CurrEO_LAN - a - c
e_f_g_h_i_k_l_m_n_o <- R2.CurrEO_NEIGH_LAN - c - b - a

e <- e_f_g_h_i_j_k_l_m_n_o - f_g_h_i_j_k_l_m_n_o
f <- e_f_g_h_i_j_k_l_m_n_o - e_g_h_i_j_k_l_m_n_o
g <- e_f_g_h_i_j_k_l_m_n_o - e_f_h_i_j_k_l_m_n_o
h <- e_f_g_h_i_j_k_l_m_n_o - e_f_g_i_j_k_l_m_n_o
i <- e_f_g_h_i_j_k_l_m_n_o - e_f_g_h_j_k_l_m_n_o
j <- e_f_g_h_i_j_k_l_m_n_o - e_f_g_h_i_k_l_m_n_o

k <- R2.CurrEO_NEIGH_LAN - R2.CurrEO - b - e - a - h
l <- R2.CurrEO_HistEO_LAN - R2.HistEO - f - b - e - a - g
m <- R2.CurrEO_HistEO_LAN - R2.LAN - f - c - j - d - i
n <- R2.CurrEO_HistEO_NEIGH - R2.NEIGH - g - c- j - d - h
o <- e_f_g_h_i_j_k_l_m_n_o - e - f - g - h - i - j - k - l - m - n

p <- 1 - R2.CurrEO_HistEO_NEIGH_LAN 

resu_par <- c(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)
show2 <- round(resu_par, 4) * 100
show2[show2 < 0] <- 0
show2 <- paste(show2, " %", sep = "")


tiff("Figures/Figure3_4vars.tiff", 20, 30, units = "cm", res = 300)
par(mfrow = c(2, 1),
    oma = c(.2, .2, 0.2, 0.2),
    mar=c(1, 1, 1, 1))
showvarparts(
  4,
  labels = show,
  Xnames = c("Vert", "Hor", "CEO", "HEO"),
  bg = c("darkorange3", "red", "forestgreen",'lightblue'),
)
title("(a) With dispersal constrains", line = .3, adj = 0)
#legend("topleft", as.expression(bquote(bold("(a) With dispersal constrains"))),
#       bty = "n")

showvarparts(
  4,
  labels = show2,
  Xnames = c("Vert", "Hor", "CEO", "HEO"),
  bg = c("darkorange3", "red", "forestgreen",'lightblue')
)
title("(b) Without dispersal constrains", line = .3, adj = 0)
dev.off()

## TABLE 1 
# Table 1 model
data_mod <- mydata %>%
  select(
    farming01,
    rich_plot,
    MammalDiv,
    VascPlant_CoKrig_KreftJetz2007,
    neighborhoodScore,
    Language_family
  )
data_mod[,-ncol(data_mod)] <- apply(data_mod[,-ncol(data_mod)], 2, scale)

full_mix <- lmer(
  farming01 ~
    rich_plot +
    MammalDiv +
    VascPlant_CoKrig_KreftJetz2007 +
    neighborhoodScore +
    (1|Language_family),
  data = data_mod
)
suma <- summary(full_mix)
tabela1 <- suma$coefficients
rownames(tabela1) <- c(
  "Intercept",
  "Potential number of domesticated species", 
  "Mammal diversity",
  "Vascular plants diversity",
  "Horizontal transmission (neighborhood effect)")
write.csv(round(tabela1, 3), file = "tabela1.csv")

data_mod <- mydata %>%
  select(
    farming01,
    rich.nobuf,
    MammalDiv,
    VascPlant_CoKrig_KreftJetz2007,
    neighborhoodScore,
    Language_family
  )
data_mod[,-ncol(data_mod)] <- apply(data_mod[,-ncol(data_mod)], 2, scale)

full_mix <- lmer(
  farming01 ~
    rich.nobuf +
    MammalDiv +
    VascPlant_CoKrig_KreftJetz2007 +
    neighborhoodScore +
    (1|Language_family),
  data = data_mod
)
suma <- summary(full_mix)
tabela1 <- suma$coefficients
rownames(tabela1) <- c(
  "Intercept",
  "Potential number of domesticated species", 
  "Mammal diversity",
  "Vascular plants diversity",
  "Horizontal transmission (neighborhood effect)")
write.csv(round(tabela1, 3), file = "tabela1_2.csv")
