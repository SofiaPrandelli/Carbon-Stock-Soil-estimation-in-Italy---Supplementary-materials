###################### CARBON STOCK IN SOIL IN FUNCTION OF ESA LAND USE CLASSES - Italian provinces level
# Carbon in soil is considered divided into principal components Soil organic carbon, above ground carbon and below ground carbon 

# Map from UNEP of Above and below ground carbon biomass storage 
# SOC estimates starting from SOC values 0-30 cm depth in Emilia-Romagna --> prediction of SOC in italy with RF 

######## setup ##########
setwd("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4")

#setwd("C:/Users/Sofia/PhD/Analisi/c_stock/provinces_level/") 
#setwd("~/Desktop/analisi/c_stock/")
library(geodata)
library(osmdata)
library(osmextract)
library(terra)
library(sf)
library(exactextractr)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(tidyr)
library(rgbif)
library(maps)
library(tidyverse)
library(readxl)
library(openxlsx)
library(TNRS)
library(readr)
library(geodata)
library(randomForest)
library(viridis)
library(ranger)
library(patchwork)

#### parallelizze ####
library(future) # oppure library(doParallel) o library(future.apply)
library(furrr)
# parallel computing 
plan(multisession, workers=4) # per non occupare tutti i core su server, noi abbiamo 16 core

# per disattivare il parallelismo:
# plan(sequential)  # esecuzione su un solo core

######## WCMC raster #########
# Filtering the dataset in gee: 
wcmc_rast <- rast("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_italy.tif")
# wcmc_df <- as.data.frame(wcmc_rast, xy = TRUE) #convert raster to spat points df
# coordinates(wcmc_df) <- c("x", "y") #convert spatial points df to spat object
# wcmc_sdf <- st_as_sf(wcmc_df) #convert to sf object
wcmc_rast

wcmc_df <- as.data.frame(wcmc_rast, xy = TRUE, na.rm = TRUE)

ggplot()+
  geom_tile(data=wcmc_df, aes(x=x, y=y, fill=carbon_tonnes_per_ha)) +
  scale_fill_viridis_c(name = "C (T/Ha)", na.value = "transparent") +
  coord_fixed() +
  labs(title = "Above- and Below-ground Carbon Storage in the Italian Territory",
       subtitle="Tonnes of Carbon storage per Hectare")+
  theme_minimal()
ggsave("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/WCMC_italy.png", width=12, height=7)

######## Italian regions and borders #########
# https://www.paulamoraga.com/tutorial-open-spatial-data/#:~:text=For%20example%2C%20the%20worldclim_country(),codes%20of%20the%20world%20countries
#borders <- ne_countries(type = "countries", country = "Italy", scale = "medium", returnclass = "sf")
regions <- rnaturalearth::ne_states("Italy", returnclass = "sf")

er <- regions %>% filter(region=="Emilia-Romagna")
ext_er <- ext(er)

####### ESA LC 2009 300 m res ###### 
# GlobCover 2009 is a global land cover map based on ENVISAT's Medium Resolution Imaging Spectrometer (MERIS) Level 1B data acquired in full resolution mode with a spatial resolution of approximately 300 meters.
# cropped from gee, resolution of 300 m. We'll create a stack with wcmc_rast in order to obtain backwards the landcover types associated with the carbon stock value in a given pixel
esalc <- rast("/media/r_projects/sofia.prandelli2/esa_LC_300m.tif")
res(esalc)
crs(esalc)
freq(esalc)

####### mask to Italy ##########
esalc <- mask(crop(esalc, regions), regions)

esalc_df <- as.data.frame(esalc, xy=T) 
esalc_df$landcover <- as.numeric(esalc_df$landcover)
length(unique(esalc_df$landcover)) # 18 out of 22

####### mask to ER region #######
esalc_er <- mask(crop(esalc, er), er)
freq(esalc_er)
esalc_er_df <- as.data.frame(esalc_er, xy=T) 
esalc_er_df$landcover <- as.numeric(esalc_er_df$landcover)
length(unique(esalc_er_df$landcover)) # 15 out of 22

####### DEM ######
# Scarica DEM per l'Italia con ris 7.5 arc-seconds (circa 250m)
elev <- elevation_30s(country = "Italy", path = tempdir())  # oppure usa "0.002694946" gradi (~250m)
crs(elev)
res(elev)
# Downscaling a circa 300 m (0.002°) per allinearlo agli altri raster
elev_300 <- resample(elev, esalc, method = "bilinear") # per variabili continue come dati bioclimatici o dem, method = "bilinear" o "cubic" è appropriato
res(elev_300)

######## WORLDCLIM VARIABLES #########
# Tramite Geodata package, res 1 km (0.08333333 minutes of a degree)
worldclim_bio_ita <- worldclim_country("Italy", res=0.002694946, var="bio", path=tempdir(), version="2.1")
res(worldclim_bio_ita)
# downscaling
worldclim_bio_ita <- resample(worldclim_bio_ita, esalc, method = "bilinear")

bio_var_df <- as.data.frame(worldclim_bio_ita, xy=T)

# considero solo variabili più utilizzate e significative nella modellazione ecologica e pedologica, che rappresentano bene condizioni climatiche medie + stagionalità
# evita collinearità tra variabili
bio_selected <- subset(worldclim_bio_ita, c("wc2.1_30s_bio_1", "wc2.1_30s_bio_4", "wc2.1_30s_bio_10", "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_15"))

names(bio_selected) <- c("AnnualMeanTemperature", "TemperatureSeasonality", "MeanTemperatureofWarmestQuarter",
                         "MeanTemperatureofColdestQuarter", "AnnualPrecipitation", "PrecipitationSeasonality")

######## Info worldclim
ext(bio_selected) 
if (crs(bio_selected) == crs(clc_rast)) {
  print("same CRS")} else {print("different CRS")} 

# rast dimensions
nrows <- nrow(worldclim_bio_ita)
ncols <- ncol(worldclim_bio_ita)

####### SOC emilia-romagna #####
# Valori stimati di Carbonio Organico stoccato nel suolo da 0 a 30 cm 
# https://ambiente.regione.emilia-romagna.it/it/geologia/suoli/proprieta-e-qualita-dei-suoli/carbonio-organico-immagazzinato-nei-suoli
# per scaricamento: https://datacatalog.regione.emilia-romagna.it/catalogCTA/dataset/r_emiro_2023-08-09t162508 
# unità di misura: Mg/Ha cioè tonnes/Ha
soc_er <- rast("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/C_org_stock_0_30_cm_EmiliaRomagna/SOC_STOCK_0_30.tif")
soc_er <- project(soc_er, "epsg:4326")
crs(soc_er)
res(soc_er) # initial resolution of soc_er = 0.001x0.001 --> we resampled to.. 100m? NO! usiamo la res degli altri raster, circa 250/300 m

esalc_er <- mask(crop(esalc, er_vect), er_vect)
soc_er <- resample(soc_er, esalc_er, method="bilinear") # weighted average of 4 pixels in the original image nearest to the new pixel location

soc_er_df <- as.data.frame(soc_er, xy=T,na.rm=T)

names(soc_er_df)[3] <- "SOC"
ggplot()+
  geom_raster(data=soc_er_df, aes(x=x, y=y, fill=SOC)) + 
  scale_fill_viridis_c(name = "SOC (Mg/Ha)",
                       limits = c(0, 500),
                       breaks = seq(0, 500, by = 100),
                       na.value = "transparent") +
  coord_fixed() +
  labs(title = "Estimation of Soil Organic Carbon Stock. Resampled map (from Emilia-Romagna Region Map of SOC, 2023)",
       subtitle="Tonnes of SOC per Ha at 0-30 cm depth")+
  theme_minimal()
ggsave("SOC_stock_ER_resampled.png", width=12, height=7)

######## Resample for rasters
soc_er_cut <- resample(soc_er_cut, esalc_er)
elev_cut <- resample(elev_cut, esalc_er)
bio_cut  <- resample(bio_cut, esalc_er)

ext(bio_cut)
crs()
res()

##### Stack raster in er: final raster "stack_train" ####
stack_train <- c(soc_er_cut, esalc_er, elev_cut, bio_cut)
names(stack_train) <- c("SOC", "ESAlc", "Elevation", "AnnualMeanTemperature", "TemperatureSeasonality", "MeanTemperatureofWarmestQuarter", "MeanTemperatureofColdestQuarter", "AnnualPrecipitation", "PrecipitationSeasonality") # bioclim: 1 4 10 11 12 15

freq(stack_train$ESAlc)
plot(stack_train)

###### Save and import final stack stack_train (solo su er) #####
#writeRaster(stack_train, "/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/stack_train_ESA.tif", overwrite=T)
#stack_train <- rast("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/stack_train_ESA.tif")

#stack_train_300m_res <- rast("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/stack_train_300m_res.tif") # QUA CON CORINE LC MA SENZA ESA LC
# mentre il tif "/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/stack_train.tif" ha la risoluzione iniziale (di 0.008°)

######### aggiungi a stack_train_300m_res il raster ESAlc
#esalc_er <- resample(esalc, stack_train_300m_res)
#stack_train_300m_res <- c(stack_train_300m_res, esalc_er)
#names(stack_train_300m_res)[10] <- "ESAlc"
#writeRaster(stack_train_300m_res, "/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/stack_train_300m_res.tif", overwrite=T)
#stack_train_300m_res <- stack_train_300m_res[[-2]] # tolgo corine lc

#######################
###### Predizione valori SOC tramite Random Forest ######
# robusto, non lineare, ideale per dati spaziali rumorosi
library(randomForest)

#### TRAINING Data ####
# rimuovi classi 190,200,210,220
# escludiamo dal training classi di lc dove SOC non è definito o non è significativo (acque, ghiacci, superfici completamente impermeabili, aree “bare” senza suolo)
remove <- c(190, 200, 210, 220)
# Converti in df
df_emr <- as.data.frame(stack_train, xy=F, na.rm=T) %>%
  filter(!is.na(ESAlc), !is.na(SOC), !is.na(Elevation)) %>% 
  filter(!(ESAlc %in% remove)) %>%
  mutate(ESAlc = droplevels(as.factor(ESAlc)))

# Fattorizza ESA
df_emr$"ESAlc" <- as.factor(df_emr$"ESAlc")

#### Model RF res 0.002° #####
# per parallelizzazione
#install.packages("ranger")

set.seed(42)
rf_model_4_permutation_ESAlc <- ranger::ranger(
  formula = SOC ~ .,
  data = df_emr,
  num.trees = 500,
  importance = "permutation",  # "impurity" o "permutation" se serve un confronto accurato dell'importanza tra variabili, ma permutation è più lenta
  # quindi --> Alleniamo il modello con importance = "impurity" per testare e ottimizzare
  # Quando hai il modello finale, rifacciamo il training con importance = "permutation" per una stima più accurata delle variabili importanti
  keep.inbag = TRUE,
  num.threads = 4  # quanti core usare
)

#### save rf model ######
saveRDS(rf_model_4_permutation_ESAlc, file = "rf_model_4_permutation_ESAlc.rds")
rf_model_4_permutation_ESAlc <- readRDS("rf_model_4_permutation_ESAlc.rds")

##### Model Results ######
# https://scikit-learn.org/stable/modules/permutation_importance.html 
print(rf_model_4_permutation_ESAlc)
# Type:                             Regression 
# Number of trees:                  500 
# Sample size:                      319826 
# Number of independent variables:  8 
# Mtry:                             2 
# Target node size:                 5 
# Variable importance mode:         permutation 
# Splitrule:                        variance 
# OOB prediction error (MSE):       61.83894 
# R squared (OOB):                  0.9143552 

print(importance(rf_model_4_permutation_ESAlc))
# ESAlc                       Elevation           AnnualMeanTemperature 
# 53.16251                       632.21984                       356.03759 
# TemperatureSeasonality MeanTemperatureofWarmestQuarter MeanTemperatureofColdestQuarter 
# 355.30404                       416.07269                       323.07161 
# AnnualPrecipitation        PrecipitationSeasonality 
# 315.66762                       227.99046 
var_imp_perm <- importance(rf_model_4_permutation_ESAlc)

var_imp_perm_df <- data.frame(
  Variable = names(var_imp_perm),
  Importance = as.numeric(var_imp_perm)
)

var_imp_perm_df$Percentage <- 100 * var_imp_perm_df$Importance / sum(var_imp_perm_df$Importance)
var_imp_perm_df <- var_imp_perm_df[order(var_imp_perm_df$Percentage, decreasing=T), ] # Ordina per impo

ggplot(var_imp_perm_df, aes(x = reorder(Variable, Percentage), y = Percentage, fill = Percentage)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis(option = "D", direction = -1) +
  labs(title = "Variable importance (Impurity)",
       x = "Variabile",
       y = "Percentage") +
  theme_minimal()
ggsave("var_imp_mod4_perm_ESAlc.png", width=12, height=7)

##### Visualizza relazione tra SOC in ER e tutte le variabili considerate - MODELLO 4 ######
# Separa variabili numeriche da variabile fattore(ESA)
#### Plot 1: var numeriche ####
num_vars <- c("Elevation", "AnnualMeanTemperature", "TemperatureSeasonality", "MeanTemperatureofWarmestQuarter",
              "MeanTemperatureofColdestQuarter", "AnnualPrecipitation", "PrecipitationSeasonality")

# Plot 1: Curva di SOC rispetto alle variabili numeriche
df_num <- df_emr %>%
  select(SOC, all_of(num_vars)) %>%
  pivot_longer(cols = all_of(num_vars), names_to = "variable", values_to = "value")

ggplot(df_num, aes(x = value, y = SOC)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  facet_wrap(~variable, scales = "free_x", ncol = 2) +
  theme_minimal() +
  labs(title = "Relazione tra SOC e variabili climatiche/topografiche",
       x = "Valore predittore",
       y = "SOC (Mg/ha)") +
  theme(strip.text = element_text(face = "bold"))
ggsave("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_numeric_variables_modelESA.png", width=12, height=7)

#### Plot 2: SOC + CLC (categorico) ####
#Calcola la mediana SOC per ciascuna classe
soc_medians <- df_emr %>%
  group_by(ESAlc) %>%
  summarise(median_SOC = median(SOC, na.rm = TRUE))
# Aggiungi valori mediani al df
df_emr <- df_emr %>%
  left_join(soc_medians, by = "ESAlc")

######## rimuovi Water bodies (class 210) and Permanent snow and ice (classes 220), Artificial surfaces and associated areas (class 190), Bare areas (class 200)
#df_emr <- subset(df_emr, !(ESAlc %in% c(190,200,210,220)))

ggplot(df_emr, aes(x =ESAlc, y = SOC, fill = median_SOC)) +
  #geom_violin(trim = FALSE, alpha = 0.6, fill="skyblue") +
  #geom_jitter(width = 0.1, alpha=0.1) +
  geom_boxplot(outlier.size = 0.1) + #, fill = "skyblue") +
  scale_y_continuous(n.breaks = 6)+
  theme_minimal() +
  scale_fill_viridis_c(option = "D", name = "SOC values", direction = 1) +
  labs(title = "SOC values distributions per ESA Land Cover in Emilia-Romagna Region",
       x = "ESA Land Cover",
       y = "SOC (Mg/Ha)")
ggsave("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_ESAlc.png", width=12, height=7)



###### Predizione su tutta italia ####
#remove <- c(190, 200, 210, 220)
train_lvls <- sort(unique(as.integer(as.character(levels(df_emr$ESAlc)))))  # livelli visti nel training
## stack nazionale
stack_national <- c(esalc, elev_300, bio_selected)

## maschera: NA ovunque esalc è da escludere o non visto nel training
mask <- (stack_national$ESAlc %in% remove) | !(stack_national$ESAlc %in% train_lvls)
# applica maschera a tutti i layer (così i NA sono sincronizzati)
for (i in 1:nlyr(stack_national)) {
  stack_national[[i]] <- mask(stack_national[[i]], mask, maskvalues = 1)
}

predictors_df <- as.data.frame(stack_national, na.rm = F)

## fattorizza ESAlc con i livelli del training
predictors_df$ESAlc <- factor(predictors_df$ESAlc, levels = train_lvls)

pred_cols <- c("ESAlc","Elevation",
               "AnnualMeanTemperature","TemperatureSeasonality",
               "MeanTemperatureofWarmestQuarter","MeanTemperatureofColdestQuarter",
               "AnnualPrecipitation","PrecipitationSeasonality")

## predici solo su righe complete (nessun NA in nessuna colonna)
ok <- stats::complete.cases(predictors_df[, pred_cols])

SOC_vec <- rep(NA_real_, nrow(predictors_df))
if (any(ok)) {
  SOC_vec[ok] <- predict(rf_model_4_permutation_ESAlc,
                         data = predictors_df[ok, pred_cols])$predictions
}
## ricomponi raster di output
SOC_pred_rast <- rast(stack_national[[1]])
values(SOC_pred_rast) <- SOC_vec

writeRaster(SOC_pred_rast, "/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_pred_ita_model4_ESAlc_new.tif")
SOC_pred_rast <- rast("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_pred_ita_model4_ESAlc_new.tif")


###### Prediction on Italy
predictors_df <- as.data.frame(stack_national, na.rm=F)
SOC_pred <- predict(rf_model_4_permutation_ESAlc, data = predictors_df)$predictions

# Crea nuovo raster con predizioni
SOC_rast <- rast(stack_national[[1]])
values(SOC_rast) <- NA  # inizializza
# inserisce predizioni dove c’erano dati validi
valid_cells <- !is.na(values(stack_national[[1]]))
values(SOC_rast)[valid_cells] <- SOC_pred
SOC_rast <- predict(stack_national, rf_model_4_permutation_ESAlc, na.rm = TRUE)

# writeRaster(SOC_rast, "/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_pred_ita_model4_ESAlc.tif") #, overwrite = TRUE)
SOC_rast <- rast("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_pred_ita_model4_ESAlc.tif")

####### plot predizione ######
SOC_rast_df <- as.data.frame(SOC_pred_rast, xy = TRUE, na.rm = TRUE)
names(SOC_rast_df)[3] <- "SOC"

ggplot()+
  geom_raster(data=SOC_rast_df, aes(x=x, y=y, fill=SOC)) +
  scale_fill_viridis_c(name = "SOC (Mg/Ha)", na.value = "transparent") +
  coord_fixed() +
  labs(title = "Estimation of Soil Organic Carbon in the Italian Territory",
       subtitle="Tonnes of SOC per Ha at 0-30 cm depth")+
  theme_minimal()
ggsave("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_pred_ita_model4_ESAlc2.png", width=12, height=7)

######## plot solo su er per confrontare raster originale con raster predetto #####
SOC_rast_er <- crop(SOC_pred_rast, vect(er)) |> mask(vect(er))
SOC_rast_er_df <- as.data.frame(SOC_rast_er, xy=T, na.rm=T)
names(SOC_rast_er_df)[3] <- "SOC"
p1 <- ggplot()+
  geom_raster(data=soc_er_df, aes(x=x, y=y, fill=SOC)) + 
  scale_fill_viridis_c(name = "SOC (Mg/Ha)",
                       limits = c(0, 500),
                       breaks = seq(0, 500, by = 100),
                       na.value = "transparent") +
  coord_fixed() +
  labs(title = "A) ")+
  theme_minimal()

p2 <- ggplot() +
  geom_raster(data = SOC_rast_er_df, aes(x = x, y = y, fill = SOC)) +
  scale_fill_viridis_c(name = "SOC (Mg/Ha)",
                        limits = c(0, 500),
                        breaks = seq(0, 500, by = 100),
                        na.value = "transparent") +
#  coord_sf(crs = st_crs(er)) +
  coord_fixed() +
  labs(title = "B) ") +
  theme_minimal()

p1 + p2 +
  # theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
  #       axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  plot_annotation(title = "Resampled map of SOC Stock from Emilia-Romagna Region, 2023 vs Predicted map with RF model",
                  subtitle="Tonnes of SOC per Ha at 0-30 cm depth") + 
  plot_layout(guides = "collect") 

ggsave("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/SOC_pred_ER_model4_confronto.png", width = 14, height = 7)

##### Valutazione quantitativa del modello ######

# df con x,y,SOC dai raster di training ER
df_emr <- as.data.frame(stack_train, xy = T, na.rm = T)  # stack_train = ER (SOC + predittori)
# crea punti dalle coordinate del df (CRS di stack_train)
pts_emr <- vect(df_emr, geom = c("x","y"), crs = crs(stack_train))

# estrai predizioni al centro cella dei punti ER
pred_vec <- terra::extract(SOC_pred_rast, pts_emr)[, 2]  # colonna 2 = valori del raster estratto
# osservazioni dal df (SOC target usato nel training)
obs_vec <- df_emr$SOC
# allinea ed elimina NA
valid <- !is.na(obs_vec) & !is.na(pred_vec)
obs <- obs_vec[valid]
pred <- pred_vec[valid]

## metriche quantitative
R2    <- cor(obs, pred)^2
RMSE  <- sqrt(mean((pred - obs)^2))
MAE   <- mean(abs(pred - obs))
BIAS  <- mean(pred - obs)
cat(sprintf("R² = %.3f\nRMSE = %.3f\nMAE = %.3f\nBias = %.3f\n", R2, RMSE, MAE, BIAS))
# R² = 0.983
# RMSE = 3.589
# MAE = 2.375
# Bias = -0.010

# Scatterplot osservato vs predetto
df_val <- data.frame(Observed = obs, Predicted = pred)
lm_fit <- lm(Predicted ~ Observed, data = df_val)
r2_val <- summary(lm_fit)$r.squared
intercept_val <- coef(lm_fit)[1]
slope_val <- coef(lm_fit)[2]

a <- ggplot(df_val, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "darkred", linetype = "dashed") +  # 1:1 line
  geom_smooth(method = "lm", color = "darkblue", se = F) + # Regression line
  annotate("text", x = min(df_val$Observed), y = max(df_val$Predicted), 
           label = paste0("R² = ", round(r2_val, 3), 
                          "\nIntercept = ", round(intercept_val, 2), 
                          "\nSlope = ", round(slope_val, 2)),
           hjust = 0, vjust = 1, size = 4) +
  labs(#title = "Observed vs Predicted SOC (in Emilia-Romagna Region)",
       x = "Observed SOC",
       y = "Predicted SOC") +
  #theme_minimal()
  theme_bw()

b <- ggplot(data=df_val, aes(x=Observed, y=Predicted) ) +
  geom_hex(bins = 80) +
  ylab("Predicted SOC") +
  xlab("Observed SOC") +
  #ggtitle("Observed vs Predicted SOC (in Emilia-Romagna Region)") +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

a + b +
  plot_annotation(title = "Observed vs Predicted SOC (in Emilia-Romagna Region)") 
ggsave("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/scatterplot_obs_pred_hexbins.png", width=16, height=7)

##########################
###### MODELLO LINEARE WCMC STOCK PER VALUTARE IN FN DI VARIABILI BIOCLIM, ALTITUD, ESALC ######
library(biglm)
# 1. Stack dei raster (assicurati che abbiano stesso extent e risoluzione)
# wcmc_rast: carbon stock
# esalc: land cover classes
# elev_300: elevation raster
# bio_selected: stack con le variabili bioclimatiche (es bio1, bio12)
stack_wcmc <- c(wcmc_rast, esalc, elev_300, worldclim_bio_ita)

names(stack_wcmc) <- c("C_stock", "LandCover", "Elevation",
                       "bio1", "bio2","bio3","bio4","bio5","bio6",
                       "bio7","bio8","bio9","bio10", "bio11", "bio12", 
                       "bio13","bio14","bio15","bio16","bio17","bio18",
                       "bio19")
df_full <- as.data.frame(stack_wcmc, na.rm = TRUE)
print(dim(df_full))
# 3. Imposta LandCover come fattore
df_full$LandCover <- as.factor(df_full$LandCover)


model_full <- lm(C_stock ~ ., data = df_full)
summary(model_full)
coef(model_full)
# analisi significatività di variabili
anova(model_full)


##### stack tra WCMC ABGCB ed esalc ######
SOC_pred_rast
wcmc_rast <- mask(crop(wcmc_rast, regions), regions)

stack <- c(SOC_pred_rast, wcmc_rast, esalc)
names(stack) <- c("SOC", "ABGC", "esalc")

## maschera: NA ovunque esalc è da escludere o non visto nel training
mask <- (stack$esalc %in% remove) | !(stack$esalc %in% train_lvls)
# applica maschera a tutti i layer
for (i in 1:nlyr(stack)) {
  stack[[i]] <- mask(stack[[i]], mask, maskvalues = 1)
}
stack_df <- as.data.frame(stack,  xy=T, na.rm = F)


######## Extents at provinces level ########
# Creating extents for every province
# Create a list to store SpatExtent objects for each province
province_extent_list <- list()
# Loop through each province
for (i in 1:nrow(regions)) {
  province <- regions[i, ] # selecting rows of regions df
  province_geometry <- st_geometry(province) #extract geometry
  province_bbox <- st_bbox(province_geometry) #return bounding of sf
  province_extent <- ext(province_bbox) #create SpatExtent object for the province
  province_extent_list[[i]] <- province_extent
}

# Assign region names to every SpatRaster object
for (i in 1:nrow(regions)) {
  names(province_extent_list)[i] <- regions$name[i]
}

######## Create a list of cropped rasters_stack: cropped_raster_list #######
cropped_raster_list <- list()
# Loop through each SpatExtent object
for (i in 1:length(province_extent_list)) {
  cropped_raster <- crop(stack, province_extent_list[[i]]) 
  cropped_raster_list[[i]] <- cropped_raster
}

# assign name regions to the list
for (i in 1:nrow(regions)) {
  names(cropped_raster_list)[i] <- regions$name[i]
}

### remove null values in cropped_raster_list
cropped_raster_list <- lapply(cropped_raster_list, function(r) {
  df <- as.data.frame(r, xy=F)
  df <- df[!is.na(df$SOC) & !is.na(df$ABGC),]
  return(df)
})

sum(is.na(cropped_raster_list$Bari$ABGC))
view(cropped_raster_list$Bari)


##### Plot boxplots for esalc for each province ####
out_dir <- "/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/plots_c_stock_per_prov/boxplot_provinces"

# unisci i df in uno solo, aggiungendo colonna Province
df_long <- bind_rows(
  lapply(seq_along(cropped_raster_list), function(i) {
    df <- cropped_raster_list[[i]]
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
    prov_name <- names(cropped_raster_list)[i]
    if (is.null(prov_name) || is.na(prov_name) || prov_name == "") {
      prov_name <- paste0("prov_", i)
    }
    # tieni solo colonne necessarie e aggiungi Province
    needed <- c("SOC", "ABGC", "esalc")
    df %>%
      transmute(
        Province = prov_name,
        esalc    = esalc,
        SOC      = SOC,
        ABGC     = ABGC
      ) %>%
      filter(!is.na(esalc), (!is.na(SOC) | !is.na(ABGC)))
  })
)

# Livelli coerenti di esalc per tutte le province
esa_levels <- sort(unique(df_long$esalc))
df_long$esalc <- factor(df_long$esalc, levels = esa_levels)

# Long format SOC/ABGC
df_long_all <- df_long %>%
  pivot_longer(cols = c(SOC, ABGC), names_to = "Variable", values_to = "Value") %>%
  filter(!is.na(Value))

# Split per provincia
by_prov <- split(df_long_all, df_long_all$Province)

for (prov in names(by_prov)) {
  d <- by_prov[[prov]]
  # se vuoi ordinare l’asse X per livello numerico di esalc:
  d$esalc <- factor(d$esalc, levels = esa_levels)
  p <- ggplot(d, aes(x = esalc, y = Value, fill = Variable)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8,
                 position = position_dodge(width = 0.75)) +
    labs( title = paste("SOC & ABGC per ESA lc —", prov),
      x = "Esa lc",
      y = "Value",
      fill = "Variabile" ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  outfile <- file.path(out_dir, paste0("boxplot_", str_replace_all(prov, "[^A-Za-z0-9_-]", "_"), ".png"))
  ggsave(outfile, plot = p, width = 12, height = 8, dpi = 300)
}



######## Create a list of df with Mean C values for Corine LC for provinces ########
# Aggiunta di errore standard (meglio rispetto a deviazione standard)
mean_c_values_esalc <- list()

sem <- function(x) {
  sd(x, na.rm=T)/sqrt(length(x[!is.na(x)]))
}

for (i in seq_along(cropped_raster_list)) {
  df <- cropped_raster_list[[i]]
  mean_ABGC <- tapply(df$ABGC, df$esalc, mean, na.rm = T) #calcola media di abgc carbon_tonnes_per_ha
  sem_ABGC <- tapply(df$ABGC, df$esalc, sem)
  mean_SOC <- tapply(df$SOC, df$esalc, mean, na.rm = T) #calcola media di soc
  sem_SOC <- tapply(df$SOC, df$esalc, sem)
  result <- data.frame(
    esalc = as.numeric(names(mean_ABGC)),
    mean_ABGC = mean_ABGC,
    sem_ABGC = sem_ABGC,
    mean_SOC = mean_SOC,
    sem_SOC = sem_SOC
  )
  mean_c_values_esalc[[i]] <- result
}

# assign provinces names
for (i in 1:nrow(regions)) {
  names(mean_c_values_esalc)[i] <- regions$name[i]
}

view(mean_c_values_esalc$Roma)

##### aggiungi classi con valori nulli ####
add_missing_esalc <- function(df, classes_to_add = c(190, 200)) {
  # Normalizza ESAlc a numerico (senza perdere info)
  esalc_num <- df$esalc
  if (is.factor(esalc_num)) esalc_num <- as.numeric(as.character(esalc_num))
  if (is.character(esalc_num)) esalc_num <- as.numeric(esalc_num)
  # Classi mancanti
  missing <- setdiff(classes_to_add, unique(esalc_num))
  if (length(missing) == 0) return(df)  # niente da aggiungere
  # Costruisci righe "scheletro" con stesse colonne del df, tutte NA
  skel <- df[0, , drop = FALSE]
  add  <- skel[rep(1, length(missing)), , drop = FALSE]
  for (nm in names(add)) add[[nm]] <- NA
  # Imposta esalc nelle righe aggiunte, rispettando il tipo originale
  if (is.factor(df$esalc)) {
    levs <- union(levels(df$esalc), as.character(classes_to_add))
    add$esalc <- factor(as.character(missing), levels = levs)
    df$esalc  <- factor(as.character(df$esalc), levels = levs)
  } else if (is.character(df$esalc)) {
    add$esalc <- as.character(missing)
  } else {
    add$esalc <- as.numeric(missing)
  }
  # Unisci
  out <- bind_rows(df, add)
}

mean_c_values_esalc <- lapply(mean_c_values_esalc, add_missing_esalc)

######## Add description of ESA LC classes (23) ####
lc_descriptions <- data.frame(
  esalc=c(11,14,20,30,40,50,60,70,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230),
  description=c(
    "Post-flooding or irrigated croplands",
    "Rainfed croplands",
    "Mosaic cropland (50-70%)/ vegetation (grassland, shrubland, forest) (20-50%)",
    "Mosaic vegetation (grassland, shrubland, forest) (50-70%) / cropland (20-50%)",
    "Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m)",
    "Closed (>40%) broadleaved deciduous forest (>5m)",
    "Open (15-40%) broadleaved deciduous forest (>5m)",
    "Closed (>40%) needleleaved evergreen forest (>5m)",
    "Open (15-40%) needleleaved deciduous or evergreen forest (>5m)",
    "Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)",
    "Mosaic forest-shrubland (50-70%) / grassland (20-50%)",
    "Mosaic grassland (50-70%) / forest-shrubland (20-50%)",
    "Closed to open (>15%) shrubland (<5m)",
    "Closed to open (>15%) grassland",
    "Sparse (>15%) vegetation (woody vegetation, shrubs, grassland)",
    "Closed (>40%) broadleaved forest regularly flooded - Fresh water",
    "Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded - saline water",
    "Closed to open (>15%) vegetation (grassland, shrubland, woody vegetation) on regularly flooded or waterlogged soil - fresh, brackish or saline water",
    "Artificial surfaces and associated areas (urban areas >50%)",
    "Bare areas",
    "Water bodies",
    "Permanent snow and ice",
    "Unclassified"
  ),
  category=c(
    "cropland",
    "cropland",
    "cropland",
    "mosaic vegetation",
    "closed forest",
    "closed forest",
    "open forest",
    "closed forest",
    "open forest",
    "closed to open forest",
    "mosaic vegetation",
    "mosaic vegetation",
    "shrubland",
    "grassland",
    "sparse vegetation",
    "closed forest",
    "closed forest",
    "closed to open vegetation",
    "artificial surfaces",
    "bare areas",
    "water bodies",
    "permanent snow and ice",
    "unclassified"
  )
)

# fn to add descriptions to df
add_descriptions <- function(df) {
  df$description <- lc_descriptions$description[match(df$esalc, lc_descriptions$esalc)]
  df$category <- lc_descriptions$category[match(df$esalc, lc_descriptions$esalc)]
  return(df)
}

# loop through each df in c_values_corine list and add descriptions
for (i in 1:length(mean_c_values_esalc)) {
  mean_c_values_esalc[[i]] <- add_descriptions(mean_c_values_esalc[[i]])
}

view(mean_c_values_esalc$Bologna)

##### SAVE DF as excel file ####
library("openxlsx")
write.xlsx(mean_c_values_esalc, file="mean_c_values_esalc.xlsx")

######## PLOT c_stock in function of lc classes for provinces ##########
setwd("/media/r_projects/sofia.prandelli2/wcmc_rail/WCMC_SOC/SOC_MODEL_4/plots_c_stock_per_prov")

category_colors <- c(
  "Post-flooding or irrigated croplands"="#aaefef",
  "Rainfed croplands"="#ffff63",
  "Mosaic cropland (50-70%)/ vegetation (grassland, shrubland, forest) (20-50%)"="#dcef63",
  "Mosaic vegetation (grassland, shrubland, forest) (50-70%) / cropland (20-50%)"="#cdcd64",
  "Closed to open (>15%) broadleaved evergreen and/or semi-deciduous forest (>5m)"="#006300",
  "Closed (>40%) broadleaved deciduous forest (>5m)"="#009f00",
  "Open (15-40%) broadleaved deciduous forest (>5m)"="#aac700",
  "Closed (>40%) needleleaved evergreen forest (>5m)"="#003b00",
  "Open (15-40%) needleleaved deciduous or evergreen forest (>5m)"="#286300",
  "Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)"="#788300",
  "Mosaic forest-shrubland (50-70%) / grassland (20-50%)"="#8d9f00",
  "Mosaic grassland (50-70%) / forest-shrubland (20-50%)"="#bd9500",
  "Closed to open (>15%) shrubland (<5m)"="#956300",
  "Closed to open (>15%) grassland"="#ffb431",
  "Sparse (>15%) vegetation (woody vegetation, shrubs, grassland)"="#ffebae",
  "Closed (>40%) broadleaved forest regularly flooded - Fresh water"="#00785a",
  "Closed (>40%) broadleaved semi-deciduous and/or evergreen forest regularly flooded - saline water"="#009578",
  "Closed to open (>15%) vegetation (grassland, shrubland, woody vegetation) on regularly flooded or waterlogged soil - fresh, brackish or saline water"="#00dc83",
  "Artificial surfaces and associated areas (urban areas >50%)"="#c31300",
  "Bare areas"="#fff5d6",
  "Water bodies"="#0046c7",
  "Permanent snow and ice"="lightblue2",
  "Unclassified"="darkgrey"
)

# category_colors <- c(
# "cropland"="#ffff63",
# "mosaic vegetation"="#cdcd64",
# "closed forest"="#003b00",
# "open forest",
# "closed to open forest",
# "mosaic vegetation",
# "shrubland",
# "grassland",
# "sparse vegetation",
# "closed to open vegetation",
# "artificial surfaces",
# "bare areas",
# "water bodies",
# "permanent snow and ice",
# "unclassified")

legend_data <- data.frame(Category = names(category_colors), Color = category_colors, stringsAsFactors = FALSE)
legend_data$Category <- factor(legend_data$Category, levels = names(category_colors))

############# LOOP FOR ALL PROVINCES
# print max values
library(purrr)
highest_values <- map_dbl(mean_c_values_esalc, ~ max(.$mean_SOC))
highest_values <- as.data.frame(highest_values)
max(highest_values, na.rm=T) # 112.6569 C t/ha per mean_ABGC
                             # 127.2578  C t/ha per SOC

#library(patchwork) # per combinare più plot in unico file

for (i in seq_along(mean_c_values_esalc)) {
  province_name <- names(mean_c_values_esalc)[i]
  df_prov <- mean_c_values_esalc[[i]]
  # Spezza la descrizione su due righe (alla prima occorrenza di spazio)
  # df_prov$label_text <- str_replace(df_prov$description, " ", "\n")
  # ABGC Plot
  plot_abgc <- ggplot(df_prov, aes(x = description, y = mean_ABGC, fill = description)) +
    geom_col(alpha=0.7) +
    geom_text(aes(label = category),
              angle = 60, 
              vjust=-0.5,
              hjust = 0, 
              size = 2.5, 
              color = "black") +
    scale_fill_manual(values = category_colors) +
    scale_x_discrete(labels = df_prov$esalc) +
    labs(title = paste("Mean values of Carbon stock in different land use types in the municipality of", province_name),
         subtitle= "Above and Below Ground Carbon (T/ha)",
         x = NULL, y = "ABGC (t/ha)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, vjust = 1, size = 8),
      axis.ticks.x = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 115))
  # SOC Plot
  plot_soc <- ggplot(df_prov, aes(x = description, y = mean_SOC, fill = description)) +
    geom_col(alpha=0.7) +
    geom_text(aes(label = category),
              angle = 60, 
              vjust=-0.5,
              hjust = 0, 
              size = 2.5, 
              color = "black") +
    scale_fill_manual(values = category_colors) +
    scale_x_discrete(labels = df_prov$esalc) +
    labs(title = " ",
         subtitle = "Soil Organic Carbon (T/ha)",
         x = NULL, y = "SOC (t/ha)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, vjust = 1, size = 8),
      axis.ticks.x = element_blank()
    ) +
    scale_y_continuous(limits = c(0, 130))
  # Combine and save
  combined_plot <- plot_abgc + plot_soc
  ggsave(
    filename = paste0("carbon_stock_", province_name, ".png"),
    plot = combined_plot,
    width = 18,
    height = 10,
  )
}


######### Analysis of distribution of c stock values across provinces #########
combined_df <- bind_rows(mean_c_values_esalc, .id = "province")

abgc_means <- combined_df %>%
  group_by(province) %>%
  summarise(mean_abgc = mean(mean_ABGC, na.rm = TRUE))
combined_df <- combined_df %>%
  left_join(abgc_means, by = "province")

soc_means <- combined_df %>%
  group_by(province) %>%
  summarise(mean_soc = mean(mean_SOC, na.rm = TRUE))
combined_df <- combined_df %>%
  left_join(soc_means, by = "province")

# ABGC boxplot per provincia
box_abgc <- ggplot(combined_df, aes(x = province, y = mean_ABGC, fill = mean_abgc)) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.2) + #fill = "#56B4E9"
  scale_fill_viridis(name = "Mean ABGC", option = "C", direction = -1) +
  labs(title = "Distribution of mean Above and Below ground carbon values per provinces",
       x = "Province",
       y = "Mean ABGC (t/ha)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# SOC boxplot per provincia
box_soc <- ggplot(combined_df, aes(x = province, y = mean_SOC, fill=mean_soc)) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.2) + #fill = "#E69F00", 
  scale_fill_viridis(name = "Mean SOC", option = "C", direction = -1) +
  labs(title = "Distribution of mean Soil Organic Carbon values per provinces",
       x = "Province",
       y = "Mean SOC (t/ha)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

combined_plot <- box_abgc + box_soc
ggsave(
  filename = "boxplot_ABGC_SOC_per_province.png",
  plot = combined_plot,
  width = 30,
  height = 7,
)



