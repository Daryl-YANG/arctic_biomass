###########################################################################################
#
#  this script extracts spectral reflectance (and quality control layer if preferred) from
#  original modis brdf-corrected reflectance hdf files downloaded from DACC, convert hdf 
#  to tif and merge spectral bands. the output is a single raster file for each hdf file
#  that contains all desired spectral bands and quality control layers
#
#    --- Last updated:  2021.06.22 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2", "randomForest", "caTools", "ggpmisc", "terra", 
                      "foreach", "doParallel")  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define an output directory to store outputs
out.dir <- "/Volumes/data1/projects/biomass/kougarok/analysis/biomass_model/output"
# create output directory if not exist
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
# creat an temporary to store files temporarily generated during the course of processing
temp.dir <- file.path(paste0(out.dir, "/", 'temporary'))
if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)

# define how many models you want to build?
nmodel <- 100

# define the output resolution of biomass map
reso <- 5 # m
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define directory to biomass training files
source.dir <- "/Volumes/data1/projects/biomass/kougarok/analysis/biomass_model/input/"
data.dir <- list.files(source.dir, pattern = 'biomass_database_v2.csv',
                       full.names = TRUE)
data.org <- read.csv(data.dir)
str(data.org)
data.org$ABG_kg_m2[data.org$HAG_plotmean > 6] <- NA
data.org <- na.omit(data.org)

# load in canopy height model map
chm.dir <- list.files(source.dir, pattern = 'kougarok_chm_1m.tif',
                      full.names = TRUE)
chm.rst <- terra::rast(chm.dir)
#*****************************************************************************************#

#************************************* biomass model *************************************#
### split data into training and validation
data.train <- data.org
data.test <- data.org
# extract data for random forest regression
rf.data.train <- data.train[, c(15:20)]
# build a random forest model 
test.pred <- c()
for (i in 1:nmodel)
{
  # select a portion of training data for each model train
  train.split <- sample.split(rf.data.train, SplitRatio = 0.9) 
  data.train <- subset(rf.data.train, train.split == "TRUE") 
  # build a random forest model
  biomass.model <- randomForest(ABG_kg_m2 ~ ., data = data.train,
                                ntree = 90, mtry = 2, importance = TRUE)
  # apply model to validation data
  rf.model.test <- predict(biomass.model, newdata = data.test)
  # store prediction
  test.pred <- cbind(test.pred, rf.model.test)
}
test.pred <- data.frame(test.pred)
pred.mean <- apply(test.pred, 1, FUN = mean)
pred.sd <- apply(test.pred, 1, FUN = sd)
#*****************************************************************************************#
#*
#************************************** make a plot **************************************#
# make a plot of the predict result
plt.data <- data.frame('truth' = data.test$ABG_kg_m2,
                       'predicted' = pred.mean,
                       'unc' = pred.sd)

rmse <- sqrt(sum((plt.data$predicted-plt.data$truth)^2)/nrow(plt.data))

formula <- y ~ x
ggplot(data = plt.data, aes(x = truth, y = predicted)) + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5) +
  geom_point(size = 2, shape = 1) + 
  geom_smooth(method = 'lm', formula = formula, col = 'black') +
  geom_errorbar(aes(ymin = predicted-unc, ymax = predicted+unc)) + 
  ylim(c(0, 5)) + xlim(c(0, 5)) +
  labs(x = 'Ground Biomass (g/m2)', y = 'Predicted Biomass (g/m2)') +
  theme(legend.position = 'none') +
  theme(axis.text = element_text(size=12), axis.title=element_text(size=13)) +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x.npc = 0.65, label.y.npc = 0.33,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = 0.65, label.y.npc = 0.25,
               formula = formula, parse = TRUE, size = 4, hjust = 0) +
  annotate('text', x= 2.65, y = 0.5, label = 'RMSE = 0.32', size = 4, hjust = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pngName = paste0(out.dir, "/",'model_result_withunc.pdf')
ggsave(pngName, plot = last_plot(), width = 12, height = 12, units = 'cm')
#*****************************************************************************************#

#************************************** apply model **************************************#
# calculate window based height stats for applying biomass model
window.size <- reso/xres(chm.rst)
grid.min <- aggregate(chm.rst, window.size, min, na.rm = TRUE)
grid.max <- aggregate(chm.rst, window.size, max, na.rm = TRUE)
grid.mean <- aggregate(chm.rst, window.size, mean, na.rm = TRUE)
grid.median <- aggregate(chm.rst, window.size, median, na.rm = TRUE)
grid.90percent <- aggregate(chm.rst, window.size, 
                            fun=function(i) quantile(i, probs=0.9,na.rm = TRUE),
                            cores = 6)
# stack window stats to multiple layer raster
grid.stack <- c(grid.min, grid.max, grid.mean, grid.median, grid.90percent) 
rm(grid.min, grid.max, grid.mean, grid.median, grid.90percent)
# convert raster to dataframe for apply random forest model
stack.dataframe <- as.data.frame(grid.stack, xy=TRUE)
coords <- stack.dataframe[, c(1,2)]
predictors <- stack.dataframe[, -c(1,2)]
names(predictors) <- names(data.org[, c(16:20)])

rm(grid.stack)
### apply model on application data
#Setup backend to use many processors
totalCores = detectCores()
#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)
result <- foreach(i=1:100) %dopar% {
  train.split <- caTools::sample.split(rf.data.train, SplitRatio = 0.9)
  data.train <- subset(rf.data.train, train.split == "TRUE")
  biomass.model <- randomForest::randomForest(ABG_kg_m2 ~ ., data = data.train,
                                              ntree = 90, mtry = 2, importance = TRUE)
  biomass.pred <- predict(biomass.model, newdata = predictors)
}
pred.df <- data.frame(result)

### calculate mean biomass prediction
pred.mean <- apply(pred.df, 1, FUN = mean, na.rm = TRUE)
# convert biomass dataframe to raster
pred.mean <- cbind(coords, pred.mean)
biomass.pred.rst <- terra::rast(pred.mean, type = 'xyz')
crs(biomass.pred.rst) <- crs(chm.rst)
basename <- basename(chm.dir)
outname <- paste0(out.dir, '/', 'parallel_biomass_', basename)
writeRaster(biomass.pred.rst, outname, overwrite=TRUE)

### calculate biomass prediction uncertainty
pred.unc <- apply(pred.df, 1, FUN = sd, na.rm = TRUE)
# convert biomass dataframe to raster
pred.unc <- cbind(coords, pred.unc)
biomass.unc.rst <- terra::rast(pred.unc, type = 'xyz')
crs(biomass.unc.rst) <- crs(chm.rst)
basename <- basename(chm.dir)
outname <- paste0(out.dir, '/', 'parallel_biomass_unc_', basename)
writeRaster(biomass.unc.rst, outname, overwrite=TRUE)
#*****************************************************************************************#













