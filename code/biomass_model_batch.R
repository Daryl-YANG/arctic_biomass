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
list.of.packages <- c("ggplot2", "randomForest", "caTools", "ggpmisc")  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
# define an output directory to store outputs
our.dir <- "/Volumes/data1/projects/biomass/kougarok/analysis/biomass_model/output"
# create output directory if not exist
if (! file.exists(our.dir)) dir.create(our.dir,recursive=TRUE)
# creat an temporary to store files temporarily generated during the course of processing
temp.dir <- file.path(paste0(our.dir, "/", 'temporary'))
if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)

# define how many models you want to build?
nmodel <- 100
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define the directory to hdf files
source.dir <- "/Volumes/data1/projects/biomass/kougarok/analysis/biomass_model/input/"
data.dir <- list.files(source.dir, pattern = 'biomass_database_v2.csv',
                       full.names = TRUE)
data.org <- read.csv(data.dir)
str(data.org)
data.org$ABG_kg_m2[data.org$HAG_plotmean > 6] <- NA
data.org <- na.omit(data.org)
#*****************************************************************************************#

#************************************* biomass model *************************************#
### split data into training and validation
data.split <- sample.split(data.org, SplitRatio = 0.7) 
# get training and validation data
data.train <- subset(data.org, data.split == "TRUE") 
data.test <- subset(data.org, data.split == "FALSE") 
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
  labs(x = 'Ground Biomass (kg/m2)', y = 'Predicted Biomass (kg/m2)') +
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
  annotate('text', x= 3.33, y = 0.5, label = 'RMSE = 0.5', size = 4, hjust = 0) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pngName = paste0(our.dir, "/",'model_result_withunc.pdf')
ggsave(pngName, plot = last_plot(), width = 12, height = 12, units = 'cm')
#*****************************************************************************************#
