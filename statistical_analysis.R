# Statistical analysis

library(tidyverse)
library(readxl)
library(factoextra)
library(emmeans)

# Import and prepare data ----
df <- read_excel(path = "dataset.xlsx", col_names = T)

df$minute <- factor(df$minute, ordered = TRUE)
df$uID <- as.factor(df$uID)
df$sex <- as.factor(df$sex)
df$condition <- as.factor(df$condition)
df$substance <- as.factor(df$substance)
df$concentration <- as.factor(df$concentration)
glimpse(df)

# Shape data for analysis ----
summarydf <- df %>% 
  gather(key = "measure", 
         value = "value", 
         -top_half, # same value from 2Dside and 3D, because absent in 2Dtop
         -minute, -uID, -condition, -substance, -concentration, -sex
  ) %>% 
  separate(col = measure,
           into = c("measure", "view"),
           sep = "_",
           remove = T) %>% 
  spread(key = measure, value = value)
summarydf$view      <- as.factor(summarydf$view)
summarydf$minute    <- as.factor(summarydf$minute)
summarydf$condition <- as.factor(summarydf$condition)
summarydf <- summarydf %>% dplyr::select(uID, condition, concentration, substance, minute, view, 
                               avgSpeed, topSpeed, 
                               turnrate, topAngSpeed,
                               avgAccel, topAccel,
                               pctTimeFreez,
                               wall,
                               top_half,
                               sex)
glimpse(summarydf)




################################################################
#
#  Principal Component Analysis
#  on the observables computed from 3D only


# Subset citalopram data
dfcit <- summarydf %>% 
  filter(view == "3D") %>%
  filter(str_detect(condition, "control|citalopram" ))

# Subset ethanol data
dfeth <- summarydf %>% 
  filter(view == "3D") %>%
  filter(str_detect(condition, "control|ethanol" ))


# select dataset to do the PCA on
# data <- dfcit
data <- dfeth


# Use varimax-rotation method
ncomp <- 3
pca_res         <- prcomp(data[, 7:15], center = T, scale = T)
get_eigenvalue(pca_res)
rawLoadings     <- pca_res$rotation[,1:ncomp] %*% diag(pca_res$sdev, ncomp, ncomp)
rotatedLoadings <- varimax(rawLoadings)$loadings
rotatedLoadings[,3] <- -rotatedLoadings[,3] # Uncomment for ethanol
rotatedLoadings[,1] <- -rotatedLoadings[,1]
rotatedLoadings # Report this in the table
# To see loadings < 0.1, call "rotatedLoadings[,]"
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- scale(data[, 7:15]) %*% invLoadings
print(scores[1:5,])                   # Scores computed via rotated loadings


# Prepare for anova
pca_ind_coord <- as.data.frame(scores)

pca_for_anova <- data %>% dplyr::select(uID, condition, concentration, substance, minute, view, sex)
pca_for_anova$PC1 <- pca_ind_coord[,1]
pca_for_anova$PC2 <- pca_ind_coord[,2]
pca_for_anova$PC3 <- pca_ind_coord[,3] # Change the sign for ethanol

data <- pca_for_anova

# Uncomment observable of interest:
data$observable <- data$PC1
# data$observable <- data$PC2
# data$observable <- data$PC3

# Repeated measures anova
anov.lm <- aov(observable ~ (condition * minute * sex) +
                 Error(uID/minute),
               data = data)
summary(anov.lm)

# Post-hoc test
contr <- emmeans(anov.lm, "minute")
contr <- emmeans(anov.lm, "condition")
contr <- emmeans(anov.lm, "minute", by = "condition")
contrast(contr, "trt.vs.ctrl")




################################################################
#
#  Repeated measures ANOVAs ----

# Subset citalopram data
dfcit <- summarydf %>% 
  filter(view == "3D") %>% # uncomment for analysis of 3D only
  filter(str_detect(condition, "control|citalopram" ))
glimpse(dfcit)

# Subset ethanol data
dfeth <- summarydf %>% 
  filter(view == "3D") %>% # uncomment for analysis of 3D only
  filter(str_detect(condition, "control|ethanol" ))
glimpse(dfeth)


# Uncomment condition to analyse ----
data <- dfcit
# data <- dfeth


# Uncomment observable of interest ----
data$observable <- data$avgSpeed
# data$observable <- data$topSpeed
# data$observable <- data$turnrate
# data$observable <- data$topAngSpeed
# data$observable <- data$avgAccel
# data$observable <- data$topAccel
# data$observable <- data$pctTimeFreez
# data$observable <- data$wall

# For top_half, because computed only from 2Dside ----
# data <- data %>% filter(data$view == "2Dside")
# data$observable <- data$top_half


# Repeated measures anova on 3D data ----
anov.lm <- aov(observable ~ (condition * sex * minute) +
               Error(uID/minute),
               data = data)


# Repeated measures anova on data from all views ----
# anov.lm <- aov(observable ~ (condition * view * sex * minute) +
#                Error(uID/(view * minute)),
#                data = data)

summary(anov.lm)


# Post-hoc test depending on anova results
contr <- emmeans(anov.lm, "view")
contr <- emmeans(anov.lm, "condition")
contr <- emmeans(anov.lm, "condition", by = "view")
contr <- emmeans(anov.lm, "minute")
contr <- emmeans(anov.lm, "minute", by = "condition")
contr <- emmeans(anov.lm, "minute", by = "view")
contr <- emmeans(anov.lm, "minute", by = c("condition", "view"))
contrast(contr, "trt.vs.ctrl") 

