#Load the relevant libraries
library(reshape2)
library(tidyverse)
library(tidyr)
library(moments)
library(ggplot2)
library(caret)
library(DescTools)
library(dplyr)
library(ISLR)

# Read the Bio_measure dataset
Bio_measure <- read.csv("proportional_species_richness_V3.csv")

#explore Bio measures
head(Bio_measure)
summary(Bio_measure)
attach(Bio_measure)
str(Bio_measure)

# converting variables
Bio_measure$period <- as.factor(Bio_measure$period) #categorical variable
Bio_measure$dominantLandClass <- as.factor(Bio_measure$dominantLandClass)

# Setting column name for the Bio_measure dataset
colnames(Bio_measure) <- c("Location", "Bees", "Bird", "Bryophytes", "Butterflies", "Carabids", "Hoverflies", "Isopods", "Ladybirds", "Macromoths", "Grasshoppers_._Crickets", "Vascular_plants", "Easting", "Northing", "dominantLandClass", "ecologicalStatus", "period")

#Create BD7, BD4, and BD11 subsets from Bio_measure
BD7 <- Bio_measure %>%
  select(Bird, Bryophytes, Butterflies, Hoverflies, Isopods, Ladybirds, Macromoths)
View(BD7)
summary(BD7)
BD4 <- Bio_measure %>%
  select(Bees, Carabids, Grasshoppers_._Crickets, Vascular_plants)

BD11 <- Bio_measure %>%
  select(Bird, Bryophytes, Carabids, Hoverflies, Isopods, Grasshoppers_._Crickets, Vascular_plants, Bees, Macromoths, Butterflies, Ladybirds)

#Univariate analysis of BD7 variables:
summary(BD7)

#Creating boxplots for each variable in BD7
boxplot(BD7$Bird)
boxplot(BD7$Bryophytes)
boxplot(BD7$Butterflies)
boxplot(BD7$Hoverflies)
boxplot(BD7$Isopods)
boxplot(BD7$Ladybirds)
boxplot(BD7$Macromoths)

# Multivariate analysis for BD7 correlation
plot(BD7)
Cor(BD7)

# hypothesis test on the variables
model <- lm(Bio_measure$Bird ~ Bio_measure$Easting, data = Bio_measure)
summary(model)

model1 <- lm(Bio_measure$Macromoths ~ Bio_measure$Easting, data = Bio_measure)
summary(model1)

#Correlation between BD7 variables and other variables
BD7_with_other_variables <- cbind(BD7, Bio_measure %>%
                                    select(period, Easting, Northing, dominantLandClass))

cor_matrix <- cor(BD7_with_other_variables %>%
                    select_if(is.numeric))

#Visualize the correlation matrix using ggplot2
melted_cor_matrix <- melt(cor_matrix)
ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12)) +
  coord_fixed()


#linear regression for BD7 variable
summary(model.Bird<- lm(ecologicalStatus~Bird))
summary(model.Bryophytes<-lm(ecologicalStatus~Bryophytes))
summary(model.Butterflies <- lm(ecologicalStatus~Butterflies))
summary(model.Hoverflies<-lm(ecologicalStatus~Hoverflies))
summary(model.Isopods<-lm(ecologicalStatus~Isopods))
summary(model.Ladybirds<-lm(ecologicalStatus~Ladybirds))
summary(model.Macromoths<-lm(ecologicalStatus~Macromoths))

par(mfrow = c(1,3))
plot(Bird, ecologicalStatus, col= "red")
abline(model.Bird, lwd = 3, col = "blue")

par(mfrow = c(1,3))
plot(Bryophytes, ecologicalStatus, col= "red")
abline(model.Bryophytes, lwd = 3, col = "blue")

par(mfrow = c(1,3))
plot(Butterflies, ecologicalStatus, col= "red")
abline(model.Butterflies, lwd = 3, col = "blue")

par(mfrow = c(1,3))
plot(Hoverflies, ecologicalStatus, col= "red")
abline(model.Hoverflies, lwd = 3, col = "blue")

par(mfrow = c(1,3))
plot(Isopods, ecologicalStatus, col= "red")
abline(model.Isopods, lwd = 3, col = "blue")

par(mfrow = c(1,3))
plot(Ladybirds, ecologicalStatus, col= "red")
abline(model.Ladybirds, lwd = 3, col = "blue")

par(mfrow = c(1,3))
plot(Macromoths, ecologicalStatus, col= "red")
abline(model.Macromoths, lwd = 3, col = "blue")

#Simple linear regression between BD7 and BD11
lm_BD11 <- lm(Bees + Carabids + Grasshoppers_._Crickets + Vascular_plants ~ ., data = BD11)
summary(lm_BD11)

#Simple linear regression between BD7 and BD11 for each period
lm_BD11_period_Y70 <- lm(Bees + Carabids + Grasshoppers_._Crickets + Vascular_plants  ~ ., data = filter(Bio_measure, period == "Y70") %>% select(colnames(BD11)))
summary(lm_BD11_period_Y70)

lm_BD11_period_Y00 <- lm(Bees + Carabids + Grasshoppers_._Crickets + Vascular_plants ~ ., data = filter(Bio_measure, period == "Y00") %>% select(colnames(BD11)))
summary(lm_BD11_period_Y00)

#select
cont_vars <- Bio_measure %>%
  select(Bird, Bryophytes, Butterflies, Hoverflies, Isopods, Ladybirds, Macromoths,Northing,Easting)
names(cont_vars)
View(cont_vars)
str(cont_vars)

# Get the column indices corresponding to the columns in the BD7 data frame
# Get the column indices corresponding to the columns in the BD7 data frame
all <- 1:ncol(BD7)

# BD 7 predictors to form the eco_stat
eco_selected <- BD7[,1:7]
eco_selected <- c(3,4,5,7,8,9,10)

#predictors not selected
eco_not_selected <- all[!(all %in% eco_selected)]

#predictors in the BD7 data frame
eco_names <- colnames(BD11)

#selected predictors
eco_selected_names <- colnames(BD7)[1:7]

# Display the names of the selected predictors
eco_selected_names
eco_selected
?rowmeans

# Bio div measure over 7 taxonomic groups
mean_selected <- rowMeans(Bio_measure[ ,eco_selected],na.rm=TRUE) # mean the 7 columns
sum(is.na(mean_selected)) # check that there are no NAs in mean_selected

# Add in the biodiversity measure which is the mean over 7 taxonomic groups
Bio_measure_MA334 <- Bio_measure %>% mutate(eco_status_7 = mean_selected)
View(Bio_measure_MA334)

# the data exploration phase
table <- data.frame()
for(i in eco_selected){
  table <- rbind(table,
                 c(eco_names[i-1],
                   round(mean(Bio_measure_MA334[,i],na.rm = TRUE),digits = 2),
                   round(sd(Bio_measure_MA334[,i],na.rm = TRUE),digits = 2),
                   round(skewness(Bio_measure_MA334[,i],na.rm = TRUE),digits = 2)
                 ))}
colnames(table) <- c("taxi_group","mean","sd","skewness")
table%>%arrange(sd,skewness)

plot(cont_vars$Northing~cont_vars$Easting) 

# eastings and northings as predictors
plot(Bio_measure_MA334$eco_status_7~Bio_measure_MA334$Easting)
cor(Bio_measure_MA334$eco_status_7,Bio_measure_MA334$Easting)
plot(Bio_measure_MA334$eco_status_7~Bio_measure_MA334$Northing)  # for BD7
cor(Bio_measure_MA334$eco_status_7,Bio_measure_MA334$Northing)

# linear regression with only Northing as a predictor
lin_mod <- lm(Bio_measure_MA334$eco_status_7~Bio_measure$Northing)
summary(lin_mod)
abline(lin_mod,col="green")
plot(jitter(fitted(lin_mod)),residuals(lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0)
qqnorm(lin_mod$residuals)
qqline(lin_mod$residuals,col="red")

# box plot comparisons for the two periods ignoring all other varaibles
eco_status <- Bio_measure_MA334%>%pull(eco_status_7)
eco_period <- Bio_measure_MA334%>%pull(period)
plot(eco_status~eco_period)

#BD7 distribution
names(Bio_measure_MA334)
Bio_measure_MA334_period <- Bio_measure_MA334%>%select(Location,period,eco_status_7)
Bio_measure_MA334_split <- Bio_measure_MA334_period%>%pivot_wider(names_from =period,values_from=eco_status_7)
Bio_measure_MA334_split <- Bio_measure_MA334_split%>%mutate(BD7_change=Y00-Y70)
head(Bio_measure_MA334_split)
hist(Bio_measure_MA334_split$BD7_change)  # the distribution of the BD7 change

BD7_change <- Bio_measure_MA334_split%>%pull(BD7_change)
t.test(BD7_change,mu=0)  # t test with H0: mu=0

# explicit calculation for the same t test 
s_mean <- mean(BD7_change)
s_sd <- sd(BD7_change)
sample_size <- length(BD7_change)
t_value <- s_mean/(s_sd/sqrt(sample_size )) # calculate the t statistic
t_value
2*pt(-abs(t_value),sample_size-1) # calculate the two tail probability for t_value

# comparing the two distributions of bio div based on 7 and 11 taxonomic groups
par(mfrow=c(1, 1))  # divide graph area in 1 columns
qqplot(Bio_measure_MA334$eco_status_7,Bio_measure_MA334$ecologicalStatus)
abline(0,1,col="red")
# both cdfs together  and do a kolmogorov test H0: distributions are the same
BD7_cdf <- ecdf(Bio_measure_MA334$eco_status_7)
BD11_cdf <- ecdf(Bio_measure_MA334$ecologicalStatus)
plot(BD11_cdf,col="red")
lines(BD7_cdf,col="green")
ks.test(Bio_measure_MA334$eco_status_7,Bio_measure_MA334$ecologicalStatus)

#Multiple linear regression
BD4_mean <- rowMeans(BD4)
BD7_mean <- rowMeans(BD7)

# Multiple linear regression with all predictors
lm_BD4_BD7 <- lm(BD4_mean ~ ., data = BD7)
summary(lm_BD4_BD7)

lm_BD7_pop <- lm(BD7_mean ~ ., data = BD7)
summary(lm_BD7_pop)

#Assessment of the fitted regression model
BD7$residuals <- residuals(lm_BD4_BD7)
BD7$predicted <- predict(lm_BD4_BD7)

BD7$residuals <- residuals(lm_BD7_pop)
BD7$predicted <- predict(lm_BD7_pop)

# Feature selection
step_model <- step(lm_BD4_BD7, trace = 0)
summary(step_model)

# Splitting into training and test sets
set.seed(123)
trainIndex <- createDataPartition(BD4_mean, p = 0.8, list = FALSE)
train_data <- BD7[trainIndex, ]
test_data <- BD7[-trainIndex, ]
train_data_BD4_mean <- BD4_mean[trainIndex]
test_data_BD4_mean <- BD4_mean[-trainIndex]

# Build training data model
lm_BD4_BD7_train <- lm(train_data_BD4_mean ~ ., data = train_data)
summary(lm_BD4_BD7_train)

# Testing model on test data
predictions <- predict(lm_BD4_BD7_train, test_data)

# Calculate residuals
residuals <- test_data_BD4_mean - predictions

# Check normality of residuals using a histogram
hist(residuals, main = "Histogram of Residuals", xlab = "Residuals", col = "lightblue", border = "black")


