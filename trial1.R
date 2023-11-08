
library(readxl)
library(ggplot2)
library(dplyr)
library(reshape2)
biomarkers_filepath<-'C:/Users/olive/Desktop/MSc Ed/Introduction to Probability and Statistics/Assignment/biomarkers.xlsx'
covariates_filepath<-'C:/Users/olive/Desktop/MSc Ed/Introduction to Probability and Statistics/Assignment/covariates.xlsx'
biomarkers_data <- read_excel(biomarkers_filepath)
covariates_data <-read_excel(covariates_filepath)
biomarkers_data <- as.data.frame(biomarkers_data)
covariates_data <- as.data.frame(covariates_data)


#biomarkers_data
#covariates_data



## HYPOTHESIS TEST

# considering whether biomarkers vary between males and females at inclusion
# mu_m: biomarkers level male
# mu_f biomarkers level female
# H_0: mu_m = mu_f
# H_a: mu_m != mu_f


## in these lines we combine the biomarkers and covariates data by sorting by patient number
## we can clean up the data by looking at biomarkers at inclusion only

patient_number <- as.numeric(sapply(strsplit(biomarkers_data$Biomarker, "-"), '[', 1))                    
biomarkers_data <- biomarkers_data[order(patient_number),]                                                       
filtered_data <- biomarkers_data[grep("-0weeks$", biomarkers_data$Biomarker),]                               
filtered_data$patient_number <- as.numeric(sapply(strsplit(filtered_data$Biomarker, "-"), '[', 1))            
master_data <- merge(covariates_data, filtered_data, by.x = "PatientID", by.y = "patient_number", all.x = TRUE)       

combined_data <- master_data[master_data$PatientID != 40, ]   ## removing patient 40 since they have no biomarkers data
biomarkers <- c("IL-8", "VEGF-A", "OPG", "TGF-beta-1", "IL-6", "CXCL9", "CXCL1", "IL-18", "CSF-1")

male_data <- combined_data[combined_data$`Sex (1=male, 2=female)` == 1, ]   ## define new dataframes for male and female
female_data <- combined_data[combined_data$`Sex (1=male, 2=female)` == 2, ]

gender_matrix <- matrix(0, nrow=3, ncol=length(biomarkers))
rownames(gender_matrix) <- c("t value", "df", "p-value")
colnames(gender_matrix) <- biomarkers


for (biomarker in biomarkers) {                                                 ## Loop through biomarkers and conduct t-tests
  test_result <- t.test(male_data[[biomarker]], female_data[[biomarker]])
  
  gender_matrix["t value", biomarker] <- test_result$statistic                  ## Populate the results matrix
  gender_matrix["df", biomarker] <- test_result$parameter
  gender_matrix["p-value", biomarker] <- test_result$p.value
}


gender_df <- as.data.frame(gender_matrix)
## print(gender_df)



## find QQ plots to judge how well data follows normal distribution

biomarker_vars <- c("IL-8", "VEGF-A", "OPG", "TGF-beta-1", "IL-6", "CXCL9", "CXCL1", "IL-18", "CSF-1")


male_biomarkers <- male_data %>% select(one_of(biomarker_vars))                 ## select only biomarkers columns               
male_long <- melt(male_biomarkers, id.vars = NULL)                              ## separate biomarkers into rows
male_long$sex <- 'Male'                                                         ## Sex column for plotting
male_long

female_biomarkers <- female_data %>% select(one_of(biomarker_vars))             ## same for female data
female_long <- melt(female_biomarkers, id.vars = NULL)
female_long$sex <- 'Female'


combined_data <- rbind(male_long, female_long)                                  ## Combine the datasets


ggplot(combined_data, aes(sample = value, color = sex)) +                       ## plot combined data
  geom_qq() +
  geom_qq_line(color = "black") +
  facet_wrap(~ variable, scales = "free") +
  scale_color_manual(values = c("Male" = "red", "Female" = "blue")) +
  theme_bw() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(legend.position = "none")


male_count <- male_data %>%                                                     ## check male and female datapoints
  summarise(count = n()) %>%
  pull(count)
##print(paste("Number of male data points for IL-8:", male_count))

# Count female data points for IL-8
female_count <- female_data %>% 
  summarise(count = n()) %>%
  pull(count)
##print(paste("Number of female data points for IL-8:", female_count))

hist(male_data$`CXCL1`, main="Histogram for IL-8", xlab="IL-8 values", ylab="Frequency")









##################################################################
## regression, with vas 12 months as response variable, and biomarker level + covariates as explanatory variables

## first, get data into relevant format, starting from combined_data

combined_data <- combined_data[!is.na(combined_data$`Vas-12months`), ]          ## remove rows where 'Vas-12months' is NA

set.seed(123)                                                                   ## to reproduce results

index <- sample(1:nrow(combined_data), size = 0.8 * nrow(combined_data))        ## row indices for 80% sample

combined_data_80 <- combined_data[index, ]                                      ## initialize datasets
combined_data_20 <- combined_data[-index, ]

## run regression model
model <- lm(`Vas-12months` ~ Age + `Sex (1=male, 2=female)` + `Smoker (1=yes, 2=no)` + `IL-8` + `VEGF-A` + `OPG` + `TGF-beta-1` + `IL-6` + `CXCL9` + `CXCL1` + `IL-18` + `CSF-1`, data=combined_data_80)
##summary(model)


plot(model$residuals, ylab='Residuals', xlab='Fitted Values', main='Residual Plot')  ## plot residuals
abline(h=0, col='red')

hist(model$residuals, breaks=20, main="Histogram of Residuals", xlab="Residuals", col="blue")  ## histogram of residuals




#################################
## predict values for last 20%


predictions <- predict(model, newdata = combined_data_20)                       ## make predictions

combined_data_20$predicted_Vas_12months <- predictions                          ## add to dataframe


actual <- combined_data_20$`Vas-12months`                                       ## actual values
predicted <- combined_data_20$predicted_Vas_12months                            ## predicted values


difference_vals <- abs(actual-predicted)


mae <- mean(abs(actual - predicted))                                            ## mean absolute error
rmse <- sqrt(mean((actual - predicted)^2))                                      ## rms error

percentage_difference <- (difference_vals / actual) * 100                       ## finding percentage difference
non_zero_actuals <- combined_data_20$`Vas-12months` != 0                        ## remove zero values
results_dataframe <- data.frame(                                                ## add results to dataframe
  Actual = actual[non_zero_actuals],
  Predicted = predicted[non_zero_actuals],
  AbsoluteDifference = difference_vals[non_zero_actuals]
)

results_dataframe$PercentageDifference <- (results_dataframe$AbsoluteDifference / results_dataframe$Actual) * 100



mpe <- mean(results_dataframe$PercentageDifference)
mape <- mean(abs((actual - predicted) / actual)) * 100
r_squared <- summary(model)$r.squared
adjusted_r_squared <- summary(model)$adj.r.squared



stats_df <- data.frame(                                                         ## add all results to dataframe
  MAE = mae,
  RMSE = rmse,
  MPE = mpe,
  MAPE = mape,
  R_Squared = r_squared,
  Adjusted_R_Squared = adjusted_r_squared,
)

## print(stats_df)


plot(predicted, actual, xlab="Predicted Values", ylab="Actual Values", pch=19)  ## plot predicted and actual values
abline(0, 1, col="red")
