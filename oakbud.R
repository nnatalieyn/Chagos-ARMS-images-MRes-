#------------data----------------------

tree <- read.csv("trees.csv")
girth <- read.csv("girth.csv")
pheno <- read.csv("phenology.csv")
weather <- read.csv("SilwoodWeatherDaily.csv")

#------------packages---------------------------

library(dplyr)
library(ggplot2)
library(scales) #for plot 1 and 2 timescale
library(lubridate) #convert to julian date
library(pals) #>8 color colorblind friendly palette
library(ggeffects) #aid in drawing mixed model plots
library(lme4)
library(lmtest) #likelihood test
library(lmerTest) #include p-value to lmm
library(gtsummary) #generate lmm table for publishing
library(broom.mixed) #necessary for lmm table

#-------TYPE of data + TRANSFORM data --------------

str(weather)
str(pheno)

weather$TIMESTAMP <- as.Date(substr(weather$TIMESTAMP, 1, 10), format = "%d/%m/%Y")

pheno$Score <- as.numeric(gsub("[><]", "", pheno$Score))
pheno <- pheno[pheno$Score <= 6, ]
names(pheno)[names(pheno) == "Date"] <- "TIMESTAMP"
pheno$TIMESTAMP <- as.Date(pheno$TIMESTAMP, format = "%d/%m/%Y")

summary(pheno$TIMESTAMP) #RANGE of dates
summary(weather$TIMESTAMP) #RANGE of dates

#-----------------MERGE data-------------------

treelocationpheno <- merge(tree, pheno, by = "TreeID") #(if tree location is necessary)

treephenoweather <- merge(pheno, weather, by = "TIMESTAMP")
###assuming that non-oak species are omited

#----------------SUBSET spring time-----------

summary(treephenoweather$TIMESTAMP) #RANGE of dates to analyse

spring_data <- treephenoweather %>%
  filter(format(TIMESTAMP, "%m-%d") >= "03-01" & format(TIMESTAMP, "%m-%d") <= "05-31")

ggplot(spring_data, aes(x = TIMESTAMP)) +
  geom_bar(stat = "count") +
  labs(title = "Distribution of Dates in Spring Data",
       x = "Date",
       y = "Frequency") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month") +  # Add year-month format
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate axis labels for readability

#-----------Summary dataset (calculate mean for new dataset)------------------

data <- spring_data %>%
  group_by(TIMESTAMP) %>%
  summarise(
    score_mean = mean(Score, na.rm = TRUE),
    air_temp_max_mean = mean(Air_Temp_Max...Deg.C...Max., na.rm = TRUE),
    air_temp_min_mean = mean(Air_Temp_Min...Deg.C...Min., na.rm = TRUE),
    .groups = "drop" #drop grouping after summarizing
  )

#-----------Year, Accumulated/julian dates, GDD variables---------

#SEPERATE timestamp into year and days accumulated columns
data$Year <- format(data$TIMESTAMP, "%Y")

data$Accum_days <- format(data$TIMESTAMP, "%m-%d")
data$Accum_days <- as.character(data$TIMESTAMP - as.Date(paste0(data$Year, "-03-01")))

data <- data %>%
  mutate(GDD = if_else(((air_temp_max_mean + air_temp_min_mean) / 2) < 5, 
                       0, 
                       ((air_temp_max_mean + air_temp_min_mean) / 2) - 5))

#-----------------------Explore data-----------------------
summary(data)
plot(data)
par(mfrow = c(1,2))
hist(data$score_mean)
hist(data$GDD)

#z-transform explanatory var - wtf i don't get it...
#data$GDDz <- scale(data$GDD, center = TRUE, scale = TRUE)
#hist(data$GDDz)

dev.off()
#par(mfrow = c(1,2))

plot(data$score_mean, data$GDD)
cor.test(data$score_mean, data$GDD)

#data <- data %>%
#  mutate(logscore_mean = log(score_mean + 1))
#hist(data$logscore_mean)

#plot(data$logscore_mean, data$GDD)
#cor.test(data$logscore_mean, data$GDD)

#-----------------------Model-----------------------

# use linear model - no
model1 <- lm(score_mean ~ GDD + Year , data = data)
summary(model1)

par(mfrow = c(2,2))
plot(model1)

#take off outliers - don't change plot
#rows_removed <- c(60, 235, 339)
#data_noutlier1 <- data[-rows_removed, ]
#par(mfrow = c(2,4))
#model2 <- lm(score_mean ~ GDD + Year, data = data_noutlier1)
#summary(model2)
#plot(model2)

# use linear mixed model - yes
model3 <- lmer(score_mean ~ GDD * Year + (1 | Accum_days), data = data) #interaction between GDD and score_mean varies by year (as a fixed effect)

#model4 <- lmer(score_mean ~ GDD + (1 + GDD | Year) + (1 | Accum_days), data = data) #year as source of variability rather than fixed interaction
#summary(model4)

#Model 3 - Fixed Interaction (GDD * Year): In model3, Year is treated as a fixed effect, which means you explicitly model the interaction between GDD and each Year. This allows you to test and interpret the effect of GDD within each specific year, but you must estimate separate coefficients for each year's interaction with GDD.
#Interpretation: The impact of GDD on score_mean can be compared directly across years with specific coefficients for each year's interaction term.
#Model 4 - Random Slope and Intercept for Year: In model4, Year is a random effect, meaning you assume that each year has a unique response to GDD but you do not estimate individual coefficients for each year. Instead, the model allows the intercepts and slopes to vary across years according to a distribution.
#Interpretation: This is more flexible than model3 because it accounts for variation in GDD's effect without having to fit a separate parameter for each year. You interpret it as allowing for year-to-year variability without needing to directly estimate each year's effect.
#When to Use Each Model:
#Model 3 (Fixed Year):Use this when you want to explicitly model and interpret the interaction between GDD and each specific Year. It is helpful when you have enough data to estimate separate effects for each year. Suitable for understanding if and how the impact of GDD differs by each year with a clear output for each one.
#Model 4 (Random Year):Use this when you want to account for variability between years but do not need or want to estimate specific effects for each year. This is beneficial when Year should be treated as a random sample from a larger population of possible years.Appropriate when you expect the effect of GDD to vary by year, but you’re interested in the overall variability rather than individual year-by-year coefficients.
#Summary:
#Model 3 gives you specific, interpretable estimates for each year but is more parameter-heavy.
#Model 4 captures the variation in GDD's effect across years in a more general way without estimating each year separately, leading to a more parsimonious model if you view Year as a source of random variation.

summary(model3) #interaction GDD * year help identify if response of score_mean to GDD change year on year. GDD have a stronger/weaker influence in some years compared to others
print(model3, correlation = T)

plot(model3) #no patterns evident, good
plot(residuals(model3))
hist(residuals(model3)) #symmetrically distributed around 0
qqnorm(resid(model3))
qqline(resid(model3))

#variance explained
2.84/(2.84+0.31) #accum_days explains ~90.2% of variance in data

#likelihood test
model3_reduced <- lmer(score_mean ~ GDD + (1 | Accum_days), data = data)
lrtest(model3, model3_reduced) #determine if adding or removing Year significantly improves nested model fit

class(model3)

#summary table with AIC, BIC and logLik
table1 <-
  tbl_regression(model3, exponentiate = TRUE) %>%
  add_glance_table(include = c(AIC, BIC, logLik))
glance(model3)
print(table1)

#------------------------Plots---------------------------------------

##plot 1 - mean temp variability throughout sample period

# create Julian dates
start_date <- as.Date("2010-04-07")
data_scaled$Julian <- as.numeric(data_scaled$TIMESTAMP - start_date + 1)

#create axis labels
years <- unique(year(data_scaled$TIMESTAMP)) #create vector of unique years in timestamp column

#create breaks at julian date of every year (axis breaks)
breaks_for_each_year <- as.numeric(as.Date(paste0(years, "-03-01")) - as.Date("2010-03-01") + 1) # Create a vector of Julian dates for March 1st of each year

# Create a plot with maximum and minimum air temperature plotted separately with fitted lines
ggplot(data_scaled, aes(x = Julian)) +
  geom_point(aes(y = air_temp_max_mean, color = "Max Air Temp"), alpha = 0.5) +
  geom_point(aes(y = air_temp_min_mean, color = "Min Air Temp"), alpha = 0.5) +
  geom_smooth(aes(y = air_temp_max_mean, color = "Max Air Temp"), method = "lm", se = TRUE) +
  geom_smooth(aes(y = air_temp_min_mean, color = "Min Air Temp"), method = "lm", se = TRUE) +
  scale_x_continuous(breaks = breaks_for_each_year, labels = years) +
  scale_color_manual(values = c("Max Air Temp" = "#DC3220", "Min Air Temp" = "#005AB5")) +
  labs(x = "Sample Years (March - May)", 
       y = "Temperature (°C)", 
       title = "Temperature Variation of Spring 2010-2023") +
  theme_classic()

##plot 2 - GDD accumulation throughout sample period

# create new column for accumulated GDD separately for each year
data_scaled <- data %>%
  group_by(Year) %>%
  mutate(accumulated_GDD = cumsum(GDD)) %>%
  ungroup()

data_scaled$Accum_days <- as.numeric(data_scaled$Accum_days)

pal <- c("#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff",
         "#920000","#db6d00","#24ff24","#ffff6d")

ggplot(data_scaled, aes(x = Accum_days, y = accumulated_GDD, group = Year)) +
  geom_line(aes(color = factor(Year)), size = 1) +  # Line plot with color-coded years
  xlim(c(20, 95)) +
  scale_x_continuous(breaks = c(21, 25, 35, 45, 55, 65, 75, 85, 91),
                     labels = c("22/3", "30/3", "5/4", "15/4", "25/4", "5/5", "15/5", "25/5", "31/05" )) +
  scale_color_manual(values = pal) +
  labs(x = "Days in Spring (DD/MM)", 
       y = "Accumulated GDD5",
       title = "Accumulated GDD from March to May (2010-2023)") +
  theme_classic() +
  theme(legend.position = "right")

##plot 3 - effect of GDD on mean_score by plotting data directly + trendline for predicted values

pred1 <- ggpredict(model3, terms = c("GDD"), type = "fe")

(ggplot(pred1) + 
    geom_line(aes(x = x, y = predicted), size = 1) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = data,                      # adding the raw data (scaled values)
               aes(x = GDD, y = score_mean, colour = Year)) + 
    ylim(c(0, 6)) +
    scale_color_manual(values = pal) +
    labs(x = "GDD5", y = "Bud Burst Score", 
         title = "Effect of GDD on Bud Burst") + 
    theme_minimal()
)

##plot 4 - effect of GDD on mean_score as estimated by model3

#create a data frame with predicted values
#data$predicted_values <- predict(model3, re.form = NA)  #re.form = Only use fixed effects

#ggplot(data, aes(x = GDD, y = score_mean, color = factor(Year))) +
#  geom_point(alpha = 0.6) +  # Scatter plot for actual data points
#  geom_line(aes(y = predicted_values), size = 1) +  # Line for predicted values
#  labs(
#    title = "Mean Bud Burst Score vs. GDD Across Years",
#    x = "Growing Degree Days (GDD)",
#    y = "Mean Bud Burst Score",
#    color = "Year"
#  ) +
#  theme_minimal() +
#  theme(
#    legend.position = "bottom",
#    plot.title = element_text(hjust = 0.5)
#  )

#- OR -

pred2 <- ggpredict(model3, terms = c("GDD", "Year"), type = "fe")
#plot(pred2) +
  #scale_color_manual(values = scales::hue_pal()(length(unique(pred2$group)))) +  # Dynamically generate colors
  #labs(x = "GDDz", y = "Mean Bud Burst Score", title = "Effect of GDD on Bud Burst Score") + 
  #theme_classic()

##plot 5 - plot 4 but seperated by year

pred2_df <- as.data.frame(pred2) #turn pred2 into data frame so easier for plotting

ggplot() +
  geom_line(data = pred2_df, aes(x = x, y = predicted, color = factor(group)), size = 1) +  # Plot predicted trendlines
  geom_ribbon(data = pred2_df, aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error, fill = factor(group)), alpha = 0.4) +  # Set ribbon colors to factor(group)
  facet_wrap(~ group, nrow = 3, scales = "free") +  # Create a panel for each year, assuming 'group' is the Year column
  labs(x = "GDD5", y = "Bud Burst Score", title = "Effect of GDD on Bud Burst Score") +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.position = "none", panel.spacing = unit(2, "lines")) +   # adding space between panels
  coord_cartesian(xlim = c(min(pred2_df$x), max(pred2_df$x)), ylim = c(min(pred2_df$predicted - pred2_df$std.error), max(pred2_df$predicted + pred2_df$std.error)))  #axis limits

#----------------------Data analysis---------------------------------------

# Assuming `data` is your dataset with columns Year, Accum_days, and score_mean

# Set a threshold for the bud burst score (e.g., 3.0)
threshold <- 3.0

# Create a function to find the first day when the score_mean reaches or exceeds the threshold
find_first_day <- function(year_data) {
  result <- year_data[year_data$score_mean >= threshold, ]
  if (nrow(result) > 0) {
    return(min(result$Accum_days))
  } else {
    return(NA)  # Return NA if the threshold isn't met
  }
}

# Apply this function for each year

threshold_days <- data_scaled %>%
  group_by(Year) %>%
  summarize(first_day = find_first_day(cur_data()))

print(threshold_days)

# Plot accumulated GDD over days for a sample year (e.g., 2017)
ggplot(data_scaled %>% filter(Year == 2021), aes(x = Accum_days, y = accumulated_GDD, group = 1)) +
  geom_line() +
  labs(title = "Accumulated GDD for 2021", x = "Day of the Year", y = "Accumulated GDD")
