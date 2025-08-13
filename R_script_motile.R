# MOTILE-0 LOAD DIRECTORY --------------------------------------------------------

#0.1 packages
library(ggplot2)        #load in this order for code to run
library(dplyr)          #data wrangling
library(webr)           #for piedonut
library(lme4)           #for GLMM
library(car)            #type II/III test
library(ggeffects)      #visualise model predictions
library(tidyverse)      #pivot dataframes
library(vegan)          #for NMDS
library(emmeans)        #for GLLM tukey
library(lmerTest)       #lmer p-values
library(tidyr)          #pivot data for UpSet
library(ComplexUpset)   #UpSet plot
library(betapart)       #partition b-diversity
library(rcartocolor)    #color palettes
library(ggtext)         #custom text
library(stringr)        #help UpSet plot

#0.2 data

  #A motile abundance
motile_count <- read.csv("assign_morpho.csv", stringsAsFactors = TRUE)

motile_count <- motile_count %>%
  mutate(
    Phyla = as.character(Phyla),
    General_name = as.character(General_name),
    ARMS_number = as.character(ARMS_number),
    Year = as.character(Year),
  )
motile_count <- motile_count %>%
  mutate(Abundance_in_image = as.numeric(as.character(Abundance_in_image)))

sapply(motile_count, class)

  #B motile size
motile_size <- read.csv("size.csv", stringsAsFactors = TRUE)

#
#
#

# MOTILE-1 SUMMARY --------------------------------------------------------

#M1.1 donut pie chart of number of morphospecies

  #summary dataset and remove empty values
motile_count %>% 
  filter(
    Phyla != "", !is.na(Phyla),
    General_name != "", General_name != "???", !is.na(General_name),
    !is.na(MorphoID)
  ) %>%
  group_by(Phyla, General_name) %>%
  summarise(n_morphosp = n_distinct(MorphoID),
            .groups = 'drop') -> PD1

  # Create chart
webr::PieDonut(PD1, 
               aes(Phyla, General_name, count = n_morphosp), 
               ratioByGroup = FALSE, #donut slices based on total number of specimens, not relative to phyla
               showRatioThreshold = 1, #only show labels for segments >XX% of total
               r0=0.2, 
               r1=0.6,
               donutLabelSize = 1,
               pieLabelSize = 1,
               labelpositionThreshold = 0.02, #if slice is >XX% of total, label will go inside slice
               showPieName = FALSE)

#
#
#

#M1.2 donut pie chart of number of specimens

  #summary dataset and remove empty values
motile_count %>% 
  filter(
    Phyla != "", !is.na(Phyla),
    General_name != "", General_name != "???", !is.na(General_name),
    !is.na(Abundance_in_image)
  ) %>%
  group_by(Phyla, General_name) %>%
  summarise(total_specimens = sum(Abundance_in_image),
            .groups = 'drop') -> PD2

  # Create chart
webr::PieDonut(PD2, 
               aes(Phyla, General_name, count = total_specimens), 
               ratioByGroup = FALSE, #donut slices based on total number of specimens, not relative to phyla
               showRatioThreshold = 1, #only show labels for segments >XX% of total
               r0=0.2, 
               r1=0.6,
               donutLabelSize = 0,
               pieLabelSize = 0,
               labelpositionThreshold = 0.02, #if slice is >XX% of total, label will go inside slice
               showPieName = FALSE)
#
#
#

#M1.3 summary stats to report
x <- sum(PD1$n_morphosp)   #total number of morphosp considered
PD1$percent_cover <- (PD1$n_morphosp/x) * 100   #% cover

y <- sum(PD2$total_specimens)   #total number of specimens
PD2$percent_cover <- (PD2$total_specimens/y) * 100   #% cover 

motile_count %>% 
  filter(
    Phyla != "", !is.na(Phyla),
    General_name != "", General_name != "???", !is.na(General_name),
    !is.na(MorphoID)
  ) %>%
  group_by(General_name, MorphoID) %>%
  summarise(n_morphosp = sum(Abundance_in_image),
            .groups = 'drop') -> abundant

#
#
#

# MOTILE-2 ABUNDANCE ------------------------------------------

#M2.1 total abundance per site and year, averaged across ARMS
abundance_summary1 <- motile_count %>%
  filter(!is.na(Abundance_in_image), MorphoID != "") %>% #filter out specimens missing images or unknown abundance
  group_by(Year, Site, ARMS_number) %>%
  summarise(n_specimens = sum(Abundance_in_image, na.rm = TRUE),
            .groups = "drop"
  )
summ1 <- abundance_summary1 %>%
  group_by(Year, Site) %>%
  summarise (
    mean = mean(n_specimens),
    sd = sd(n_specimens),
    .groups = "drop"
  )
 
#
#
#

# M2.2 GLMM of total abundance

  #total abundance including ARMS unit 
abundance_summary2 <- motile_count %>%
  filter(!is.na(Abundance_in_image), MorphoID != "") %>% #filter out specimens missing images or unknown abundance
  group_by(Year, Site, ARMS_number) %>%
  summarise(n_specimens = sum(Abundance_in_image, na.rm = TRUE),
            .groups = "drop"
  )

  #explore data structure
summary(abundance_summary2)
hist(abundance_summary2$n_specimens)
  #poisson assume mean = var & no zero-inflation
mean(abundance_summary2$n_specimens)
var(abundance_summary2$n_specimens)

  #create model
Mmodel1 <- glmer(
  n_specimens ~ Site * Year + (1 | ARMS_number),
  family = poisson(link = "log"),
  data = abundance_summary2
  )

  #model fitting
print(Mmodel1, correlation = T) #high correlation = unstable estimates
plot(Mmodel1) #no pattern = model fit well
plot(residuals(Mmodel1)) #no skew/outliers = fit well
hist(residuals(Mmodel1)) #symmetrically distributed around 0
qqnorm(resid(Mmodel1)) 
qqline(resid(Mmodel1))

summary(Mmodel1) #"dropping 2 columns" = because Court have no data for 2019 and 2021

  #visualise model prediction
pred <- ggpredict(Mmodel1, terms = c("Site", "Year"))
plot(pred) + theme_minimal()

  #likelihood ratio test (whether random effect improve model)
Mmodel0 <- glm(
  n_specimens ~ Site * Year,
  family = poisson(link = "log"),
  data = abundance_summary2
)
anova(Mmodel0, Mmodel1, test = "LRT") #p<0.05, yes significantly

#
#
#

#M2.3 ANODEV (effect of site and year)
Anova(Mmodel1, type = 3) #Type III - test each variable and interaction

#
#
#

#M2.4 Tukey
emmeans(Mmodel1, pairwise ~ Site | Year, adjust = "tukey") #effect of site within each year
emmeans(Mmodel1, pairwise ~ Year | Site, adjust = "tukey") #effect of year within each site

#
#
#

#M2.5 abundance boxplot

  #sum abundance for each replicate
narrow_motile_count <- abundance_summary2 %>%
  #create unique year-site combo
  mutate(SiteYear = paste(Site, Year, sep = "_"))

  #order data
narrow_motile_count$SiteYear <- factor(narrow_motile_count$SiteYear,
                                       levels = c("Anglaise_2019", "Anglaise_2021", "Anglaise_2022", "Moresby_2019", "Moresby_2021", "Moresby_2022", "Coin_2019", "Coin_2021", "Coin_2022", "Court_2022"))
narrow_motile_count <- narrow_motile_count %>%
  arrange(SiteYear)

  #plot
colvec2 <- c("2019" = "#60FF00",
             "2021" = "#9C9E4C", 
             "2022" = "#2ACDCC")

ggplot(narrow_motile_count, aes (x = SiteYear, y = n_specimens, fill = as.factor(Year))) +
  geom_boxplot(size = 1)+
  labs(x = "Site-Year", y = "Abundance (n)") +
  theme_minimal()+
  scale_fill_manual(values = colvec2)

#
#
#

# MOTILE-3 BIOMASS ----------------------------------------------------------

#M3.1 sums
length(unique(motile_size$MorphoID)) #n of morphosp suitable to measure
sum(!is.na(motile_size$L_mm) & motile_size$L_mm != "") #n of specimens measured
sum(!is.na(motile_size$W_g) & motile_size$W_g != "")#n of specimens weighed

#
#
#

#M3.2 min, max, mean, stddev and specimens measured by phyla
size_summary1 <- motile_size %>%
  filter(!is.na(L_mm), !is.na(W_g)) %>%
  group_by(Phyla) %>%
  summarise(
    n_specimens = n(),
    
    max_length = max(L_mm),
    min_length = min(L_mm),
    mean_length = mean(L_mm),
    sd_length = sd(L_mm),
    
    max_biomass = max(W_g),
    min_biomass = min(W_g),
    mean_biomass = mean(W_g),
    sd_biomass = sd(W_g),
    
    .groups = "drop"
  )

#M3.3 min, max, mean, stddev and specimens measured by morphoID
size_summary2 <- motile_size %>%
  filter(!is.na(L_mm), !is.na(W_g)) %>%
  group_by(MorphoID) %>%
  summarise(
    n_specimens = n(),
    
    max_length = max(L_mm),
    min_length = min(L_mm),
    mean_length = mean(L_mm),
    sd_length = sd(L_mm),
    
    max_biomass = max(W_g),
    min_biomass = min(W_g),
    mean_biomass = mean(W_g),
    sd_biomass = sd(W_g),
    
    .groups = "drop"
  )

#
#
#

# M3.4 piechart of biomass
motile_size %>% 
  filter(
    W_g != "", !is.na(W_g),
  ) %>%
  group_by(Phyla, General_name) %>%
  summarise(total_weight = sum(W_g),
            .groups = 'drop') %>%
  mutate(total_weight_4th = total_weight^(1/4)) -> PD3

PD3_clean <- PD3 %>%
  mutate(
    Phyla = as.character(Phyla),
    General_name = as.character(General_name),
    total_weight_4th = as.numeric(total_weight_4th)
  )

webr::PieDonut(PD3_clean, 
               aes(Phyla, General_name, count = total_weight_4th), 
               ratioByGroup = FALSE, #donut slices based on total number of specimens, not relative to phyla
               showRatioThreshold = 1, #only show labels for segments >XX% of total
               r0=0.2, 
               r1=0.6,
               donutLabelSize = 1,
               pieLabelSize = 1,
               labelpositionThreshold = 0, #if slice is >XX% of total, label will go inside slice
               showPieName = FALSE)

#
#
#

# M3.5 LMM of total biomass

  #total biomass including metadata
size_summary3 <- motile_size %>%
  filter(!is.na(W_g), MorphoID != "") %>% #filter out specimens missing images or unknown abundance
  group_by(Year, Site, ARMS_number) %>%
  summarise(total_biomass = sum(W_g, na.rm = TRUE),
            .groups = "drop"
  ) %>%
  mutate(total_biomass_4th = total_biomass^(1/4)) #4th root transformation to remove right skew and align with PERMANOVA/NMDS
summ2 <- size_summary3 %>%
  group_by(Year, Site) %>%
  summarise (
    mean = mean(total_biomass),
    sd = sd(total_biomass),
    .groups = "drop"
  )

size_summary3$Year <- as.factor(size_summary3$Year)

  #explore data structure
summary(size_summary3)
hist(size_summary3$total_biomass_4th) #normal distribution

  #create model
Mmodel2 <- lmer(
  total_biomass_4th ~ Site * Year + (1 | ARMS_number),
  data = size_summary3
)

  #model fitting
print(Mmodel2, correlation = T) #high correlation = unstable estimates
plot(Mmodel2) #no pattern = model fit well
plot(residuals(Mmodel2)) #no skew/outliers = fit well
hist(residuals(Mmodel2)) #symmetrically distributed around 0
qqnorm(resid(Mmodel2)) 
qqline(resid(Mmodel2))

summary(Mmodel2) #"dropping 2 column" because no Court 2019 and 2021

  #visualise model prediction
pred <- ggpredict(Mmodel2, terms = c("Site", "Year"))
plot(pred) + theme_minimal()

  #test significance of terms (remove Year:Site interaction)
Mmodel0 <- lmer(
  total_biomass_4th ~ Site + Year + (1 | ARMS_number),
  data = size_summary3
)
anova(Mmodel0, Mmodel2, test = "LRT") #p<0.05, full model better explain my data

#
#
#

#M3.6 ANOVA 
anova(Mmodel2, type = 3) #Type III - satterthwaite
anova_model <- aov(total_biomass_4th ~ Site * Year, data = size_summary3)
summary(anova_model)


#
#
#

#M3.7 Tukey
emmeans(Mmodel2, pairwise ~ Site | Year, adjust = "tukey") #give me the effect of Site within each year
emmeans(Mmodel2, pairwise ~ Year | Site, adjust = "tukey") #can't do

#
#
#

#M3.8 boxplot

  #sum abundance for each replicate
narrow_motile_size <- size_summary3 %>%
  #create unique year-site combo
  mutate(SiteYear = paste(Site, Year, sep = "_"))

#order data
narrow_motile_size$SiteYear <- factor(narrow_motile_size$SiteYear,
                                       levels = c("Anglaise_2019", "Anglaise_2021", "Anglaise_2022", "Moresby_2019", "Moresby_2021", "Moresby_2022", "Coin_2019", "Coin_2021", "Coin_2022", "Court_2022"))
narrow_motile_count <- narrow_motile_size %>%
  arrange(SiteYear)

#plot
ggplot(narrow_motile_size, aes (x = SiteYear, y = total_biomass_4th, fill = as.factor(Year))) +
  geom_boxplot(size = 1)+
  labs(x = "Site-Year", y = "Biomass (g1/4)") +
  theme_minimal()+
  scale_fill_manual(values = colvec2)

#
#
#

# MOTILE-4 ABUNDANCE TURNOVER ---------------------------------------------

#M4.1 - PERMANOVA

  #pivot and filter necessary data
wide_motile_count <- motile_count %>%
  filter(!is.na(MorphoID) & MorphoID != "") %>%
  select(Year, Site, ARMS_number, Abundance_in_image, MorphoID) %>%
  group_by(Year, Site, ARMS_number, MorphoID) %>%
  summarise(Total_abundance = sum(Abundance_in_image, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = MorphoID, 
    values_from = Total_abundance, 
    values_fill = 0)

metadata <- wide_motile_count %>% select(Year, Site, ARMS_number)
abundance <- wide_motile_count %>% select(-Year, -Site, -ARMS_number)

abundance_4th <- abundance^(1/4)
metadata$Year <- factor(metadata$Year)

  #PERMANOVA assumption - homogeneity of mulvariate dispersions
dist_mat1 <- vegdist(abundance_4th, method = "bray")
bd <- betadisper(dist_mat1, group = metadata$Year)
permutest(bd) #check if Site and Year are p >0.05 = good

  #PERMANOVA
adonis2(dist_mat1 ~ Year, data = metadata, permutations = 999, method = "bray") #year only
adonis2(dist_mat1 ~ Site, data = metadata, permutations = 999, method = "bray") #site only
adonis2(dist_mat1 ~ Year * Site, data = metadata, permutations = 999, method = "bray") #year and site


#
#
#

#M4.2 - NMDS

set.seed(10) #reproducible results
abundance_4th_matrix <- as.matrix(abundance_4th) #convert to matrix

nmds = metaMDS(abundance_4th_matrix, k=2, trymax = 100, distance = 'bray')  #calculate Bray-Curtis distances, lower stress = greater fit
stressplot(nmds) #little scatter around line = good

site_scores <- as.data.frame(scores(nmds, display = "sites"))
ordiplot_data <- cbind(metadata, site_scores) 

  #NMDs of sites
hull_dataS <- ordiplot_data %>%
  group_by(Site) %>%
  slice(chull(NMDS1, NMDS2))

colvec1 <- c("Anglaise" = "#D81B60", 
             "Moresby" = "#1E88E5", 
             "Coin" = "#FFC107", 
             "Court" = "#004D40")
pchvec1 <- c("2019" = 21,
             "2021" = 22,
             "2022" = 24)

ordiplot_data$Site <- factor(ordiplot_data$Site, levels = names(colvec1))
unique(ordiplot_data$Site)
names(colvec1)
levels(ordiplot_data$Site)

ggplot(ordiplot_data, aes(x = NMDS1, y = NMDS2, color = Site, fill = Site, shape = as.factor(Year))) +
  stat_ellipse(aes(group = Site, fill = Site),
               geom = "polygon", type = "norm", level = 0.95, alpha = 0.3, color = NA) +
  geom_polygon(data = hull_dataS, aes(group = Site), 
               fill = NA, color = "black", alpha = 1, linetype = "dotted") +
  geom_point(size = 3, stroke = 1, color = "black") +
  scale_shape_manual(values = pchvec1) +
  scale_color_manual(values = colvec1) +
  scale_fill_manual(values = colvec1) +
  theme_minimal()

  #NMDs of years
hull_dataY <- ordiplot_data %>%
  group_by(Year) %>%
  slice(chull(NMDS1, NMDS2))

colvec2 <- c("2019" = "#60FF00",
             "2021" = "#9C9E4C", 
             "2022" = "#2ACDCC")
pchvec2 <- c("Anglaise" = 21, 
             "Moresby" = 22, 
             "Coin" = 24, 
             "Court" = 23)
ordiplot_data$Year <- factor(ordiplot_data$Year, levels = names(colvec2))

ggplot(ordiplot_data, aes(x = NMDS1, y = NMDS2, color = Year, fill = Year, shape = Site)) +
  stat_ellipse(aes(group = Year, fill = Year),
               geom = "polygon", type = "norm", level = 0.95, alpha = 0.3, color = NA) +
  geom_polygon(data = hull_dataY, aes(group = Year), 
               fill = NA, color = "black", alpha = 1, linetype = "dotted") +
  geom_point(size = 3, stroke = 1, color = "black") +
  scale_shape_manual(values = pchvec2) +
  scale_color_manual(values = colvec2) +
  scale_fill_manual(values = colvec2) +
  theme_minimal()

#
#
#

#M4.3 SIMPER

  #create new dataframe sorting specimens by general_name
wide_motile_count2 <- motile_count %>%
  filter(!is.na(MorphoID) & MorphoID != "") %>%
  select(Year, Site, ARMS_number, Abundance_in_image, General_name) %>%
  group_by(Year, Site, ARMS_number, General_name) %>%
  summarise(Total_abundance = sum(Abundance_in_image, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = General_name, 
    values_from = Total_abundance, 
    values_fill = 0)

metadata <- wide_motile_count2 %>% select(Year, Site, ARMS_number)
abundance <- wide_motile_count2 %>% select(-Year, -Site, -ARMS_number)

year <- metadata$Year
site <- metadata$Site

simper1 <- simper(abundance, year, permutations = 999)
summary(simper)

  #stacked barplot of significant groups
  #extract pairwise comparisons
simper_df1 <- bind_rows(
  lapply(names(simper1), function(group) {
    res <- as.data.frame(simper1[[group]])
    res$Comparison <- group
    res$Taxon <- rownames(res)
    rownames(res) <- NULL
    return(res)
  })
)
  #filter for p <0.05
simper_df1 <- simper_df1 %>%
  filter(p < 0.05) %>%
  #create data to plot
  group_by(Comparison, Taxon) %>%
  #"average" = how much each taxon contribute to average dissim between groups
  summarise(Dissimilarity = sum(average), .groups = "drop") %>%
  group_by(Comparison) %>%
  #convert similarity contributions to % within each comparison = sum to 100%
  mutate(Dissimilarity_percent = 100 * Dissimilarity/sum(Dissimilarity)) %>%
  ungroup()


simper2 <- simper(abundance, site, permutations = 999)

  #extract pairwise comparisons
simper_df2 <- bind_rows(
  lapply(names(simper2), function(group) {
    res <- as.data.frame(simper2[[group]])
    res$Comparison <- group
    res$Taxon <- rownames(res)
    rownames(res) <- NULL
    return(res)
  })
)
#filter for p <0.05
simper_df2 <- simper_df2 %>%
  filter(p < 0.05) %>%
  #create data to plot
  group_by(Comparison, Taxon) %>%
  #"average" = how much each taxon contribute to average dissim between groups
  summarise(Dissimilarity = sum(average), .groups = "drop") %>%
  group_by(Comparison) %>%
  #convert similarity contributions to % within each comparison = sum to 100%
  mutate(Dissimilarity_percent = 100 * Dissimilarity/sum(Dissimilarity)) %>%
  ungroup()

write.csv(simper_df)

  #plot SIMPER results as stacked barplot
simper_master <- bind_rows(
  simper_df1 %>% mutate(plot_id = "Plot1"),
  simper_df2 %>% mutate(plot_id = "Plot2")
)

order <- unique(simper_master$Taxon)
base_colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
shared_colors <- colorRampPalette(base_colors)(length(order))
shared_colors <- setNames(shared_colors, order)

  #plot
ggplot(simper_df2, aes(x = Comparison, y = Dissimilarity_percent, fill = Taxon)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = shared_colors, drop = FALSE) +
  theme_minimal()

#
#
#

# MOTILE-5 BIOMASS TURNOVER -----------------------------------------------

# M5.1 PERMANOVA

  #pivot and filter necessary data
wide_motile_size <- motile_size %>%
  filter(!is.na(L_mm), !is.na(W_g)) %>%
  select(Year, Site, ARMS_number, MorphoID, W_g) %>%
  group_by(Year, Site, ARMS_number, MorphoID) %>%
  summarise(Total_biomass = sum(W_g, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = MorphoID, 
    values_from = Total_biomass, 
    values_fill = 0)

metadata <- wide_motile_size %>% select(Year, Site, ARMS_number)
biomass <- wide_motile_size %>% select(-Year, -Site, -ARMS_number)

biomass_4th <- biomass^(1/4)
metadata$Year <- factor(metadata$Year)

  #PERMANOVA assumption - homogeneity of mulvariate dispersions
dist_mat2 <- vegdist(biomass_4th, method = "bray")
bd <- betadisper(dist_mat2, group = metadata$Site)
permutest(bd) #check if Site and Year are p >0.05 = good

  #PERMANOVA
adonis2(dist_mat2 ~ Year, data = metadata, permutations = 999, method = "bray") #year only
adonis2(dist_mat2 ~ Site, data = metadata, permutations = 999, method = "bray") #site only
adonis2(dist_mat2 ~ Year*Site, data = metadata, permutations = 999, method = "bray") #year and site

#
#
#

# M5.2 NMDS

set.seed(20) #reproducible results
biomass_4th_matrix <- as.matrix(biomass_4th)

nmds = metaMDS(biomass_4th_matrix, k=2, trymax = 100, distance = 'bray')  #calculate Bray-Curtis distances, lower stress = greater fit
stressplot(nmds) #little scatter around line = good

site_scores <- as.data.frame(scores(nmds, display = "sites"))
ordiplot_data <- cbind(metadata, site_scores) 

  #NMDS of sites
hull_dataS <- ordiplot_data %>%
  group_by(Site) %>%
  slice(chull(NMDS1, NMDS2))

ordiplot_data$Site <- factor(ordiplot_data$Site, levels = names(colvec1))
unique(ordiplot_data$Site)
names(colvec1)
levels(ordiplot_data$Site)

ggplot(ordiplot_data, aes(x = NMDS1, y = NMDS2, color = Site, fill = Site, shape = as.factor(Year))) +
  stat_ellipse(aes(group = Site, fill = Site),
                   geom = "polygon", type = "norm", level = 0.95, alpha = 0.3, color = NA) +
  geom_polygon(data = hull_dataS, aes(group = Site), 
               fill = NA, color = "black", alpha = 1, linetype = "dotted") +
  geom_point(size = 3, stroke = 1, color = "black") +
  scale_shape_manual(values = pchvec1) +
  scale_color_manual(values = colvec1) +
  scale_fill_manual(values = colvec1) +
  theme_minimal()

  #NMDs of years
hull_dataY <- ordiplot_data %>%
  group_by(Year) %>%
  slice(chull(NMDS1, NMDS2))

ordiplot_data$Year <- factor(ordiplot_data$Year, levels = names(colvec2))

ggplot(ordiplot_data, aes(x = NMDS1, y = NMDS2, color = Year, fill = Year, shape = Site)) +
  stat_ellipse(aes(group = Year, fill = Year),
                   geom = "polygon", type = "norm", level = 0.95, alpha = 0.3, color = NA) +
  geom_polygon(data = hull_dataY, aes(group = Year), 
               fill = NA, color = "black", alpha = 1, linetype = "dotted") +
  geom_point(size = 3, stroke = 1, color = "black") +
  scale_shape_manual(values = pchvec2) +
  scale_color_manual(values = colvec2) +
  scale_fill_manual(values = colvec2) +
  theme_minimal()

#
#
#

# # MOTILE-6 BETA DIVERSITY PARTITIONING --------------------------------------------

#M6.1 - partition b-diversity into turnover and nestedness

  #load wide_motile_count, abundance and metadata

  #convert abundance into P/A matrix
pa <- abundance
pa[pa > 0] <- 1

  #split data by site
sites <- split(pa, metadata$Site)

  #calculate temporal b-diversity (partition across years for each site)   
temporal <- lapply(names(sites), function(site_name) {
  subdata <- sites[[site_name]]
  
  comm <- subdata[, !colnames(subdata) %in% c("Site", "Year")]
  beta <- beta.multi(comm, index.family = "jaccard")
  
  data.frame(
    Site = site_name,
    turnover = beta$beta.JTU,
    nestedness = beta$beta.JNE,
    total = beta$beta.JAC
  )
})

temporal <- bind_rows(temporal)

  #split data by year
years <- split(pa, metadata$Year)

  #calculate spatial b-diversity (partition across sites within each year)   
spatial <- lapply(names(years), function(year_name) {
  subdata <- years[[year_name]]
  
  comm <- subdata[, !colnames(subdata) %in% c("Site", "Year")]
  beta <- beta.multi(comm, index.family = "jaccard")
  
  data.frame(
    Site = year_name,
    turnover = beta$beta.JTU,
    nestedness = beta$beta.JNE,
    total = beta$beta.JAC
  )
})

spatial <- bind_rows(spatial)

#
#
#

# M6.2 - UpSet plot - unique species per site/year

  # Example site-year column
motile_count <- motile_count %>%
  mutate(SiteYear = paste(Site, Year, sep = "_"))

  # Build presence/absence matrix
binary_matrix <- motile_count %>%
  group_by(MorphoID, SiteYear) %>%
  summarise(Present = ifelse(sum(Abundance_in_image) > 0, 1, 0), .groups = "drop") %>%
  pivot_wider(names_from = SiteYear, values_from = Present, values_fill = 0)

  #reorder site_years
ordered_SiteYear <- rev(c("Anglaise_2019", "Anglaise_2021", "Anglaise_2022", "Moresby_2019", "Moresby_2021", "Moresby_2022", "Coin_2019", "Coin_2021", "Coin_2022", "Court_2022"))

  #convert PA to logical  
binary_matrix[-1] <- lapply(binary_matrix[-1], as.logical)
binary_matrix[is.na(binary_matrix)] <- FALSE

  #create new column to represent each intersection
binary_matrix <- binary_matrix %>%
  mutate(
    intersection_pattern = apply(select(., -MorphoID), 1, function(row) paste(which(row), collapse = "-"))
  )

  # count how many species per unique combination of SiteYears
intersection_counts <- binary_matrix %>%
  count(intersection_pattern, name = "pattern_count")
  # join this count back to the data
binary_matrix <- binary_matrix %>%
  left_join(intersection_counts, by = "intersection_pattern")
  # filter out patterns with only 1 species
binary_matrix_filtered <- binary_matrix %>%
  filter(pattern_count > 1) %>%
  select(-intersection_pattern, -pattern_count)

  # extract the ordered SiteYear names (excluding MorphoID column)
siteyear <- setdiff(names(binary_matrix_filtered), "MorphoID")

  #classify morphosp as unique or repeated
binary_matrix_filtered <- binary_matrix_filtered %>%
  mutate(
    species_type = ifelse(rowSums(select(., -MorphoID)) == 1, "Unique", "Repeated")
  )

  #reorder intersection columns based on years of site
binary_matrix <- binary_matrix[, c("MorphoID", ordered_SiteYear)]


# #reporting numbers before plotting: ----------------------------------------------------------
  #total number of site-year intersections with only 1 species (i filtered out)
a <- intersection_counts %>%
  filter(pattern_count == 1)
nrow(a)

  #morphosp that appear most frequently across intersections (number of TRUE for each morphosp)
b <- binary_matrix_filtered %>%
  rowwise() %>%
  mutate(n = sum(c_across(-MorphoID))) %>%
  ungroup() %>%
  #arrange(desc(n))
  filter(n == 1)

  #the taxonomic composition of these unique sp
c <- b$MorphoID

filter <- motile_count %>%
  distinct(MorphoID, .keep_all = TRUE)

filter <- filter %>%
  filter(MorphoID %in% c)

filter %>% 
  group_by(Phyla, General_name) %>%
  summarise(n_morphosp = n_distinct(MorphoID),
            .groups = 'drop') -> PD4

webr::PieDonut(PD4, 
               aes(Phyla, General_name, count = n_morphosp), 
               ratioByGroup = FALSE, #donut slices based on total number of specimens, not relative to phyla
               showRatioThreshold = 0, #only show labels for segments >XX% of total
               r0=0.2, 
               r1=0.6,
               donutLabelSize = 3,
               pieLabelSize = 3,
               labelpositionThreshold = 0.02, #if slice is >XX% of total, label will go inside slice
               showPieName = FALSE)

  #morphosp only in one site-year combo + have 1 specimen
species_presence_counts <- rowSums(binary_matrix[, -1])  # Assuming first column is species ID
filtered_matrix <- binary_matrix[species_presence_counts == 1, ] # Keep only species present in more than one site-year
filtered_count <- motile_count %>%
  group_by(MorphoID,SiteYear)%>%
  summarise(total = sum(Abundance_in_image), .groups = "drop") %>%
  pivot_wider(names_from = SiteYear, values_from = total, values_fill = 0) %>%
  rowwise() %>%
  filter(
    sum(c_across(-MorphoID) > 0) == 1, #only 1 non-zero entry
    sum(c_across(-MorphoID) == 1) == 1 #entry is exactly 1
  ) %>%
  ungroup()

# #plot repeated species as UpSet plot seperately ------------------------------

binary_matrix_unique <- binary_matrix_filtered %>%
  filter(species_type == "Unique")

unique_species_count <- binary_matrix_unique %>%
  select(-species_type) %>%
  mutate_all(as.numeric) %>%
  colSums()
unique_species_count

binary_matrix_shared <- binary_matrix_filtered %>%
  filter(species_type == "Repeated")
shared_species_count <- binary_matrix_shared %>%
  select(-species_type) %>%
  mutate_all(as.numeric) %>%
  colSums()
shared_species_count

total_species_count <- unique_species_count + shared_species_count
total_species_count

  # Plot
upset(
  binary_matrix_shared,
  intersect = ordered_SiteYear,
  name = "Species Presence",
  base_annotations = list(
    'Intersection size' = intersection_size() +
      aes(fill = species_type) +
      scale_fill_manual(values = c("Repeated" = "grey50"))
  ),
  set_sizes = upset_set_size() +
    geom_bar(aes(fill = species_type)) +
    scale_fill_manual(values = c("Repeated" = "grey50")) +
    labs(fill = "Species Type"),
  sort_sets = FALSE,
  sort_intersections_by = NULL
)

# ------------------------------------------------------------------------------

# MOTILE&SESSILE MANTEL (load M4.1 wide_motile_count/dist_mat1)

  #lower sampling resolution for sessiles/aggregate data so 1 ARMS unit = 1 community
sessile_percent_agg <- sessile_percent %>%
  group_by(Date, Site, ARMS_number) %>%
  summarise(across(all_of(sessile), sum, na.rm = TRUE), .groups = "drop")

  #turn into a matrix of abundances only
sessile_percent_agg <- sessile_percent_agg %>% select(-Date, -Site, -ARMS_number)

dist_mat4 <- vegdist(sessile_percent_agg, method = "bray")
mantel(dist_mat1, dist_mat4, method = "spearman")
