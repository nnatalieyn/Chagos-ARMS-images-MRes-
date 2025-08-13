# SESSILE-0 LOAD DIRECTORY --------------------------------------------------------

#0.1 packages (best to load packages from R_script_motile.R too)
library(dplyr)          #data wrangling
library(tidyr)          #data wrangling
library(terra)          #to input CRS
library(sf)             #manipulate geojson
library(lubridate)      #convert dates
library(webr)           #for piedonut
library(ggplot2)        #load in this order for code to run
library(lme4)           #for GLMM
library(ggeffects)      #visualise model predictions
library(car)            #type II/III test
library(emmeans)        #for GLLM tukey
library(purrr)          #for pearson p-value
library(RColorBrewer)   #color palette
library(glmmTMB)        #for nbinom2

#0.2 data

  #A freq counts of all CN labels (total_count)
pointcount <- read.csv("CN_pointcount_accuracy.csv", stringsAsFactors = TRUE)
pointcount$Date <- year(dmy(pointcount$Date))
pointcount <- pointcount %>%
  mutate(
    Name = as.character(Name),
    Date = as.factor(Date),
    ARMS_number = as.character(ARMS_number),
    Aux_5 = as.character(Aux_5)
  )

total_count <- pointcount[, c(1:7,12)]
total_count <- total_count %>%
  group_by(Name, Date, Site, ARMS_number, top.bottom, open.closed, Aux_5, Label.code) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = Label.code,
    values_from = count,
    values_fill = 0 #fill missing label counts with zero
  )
total_count$available_count <- 225 - total_count$`_UNAV`
total_count$recruited_count <- 225 - total_count$`_UNAV` - total_count$`_NR`
total_count$rec_on_ava_percent <- (total_count$recruited_count/total_count$available_count) * 100

  #B transformed % cover of sessile CN labels on ava space (sessile_percent)
sessile <- c("_BI", "_BREN", "_BRY", "_CAWT", "_CCA", "_CO", "_GREN", "_RDEN", "_SOWT", "_SP", "Turf", "_FORM", "_GRUP", "_RDUP", "_TUNC", "_BRUP", "Tun_sol", "_BGMA", "Anemone", "_OCTO", "_HYD", "_CMOR", "_ZO")

sessile_percent <- total_count %>%
  rowwise() %>%
  mutate(
    across(all_of(sessile), ~ (.x / available_count) *100)
  )
sessile_percent <- sessile_percent[, c(1:14, 16:18, 20, 24:27, 30, 32, 33, 35:37, 39:40)]

  #C georeferenced plate images (.tiff) (for CN-T to export segmentation data)
input <- "plate_tiff_raw" #locate folder of ungeoreferenced plate images
output <- "plate_tiff_crs" #folder for georeferenced plate images

files <- list.files(input, pattern = "\\.tif$", full.names = TRUE)

for (input_path in files) { #input CRS info
  r <- rast(input_path) #load img
  ext(r) <- c(0, 225, 0, 225) #set extent, each image  covers 0.225m x 0.225m
  crs(r) <- "EPSG:3857" #set CRS: Web Mercator is a dummy meter CRS
  output_path <- file.path(output, basename(input_path))
  writeRaster(r, output_path, overwrite = TRUE)
  cat("Georeferenced:", basename(input_path), "\n")
}

  #D segmentation polygons (.geojson)
segments <- st_read("CNT_area.geojson")
print(st_crs(segments))  #check if CRS correct
print(st_bbox(segments)) #check if between 0-0.22 m 
st_crs(segments) <- NA

segments$valid <- st_is_valid(segments) #check if polygons have geometry issues
table(segments$valid)
segments <- segments %>%
  mutate(geometry = st_buffer(geometry, 0))

  #E CN accuracy  
CNaccuracy <- pointcount[, c(1:7,12, 16:25)]

  #F CN-T accuracy (.csv)
CNTaccuracy <- read.csv("CNT_accuracy.csv")

#
#
#

# SESSILE-1 SUMMARY -------------------------------------------------------

#S1.1 distribution of all point annotations
totalpercent <- pointcount %>%
  count(Label.code) %>%
  mutate(percentage = (n/114300)*100)

#S1.2 distribution of live sessiles on avaliable
sessiledis <- colSums(sessile_percent[, sessile], na.rm = TRUE)
totalpercent2 <- (sessiledis / sum(sessiledis)) * 100
print(totalpercent2)

#S1.3 pie chart of % cover of live sessile CN labels
  #convert to dataframe
totalpercent2 <- data.frame(
  Label = names(totalpercent2),
  TotalPercent = as.numeric(totalpercent2)
)
  #add phylum info
totalpercent2 <- totalpercent2 %>%
  mutate(Phylum = case_when(
    Label %in% c("_SP") ~ "Porifera",
    Label %in% c("_BI") ~ "Mollusca",
    Label %in% c("_BRY") ~ "Bryozoa",
    Label %in% c("_CAWT", "_SOWT") ~ "Annelida",
    Label %in% c("_CMOR", "_CO", "Anemone", "_OCTO", "_HYD", "_ZO") ~ "Cnidaria",
    Label %in% c("_FORM") ~ "Retaria",
    Label %in% c("_TUNC", "Tun_sol") ~ "Chordata",
    Label %in% c("_GREN", "_GRUP") ~ "Chlorophyta",
    Label %in% c("_CCA", "_RDEN", "_RDUP") ~ "Rhodophyta",
    Label %in% c("_BREN", "_BRUP") ~ "Ochrophyta",
    Label %in% c("_BGMA") ~ "Cyanobacteria", 
    Label %in% c("Turf") ~ "Turf Algae"
  ))
totalpercent2$Phylum <- factor(totalpercent2$Phylum,
                               levels = c("Porifera", "Mollusca", "Bryozoa", "Annelida", "Cnidaria", "Chordata", "Retaria", "Chlorophyta", "Rhodophyta", "Ochrophyta", "Cyanobacteria", "Turf Algae"))

  #plot
webr::PieDonut(totalpercent2,
               aes(Phylum, Label, count = TotalPercent),
               ratioByGroup = FALSE, #donut slices based on total number of specimens, not relative to phyla
               showRatioThreshold = 0, #only show labels for segments >XX% of total
               r0=0.2, 
               r1=0.6,
               donutLabelSize = 3,
               pieLabelSize = 3,
               labelpositionThreshold = 0.01, #if slice is >XX% of total, label will go inside slice
               showPieName = FALSE)

#
#
#

# SESSILE-2 ABUNDANCE (CN) ------------------------------------------------

#S2.1 total live sessile cover per site and year and face

abundance_summary3 <- sessile_percent %>%
  #sum all the sessile percentage covers for each image
  rowwise() %>%
  mutate(totalpercent3 = sum(c_across(all_of(sessile)), na.rm = TRUE)) %>%
  ungroup() %>%
  #create new summary dataframe
  select(-all_of(sessile))

#
#
#

#S.2 GLM of total abundance

  #explore data structure
summary(abundance_summary3)
hist(abundance_summary3$totalpercent3)
  #overdispersed = var > mean
mean(abundance_summary3$totalpercent3)
var(abundance_summary3$totalpercent3)

  #create model
Smodel1 <- glmmTMB(
  totalpercent3 ~ Site * Date + top.bottom + open.closed + (1 | Aux_5) + (1 | ARMS_number),
  family = nbinom2,
  data = abundance_summary3
)

  #"4 not defined because of singularities" Courts Knoll
table(abundance_summary3$Site, abundance_summary3$Date, abundance_summary3$top.bottom)

  #model fitting
print(Smodel1, correlation = T) #high correlation = unstable estimates
plot(resid(Smodel1), fitted(Smodel1)) #no skew/outliers = fit well
hist(residuals(Smodel1)) #symmetrically distributed around 0
qqnorm(resid(Smodel1)) 
qqline(resid(Smodel1))

summary(Smodel1)

  #visualise model prediction
pred <- ggpredict(Smodel1, terms = c("Site", "Date", "top.bottom", "open.closed"))
plot(pred) + theme_minimal()

  #likelihood ratio test (whether random effect improve model)
Smodel0 <- glmmTMB(
  totalpercent3 ~ Site * Date + top.bottom + open.closed,
  family = nbinom2,
  data = abundance_summary3
)
anova(Smodel0, Smodel1, test = "LRT") #p>0.05, NOT sig

#
#
#

#S2.3 ANODEV (effect of site and year)
  
  #can't run anodev with undefined var
subset_abundance_summary3 <- abundance_summary3[abundance_summary3$Site != "Courts Knoll",]
subset_abundance_summary3 <- droplevels(subset_abundance_summary3)
table(subset_abundance_summary3$Site, subset_abundance_summary3$Date, subset_abundance_summary3$top.bottom)

Smodel2 <- glmmTMB(
  totalpercent3 ~ Site * Date + top.bottom + open.closed + (1 | ARMS_number) + (1 | Aux_5),
  family = nbinom2,
  data = subset_abundance_summary3
)

Anova(Smodel2, type = 3) #Type III - test each variable and interaction

#
#
#

#S2.4 Tukey
emmeans(Smodel1, pairwise ~ Site | Date, adjust = "tukey") #effect of site within each year
emmeans(Smodel1, pairwise ~ Date | Site, adjust = "tukey") #effect of year within each site

#
#
#

#S2.5 barplot

  #site-year combo
narrow_sessile_percent <- sessile_percent %>%
  mutate(SiteYear = paste(Site, Date, sep = "_")) %>%
  pivot_longer(
    cols = all_of(sessile),
    names_to = "Label",
    values_to = "Abundance"
  )
narrow_sessile_total <- narrow_sessile_percent %>%
  group_by(SiteYear, Site, Date, ARMS_number, Name) %>%
  summarise(
    total = sum(Abundance),
    .groups = ("drop")
  ) %>%
  ungroup() %>%
  group_by(SiteYear, Site, Date, ARMS_number) %>%
  summarise(
    total = mean(total),
    .groups = ("drop")
  )

  #order data
narrow_sessile_total$SiteYear <- factor(narrow_sessile_total$SiteYear,
                                          levels = c("Ile Anglaise_2019", "Ile Anglaise_2021", "Ile Anglaise_2022", "Moresby_2019", "Moresby_2021", "Moresby_2022", "Ile du Coin_2019", "Ile du Coin_2021", "Ile du Coin_2022", "Courts Knoll_2022"))
narrow_sessile_total <- narrow_sessile_total %>%
  arrange(SiteYear)

  #plot
ggplot(narrow_sessile_total, aes (x = SiteYear, y = total, fill = as.factor(Date))) +
  geom_boxplot(size = 1)+
  scale_y_continuous(limits = c(40,100)) + 
  labs(x = "Site-Year", y = "Average percentage cover (%)") +
  theme_minimal()+
  scale_fill_manual(values = colvec2)

#
#
#

# SESSILE-3 ABUNDANCE + SURFACE AREA (CN-T) ----------------------------------------------

# S3.1 - area per group
  #calculate area in mm^2
segments <- segments %>%
  mutate(area_mm2 = st_area(.) %>% as.numeric())

  #summarise total area per class
area <- segments %>%
  group_by(source_image, label) %>%
  summarise(total_area_mm2 = sum(area_mm2, na.rm = TRUE))

  #double check units are correctly scaled, 48400mm2 = 220mm x 220mm
  #check how much space is annotated 
TOTAL_area <- segments %>%
  group_by(source_image) %>%
  summarise(total_area_mm2 = sum(area_mm2, na.rm = TRUE))
sum(TOTAL_area$total_area_mm2)

#
#
#

#S3.2 - percentage cover per sessile group

  #standardise to available area
NR_area <- area %>%
  filter(label == "_NR") %>%
  select(source_image, NR_area_mm2 = total_area_mm2)

UNAV_area <- area %>%
  filter(label == "_UNAV") %>%
  select(source_image, UNAV_area_mm2 = total_area_mm2)

TOTAL_area <- TOTAL_area %>%
  left_join(st_drop_geometry(UNAV_area), by = "source_image") %>%
  left_join(st_drop_geometry(NR_area), by = "source_image") %>%
  mutate(
    available_area = total_area_mm2 - UNAV_area_mm2,
    recruited_area = total_area_mm2 - UNAV_area_mm2 - NR_area_mm2,
    rec_on_ava_area = (recruited_area/available_area) * 100
    )

  #calculate percentage cover
area <- area %>%
  left_join(st_drop_geometry(TOTAL_area %>% select(source_image, available_area)),
            by = "source_image")

sessile_percent2 <- area %>%
  filter(label %in% sessile) %>%
  mutate(percent = (total_area_mm2 / available_area) *100) %>%
  select(source_image, label, percent) %>%
  st_drop_geometry

#
#
#

# SESSILE-4 ABUNDANCE TURNOVER --------------------------------------------

#S4.1 - PERMANOVA

  #filter necessary data (group composition by ARMS)
metadata <- sessile_percent %>% select(Name, Date, Site, ARMS_number, top.bottom, open.closed, Aux_5)
abundance <- sessile_percent %>% select(-Name, -Date, -Site, -ARMS_number, -top.bottom, -open.closed, -Aux_5)

abundance_4th <- abundance^(1/4)
metadata$Year <- factor(metadata$Date)

  #PERMANOVA assumption - homogeneity of multivariate dispersions
dist_mat3 <- vegdist(abundance_4th, method = "bray")
bd <- betadisper(dist_mat3, group = metadata$Year)
permutest(bd) #check if Site and Year are p >0.05 = good

  #PERMANOVA
adonis2(dist_mat3 ~ Date, data = metadata, permutations = 999, method = "bray") #year only
adonis2(dist_mat3 ~ Site, data = metadata, permutations = 999, method = "bray") #site only
adonis2(dist_mat3 ~ top.bottom, data = metadata, permutations = 999, method = "bray") 
adonis2(dist_mat3 ~ open.closed, data = metadata, permutations = 999, method = "bray")
adonis2(dist_mat3 ~ Date * Site * top.bottom * open.closed, data = metadata, permutations = 999, method = "bray") 

#
#
#

#S4.2 - NMDS
abundance_summary4 <- sessile_percent %>%
  group_by(Site, Date, ARMS_number) %>%
  summarise(across(all_of(sessile), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

metadata <- abundance_summary4 %>% select(Date, Site, ARMS_number)
abundance <- abundance_summary4 %>% select(-Date, -Site, -ARMS_number)

abundance_4th <- abundance^(1/4)
metadata$Year <- factor(metadata$Date)

set.seed(30) #reproducible results
abundance_4th_matrix <- as.matrix(abundance_4th) #convert to matrix

nmds = metaMDS(abundance_4th_matrix, k=2, trymax = 100, distance = 'bray')  #calculate Bray-Curtis distances, lower stress = greater fit
stressplot(nmds) #little scatter around line = good

site_scores <- as.data.frame(scores(nmds, display = "sites"))
ordiplot_data <- cbind(metadata, site_scores) 

  #NMDS of sites
hull_dataS <- ordiplot_data %>%
  group_by(Site) %>%
  slice(chull(NMDS1, NMDS2))

colvec1 <- c("Ile Anglaise" = "#D81B60", 
             "Moresby" = "#1E88E5", 
             "Ile du Coin" = "#FFC107", 
             "Courts Knoll" = "#004D40")  
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
  group_by(Date) %>%
  slice(chull(NMDS1, NMDS2))

pchvec2 <- c("Ile Anglaise" = 21, 
             "Moresby" = 22, 
             "Ile du Coin" = 24, 
             "Courts Knoll" = 23)  
ordiplot_data$Date <- factor(ordiplot_data$Date, levels = names(colvec2))

ggplot(ordiplot_data, aes(x = NMDS1, y = NMDS2, color = Date, fill = Date, shape = Site)) +
  stat_ellipse(aes(group = Date, fill = Date),
               geom = "polygon", type = "norm", level = 0.95, alpha = 0.3, color = NA) +
  geom_polygon(data = hull_dataY, aes(group = Date), 
               fill = NA, color = "black", alpha = 1, linetype = "dotted") +
  geom_point(size = 3, stroke = 1, color = "black") +
  scale_shape_manual(values = pchvec2) +
  scale_color_manual(values = colvec2) +
  scale_fill_manual(values = colvec2) +
  theme_minimal()

#
#
#

#S4.3 SIMPER
year <- metadata$Date
site <- metadata$Site

simper <- simper(abundance, year, permutations = 999)
summary(simper)

  #stacked barplot of significant groups

  #extract pairwise comparisons
simper_df <- bind_rows(
  lapply(names(simper), function(group) {
    res <- as.data.frame(simper[[group]])
    res$Comparison <- group
    res$Taxon <- rownames(res)
    rownames(res) <- NULL
    return(res)
  })
)
  #filter for p <0.05
simper_df <- simper_df %>%
  filter(p < 0.05) %>%
  #create data to plot
  group_by(Comparison, Taxon) %>%
  summarise(Dissimilarity = sum(average), .groups = "drop") %>%
  group_by(Comparison) %>%
  mutate(Dissimilarity_percent = 100 * Dissimilarity/sum(Dissimilarity)) %>%
  ungroup()

write.csv(simper_df)

  #use the same color palette as cn/cnt barplot

order <- c("_SP", "_BI", "_BRY", "Turf", "_BGMA", "_BRUP", "_BREN", "_CCA", "_RDEN", "_RDUP", "_FORM", "_CMOR", "_CO", "Anemone", "_OCTO", "_HYD", "_ZO", "_CAWT", "_SOWT", "_TUNC", "Tun_sol", "_GREN", "_GRUP")
base_colors <- c(rcartocolor::carto_pal(12, "Safe"), rcartocolor::carto_pal(11, "Antique"))
palette_colors <- colorRampPalette(base_colors)(length(order))
# Create named palette
label_palette <- setNames(palette_colors, order)
  
  #plot
ggplot(simper_df, aes(x = Comparison, y = Dissimilarity_percent, fill = Taxon)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = label_palette)+
  theme_minimal()

#
#
#

# SESSILE-5 ABUNDANCE + ACCURACY COMPARISON (CN vs CNT) -------------------------------------------

# S5.0 create abundance combo dataset for CN and CNT
sessile_percent3 <- subset(narrow_sessile_percent, Name %in% c( 
                           "Steyaert_ANG_ARMS2_180301_210419_3B.JPG", 
                           "Steyaert_ANG_ARMS2_180301_210419_7B.JPG",
                           "Steyaert_ANG_ARMS4_180301_190310_9T.JPG",
                           "Steyaert_COIN_ARMS2_180301_190305_4B.JPG",
                           "Steyaert_COIN_ARMS7_180301_190321_1B.JPG",
                           "Steyaert_COIN_ARMS7_180301_210429_2B.JPG",
                           "Steyaert_COIN_ARMS8_180301_221019_6T.JPG",
                           "Steyaert_MOR_ARMS4_180301_210425_3B.JPG",
                           "Steyaert_MOR_ARMS6_180301_221029_7B.JPG"))

  #remove image extensions so image names are the same
sessile_percent3 <- sessile_percent3 %>%
  mutate(Name = gsub("\\.JPG$", "", Name, ignore.case = TRUE))
sessile_percent2 <- sessile_percent2 %>%
  mutate(source_image = gsub("\\.tif$", "", source_image, ignore.case = TRUE))

sessile_percent3 <- sessile_percent3 %>%
  left_join(sessile_percent2 %>% select(source_image, label, percent),
            by = c("Name" = "source_image", "Label" = "label")
            )
sessile_percent3[is.na(sessile_percent3)] <- 0 #replace labels with NA with 0

#
#
#

# S5.1 pearson correlation

sessile_percent3 %>%
  group_by(Label) %>%
  #filter groups with 0 annotations across both methods 
  filter(sd(Abundance, na.rm = TRUE) > 0,
         sd(percent, na.rm = TRUE)> 0) %>%
  nest() %>%
  mutate(
    pearson = map(data, ~cor.test(~ Abundance + percent, data = ., method = "pearson")),
    r = map_dbl(pearson, ~ .x$estimate),
    p = map_dbl(pearson, ~ .x$p.value)
  ) %>%
  select(Label, r, p)

#
#
#

# S5.2 create accuracy combo dataset for CN and CNT
CNaccuracy <- subset(CNaccuracy, Name %in% c( 
  "Steyaert_ANG_ARMS2_180301_210419_3B.JPG", 
  "Steyaert_ANG_ARMS2_180301_210419_7B.JPG",
  "Steyaert_ANG_ARMS4_180301_190310_9T.JPG",
  "Steyaert_COIN_ARMS2_180301_190305_4B.JPG",
  "Steyaert_COIN_ARMS7_180301_190321_1B.JPG",
  "Steyaert_COIN_ARMS7_180301_210429_2B.JPG",
  "Steyaert_COIN_ARMS8_180301_221019_6T.JPG",
  "Steyaert_MOR_ARMS4_180301_210425_3B.JPG",
  "Steyaert_MOR_ARMS6_180301_221029_7B.JPG"))

CNaccuracy <- CNaccuracy %>%
  mutate(across(starts_with("Machine."), as.character))

narrow_CNaccuracy <- bind_rows(
  CNaccuracy %>% transmute(Name, rank = 1,
                           confidence = `Machine.confidence.1`,
                           suggestion = `Machine.suggestion.1`),
  CNaccuracy %>% transmute(Name, rank = 2,
                           confidence = `Machine.confidence.2`,
                           suggestion = `Machine.suggestion.2`),
  CNaccuracy %>% transmute(Name, rank = 3,
                           confidence = `Machine.confidence.3`,
                           suggestion = `Machine.suggestion.3`),
  CNaccuracy %>% transmute(Name, rank = 4,
                           confidence = `Machine.confidence.4`,
                           suggestion = `Machine.suggestion.4`),
  CNaccuracy %>% transmute(Name, rank = 5,
                           confidence = `Machine.confidence.5`,
                           suggestion = `Machine.suggestion.5`)
) %>%
  mutate(
    rank = as.integer(rank),
    confidence = as.numeric(confidence),
    suggestion = trimws(suggestion)
  ) %>%
  filter(!is.na(confidence)) %>% #exclude NAs
  mutate(method = "CN") # add method

CNTaccuracy <- CNTaccuracy %>%
  mutate(across(starts_with("Machine."), as.character))

narrow_CNTaccuracy <- bind_rows(
  CNTaccuracy %>% transmute(Name, rank = 1,
                        confidence = `Machine.confidence.1`,
                        suggestion = `Machine.suggestion.1`),
  CNTaccuracy %>% transmute(Name, rank = 2,
                        confidence = `Machine.confidence.2`,
                        suggestion = `Machine.suggestion.2`),
  CNTaccuracy %>% transmute(Name, rank = 3,
                        confidence = `Machine.confidence.3`,
                        suggestion = `Machine.suggestion.3`),
  CNTaccuracy %>% transmute(Name, rank = 4,
                        confidence = `Machine.confidence.4`,
                        suggestion = `Machine.suggestion.4`),
  CNTaccuracy %>% transmute(Name, rank = 5,
                        confidence = `Machine.confidence.5`,
                        suggestion = `Machine.suggestion.5`)
) %>%
  mutate(
    rank = as.integer(rank),
    confidence = as.numeric(confidence),
    suggestion = trimws(suggestion)
  ) %>%
  filter(!is.na(confidence)) %>% # exclude all NAs = manually drawn with SAM
  mutate(confidence = confidence * 100) %>% # transform %%
  mutate(method = "CNT") # add method

combined_accuracy <- bind_rows(narrow_CNaccuracy, narrow_CNTaccuracy)
combined_accuracy <- combined_accuracy %>% filter(rank == 1) #exclude ALL 2th-5th suggestions

#
#
#

# S5.3 Wilcoxon rank-sum test
wilcox <- combined_accuracy %>%
  group_by(suggestion) %>%
  filter(n_distinct(method) == 2) %>% #keep only labels seen in both methods
  nest() %>%
  mutate(
    wilcox_test = map(data, ~ wilcox.test(confidence ~ method, data = .x)),
    p_value = map_dbl(wilcox_test, ~ .x$p.value),
    CN_mean = map_dbl(data, ~ mean(.x$confidence[.x$method == "CN"], na.rm = TRUE)),
    CNT_mean = map_dbl(data, ~ mean(.x$confidence[.x$method == "CNT"], na.rm = TRUE))
  ) %>%
  select(suggestion, CN_mean, CNT_mean, p_value)
wilcox

#
#
#

#S5.4 stacked bar plot of CN and CN-T

sessile_percent3 <- sessile_percent3 %>%
  mutate(Name = case_when(
    Name == "Steyaert_ANG_ARMS2_180301_210419_3B" ~ "1",
    Name == "Steyaert_ANG_ARMS2_180301_210419_7B" ~ "2",
    Name == "Steyaert_ANG_ARMS4_180301_190310_9T" ~ "3",
    Name == "Steyaert_COIN_ARMS2_180301_190305_4B" ~ "4",
    Name == "Steyaert_COIN_ARMS7_180301_190321_1B" ~ "5",
    Name == "Steyaert_COIN_ARMS7_180301_210429_2B" ~ "6",
    Name == "Steyaert_COIN_ARMS8_180301_221019_6T" ~ "7",
    Name == "Steyaert_MOR_ARMS4_180301_210425_3B" ~ "8",
    Name == "Steyaert_MOR_ARMS6_180301_221029_7B" ~ "9"
    ))

errors4 <- sessile_percent3 %>%
  group_by(Name) %>%
  summarise(
    meanAb = mean(Abundance, na.rm = TRUE),
    seAb = sd(Abundance, na.rm = TRUE) / sqrt(n()),
    meanPe = mean(percent, na.rm = TRUE),
    sePe = sd(percent, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

  #mean abundance per taxo group per photo
bardata4 <- sessile_percent3 %>%
  group_by(Name, Label) %>%
  summarise(
    Abundance = mean(Abundance, na.rm = TRUE), 
    Percent = mean(percent, na.rm = TRUE),
    .groups = "drop"
    )

bar_heights4 <- bardata4 %>%
  group_by(Name) %>%
  summarise(
    total_meanAb = sum(Abundance),
    total_meanPe = sum(Percent),
    .groups = "drop")

errors4 <- errors4 %>%
  left_join(bar_heights4, by = "Name")

  #convert data to long format
bardata4 <- bardata4 %>%
  pivot_longer(cols = c(Abundance, Percent),
               names_to = "Measure",
               values_to = "Mean")
errors4 <- errors4 %>%
  pivot_longer(
    cols = c(meanAb, seAb, meanPe, sePe, total_meanAb, total_meanPe),
    names_to = c(".value", "Measure"),
    names_pattern = "(.*)(Ab|Pe)"
  ) %>%
  mutate(
    Measure = case_when(
      Measure == "Ab" ~ "Abundance",
      Measure == "Pe" ~ "Percent",
      TRUE ~ Measure
    )
  ) %>%
  mutate(Photo_Measure = paste(Name, Measure, sep = "_"))
    

bardata4 <- bardata4 %>%
  mutate(Photo_Measure = paste(Name, Measure, sep = "_"))

  #define order of specimens
order <- c("_SP", "_BI", "_BRY", "Turf", "_BGMA", "_BRUP", "_BREN", "_CCA", "_RDEN", "_RDUP", "_FORM", "_CMOR", "_CO", "Anemone", "_OCTO", "_HYD", "_ZO", "_CAWT", "_SOWT", "_TUNC", "Tun_sol", "_GREN", "_GRUP")
bardata4$Label <- factor(bardata4$Label, levels = order)

  # Choose a palette with enough colors or expand manually
#n <- length(unique(bardata4$Label))
#label_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(n)

base_colors <- c(rcartocolor::carto_pal(12, "Safe"), rcartocolor::carto_pal(11, "Antique"))
palette_colors <- colorRampPalette(base_colors)(length(order))
# Create named palette
label_palette <- setNames(palette_colors, order)

  #plot CN AND CN-T barplot
ggplot(data = bardata4, mapping = aes(x = Photo_Measure, y = Mean, fill = Label)) +
  geom_bar(stat = "identity") +
  geom_errorbar(
    data = errors4,
    aes(x = Photo_Measure, ymin = total_mean - se, ymax = total_mean + se),
    inherit.aes = FALSE,
    width = 0.4,
    color = "black"
  ) +
  scale_fill_manual(values = label_palette) +
  theme_minimal()

#
#
#

# -------------------------------------------------------------------------

# MOTILE&SESSILE MANTEL

  #lower sampling resolution for sessiles/aggregate data so 1 ARMS unit = 1 community
sessile_percent_agg <- sessile_percent %>%
  group_by(Date, Site, ARMS_number) %>%
  summarise(across(all_of(sessile), sum, na.rm = TRUE), .groups = "drop")

  #turn into a matrix of abundances only
sessile_percent_agg <- sessile_percent_agg %>% select(-Date, -Site, -ARMS_number)

dist_mat4 <- vegdist(sessile_percent_agg, method = "bray")
mantel(dist_mat1, dist_mat4, method = "spearman")
