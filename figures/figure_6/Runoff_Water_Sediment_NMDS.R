library(phyloseq)
library(tidyverse)
library(RColorBrewer)

theme_set(theme_bw())

phy <- readRDS("Armstrong_Rainfall_Phy.RDS")

runoff_water_sediment_manure <- subset_samples(phy, matrix != "soil")

all_data <- data.frame(sample_data(runoff_water_sediment_manure)) %>%
  mutate(Matrix = case_when(matrix %in% c("manure") ~ "Manure",
                                 matrix %in% c("runoff_water") ~ "Runoff Water",
                                 matrix %in% c("runoff_sediment") ~ "Runoff Sediment")) %>%
  mutate(Treatment = case_when(treatment %in% c("ACM") ~ "No Strip + Manure",
                             treatment %in% c("ACS") ~ "Strip + No Manure",
                             treatment %in% c("ACSM") ~ "Strip + Manure",
                             treatment %in% c("AM") ~ "Manure"))

sample_names(runoff_water_sediment_manure)

rownames(all_data) <- all_data$sample_name

sample_data(runoff_water_sediment_manure) <- all_data

GP.ord.test <- ordinate(runoff_water_sediment_manure, "NMDS", "bray")

runoff_water_sediment_sample_plot = plot_ordination(runoff_water_sediment_manure, GP.ord.test, 
                                                 type = "samples", color = "Treatment",
                                                 shape = "Matrix") +
  guides(shape = guide_legend(order = 1, reverse = TRUE), color = guide_legend(order = 2, reverse = TRUE))

runoff_water_sediment_NMDS <- runoff_water_sediment_sample_plot + geom_point(size = 2) + stat_ellipse(type = "norm", level = 0.95)  +
  scale_color_manual(values = c("#666666", "#D95F02", "#7570B3", "#E6AB02")) +
  theme(
    axis.text.x = element_text(family = "Times New Roman"),
    axis.title.x = element_text(family = "Times New Roman"),
    axis.text.y = element_text(family = "Times New Roman"),
    axis.title.y = element_text(family = "Times New Roman"),
    legend.title = element_text(family = "Times New Roman"),
    legend.text = element_text(family = "Times New Roman"))  

runoff_water_sediment_NMDS

ggsave("figure_6.tiff", plot = runoff_water_sediment_NMDS,
       width = 5, height = 5, dpi = 500, units = c("in"))

