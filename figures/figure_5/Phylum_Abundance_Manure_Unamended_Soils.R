library(tidyverse)
library(phyloseq)
library(phylosmith)
library(writexl)
library(readxl)
library(VennDiagram)
library(cowplot)

phy <- readRDS("Armstrong_Rainfall_Phy.RDS")

## BAR PLOT

# Subset to Manure and Soil

soil_manure <- subset_samples(phy, matrix == "manure" | matrix == "soil")

# Subset to Manure and Non-manured Soil

unamended_soil_manure <- subset_samples(soil_manure, manure_amendment == "N" | manure_amendment == "manure")

# Calculate Relative Abundances

unamended_soil_manure_rel <- unamended_soil_manure %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>% 
  transform_sample_counts(function(x) x / sum(x))

# Average across samples and phyla

conglom_unamended_soil_manure_rel <- conglomerate_samples(unamended_soil_manure_rel, c("soil_type"))
conglom_unamended_soil_manure_rel <- conglomerate_taxa(conglom_unamended_soil_manure_rel, "Phylum", FALSE)

melted <- psmelt(conglom_unamended_soil_manure_rel)

write_xlsx(x = melted, 
           path = "unamended_manure_rel.xlsx")

# Use unamended_manure_rel.xlsx spreadsheet to create graph_data spreadsheet

graph_data <- read_excel("unamended_manure_rel_editted.xlsx", sheet = "Combined")



levels = ordered(c("Manure", "Crop", "Interface", 
                   "Strip"))

phylum_abundance_plot <- ggplot(graph_data, aes(ordered(Label, levels = levels), y = Abundance)) +
  geom_bar(aes(fill = Phylum), stat = "identity") +
  labs(y = "Average Relative Abundance") +
  scale_fill_manual(values = c("999933", "#f97b72", "#3969AC", "#CC6677", 
                               "#F2B701", "darkgreen", "#AA4499", "#44AA99", 
                               "#88CCEE", "#DDCC77", "#4b4b8f", "#6699CC", 
                               "darkorange4", "royalblue4", "#7F3C8D")) +
  scale_x_discrete(labels = c("Manure", "Control\nCrop\nSoil", "Control\nInterface\nSoil",
                              "Control\nStrip\nSoil")) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(family = "Times New Roman", size = 12, color = "black"),
        axis.title.y = element_text(family = "Times New Roman", size = 12, face = "bold"),
        axis.text.y = element_text(family = "Times New Roman", size = 12, color = "black"),
        legend.title = element_text(family = "Times New Roman", size = 12, face = "bold"),
        legend.text = element_text(family = "Times New Roman", size = 12),
        panel.background = element_rect(colour = "black", linetype = "solid", size = 0.75))

phylum_abundance_plot

## VENN DIAGRAM

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

manure_ASVs <- subset_samples(unamended_soil_manure, matrix == "manure") %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>%
  otu_table() %>%
  rownames()

unamended_crop_soil_ASVs <- subset_samples(unamended_soil_manure, matrix == "soil" & soil_type == "Crop") %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>%
  otu_table() %>%
  rownames()

unamended_strip_soil_ASVs <- subset_samples(unamended_soil_manure, matrix == "soil" & soil_type == "Strip") %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>%
  otu_table() %>%
  rownames()

unamended_interface_soil_ASVs <- subset_samples(unamended_soil_manure, matrix == "soil" & soil_type == "Interface") %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>%
  otu_table() %>%
  rownames()

display_venn(x = list(manure_ASVs, unamended_crop_soil_ASVs, unamended_interface_soil_ASVs, unamended_strip_soil_ASVs),
             category.names = c("Manure\nASVs", "Control Crop\nSoil ASVs", 
                                "Control Interface\nSoil ASVs", "Control Strip\nSoil ASVs"),
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             cex = 1.1,
             cat.cex = 1.1)

venn_diagram <- venn.diagram(x = list(manure_ASVs, unamended_crop_soil_ASVs, unamended_interface_soil_ASVs, unamended_strip_soil_ASVs),
             category.names = c("Manure\nASVs", "Control Crop\nSoil ASVs", 
                                "Control Interface\nSoil ASVs", "Control Strip\nSoil ASVs"),
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             cex = 1.1,
             cat.cex = 1.1,
             filename = NULL)
  
combined_plot <- ggdraw() +
  draw_plot(phylum_abundance_plot, x = 0.005, y = 0, width = 0.5, height = 1) +
  draw_plot(venn_diagram, x = 0.5, y = 0, width = 0.495, height = 1) +
  draw_plot_label(label = c("a", "b"), size = 22, label_fontfamily = "Times New Roman", x = c(0, 0.5), y = c(1, 1))

combined_plot

ggsave("figure_5.tiff", bg = "white",
       plot = combined_plot, width = 12, height = 6, dpi = 500, units = c("in"))

