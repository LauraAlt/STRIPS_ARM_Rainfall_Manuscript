library(dplyr)
library(ggplot2)
library(phyloseq)
library(phylosmith)
library(writexl)
library(readxl)
library(VennDiagram)
theme_set(theme_bw())

phy <- readRDS("Armstrong_Rainfall_Phy.RDS")

## BAR PLOT

# Subset to Manure and Soil

soil_manure <- subset_samples(phy, matrix == "manure" | matrix == "soil")

# Subset to Manure and Unamended Soil

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



levels = ordered(c('Manure', 'Crop', 'Interface', 
                   'Strip'))

plot <- ggplot(graph_data, aes(ordered(Label, levels = levels), y = Abundance)) +
  geom_bar(aes(fill = Phylum), stat = "identity") +
  labs(y = "Average Relative Abundance") +
  scale_fill_manual(values = c("999933", "#f97b72", "#3969AC", "#CC6677", 
                               "#F2B701", "darkgreen", "#AA4499", "#44AA99", 
                               "#88CCEE", "#DDCC77", "#4b4b8f", "#6699CC", 
                               "darkorange4", "royalblue4", "#7F3C8D")) +
  scale_x_discrete(labels = c("Manure", "Control\nCrop\nSoil", "Control\nInterface\nSoil",
                              "Control\nStrip\nSoil")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(family = "Times New Roman", size = 20),
        axis.title.y = element_text(family = "Times New Roman", size = 24, face = "bold"),
        axis.text.y = element_text(family = "Times New Roman", size = 20),
        legend.title = element_text(family = "Times New Roman", size = 24, face = "bold"),
        legend.text = element_text(family = "Times New Roman", size = 20))

plot

ggsave("Phylum_Abundance_Manure_Unamended.tiff", 
       plot = plot, width = 10, height = 10, dpi = 300, units = c("in"))

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

venn.diagram(x = list(manure_ASVs, unamended_crop_soil_ASVs, unamended_interface_soil_ASVs, unamended_strip_soil_ASVs),
             category.names = c("Manure\nASVs", "Control Crop\nSoil ASVs", 
                                "Control Interface\nSoil ASVs", "Control Strip\nSoil ASVs"),
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
             cex = 1.1,
             cat.cex = 1.1,
             filename = "Manure_Associated_Bacteria_Venn_Diagram.png")
