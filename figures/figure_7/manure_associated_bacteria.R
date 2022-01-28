library(dplyr)
library(ggplot2)
library(phyloseq)
library(phylosmith)
library(writexl)
theme_set(theme_bw())

phy <- readRDS("Armstrong_Rainfall_Phy.RDS")

# Add labels for final plot
all_data <- data.frame(sample_data(phy))

all_data <- all_data %>%
  mutate(day = case_when(sample_day %in% c("Baseline") ~ "Baseline",
                         sample_day %in% c("T000") ~ "Day 0",
                         sample_day %in% c("T002") ~ "Day 2",
                         sample_day %in% c("T014") ~ "Day 14",
                         sample_day %in% c("T153") ~ "Day 153")) %>%
  mutate(description = case_when(treatment %in% c("ACM") ~ "No Strip + Manure",
                                 treatment %in% c("ACS") ~ "Strip + No Manure",
                                 treatment %in% c("ACSM") ~ "Strip + Manure"))

sample_names(phy)

rownames(all_data) <- all_data$sample_name

sample_data(phy) <- all_data

# Identify Manure Associated Bacteria
soil_manure <- subset_samples(phy, matrix == "soil" | matrix == "manure")
unique(sample_data(soil_manure)$matrix)
unique(sample_data(soil_manure)$manure_amendment)

unamended_soil_manure <- subset_samples(soil_manure, manure_amendment == "N" | manure_amendment == "manure")

unamended_soil_manure0 <- prune_taxa(taxa_sums(unamended_soil_manure) > 0, unamended_soil_manure)
any(taxa_sums(unamended_soil_manure0) == 0)

unique_taxa <- unique_taxa(unamended_soil_manure0, treatment = c("matrix"))
sapply(unique_taxa, length)

# Strip + Manure Data Processing
ACSM_soil_rel <- subset_samples(phy, matrix == "soil" & treatment == "ACSM" & day != "Baseline") %>%
  filter_taxa(function(x) sum(x) >= 1, T) %>% 
  transform_sample_counts(function(x) x / sum(x))

persistors <- prune_taxa(unique_taxa$manure, ACSM_soil_rel)

ntaxa(persistors)

conglom_persistors <- conglomerate_samples(persistors, c("sample_num", "day", "depth"))

conglom_persistors_phyla <- conglomerate_taxa(conglom_persistors, "Phylum")
conglom_persistors_genus <- conglomerate_taxa(conglom_persistors, "Genus")

melted_phyla <- psmelt(conglom_persistors_phyla)
melted_genus <- psmelt(conglom_persistors_genus)

write_xlsx(x = melted_phyla, 
           path = "MAB_Phyla_ACSM.xlsx")
write_xlsx(x = melted_genus, 
           path = "MAB_Genus_ACSM.xlsx")

levels = ordered(c('Day 0', 'Day 2', 'Day 14', 'Day 153'))
levels1 = ordered(c('Depth 1', 'Depth 2'))

# ACSM Graph
g_ACSM <- ggplot(melted_phyla, aes(x = as.factor(sample_num), y = Abundance*100))
plot_ACSM <- g_ACSM + 
  geom_bar(aes(fill = Phylum), stat = "identity") +
  facet_grid(ordered(depth, levels = levels1) ~ ordered(day, levels = levels)) +
  labs(x = "Plot Sample Location", y = "Percentage of Bacterial Community\nOriginating from Manure") +
  ylim(0, 27) +
  scale_fill_manual(values = c("#f97b72", "#3969AC", "#F2B701", "#E73F74",
                               "#E68310", "darkgreen", "#CF1C90", 
                               "#11A579", "#4b4b8f", "darkorange4", "cyan1", "royalblue4",
                               "#A5AA99", "#7F3C8D")) +
  theme(panel.background = element_rect(color = 'black', size = 1.4),
        strip.text.x = element_text(family = "Times New Roman", size = 20, face = 'bold'),
        strip.background = element_rect(colour = 'black', size = 1.4, fill = 'white'),
        strip.text.y = element_text(family = "Times New Roman", size = 20, face = 'bold'),
        legend.title = element_text(family = "Times New Roman", size = 20, face = 'bold'),
        legend.text = element_text(family = "Times New Roman", size = 20),
        axis.title.x = element_text(family = "Times New Roman", size = 20, face = 'bold'),
        axis.text.y = element_text(family = "Times New Roman", size = 20),
        axis.title.y = element_text(family = "Times New Roman", size = 20, face = 'bold'),
        axis.text.x = element_text(family = "Times New Roman", size = 20)
  )

plot_ACSM

ggsave("manure_associated_bacteria_ACSM.tiff", 
       plot = plot_ACSM, width = 20, height = 10, dpi = 100, units = c("in"))
