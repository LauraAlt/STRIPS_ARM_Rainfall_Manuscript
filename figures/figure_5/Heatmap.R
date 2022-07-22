library("readxl")
library("tidyverse")
library("ggh4x")
library("viridis")

raw_data <- read_excel("manure_water_sediment_combined.xlsx", sheet = "Data")

manure_ARGs <- c("tetM", "tetQ", "tet44", "tnpA6", "tetT", "tetW", "ermF", "tet36", "aph3",
                 "ermQ", "aadE", "ermB", "tetO", "sat4", "erm35", "aadD", "tetL", "tetX",
                 "sul2", "tnpA5", "intl2", "ermT", "intl3", "sul1",
                 "intl1")

raw_data1 <- subset(raw_data, Primer %in% manure_ARGs)

output_detects_count <- raw_data1 %>%
  group_by(Plot, Treatment, Matrix, Primer) %>%
  summarise(n = n(),
            n_detects = sum(Relative_Abundance > 0))

output_detects_median <- raw_data1 %>%
  group_by(Plot, Treatment, Matrix, Primer) %>%
  filter(Relative_Abundance > 0) %>%
  summarise(median = median(Relative_Abundance))

output_detects <- merge(output_detects_count, output_detects_median, by = c("Plot", "Treatment", "Matrix", "Primer"))

output_non_detects <- raw_data1 %>%
  group_by(Plot, Treatment, Matrix, Primer) %>%
  summarise(Sum_Relative_Abundance = sum(Relative_Abundance)) %>%
  ungroup() %>%
  filter(Sum_Relative_Abundance == 0) %>%
  mutate(n = 6) %>%
  mutate(n_detects = case_when(Primer == 0 ~ "Fill")) %>%
  mutate(median = 0) %>%
  within(rm(Sum_Relative_Abundance))

output <- rbind(output_detects, output_non_detects)

#Manure heatmap

manure_data <- subset(output, Matrix == "Manure")

manure_data$log.abundance <- log(manure_data$median)

manure_data <- manure_data %>%
  mutate(gene_class = case_when(Primer %in% c("aaC3", "aadD", "aadE", "aph3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm35") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfa",
                                Primer %in% c("tetL", "tetM", "tetO", "tetQ", "tetT", "tetW", "tetX",
                                              "tet33", "tet36", "tet44") ~ "Tetracycline",
                                Primer %in% c("intl1", "intl2", "intl3", "IS1247", "IS6100", "tnpA2",
                                              "tnpA5", "tnpA6") ~ "MGE"))

manure_data$gene_class_f = factor(manure_data$gene_class, 
                                levels = c('Aminoglycoside', 'MLSB', 'Sulfa', "Tetracycline", "MGE"))

levels = ordered(c('tnpA6', 'tnpA5', 'intl3',
                   'intl2', 'intl1', 'tet44', 'tet36', 'tetX', 'tetW', 'tetT', 'tetQ', 'tetO',
                   'tetM', 'tetL', 'sul2', 'sul1', 'erm35', 'ermT', 'ermQ', 'ermF', 'ermB', 'sat4', 'aph3', 
                   'aadE', 'aadD'))

manure_heatmap <- ggplot(manure_data, aes(as.factor(Plot), ordered(Primer, levels = levels))) +
  geom_tile(aes(fill = log.abundance), colour = "white", size = 0.1) +
  geom_text(aes(label = n_detects)) +
  scale_fill_viridis_c(limits = c(-11.5, -1),
                       name = "Natural Log of\nRelative Gene\nAbundance", option = "plasma")  +
  labs(y = "Gene") +
  facet_nested(gene_class_f ~ Matrix + Treatment, scales = "free", space = "free") + theme_classic() +
  theme(
    axis.text.x = element_text(
      family = "Times New Roman",
      angle = 90,
      size = 14
    ),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0.95, size = 14, family = "Times New Roman", face = 'italic'),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_rect(color = 'black', size = 1.4),
    strip.text.x = element_text(size = 12, family = "Times New Roman", face = 'bold'),
    strip.text.y = element_blank(),
    strip.background = element_rect(colour = 'black', size = 1.4)
  )  

manure_heatmap

ggsave("manure_heatmap.tiff", plot = manure_heatmap,
       width = 2.5, height = 10, dpi = 300, units = c("in"))

#Water heatmap

water <- subset(output, Matrix == "Runoff Water")

water$log.abundance <- log(water$median)

water$Treatment_f = factor(water$Treatment, 
                               levels = c('Strip + No Manure', 'Strip + Manure', 'No Strip + Manure'))

water <- water %>%
  mutate(gene_class = case_when(Primer %in% c("aaC3", "aadD", "aadE", "aph3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm35") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfa",
                                Primer %in% c("tetL", "tetM", "tetO", "tetQ", "tetT", "tetW", "tetX",
                                              "tet33", "tet36", "tet44") ~ "Tetracycline",
                                Primer %in% c("intl1", "intl2", "intl3", "IS1247", "IS6100", "tnpA2",
                                              "tnpA5", "tnpA6") ~ "MGE"))

water$gene_class_f = factor(water$gene_class, 
                             levels = c('Aminoglycoside', 'MLSB', 'Sulfa', "Tetracycline", "MGE"))

levels = ordered(c('tnpA6', 'tnpA5', 'intl3',
                   'intl2', 'intl1', 'tet44', 'tet36', 'tetX', 'tetW', 'tetT', 'tetQ', 'tetO',
                   'tetM', 'tetL', 'sul2', 'sul1', 'erm35', 'ermT', 'ermQ', 'ermF', 'ermB', 'sat4', 'aph3', 
                   'aadE', 'aadD'))

water_heatmap <- ggplot(water, aes(as.factor(Plot), ordered(Primer, levels = levels))) +
  geom_tile(aes(fill = log.abundance), colour = "white", size = 0.1) +
  geom_text(aes(label = n_detects)) +
  scale_fill_viridis_c(limits = c(-11.5, -1),
                       name = "Natural Log of\nRelative Gene\nAbundance", option = "plasma")  +
  labs(y = "Gene") +
  facet_nested(gene_class_f ~ Matrix + Treatment_f, scales = "free", space = "free") + theme_classic() +
  theme(
    axis.text.x = element_text(
      family = "Times New Roman",
      angle = 90,
      size = 14
    ),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.background = element_rect(color = 'black', size = 1.4),
    strip.text.x = element_text(size = 14, family = "Times New Roman", face = 'bold'),
    strip.text.y = element_blank(),
    strip.background = element_rect(colour = 'black', size = 1.4)
  )  

water_heatmap

ggsave("water_heatmap.tiff", plot = water_heatmap,
       width = 6, height = 10, dpi = 300, units = c("in"))

#Sediment heatmap

sediment <- subset(output, Matrix == "Runoff Sediment")

sediment$log.abundance <- log(sediment$median)

sediment$Treatment_f = factor(sediment$Treatment, 
                           levels = c('Strip + No Manure', 'Strip + Manure', 'No Strip + Manure'))

sediment <- sediment %>%
  mutate(gene_class = case_when(Primer %in% c("aaC3", "aadD", "aadE", "aph3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm35") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfa",
                                Primer %in% c("tetL", "tetM", "tetO", "tetQ", "tetT", "tetW", "tetX",
                                              "tet33", "tet36", "tet44") ~ "Tetracycline",
                                Primer %in% c("intl1", "intl2", "intl3", "IS1247", "IS6100", "tnpA2",
                                              "tnpA5", "tnpA6") ~ "MGE"))

sediment$gene_class_f = factor(sediment$gene_class, 
                            levels = c('Aminoglycoside', 'MLSB', 'Sulfa', "Tetracycline", "MGE"))

levels = ordered(c('tnpA6', 'tnpA5', 'intl3',
                   'intl2', 'intl1', 'tet44', 'tet36', 'tetX', 'tetW', 'tetT', 'tetQ', 'tetO',
                   'tetM', 'tetL', 'sul2', 'sul1', 'erm35', 'ermT', 'ermQ', 'ermF', 'ermB', 'sat4', 'aph3', 
                   'aadE', 'aadD'))

sediment_heatmap <- ggplot(sediment, aes(as.factor(Plot), ordered(Primer, levels = levels))) +
  geom_tile(aes(fill = log.abundance), colour = "white", size = 0.1) +
  geom_text(aes(label = n_detects)) +
  scale_fill_viridis_c(limits = c(-11.5, -1),
                       name = "Natural Log of\nMedian Relative\nAbundance", option = "plasma")  +
  labs(y = "Gene") +
  facet_nested(gene_class_f ~ Matrix + Treatment_f, scales = "free", space = "free") + theme_classic() +
  theme(
    axis.text.x = element_text(
      family = "Times New Roman",
      angle = 90,
      size = 14
    ),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(color = 'black', size = 1.4),
    strip.text.x = element_text(size = 14, family = "Times New Roman", face = 'bold'),
    strip.text.y = element_text(size = 12, family = "Times New Roman"),
    strip.background = element_rect(colour = 'black', size = 1.4),
    legend.title = element_text(size = 14, family = "Times New Roman"),
    legend.text = element_text(size = 14, family = "Times New Roman")
  )  

sediment_heatmap

ggsave("sediment_heatmap.tiff", plot = sediment_heatmap,
       width = 8, height = 10, dpi = 300, units = c("in"))
