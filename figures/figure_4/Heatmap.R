library(readxl)
library(tidyverse)
library(ggh4x)
library(viridis)

data <- read_excel("Armstrong_qPCR_processed.xlsx", sheet = "Sheet1")

output_detects_count <- data %>%
  group_by(Plot, Treatment, Matrix, Primer) %>%
  summarise(n = n(),
            n_detects = sum(Relative_Abundance > 0))

output_detects_median <- data %>%
  group_by(Plot, Treatment, Matrix, Primer) %>%
  filter(Relative_Abundance > 0) %>%
  summarise(median = median(Relative_Abundance))

output_detects <- merge(output_detects_count, output_detects_median, by = c("Plot", "Treatment", "Matrix", "Primer"))

output_non_detects <- data %>%
  group_by(Plot, Treatment, Matrix, Primer) %>%
  summarise(Sum_Relative_Abundance = sum(Relative_Abundance)) %>%
  ungroup() %>%
  filter(Sum_Relative_Abundance == 0) %>%
  mutate(n = 6) %>%
  mutate(n_detects = case_when(Primer == 0 ~ "Fill")) %>%
  mutate(median = 0) %>%
  within(rm(Sum_Relative_Abundance))

output <- rbind(output_detects, output_non_detects) %>%
  mutate(primer_ed = case_when(Primer == "intl1(clinical)" ~ "*intl1(clinical)",
                               Primer == "IS1247" ~ "*IS1247",
                               Primer == "IS6100" ~ "*IS6100",
                               Primer == "tnpA-02" ~ "*tnpA-02",
                               TRUE ~ Primer)) %>%
  mutate(matrix_ed = case_when(Matrix == "Water" ~ "Runoff Water",
                               Matrix == "Sediment" ~ "Runoff Sediment",
                               TRUE ~ Matrix)) %>%
  mutate(gene_class = case_when(Primer %in% c("aadD", "aadE", "aphA3", "sat4") ~ "AMG",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm(35)") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfa",
                                Primer %in% c("tet(L)", "tet(M)", "tet(O)", "tet(Q)", "tet(T)", "tet(W)", "tet(X)",
                                              "tet(36)", "tet(44)") ~ "Tetracycline",
                                Primer %in% c("intl2", "tnpA-05", "tnpA-06") ~ "MGE",
                                Primer %in% c("intl1(clinical)", "IS1247", "IS6100", "tnpA-02") ~ "*MGE")) %>%
  mutate(plot_place_holder = case_when(Plot == "NA" ~ " ",
                                  TRUE ~ Plot)) %>%
  mutate(treatment_place_holder = case_when(Treatment == "Manure" ~ " ",
                                   TRUE ~ Treatment))

output$gene_class_f = factor(output$gene_class, levels = c("AMG", "MLSB", "Sulfa", "Tetracycline", "MGE", "*MGE"))

output$matrix_f = factor(output$matrix_ed, levels = c ("Manure", "Runoff Water", "Runoff Sediment"))

output$treatment_f = factor(output$treatment_place_holder, levels = c(" ", "Strip + No Manure", "Strip + Manure", "No Strip + Manure"))

output$log.abundance <- log(output$median)

levels = ordered(c("*tnpA-02", "*IS6100", "*IS1247", "*intl1(clinical)",
                   "tnpA-06", "tnpA-05", "intl2",
                   "tet(X)", "tet(W)", "tet(T)", "tet(Q)", "tet(O)", "tet(M)", "tet(L)", "tet(44)", "tet(36)",
                   "sul2", "sul1",
                   "ermT", "ermQ", "ermF", "ermB", "erm(35)",
                   "sat4", "aphA3", "aadE", "aadD"))

#Heatmap

heatmap <- ggplot(output, aes(as.factor(plot_place_holder), ordered(primer_ed, levels = levels))) +
  geom_tile(aes(fill = log.abundance), colour = "white", size = 0.1) +
  geom_text(aes(label = n_detects), size = 4, family = "Times New Roman") +
  scale_fill_viridis_c(name = "Natural Log of Relative Gene Abundance", option = "plasma")  +
  facet_nested(gene_class_f ~ matrix_f + treatment_f, scales = "free", space = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, family = "Times New Roman", angle = 90, colour = "black"),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_text(size = 11, family = "Times New Roman", face = "bold", colour = "black"),
    axis.text.y = element_text(hjust = 0.95, size = 12, family = "Times New Roman", face = "italic", colour = "black"),
    axis.title.y = element_blank(),
    strip.text.y = element_text(size = 11, family = "Times New Roman", colour = "black"),
    legend.title = element_text(size = 12, family = "Times New Roman", face = "bold", colour = "black"),
    legend.text = element_text(size = 12, family = "Times New Roman", colour = "black"),
    legend.position = "top",
    panel.background = element_rect(colour = "black", size = 1.4),
    strip.background = element_rect(colour = "black", size = 1.4)) +
  guides(fill = guide_colourbar(title.position = "left", title.vjust = 0.75))

heatmap

ggsave("figure_4.tiff", plot = heatmap,
       width = 12.5, height = 7.5, dpi = 500, units = c("in"))
