library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(rstatix)
library(FSA)

raw_data <- read_excel("water_sediment_combined.xlsx")

manure_ARGs <- c("tetM", "tetQ", "tet44", "tnpA6", "tetT", "tetW", "ermF", "tet36", "aph3",
                 "ermQ", "aadE", "ermB", "tetO", "sat4", "erm35", "aadD", "tetL", "tetX",
                 "sul2", "tnpA5", "intl2", "ermT", "intl3", "sul1",
                 "intl1")

raw_data <- raw_data %>%
  filter(Primer %in% manure_ARGs)

#Table 2
output_detects <- raw_data %>%
  group_by(Treatment, Matrix, Primer) %>%
  summarise(n = n(),
            n_detects = sum(Relative_Abundance > 0))

output_median <- raw_data %>%
  group_by(Treatment, Matrix, Primer) %>%
  filter(Relative_Abundance > 0) %>%
  summarise(median = median(Relative_Abundance))

output <- merge(output_detects, output_median, by = c("Treatment", "Matrix", "Primer"))

write_xlsx(output, "resistance_genes_individual_output.xlsx")

#Stats - Dunn Test
water <- raw_data %>%
  subset(Matrix == "Water") %>%
  pivot_wider(names_from = Primer, values_from = Relative_Abundance)

water$Treatment <- as.factor(water$Treatment)

sediment <- raw_data %>%
  subset(Matrix == "Sediment") %>%
  pivot_wider(names_from = Primer, values_from = Relative_Abundance)

sediment$Treatment <- as.factor(sediment$Treatment)

lapply(water[,c("tetM", "tetQ", "tet44", "tnpA6", "tetT", "tetW", "ermF", "tet36", "aph3",
                "ermQ", "aadE", "ermB", "tetO", "sat4", "erm35", "aadD", "tetL", "tetX",
                "sul2", "tnpA5", "intl2", "ermT", "intl3", "sul1",
                "intl1")], function(x) dunnTest(x ~ Treatment, water, method = "bonferroni"))

lapply(sediment[,c("tetM", "tetQ", "tnpA6", "tetT", "tetW", "ermF", "tet36", "aph3",
                "ermQ", "aadE", "ermB", "tetO", "sat4", "erm35", "aadD", "tetX",
                "sul2", "tnpA5", "intl2", "ermT", "intl3", "sul1",
                "intl1")], function(x) dunnTest(x ~ Treatment, sediment, method = "bonferroni"))

#Individual check
dunn_test(water, tetM ~ Treatment, p.adjust.method = "bonferroni")
dunn_test(sediment, tetM ~ Treatment, p.adjust.method = "bonferroni")
