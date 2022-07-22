library(readxl)
library(tidyverse)
library(writexl)
library(rstatix)

data <- read_excel("Armstrong_qPCR_processed.xlsx", sheet = "data_for_stats")

#Table values
output_detects <- data %>%
  group_by(Treatment, Matrix, Primer) %>%
  summarise(n = n(),
            n_detects = sum(Relative_Abundance > 0))

output_median <- data %>%
  group_by(Treatment, Matrix, Primer) %>%
  filter(Relative_Abundance > 0) %>%
  summarise(median = median(Relative_Abundance))

output <- merge(output_detects, output_median, by = c("Treatment", "Matrix", "Primer"))

write_xlsx(output, "resistance_genes_individual_output.xlsx")

#Stats - Dunn Test
water <- data %>%
  subset(Matrix == "Water") %>%
  pivot_wider(names_from = Primer, values_from = Relative_Abundance) %>%
  subset(select = -c(Plot, Time, Matrix, `intl1(clinical)`)) %>%
  rename(tnpA_02 = `tnpA-02`,
         tnpA_05 = `tnpA-05`,
         tnpA_06 = `tnpA-06`)

water_genes <-  colnames(water[2:27])

water_treatment <- names(water)[1]

lapply(water_genes, function(x) {dunn_test(water, reformulate(water_treatment, x), p.adjust.method = "bonferroni")})

sediment <- data %>%
  subset(Matrix == "Sediment") %>%
  pivot_wider(names_from = Primer, values_from = Relative_Abundance) %>%
  subset(select = -c(Plot, Time, Matrix, `tet(44)`, `tet(L)`, IS1247)) %>%
  rename(tnpA_02 = `tnpA-02`,
         tnpA_05 = `tnpA-05`,
         tnpA_06 = `tnpA-06`)

sediment_genes <-  colnames(sediment[2:25])

sediment_treatment <- names(sediment)[1]

lapply(sediment_genes, function(x) {dunn_test(sediment, reformulate(sediment_treatment, x), p.adjust.method = "bonferroni")})
