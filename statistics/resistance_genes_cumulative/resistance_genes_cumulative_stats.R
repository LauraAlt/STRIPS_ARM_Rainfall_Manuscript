library(readxl)
library(tidyverse)
library(ggpubr)
library(rstatix)

data <- read_excel("Armstrong_qPCR_processed.xlsx")

water <- subset(data, Matrix == "Water")
sediment <- subset(data, Matrix == "Sediment")

#Cumulative Abundance
cumulative_abundance_water <- water %>%
  group_by(Plot, Time, Treatment) %>% 
  summarise(total_sample_relative_abundance = sum(Relative_Abundance)) %>%
  ungroup() %>%
  reorder_levels(Treatment, order = c("Strip + No Manure", "Strip + Manure", "No Strip + Manure"))
cumulative_abundance_sediment <- sediment %>%
  group_by(Plot, Time, Treatment) %>% 
  summarise(total_sample_relative_abundance = sum(Relative_Abundance)) %>%
  ungroup() %>%
  reorder_levels(Treatment, order = c("Strip + No Manure", "Strip + Manure", "No Strip + Manure"))

#Kruskal-Wallis
kruskal.test(total_sample_relative_abundance ~ Treatment, data = cumulative_abundance_water)
kruskal.test(total_sample_relative_abundance ~ Treatment, data = cumulative_abundance_sediment)

#Dunn Test (Option 1)
#Code from https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
res.kruskal_water <- cumulative_abundance_water %>% 
  kruskal_test(total_sample_relative_abundance ~ Treatment)
res.kruskal_water

pwc_water <- cumulative_abundance_water %>% dunn_test(total_sample_relative_abundance ~ Treatment,
                                       p.adjust.method = "bonferroni") 
pwc_water

pwc_water <- pwc_water %>% add_xy_position(x = "Treatment")
plot_water <- ggboxplot(cumulative_abundance_water, x = "Treatment", y = "total_sample_relative_abundance") +
  stat_pvalue_manual(pwc_water, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal_water, detailed = TRUE),
    caption = get_pwc_label(pwc_water)
  )

plot_water

res.kruskal_sediment <- cumulative_abundance_sediment %>% 
  kruskal_test(total_sample_relative_abundance ~ Treatment)
res.kruskal_sediment

pwc_sediment <- cumulative_abundance_sediment %>% dunn_test(total_sample_relative_abundance ~ Treatment,
                                       p.adjust.method = "bonferroni") 
pwc_sediment

pwc_sediment <- pwc_sediment %>% add_xy_position(x = "Treatment")
plot_sediment <- ggboxplot(cumulative_abundance_sediment, x = "Treatment", y = "total_sample_relative_abundance") +
  stat_pvalue_manual(pwc_sediment, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal_sediment, detailed = TRUE),
    caption = get_pwc_label(pwc_sediment)
  )

plot_sediment
