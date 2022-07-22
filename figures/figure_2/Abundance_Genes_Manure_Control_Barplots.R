library("readxl")
library("ggplot2")
library("tidyverse")

source("summarySE.R")

combined_data <- read_excel("manure_water_sediment_combined.xlsx", sheet = "Data")

manure_control <- combined_data %>%
  filter(Treatment == "Manure" | Treatment == "Strip + No Manure")


manure_control <- manure_control %>%
  mutate(Label = case_when(Matrix == "Manure" ~ "Manure",
                           Matrix == "Runoff Water" ~ "Control Runoff Water",
                           Matrix == "Runoff Sediment" ~ "Control Runoff Sediment"))

SE1 <- summarySE(manure_control, measurevar="Relative_Abundance", groupvars=c("Label", "Primer"))


ARG_bottom <- manure_control %>% filter(Primer %in% c("tetM", "sat4", "intl3", "tnpA2", "IS6100",
                                           "intl1", "tet33", "aaC3", "IS1247"))

SE2 <- summarySE(ARG_bottom, measurevar="Relative_Abundance", groupvars=c("Label", "Primer"))


levels1 = ordered(c("tetM", "tetQ", "tet44", "tnpA6", "tetT", "tetW", "ermF", "tet36", "aph3",
                    "ermQ", "aadE", "ermB", "tetO", "sat4", "erm35", "aadD", "tetL", "tetX",
                    "sul2", "tnpA5", "intl2", "ermT", "intl3", "tnpA2", "IS6100", "sul1",
                    "intl1", "tet33", "aaC3", "IS1247"))

levels2 = ordered(c("tetM", "sat4", "intl3", "tnpA2", "IS6100",
                    "intl1", "tet33", "aaC3", "IS1247"))

levels3 = ordered(c("Manure", "Control Runoff Water", "Control Runoff Sediment"))

all_ARG <- ggplot(SE1, aes(ordered(Primer, levels = levels1), Relative_Abundance,
                           fill = ordered(Label, levels = levels3))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Relative_Abundance - se, ymax = Relative_Abundance + se), width = 0.2,
                position=position_dodge(0.9)) +
  labs(y = "Average Relative Abundance") +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00")) +
  theme_light() +
  theme(axis.text.x = element_text(family = "Times New Roman", face = "italic", size = 14,
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(family = "Times New Roman", size = 14),
        axis.title.y = element_text(family = "Times New Roman", face = "bold", size = 16),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman", size = 14),
        legend.position = "top")

SE2$Appended_RA <- ifelse(SE2$Relative_Abundance > 0.0025, 0.0025, SE2$Relative_Abundance)

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

labels_data <- subset(SE2, Appended_RA > 0.0024)

labels <- round_df(labels_data, 3)

bottom_ARG <- ggplot(SE2, aes(ordered(Primer, levels = levels2), Appended_RA,
                              fill = ordered(Label, levels = levels3))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Relative_Abundance - se, ymax = Relative_Abundance + se), width = 0.2,
                position=position_dodge(0.9)) +
  geom_text(data = labels,
            aes(y = 0.0026, label = paste("(", Relative_Abundance, ")", sep = "")),
            colour = "black", family = "Times New Roman", size = 4, hjust = 1) +
  coord_cartesian(ylim = c(0, 0.0026)) +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#E69F00")) +
  scale_x_discrete(labels = c("tetM" = "tetM", "sat4" = "sat4", "intl3" = "intl3", "tnpA2" = "tnpA2*", "IS6100" = "IS6100*",
                              "intl1" = "intl1", "tet33" = "tet33*", "aaC3" = "aac3*", "IS1247" = "IS1247*")) +
  theme_light() +
  theme(axis.text.x = element_text(family = "Times New Roman", face = "italic", size = 18,
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(family = "Times New Roman", size = 16),
        axis.title.y = element_blank(),
        legend.position = "none")

all_ARG
bottom_ARG

ggsave("Abundance_Genes_Manure_Control.tiff", plot = all_ARG,
       width = 10, height = 6, dpi = 300, units = c("in"))

ggsave("Abundance_Genes_Manure_Control_Bottom.tiff", plot = bottom_ARG,
       width = 8, height = 6, dpi = 300, units = c("in"))
